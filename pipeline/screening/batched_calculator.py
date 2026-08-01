"""Thread-safe dynamic batching for fairchem ASE calculations.

One model is resident on each GPU.  Independent ASE optimizers submit energy
and force requests to a small batching thread, which combines the native graph
objects into one fairchem inference call and splits the predictions back out.
"""

from __future__ import annotations

import queue
import threading
import time
from dataclasses import dataclass, field

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.stress import full_3x3_to_voigt_6_stress


@dataclass
class _Request:
    atoms: object
    event: threading.Event = field(default_factory=threading.Event)
    result: dict | None = None
    error: BaseException | None = None


class BatchedInferenceService:
    """Own a fairchem model and dynamically batch requests by total atoms."""

    def __init__(self, calculator, max_batch_atoms: int = 768,
                 batch_wait_ms: float = 3.0):
        if max_batch_atoms < 1:
            raise ValueError('max_batch_atoms must be positive')
        if batch_wait_ms < 0:
            raise ValueError('batch_wait_ms cannot be negative')
        self.calculator = calculator
        self.max_batch_atoms = max_batch_atoms
        self.batch_wait_s = batch_wait_ms / 1000.0
        self._requests: queue.Queue = queue.Queue()
        self._stop = object()
        self._thread = threading.Thread(
            target=self._serve, name='fairchem-batch-server', daemon=True)
        self._thread.start()

    def calculator_proxy(self) -> Calculator:
        return BatchedFAIRChemCalculator(self)

    def predict(self, atoms) -> dict:
        request = _Request(atoms.copy())
        self._requests.put(request)
        request.event.wait()
        if request.error is not None:
            raise request.error
        return request.result

    def close(self):
        if self._thread.is_alive():
            self._requests.put(self._stop)
            self._thread.join()

    def _serve(self):
        deferred = None
        while True:
            first = deferred if deferred is not None else self._requests.get()
            deferred = None
            if first is self._stop:
                return
            batch = [first]
            atoms_in_batch = len(first.atoms)
            deadline = time.monotonic() + self.batch_wait_s
            while atoms_in_batch < self.max_batch_atoms:
                timeout = deadline - time.monotonic()
                if timeout <= 0:
                    break
                try:
                    item = self._requests.get(timeout=timeout)
                except queue.Empty:
                    break
                if item is self._stop:
                    deferred = item
                    break
                if atoms_in_batch + len(item.atoms) > self.max_batch_atoms:
                    deferred = item
                    break
                batch.append(item)
                atoms_in_batch += len(item.atoms)
            self._run_batch(batch)

    def _run_batch(self, requests):
        try:
            from fairchem.core.datasets.atomic_data import atomicdata_list_to_batch

            data = []
            for request in requests:
                atoms = request.atoms
                if len(atoms) == 0:
                    raise ValueError('Atoms object has no atoms inside.')
                self.calculator._check_atoms_pbc(atoms)
                self.calculator.predictor.validate_atoms_data(
                    atoms, self.calculator.task_name)
                data.append(self.calculator.a2g(atoms))

            predictions = self.calculator.predictor.predict(
                atomicdata_list_to_batch(data))
            n_systems = len(requests)
            atom_offset = 0
            for index, request in enumerate(requests):
                n_atoms = len(request.atoms)
                result = {}
                for key in self.calculator.implemented_properties:
                    if key == 'free_energy':
                        continue
                    pred = predictions.get(key)
                    if pred is None:
                        continue
                    array = pred.detach().cpu().numpy()
                    if key == 'energy':
                        energy = float(array[index])
                        result['energy'] = result['free_energy'] = energy
                    elif key == 'forces':
                        result['forces'] = array[atom_offset:atom_offset + n_atoms]
                    elif key == 'stress':
                        stress = array[index] if array.shape[0] == n_systems else array
                        result['stress'] = full_3x3_to_voigt_6_stress(
                            np.asarray(stress).reshape(3, 3))
                atom_offset += n_atoms
                request.result = result
        except BaseException as exc:
            for request in requests:
                request.error = exc
        finally:
            for request in requests:
                request.event.set()


class BatchedFAIRChemCalculator(Calculator):
    """ASE calculator proxy backed by a shared :class:`BatchedInferenceService`."""

    def __init__(self, service: BatchedInferenceService):
        super().__init__()
        self.service = service
        self.implemented_properties = list(service.calculator.implemented_properties)

    def calculate(self, atoms=None, properties=('energy',),
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        self.results = self.service.predict(atoms)
