"""Build a provenance-gated candidate qubit Hamiltonian from FCIDUMP integrals."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _metadata(path: str | Path) -> tuple[Path, dict]:
    source = Path(path).expanduser().resolve()
    sidecar = source.with_suffix(source.suffix + '.json')
    if not source.is_file() or not sidecar.is_file():
        raise ValueError('FCIDUMP and FCIDUMP.json metadata sidecar are required')
    data = json.loads(sidecar.read_text())
    required = (
        'schema_version', 'candidate_id', 'geometry_id', 'geometry_path',
        'geometry_sha256',
        'electronic_structure_protocol', 'basis', 'active_orbitals',
        'active_electrons', 'charge', 'multiplicity', 'frozen_orbitals',
        'integral_source', 'fermion_to_qubit_mapping', 'fcidump_sha256',
    )
    missing = [key for key in required if data.get(key) in (None, '', [])]
    if missing:
        raise ValueError(f'candidate Hamiltonian metadata missing: {missing}')
    if data['schema_version'] != 1:
        raise ValueError('candidate Hamiltonian schema_version must be 1')
    if data['fermion_to_qubit_mapping'] != 'Jordan-Wigner':
        raise ValueError('only the declared Jordan-Wigner mapping is supported')
    if _sha256(source) != str(data['fcidump_sha256']).lower():
        raise ValueError('FCIDUMP checksum mismatch')
    geometry = Path(data['geometry_path']).expanduser()
    if not geometry.is_absolute():
        geometry = sidecar.parent / geometry
    if not geometry.is_file() or _sha256(geometry) != str(data['geometry_sha256']).lower():
        raise ValueError('geometry source is missing or its checksum mismatches')
    return source, data


def build_candidate_hamiltonian(path: str | Path) -> dict:
    """Convert sourced spatial-orbital integrals into Pauli terms.

    PySCF reads the standard FCIDUMP representation and OpenFermion performs
    the spatial-to-spin-orbital and Jordan-Wigner transformations.
    """
    source, metadata = _metadata(path)
    try:
        from pyscf import ao2mo
        from pyscf.tools import fcidump
        from openfermion.chem.molecular_data import spinorb_from_spatial
        from openfermion.ops.representations import InteractionOperator
        from openfermion.transforms import jordan_wigner
    except ImportError as exc:
        raise RuntimeError(
            'candidate Hamiltonian construction requires PySCF and OpenFermion') from exc

    values = fcidump.read(str(source), verbose=False)
    norb, nelec, ms2 = (int(values['NORB']), int(values['NELEC']),
                         int(values.get('MS2', 0)))
    if norb != int(metadata['active_orbitals']):
        raise ValueError('FCIDUMP NORB does not match declared active_orbitals')
    if nelec != int(metadata['active_electrons']):
        raise ValueError('FCIDUMP NELEC does not match declared active_electrons')
    if int(metadata['multiplicity']) != abs(ms2) + 1:
        raise ValueError('FCIDUMP MS2 does not match declared multiplicity')
    if not 0 < nelec <= 2 * norb:
        raise ValueError('active electron count is outside the spin-orbital space')

    spatial_two_body = ao2mo.restore(1, values['H2'], norb)
    one_body, two_body = spinorb_from_spatial(values['H1'], spatial_two_body)
    interaction = InteractionOperator(
        float(values.get('ECORE', 0.0)), one_body, 0.5 * two_body)
    qubit_operator = jordan_wigner(interaction)
    terms = []
    for operators, coefficient in qubit_operator.terms.items():
        if abs(complex(coefficient).imag) > 1e-10:
            raise ValueError('Jordan-Wigner Hamiltonian is not Hermitian')
        word = ['I'] * (2 * norb)
        for index, symbol in operators:
            word[index] = symbol
        terms.append((float(complex(coefficient).real), ''.join(word)))
    terms.sort(key=lambda item: item[1])
    return {
        'schema_version': 1,
        'candidate_id': metadata['candidate_id'],
        'geometry_id': metadata['geometry_id'],
        'geometry_path': metadata['geometry_path'],
        'geometry_sha256': metadata['geometry_sha256'],
        'electronic_structure_protocol': metadata['electronic_structure_protocol'],
        'basis': metadata['basis'],
        'active_orbitals': norb,
        'active_electrons': nelec,
        'charge': int(metadata['charge']),
        'multiplicity': int(metadata['multiplicity']),
        'frozen_orbitals': metadata['frozen_orbitals'],
        'integral_source': metadata['integral_source'],
        'fcidump_sha256': metadata['fcidump_sha256'],
        'fermion_to_qubit_mapping': 'Jordan-Wigner',
        'n_qubits': 2 * norb,
        'pauli_terms': terms,
        'catalyst_specific_hamiltonian': True,
        'evidence_level': 'candidate_specific_integrals',
    }


def write_candidate_hamiltonian(path: str | Path, output: str | Path) -> dict:
    result = build_candidate_hamiltonian(path)
    Path(output).write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fcidump')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = write_candidate_hamiltonian(args.fcidump, args.output)
    print(json.dumps({key: result[key] for key in (
        'candidate_id', 'n_qubits', 'active_electrons', 'evidence_level')},
        indent=2))


if __name__ == '__main__':
    main()
