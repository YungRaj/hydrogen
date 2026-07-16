#!/usr/bin/env python3
"""
Multi-Task PyTorch Surrogate Model for Catalyst Property Prediction.

Trained on MACE screening data to accelerate the genetic algorithm by
predicting catalyst descriptors ~1000× faster than full GNN evaluation.

Predicts:
  1. Validity (binary classification)
  2. dE_split (C-H splitting reaction energy)
  3. Coking resistance index
  4. Segregation / binding stability energy
  5. E_act (activation barrier)
"""

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from typing import Optional, Tuple
import logging

from pipeline.catalyst_spaces import FEATURE_DIM

logger = logging.getLogger('surrogate_model')


class CatalystSurrogate(nn.Module):
    """
    Multi-task neural network for catalyst property prediction.
    
    Architecture: Shared feature extractor → task-specific heads
    """

    def __init__(self, input_dim: int = FEATURE_DIM, hidden_dims: tuple = (512, 256, 128)):
        super().__init__()

        layers = []
        prev_dim = input_dim
        for h_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, h_dim),
                nn.BatchNorm1d(h_dim),
                nn.GELU(),
                nn.Dropout(0.1),
            ])
            prev_dim = h_dim

        self.backbone = nn.Sequential(*layers)

        # Task-specific heads
        self.head_valid = nn.Sequential(
            nn.Linear(prev_dim, 32), nn.GELU(), nn.Linear(32, 1)
        )
        self.head_de_split = nn.Sequential(
            nn.Linear(prev_dim, 32), nn.GELU(), nn.Linear(32, 1)
        )
        self.head_coking = nn.Sequential(
            nn.Linear(prev_dim, 32), nn.GELU(), nn.Linear(32, 1)
        )
        self.head_seg = nn.Sequential(
            nn.Linear(prev_dim, 32), nn.GELU(), nn.Linear(32, 1)
        )
        self.head_e_act = nn.Sequential(
            nn.Linear(prev_dim, 32), nn.GELU(), nn.Linear(32, 1)
        )

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, ...]:
        features = self.backbone(x)
        valid_logit = self.head_valid(features)
        de_split = self.head_de_split(features)
        coking = self.head_coking(features)
        seg = self.head_seg(features)
        e_act = self.head_e_act(features)
        return valid_logit, de_split, coking, seg, e_act


def train_surrogate(X: np.ndarray, y_valid: np.ndarray,
                    y_de_split: np.ndarray, y_coking: np.ndarray,
                    y_seg: np.ndarray, y_e_act: np.ndarray,
                    epochs: int = 30, batch_size: int = 2048,
                    lr: float = 0.003, device: str = 'cuda:0') -> CatalystSurrogate:
    """
    Train the surrogate model on MACE screening data.
    
    Args:
        X: Feature matrix (N, FEATURE_DIM)
        y_valid: Validity labels (N,) — 1.0 if MACE succeeded
        y_de_split: C-H splitting energy (N,)
        y_coking: Coking resistance index (N,)
        y_seg: Segregation/binding energy (N,)
        y_e_act: Activation barrier (N,)
        
    Returns:
        Trained CatalystSurrogate model
    """
    model = CatalystSurrogate(input_dim=X.shape[1]).to(device)

    # Convert to tensors
    X_t = torch.tensor(X, dtype=torch.float32).to(device)
    y_val_t = torch.tensor(y_valid, dtype=torch.float32).unsqueeze(1).to(device)
    y_de_t = torch.tensor(y_de_split, dtype=torch.float32).unsqueeze(1).to(device)
    y_cok_t = torch.tensor(y_coking, dtype=torch.float32).unsqueeze(1).to(device)
    y_seg_t = torch.tensor(y_seg, dtype=torch.float32).unsqueeze(1).to(device)
    y_act_t = torch.tensor(y_e_act, dtype=torch.float32).unsqueeze(1).to(device)

    dataset = TensorDataset(X_t, y_val_t, y_de_t, y_cok_t, y_seg_t, y_act_t)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)
    bce_loss = nn.BCEWithLogitsLoss()
    mse_loss = nn.MSELoss()

    model.train()
    logger.info(f"Training surrogate on {len(X)} samples for {epochs} epochs...")

    for epoch in range(epochs):
        total_loss = 0.0
        for batch in loader:
            bX, bV, bD, bC, bS, bA = batch
            optimizer.zero_grad()

            valid_logit, de_split, coking, seg, e_act = model(bX)

            # Classification loss for validity
            loss_v = bce_loss(valid_logit, bV)

            # Regression losses only for valid candidates
            mask = (bV > 0.5).squeeze()
            if mask.sum() > 0:
                loss_d = mse_loss(de_split[mask], bD[mask])
                loss_c = mse_loss(coking[mask], bC[mask])
                loss_s = mse_loss(seg[mask], bS[mask])
                loss_a = mse_loss(e_act[mask], bA[mask])
                loss = loss_v + 2.0 * (loss_d + loss_c + loss_s + loss_a)
            else:
                loss = loss_v

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            total_loss += loss.item()

        scheduler.step()

        if (epoch + 1) % 10 == 0 or epoch == 0:
            logger.info(f"  Epoch {epoch+1}/{epochs}: loss = {total_loss/len(loader):.4f}")

    model.eval()
    logger.info("Surrogate training complete.")
    return model


@torch.no_grad()
def predict_batch(model: CatalystSurrogate, X: np.ndarray,
                  device: str = 'cuda:0') -> dict:
    """
    Predict catalyst properties for a batch of feature vectors.
    
    Returns dict with numpy arrays for each property.
    """
    model.eval()
    X_t = torch.tensor(X, dtype=torch.float32).to(device)

    valid_logit, de_split, coking, seg, e_act = model(X_t)

    return {
        'valid_prob': torch.sigmoid(valid_logit).cpu().numpy().flatten(),
        'de_split': de_split.cpu().numpy().flatten(),
        'coking_index': coking.cpu().numpy().flatten(),
        'segregation_energy': seg.cpu().numpy().flatten(),
        'E_act': e_act.cpu().numpy().flatten(),
    }
