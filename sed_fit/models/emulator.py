"""SPS Emulator MLP -- forward model: physical params -> predicted magnitudes.

Architecture (~100K params depending on n_bands):
    Input(6)  [z, log_mass, log_age, log_Z, Av, log_tau]
    -> FC(6, 256) -> ReLU
    -> FC(256, 256) -> ReLU
    -> FC(256, 128) -> ReLU
    -> FC(128, n_bands)  # predicted apparent magnitudes
"""

import torch
import torch.nn as nn


class SPSEmulatorMLP(nn.Module):
    """
    Stellar Population Synthesis emulator.

    Maps physical parameters (redshift, stellar mass, age, metallicity,
    dust extinction, SFH timescale) to broadband magnitudes.  Trained on
    a grid of SPS models so that it can replace expensive template fitting
    with a single forward pass.

    Parameters
    ----------
    n_phys_params : int
        Number of input physical parameters (default 6).
    n_bands : int
        Number of output photometric bands.
    n_hidden : int
        Width of hidden layers (default 256).
    """

    def __init__(self, n_phys_params=6, n_bands=9, n_hidden=256):
        super().__init__()
        self.n_phys_params = n_phys_params
        self.n_bands = n_bands

        self.network = nn.Sequential(
            nn.Linear(n_phys_params, n_hidden),
            nn.ReLU(inplace=True),
            nn.Linear(n_hidden, n_hidden),
            nn.ReLU(inplace=True),
            nn.Linear(n_hidden, n_hidden // 2),
            nn.ReLU(inplace=True),
            nn.Linear(n_hidden // 2, n_bands),
        )

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight, nonlinearity='relu')
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def forward(self, x):
        """
        Parameters
        ----------
        x : torch.Tensor
            Shape (batch, n_phys_params).

        Returns
        -------
        mags : torch.Tensor
            Shape (batch, n_bands). Predicted magnitudes.
        """
        return self.network(x)

    def count_parameters(self):
        """Count trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
