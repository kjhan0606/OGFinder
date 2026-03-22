"""Inverse SED model with Mixture Density Network (MDN) uncertainty.

Architecture (~80K params):
    Input(n_bands + 1)  [magnitudes + photo_z]
    -> FC(N+1, 128) -> BN -> ReLU -> Dropout(0.1)
    -> FC(128, 128) -> BN -> ReLU -> Dropout(0.1)
    -> FC(128, 64) -> BN -> ReLU
    MDN Head (K=3 components, D=5 output dims):
      FC(64, K) -> Softmax      -> pi     (mixture weights)
      FC(64, K*D) -> reshape    -> mu     (component means)
      FC(64, K*D) -> Softplus   -> sigma  (component stds, min 1e-4)

Output dimensions D=5: log_mass, log_age, log_Z, Av, log_tau
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

from ..config import OUTPUT_PARAM_NAMES


class InverseSEDMLP(nn.Module):
    """
    Inverse SED model: observed magnitudes + photo-z -> physical params.

    Uses a Mixture Density Network head to capture multi-modal posterior
    distributions (e.g., age-metallicity degeneracy).

    Parameters
    ----------
    n_input : int
        Number of input features (n_bands + 1 for photo_z).
    n_output : int
        Number of output physical parameters (default 5).
    n_hidden : int
        Width of hidden layers (default 128).
    n_components : int
        Number of MDN mixture components (default 3).
    dropout : float
        Dropout rate (default 0.1).
    """

    def __init__(self, n_input=10, n_output=5, n_hidden=128,
                 n_components=3, dropout=0.1):
        super().__init__()
        self.n_input = n_input
        self.n_output = n_output
        self.n_hidden = n_hidden
        self.n_components = n_components
        self.K = n_components
        self.D = n_output

        # Shared trunk
        self.trunk = nn.Sequential(
            nn.Linear(n_input, n_hidden),
            nn.BatchNorm1d(n_hidden),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),

            nn.Linear(n_hidden, n_hidden),
            nn.BatchNorm1d(n_hidden),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),

            nn.Linear(n_hidden, n_hidden // 2),
            nn.BatchNorm1d(n_hidden // 2),
            nn.ReLU(inplace=True),
        )

        h = n_hidden // 2  # 64

        # MDN heads
        self.pi_head = nn.Linear(h, self.K)               # mixture weights
        self.mu_head = nn.Linear(h, self.K * self.D)       # means
        self.sigma_head = nn.Linear(h, self.K * self.D)    # stds

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight, nonlinearity='relu')
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, nn.BatchNorm1d):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x):
        """
        Parameters
        ----------
        x : torch.Tensor
            Shape (batch, n_input).  Normalized magnitudes + photo_z.

        Returns
        -------
        pi : torch.Tensor
            Shape (batch, K).  Mixture weights (sum to 1).
        mu : torch.Tensor
            Shape (batch, K, D).  Component means.
        sigma : torch.Tensor
            Shape (batch, K, D).  Component standard deviations (>= 1e-4).
        """
        h = self.trunk(x)

        # Mixture weights
        pi = F.softmax(self.pi_head(h), dim=-1)            # (batch, K)

        # Means
        mu = self.mu_head(h)                                # (batch, K*D)
        mu = mu.view(-1, self.K, self.D)                    # (batch, K, D)

        # Standard deviations (positive, floored at 1e-4)
        sigma = F.softplus(self.sigma_head(h)) + 1e-4       # (batch, K*D)
        sigma = sigma.view(-1, self.K, self.D)              # (batch, K, D)

        return pi, mu, sigma

    @staticmethod
    def mdn_loss(pi, mu, sigma, y):
        """
        Negative log-likelihood loss for the MDN.

        Parameters
        ----------
        pi : torch.Tensor
            Shape (batch, K). Mixture weights.
        mu : torch.Tensor
            Shape (batch, K, D). Component means.
        sigma : torch.Tensor
            Shape (batch, K, D). Component stds.
        y : torch.Tensor
            Shape (batch, D). Target values.

        Returns
        -------
        loss : torch.Tensor
            Scalar. Mean negative log-likelihood.
        """
        K = pi.shape[1]
        D = mu.shape[2]

        # Expand y to (batch, K, D)
        y_expanded = y.unsqueeze(1).expand_as(mu)

        # Log-probability of each Gaussian component (independent dims)
        # log N(y | mu, sigma) = -0.5 * [D*log(2pi) + sum(2*log(sigma) + ((y-mu)/sigma)^2)]
        var = sigma ** 2
        log_component = -0.5 * (
            D * torch.log(torch.tensor(2.0 * 3.141592653589793, device=mu.device))
            + torch.sum(2.0 * torch.log(sigma) + (y_expanded - mu) ** 2 / var, dim=-1)
        )  # (batch, K)

        # Log mixture probability: log(sum_k pi_k * N_k)
        log_pi = torch.log(pi + 1e-10)                     # (batch, K)
        log_mix = torch.logsumexp(log_pi + log_component, dim=-1)  # (batch,)

        return -log_mix.mean()

    def predict(self, x):
        """
        Point estimates and uncertainties from the MDN.

        Takes the component with the highest mixture weight as the point
        estimate, and reports the corresponding sigma as the error.  Also
        computes a weighted mean across all components.

        Parameters
        ----------
        x : torch.Tensor
            Shape (batch, n_input).

        Returns
        -------
        result : dict
            Keys are parameter names from OUTPUT_PARAM_NAMES.
            Each value is a dict with 'value' (best component),
            'error' (best component sigma), 'weighted_mean',
            'weighted_std'.  All numpy arrays of shape (batch,).
        """
        self.eval()
        with torch.no_grad():
            pi, mu, sigma = self.forward(x)

            # Best component (highest weight)
            best_k = pi.argmax(dim=1)                       # (batch,)
            batch_idx = torch.arange(pi.size(0), device=pi.device)

            best_mu = mu[batch_idx, best_k, :]              # (batch, D)
            best_sigma = sigma[batch_idx, best_k, :]        # (batch, D)

            # Weighted mean and std across components
            pi_expanded = pi.unsqueeze(-1)                  # (batch, K, 1)
            weighted_mean = (pi_expanded * mu).sum(dim=1)   # (batch, D)
            diff_sq = (mu - weighted_mean.unsqueeze(1)) ** 2
            weighted_var = (pi_expanded * (sigma ** 2 + diff_sq)).sum(dim=1)
            weighted_std = torch.sqrt(weighted_var)         # (batch, D)

        result = {}
        for i, name in enumerate(OUTPUT_PARAM_NAMES):
            result[name] = {
                'value': best_mu[:, i].cpu().numpy(),
                'error': best_sigma[:, i].cpu().numpy(),
                'weighted_mean': weighted_mean[:, i].cpu().numpy(),
                'weighted_std': weighted_std[:, i].cpu().numpy(),
            }

        return result

    def count_parameters(self):
        """Count trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
