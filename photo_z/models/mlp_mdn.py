"""MLP + Mixture Density Network for photometric redshift estimation (~50K params).

Optionally includes Quality Head v2 for outlier detection.
When quality_head=True, the model outputs P(reliable) — the probability
that the photo-z estimate is not an outlier (|Δz/(1+z)| < 0.15).
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np


class PhotoZMDN(nn.Module):
    """
    MLP with Mixture Density Network head for photometric redshift.

    Architecture (~29K params, ~35K with quality head):
        Input(n_features)
        -> FC(n_feat, 128) + BN + ReLU + Dropout(0.1)
        -> FC(128, 128) + BN + ReLU + Dropout(0.1)
        -> FC(128, 64) + BN + ReLU
        MDN Head (K Gaussian components):
          FC(64, K) -> Softmax   -> pi_k   (weights)
          FC(64, K)              -> mu_k   (means)
          FC(64, K) -> Softplus  -> sigma_k (stds, min 1e-4)
        Quality Head v2 (optional):
          Input: backbone(64) + pi(K) + mu(K) + sigma(K)
          FC(64+3K, 64) -> ReLU -> Dropout
          FC(64, 32) -> ReLU
          FC(32, 1) -> Sigmoid  -> P(reliable)

    Parameters
    ----------
    n_features : int
        Number of input features.
    n_hidden : int
        Hidden layer width (default 128).
    n_components : int
        Number of Gaussian mixture components K (default 5).
    dropout : float
        Dropout rate (default 0.1).
    quality_head : bool
        If True, add Quality Head v2 for outlier detection (default False).
    """

    def __init__(self, n_features=20, n_hidden=128, n_components=5,
                 dropout=0.1, quality_head=False):
        super().__init__()
        self.n_features = n_features
        self.n_hidden = n_hidden
        self.n_components = n_components
        self.has_quality_head = quality_head

        # Input batch normalization
        self.input_bn = nn.BatchNorm1d(n_features)

        # Hidden layers
        self.hidden = nn.Sequential(
            # Block 1
            nn.Linear(n_features, n_hidden),
            nn.BatchNorm1d(n_hidden),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            # Block 2
            nn.Linear(n_hidden, n_hidden),
            nn.BatchNorm1d(n_hidden),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            # Block 3
            nn.Linear(n_hidden, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(inplace=True),
        )

        # MDN head: K Gaussian components
        K = n_components
        self.fc_pi = nn.Linear(64, K)       # mixture weights
        self.fc_mu = nn.Linear(64, K)       # means
        self.fc_sigma = nn.Linear(64, K)    # standard deviations

        # Quality Head v2: backbone(64) + MDN outputs(3*K) -> P(reliable)
        if quality_head:
            q_input_dim = 64 + 3 * K
            self.quality_fc = nn.Sequential(
                nn.Linear(q_input_dim, 64),
                nn.ReLU(inplace=True),
                nn.Dropout(dropout),
                nn.Linear(64, 32),
                nn.ReLU(inplace=True),
                nn.Linear(32, 1),
            )

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
        Forward pass returning MDN parameters (and quality logit if enabled).

        Parameters
        ----------
        x : torch.Tensor
            Shape (batch, n_features).

        Returns
        -------
        pi : torch.Tensor
            Mixture weights, shape (batch, K), sum to 1.
        mu : torch.Tensor
            Component means, shape (batch, K).
        sigma : torch.Tensor
            Component stds, shape (batch, K), > 0.
        q_logit : torch.Tensor or None
            Quality logit, shape (batch,). Only if quality_head=True.
        """
        x = self.input_bn(x)
        h = self.hidden(x)

        pi = F.softmax(self.fc_pi(h), dim=1)
        mu = self.fc_mu(h)
        sigma = F.softplus(self.fc_sigma(h)) + 1e-4

        if self.has_quality_head:
            q_input = torch.cat([h, pi, mu, sigma], dim=1)
            q_logit = self.quality_fc(q_input).squeeze(-1)
            return pi, mu, sigma, q_logit

        return pi, mu, sigma

    @staticmethod
    def mdn_loss(pi, mu, sigma, y):
        """
        Negative log-likelihood of the Gaussian mixture model.

        Uses logsumexp for numerical stability.

        Parameters
        ----------
        pi : torch.Tensor
            Mixture weights, shape (batch, K).
        mu : torch.Tensor
            Component means, shape (batch, K).
        sigma : torch.Tensor
            Component stds, shape (batch, K).
        y : torch.Tensor
            Target values, shape (batch,) or (batch, 1).

        Returns
        -------
        loss : torch.Tensor
            Scalar, mean negative log-likelihood.
        """
        if y.dim() == 1:
            y = y.unsqueeze(1)  # (batch, 1)

        # Log probability of each Gaussian component
        # log N(y | mu_k, sigma_k) = -0.5*log(2*pi) - log(sigma_k) - 0.5*((y-mu_k)/sigma_k)^2
        log_normal = (-0.5 * np.log(2.0 * np.pi)
                      - torch.log(sigma)
                      - 0.5 * ((y - mu) / sigma) ** 2)

        # log(sum_k pi_k * N(y | mu_k, sigma_k))
        # = log(sum_k exp(log(pi_k) + log_normal_k))
        log_pi = torch.log(pi + 1e-10)
        log_mixture = torch.logsumexp(log_pi + log_normal, dim=1)

        return -log_mixture.mean()

    def predict(self, x):
        """
        Point predictions with uncertainty estimates.

        Parameters
        ----------
        x : torch.Tensor
            Shape (batch, n_features).

        Returns
        -------
        z_point : torch.Tensor
            Weighted mean of mixture, shape (batch,).
        z_err : torch.Tensor
            Weighted std (total uncertainty), shape (batch,).
        z_q68 : torch.Tensor
            68% credible interval width from mixture CDF, shape (batch,).
        p_reliable : torch.Tensor or None
            P(reliable) from Quality Head v2, shape (batch,).
            Only returned if quality_head=True. None otherwise.
        """
        self.eval()
        with torch.no_grad():
            fwd = self.forward(x)
            if self.has_quality_head:
                pi, mu, sigma, q_logit = fwd
                p_reliable = torch.sigmoid(q_logit)
            else:
                pi, mu, sigma = fwd
                p_reliable = None

            # Point estimate: weighted mean
            z_point = (pi * mu).sum(dim=1)

            # Total variance: law of total variance
            # Var = E[sigma^2] + E[mu^2] - (E[mu])^2
            mean_var = (pi * sigma ** 2).sum(dim=1)
            mean_sq = (pi * mu ** 2).sum(dim=1)
            total_var = mean_var + mean_sq - z_point ** 2
            z_err = torch.sqrt(torch.clamp(total_var, min=1e-8))

            # 68% credible interval from mixture CDF
            z_q68 = self._compute_q68(pi, mu, sigma, z_point)

        return z_point, z_err, z_q68, p_reliable

    def _compute_q68(self, pi, mu, sigma, z_point):
        """
        Compute 68% credible interval width by evaluating mixture CDF.

        Uses bisection-free approach: evaluate CDF at grid points around
        the point estimate and find the 16th and 84th percentiles.

        Parameters
        ----------
        pi, mu, sigma : torch.Tensor
            MDN parameters, each shape (batch, K).
        z_point : torch.Tensor
            Point estimates, shape (batch,).

        Returns
        -------
        q68 : torch.Tensor
            Width of 68% interval (z_84 - z_16), shape (batch,).
        """
        batch_size = pi.shape[0]
        device = pi.device

        # Total std for determining grid range
        mean_var = (pi * sigma ** 2).sum(dim=1)
        mean_sq = (pi * mu ** 2).sum(dim=1)
        total_std = torch.sqrt(torch.clamp(
            mean_var + mean_sq - z_point ** 2, min=1e-8))

        # Grid: z_point +/- 4*total_std, 200 points
        n_grid = 200
        z_lo = z_point - 4.0 * total_std  # (batch,)
        z_hi = z_point + 4.0 * total_std  # (batch,)
        t = torch.linspace(0, 1, n_grid, device=device).unsqueeze(0)  # (1, n_grid)
        z_grid = z_lo.unsqueeze(1) + t * (z_hi - z_lo).unsqueeze(1)  # (batch, n_grid)

        # Evaluate mixture CDF at each grid point
        # CDF(z) = sum_k pi_k * Phi((z - mu_k) / sigma_k)
        # z_grid: (batch, n_grid) -> (batch, n_grid, 1)
        # mu: (batch, K) -> (batch, 1, K)
        z_exp = z_grid.unsqueeze(2)
        mu_exp = mu.unsqueeze(1)
        sigma_exp = sigma.unsqueeze(1)
        pi_exp = pi.unsqueeze(1)

        # Standard normal CDF via erfc
        u = (z_exp - mu_exp) / (sigma_exp * np.sqrt(2.0))
        phi = 0.5 * torch.erfc(-u)  # (batch, n_grid, K)
        cdf = (pi_exp * phi).sum(dim=2)  # (batch, n_grid)

        # Find z_16 and z_84
        z_16 = self._find_quantile(z_grid, cdf, 0.16)
        z_84 = self._find_quantile(z_grid, cdf, 0.84)

        return z_84 - z_16

    @staticmethod
    def _find_quantile(z_grid, cdf, q):
        """
        Find quantile from CDF evaluated on a grid via linear interpolation.

        Parameters
        ----------
        z_grid : torch.Tensor
            Grid of z values, shape (batch, n_grid).
        cdf : torch.Tensor
            CDF values at grid points, shape (batch, n_grid).
        q : float
            Quantile to find (e.g. 0.16 or 0.84).

        Returns
        -------
        z_q : torch.Tensor
            Quantile values, shape (batch,).
        """
        batch_size, n_grid = z_grid.shape
        device = z_grid.device

        # Find the index where CDF crosses q
        above = (cdf >= q).float()
        # First index where CDF >= q
        # Use argmax on the above mask (first True)
        idx = above.argmax(dim=1)  # (batch,)
        # Clamp to valid range for interpolation
        idx = torch.clamp(idx, min=1, max=n_grid - 1)

        # Gather values for linear interpolation
        idx_lo = idx - 1
        batch_idx = torch.arange(batch_size, device=device)

        cdf_lo = cdf[batch_idx, idx_lo]
        cdf_hi = cdf[batch_idx, idx]
        z_lo = z_grid[batch_idx, idx_lo]
        z_hi = z_grid[batch_idx, idx]

        # Linear interpolation
        denom = cdf_hi - cdf_lo
        denom = torch.clamp(denom, min=1e-10)
        frac = (q - cdf_lo) / denom
        frac = torch.clamp(frac, 0.0, 1.0)

        return z_lo + frac * (z_hi - z_lo)

    def count_parameters(self):
        """Count trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
