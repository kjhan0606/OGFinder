"""MLP baseline model for merge/separate classification."""

import torch
import torch.nn as nn


class MergeMLP(nn.Module):
    """
    MLP classifier for merge/separate decision.

    Architecture:
    Input(20) → BN → [Linear(128) → BN → ReLU → Dropout(0.3)]
                    → [Linear(64)  → BN → ReLU → Dropout(0.2)]
                    → [Linear(32)  → BN → ReLU]
                    → Linear(2)

    ~12K parameters, <1ms inference on CPU.
    """

    def __init__(self, n_features=20, n_classes=2,
                 hidden_dims=None, dropout_rates=None):
        super().__init__()

        if hidden_dims is None:
            hidden_dims = [128, 64, 32]
        if dropout_rates is None:
            dropout_rates = [0.3, 0.2, 0.0]

        assert len(hidden_dims) == len(dropout_rates)

        # Input batch normalization
        self.input_bn = nn.BatchNorm1d(n_features)

        layers = []
        in_dim = n_features
        for h_dim, drop in zip(hidden_dims, dropout_rates):
            layers.append(nn.Linear(in_dim, h_dim))
            layers.append(nn.BatchNorm1d(h_dim))
            layers.append(nn.ReLU(inplace=True))
            if drop > 0:
                layers.append(nn.Dropout(drop))
            in_dim = h_dim

        self.hidden = nn.Sequential(*layers)
        self.classifier = nn.Linear(in_dim, n_classes)

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
            Shape (batch, n_features).

        Returns
        -------
        logits : torch.Tensor
            Shape (batch, n_classes). Raw logits (use CrossEntropyLoss).
        """
        x = self.input_bn(x)
        x = self.hidden(x)
        return self.classifier(x)

    def predict_proba(self, x):
        """Return softmax probabilities."""
        self.eval()
        with torch.no_grad():
            logits = self.forward(x)
            return torch.softmax(logits, dim=1)

    def count_parameters(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
