"""Star/Galaxy binary classifier CNN (~80K params)."""

import torch
import torch.nn as nn


class StarFinderCNN(nn.Module):
    """
    CNN for binary star/galaxy classification.

    Architecture (~80K params):
        Input(1, 64, 64)  (grayscale)
        -> Conv(1->16, 3x3) + BN + ReLU + MaxPool   [32x32]
        -> Conv(16->32, 3x3) + BN + ReLU + MaxPool   [16x16]
        -> Conv(32->64, 3x3) + BN + ReLU + MaxPool   [8x8]
        -> AdaptiveAvgPool(4x4)                       [4x4]
        -> Flatten                                    [64*4*4=1024]
        -> Dropout(0.3) + Linear(1024, 128) + ReLU
        -> Dropout(0.2) + Linear(128, 1) + Sigmoid
    """

    def __init__(self, in_channels=1, dropout=0.3, dropout2=0.2):
        super().__init__()
        self.in_channels = in_channels

        self.features = nn.Sequential(
            # Block 1: 1 -> 16, 64x64 -> 32x32
            nn.Conv2d(in_channels, 16, 3, padding=1, bias=False),
            nn.BatchNorm2d(16),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),

            # Block 2: 16 -> 32, 32x32 -> 16x16
            nn.Conv2d(16, 32, 3, padding=1, bias=False),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),

            # Block 3: 32 -> 64, 16x16 -> 8x8
            nn.Conv2d(32, 64, 3, padding=1, bias=False),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),

            # Adaptive pool to 4x4
            nn.AdaptiveAvgPool2d(4),
        )

        self.classifier = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(64 * 4 * 4, 128),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout2),
            nn.Linear(128, 1),
            nn.Sigmoid(),
        )

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, nonlinearity='relu')
            elif isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight, nonlinearity='relu')
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, nn.BatchNorm2d):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x):
        x = self.features(x)
        x = x.view(x.size(0), -1)
        return self.classifier(x)

    def count_parameters(self):
        """Count trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
