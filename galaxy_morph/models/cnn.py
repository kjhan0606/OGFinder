"""Galaxy morphology CNN classifier (~250K params)."""

import torch
import torch.nn as nn


class ResBlock(nn.Module):
    """Residual block with optional 1x1 projection shortcut."""

    def __init__(self, in_channels, out_channels, stride=1):
        super().__init__()
        self.conv1 = nn.Conv2d(in_channels, out_channels, 3,
                               stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(out_channels)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = nn.Conv2d(out_channels, out_channels, 3,
                               stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(out_channels)

        self.shortcut = nn.Identity()
        if stride != 1 or in_channels != out_channels:
            self.shortcut = nn.Sequential(
                nn.Conv2d(in_channels, out_channels, 1,
                          stride=stride, bias=False),
                nn.BatchNorm2d(out_channels),
            )

    def forward(self, x):
        out = self.relu(self.bn1(self.conv1(x)))
        out = self.bn2(self.conv2(out))
        out = out + self.shortcut(x)
        return self.relu(out)


class GalaxyMorphCNN(nn.Module):
    """
    CNN for Galaxy10 DECaLS 10-class morphology classification.

    Architecture (~250K params):
        Input(3, 64, 64)  (RGB for training, 1ch->3ch duplicated for FITS)
        -> Conv(3->16, 3x3) + BN + ReLU + MaxPool  [32x32]
        -> Conv(16->32, 3x3) + BN + ReLU + MaxPool  [16x16]
        -> ResBlock(32->64, stride=2)                [8x8]
        -> ResBlock(64->128, stride=2)               [4x4]
        -> GlobalAvgPool                             [128]
        -> Dropout(0.3) + Linear(128,64) + ReLU + Dropout(0.2)
        -> Linear(64, 10)
    """

    def __init__(self, n_classes=10, in_channels=3, dropout=0.3, dropout2=0.2,
                 hidden_dim=64):
        super().__init__()
        self.in_channels = in_channels

        self.features = nn.Sequential(
            # Block 1: in_channels -> 16 channels, 64x64 -> 32x32
            nn.Conv2d(in_channels, 16, 3, padding=1, bias=False),
            nn.BatchNorm2d(16),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),

            # Block 2: 16 -> 32 channels, 32x32 -> 16x16
            nn.Conv2d(16, 32, 3, padding=1, bias=False),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),

            # ResBlock: 32 -> 64, 16x16 -> 8x8
            ResBlock(32, 64, stride=2),

            # ResBlock: 64 -> 128, 8x8 -> 4x4
            ResBlock(64, 128, stride=2),

            # Global average pooling -> [batch, 128]
            nn.AdaptiveAvgPool2d(1),
        )

        self.classifier = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(128, hidden_dim),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout2),
            nn.Linear(hidden_dim, n_classes),
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

    def predict_proba(self, x):
        """Return softmax probabilities."""
        self.eval()
        with torch.no_grad():
            logits = self.forward(x)
            return torch.softmax(logits, dim=1)

    def count_parameters(self):
        """Count trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
