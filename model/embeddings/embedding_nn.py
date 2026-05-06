"""
Relevant paper on embedding a DNA sequence
https://www.sciencedirect.com/science/article/pii/S2001037023005214
"""

import torch
import torch.nn as nn


class NeuralBedsCNN(nn.Module):
    def __init__(self, embedding_dim=512):
        super(NeuralBedsCNN, self).__init__()

        # Convolutional layer: 1024 neurons, kernel size 5, ReLU activation
        self.conv1 = nn.Conv2d(1, 1024, kernel_size=5, stride=1, padding=0)
        self.relu1 = nn.ReLU()

        # Pooling layer with stride 1, padding 0
        self.pool = nn.AdaptiveAvgPool2d((8, 8))  # Adaptive pooling

        # Fully connected layer: 512 neurons, ReLU, dropout 0.2
        self.fc1 = nn.Linear(1024 * 8 * 8, 512)
        self.relu2 = nn.ReLU()
        self.dropout = nn.Dropout(0.2)

        # Output layer with tanh activation
        self.fc_out = nn.Linear(512, embedding_dim)
        self.tanh = nn.Tanh()

    def forward(self, x):
        x = self.relu1(self.conv1(x))
        x = self.pool(x)
        x = self.dropout(self.relu2(self.fc1(x)))
        x = self.tanh(self.fc_out(x))
        return x
