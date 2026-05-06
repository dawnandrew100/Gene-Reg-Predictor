import torch
import torch.nn as nn


class GeneRegModel(nn.Module):
    def __init__(self, n_embeddings=512):
        super().__init__()
        # RNN layer
        self.rnn = nn.GRU(n_embeddings, 64, batch_first=True, bidirectional=True)
        # Fully connected layer
        self.fc1 = nn.Sequential(
            nn.Linear(128, 32),
            nn.ReLU(),
            nn.Dropout(0.3),
        )
        # Attention
        self.attn = nn.Linear(32, 1)
        # Fully connected layer to output binary
        self.fc2 = nn.Sequential(
            nn.Linear(32, 16),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(16, 1),
        )

    def forward(self, x):
        # (B, S, E) -> (B, S, 32)
        x, _ = self.rnn(x)
        # (B, S, 32) -> (B, S, 16)
        x = self.fc1(x)
        # (B, S, 16) -> (B, S, 1)
        attn_weights = torch.softmax(self.attn(x), dim=1)
        # (B, S, 16) * (B, S, 1) -> (B, 16)
        context = torch.sum(x * attn_weights, dim=1)
        # (B, 16) -> (B)
        return self.fc2(context).squeeze(-1)
