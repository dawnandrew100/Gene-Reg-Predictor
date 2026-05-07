import torch
import torch.nn as nn


# This took 1 hours 20 minutes to train to 100 epochs
class GeneRegRNNModel(nn.Module):
    def __init__(self, embed_dim):
        super().__init__()
        # RNN layer
        self.rnn = nn.GRU(
            embed_dim, hidden_size=64, batch_first=True, bidirectional=True
        )
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
        # (B, 2, 512) -> (B, 2, 128)
        x, _ = self.rnn(x)
        # (B, 2, 128) -> (B, 2, 32)
        x = self.fc1(x)
        # (B, 2, 1)
        attn_weights = torch.softmax(self.attn(x), dim=1)
        # (B, 2, 32) * (B, S, 1) -> (B, 32)
        context = torch.sum(x * attn_weights, dim=1)
        # (B, 32) -> (B, 1)-> (B)
        return self.fc2(context).squeeze(-1)


# This took 6 hours 22 minutes to train 100 epochs
class GeneRegTransformerModel(nn.Module):
    def __init__(self, embed_dim, n_heads, num_layers=3):
        super().__init__()

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=embed_dim, nhead=n_heads, dim_feedforward=1024, batch_first=True
        )
        self.transformer_stack = nn.TransformerEncoder(
            encoder_layer, num_layers=num_layers
        )

        self.classifier = nn.Sequential(
            nn.Linear(embed_dim * 2, 256), nn.ReLU(), nn.Linear(256, 1)
        )

    def forward(self, x):
        # (B, 2, 512) -> (B, 2, 512)
        x = self.transformer_stack(x)
        # (B, 2, 512) -> (B, 1024)
        x = x.view(x.size(0), -1)

        # (B, 1024) -> (B, 256) -> (B, 1)
        return self.classifier(x)
