# Third party packages
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    roc_auc_score,
    RocCurveDisplay,
)
from sklearn.model_selection import train_test_split
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from tqdm import trange

# Built-ins
import random

# Local imports
from nn_model import GeneRegModel


def main():
    # Load data
    print("Loading data")
    df = pd.read_csv("../sequence_mouse_Regulator_Gene.csv")
    dna_embeddings = torch.load("./embeddings/dna_embeddings.pt")

    # Pair positive embeddings and generate negative samples
    print("Generating dataset")
    x_data, y_data = get_dataset(df, dna_embeddings)

    # Split data into training, validation, and testing
    x_train, x_test, y_train, y_test = train_test_split(
        x_data, y_data, test_size=0.2, stratify=y_data
    )

    # Create DataLoaders
    print("Creating DataLoaders")
    batch_size = 256
    train_loader = DataLoader(
        TensorDataset(x_train, y_train), batch_size=batch_size, shuffle=True
    )
    test_loader = DataLoader(
        TensorDataset(x_test, y_test), batch_size=batch_size, shuffle=False
    )

    # Invoke the model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = GeneRegModel(n_embeddings=512).to(device)

    # Train the model
    pos_weight = (1 - y_train.mean()) / y_train.mean()
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=pos_weight.to(device))
    optimizer = optim.Adam(model.parameters(), lr=1e-4, weight_decay=1e-3)

    train_losses = []
    test_losses = []
    # This took 1 hours 20 minutes to run
    for epoch in trange(100, desc="Training Epochs"):
        model.train()
        train_loss = 0
        for x, y in train_loader:
            x, y = x.to(device), y.to(device)
            optimizer.zero_grad()
            y_pred = model(x).view(-1, 1)
            loss = loss_fn(y_pred, y.unsqueeze(1))
            train_loss += loss.item() * x.shape[0]
            loss.backward()
            optimizer.step()
        train_losses.append(train_loss / len(train_loader.dataset))

        model.eval()
        test_loss = 0
        with torch.no_grad():
            for x, y in test_loader:
                x, y = x.to(device), y.to(device)
                y_pred = model(x).view(-1, 1)
                loss = loss_fn(y_pred, y.unsqueeze(1))
                test_loss += loss.item() * x.shape[0]
            test_losses.append(test_loss / len(test_loader.dataset))

    # Plot training and test losses
    fig, ax = plt.subplots()
    ax.plot(train_losses, label="train")
    ax.plot(test_losses, label="test")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")
    ax.legend()
    plt.show()

    # Evaluate the model using accuracy, precision, recall, and AUC
    model.eval()
    for name, loader in [("Train", train_loader), ("Test", test_loader)]:
        y_probs, y_preds, y_true = [], [], []
        with torch.no_grad():
            for x, y in loader:
                x = x.to(device)
                logits = model(x)
                probs = torch.sigmoid(logits).cpu().numpy().flatten()

                y_probs.extend(probs)
                y_preds.extend((probs > 0.5).astype(int))
                y_true.extend(y.numpy().flatten())

        print(f"\n--- {name} Metrics ---")
        print(f"Accuracy:  {accuracy_score(y_true, y_preds):.4f}")
        print(f"Precision: {precision_score(y_true, y_preds):.4f}")
        print(f"Recall:    {recall_score(y_true, y_preds):.4f}")
        print(f"AUC:       {roc_auc_score(y_true, y_probs):.4f}")


def get_dataset(df, embeddings):
    x = []
    y = []

    all_targets = df["Target Symbol"].unique().tolist()
    all_pairs = set(zip(df["Regulator Symbol"], df["Target Symbol"]))

    for _, row in df.iterrows():
        regulator = row["Regulator Symbol"]
        target = row["Target Symbol"]

        if regulator in embeddings and target in embeddings:
            # --- POSITIVE SAMPLE ---
            pos_seq = torch.stack([embeddings[regulator], embeddings[target]])
            x.append(pos_seq)
            y.append(1.0)

            # --- NEGATIVE SAMPLES ---
            neg_tar = random.choice(all_targets)
            while (regulator, neg_tar) in all_pairs:
                neg_tar = random.choice(all_targets)

            if neg_tar in embeddings:
                neg_seq = torch.stack([embeddings[regulator], embeddings[neg_tar]])
                x.append(neg_seq)
                y.append(0.0)

    return torch.stack(x).float(), torch.tensor(y).float()


if __name__ == "__main__":
    main()
