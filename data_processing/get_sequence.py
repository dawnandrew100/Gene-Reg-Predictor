import pandas as pd


def main():
    file_path = "../filtered_mouse_Regulator_Gene.csv"
    df = pd.read_csv(file_path)
    print(df)

    regulators = df["Regulator Symbol"].unique().tolist()
    targets = df["Target Symbol"].unique().tolist()

    # List of all of the regulator and target symbols
    symbols = list(set(regulators + targets))
    print(
        f"One of the IDs in this list is '{symbols[0]}' and the list is {len(symbols)} IDs long."
    )


if __name__ == "__main__":
    main()
