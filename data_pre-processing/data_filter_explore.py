import pandas as pd


def main():
    file_name = "mouse_Regulator_Gene"
    file_extension = ".txt"
    file_path = f"../{file_name}{file_extension}"
    headers = [
        "Regulator Symbol",
        "Regulator ID",
        "Target Symbol",
        "Target ID",
        "Regulator Type",
        "Target Type",
    ]
    # Initial data exploration
    df = pd.read_csv(file_path, names=headers, header=None, sep="\t")
    print(df)
    print(df.describe())
    """
           Regulator Symbol                Regulator ID Target Symbol   Target ID Regulator Type Target Type
    count           2011023                     2011023       2011023     2011023        2011023     2011023
    unique            16938                       17768         24399       21967              4           1
    top             Gm25239  Ensembl:ENSMUSG00000092656        Elavl1  NCBI:15568          miRNA        Gene
    freq              16410                       16410          2583        2583        1203978     2011023
    """
    print(df["Regulator Type"].value_counts())
    """
    Regulator Type
    miRNA      1203978
    lncRNA      677881
    TF          128552
    circRNA        612

    """

    # Filtering out TF and circRNA regulators
    selected_regulators = ["miRNA", "lncRNA"]
    # Removing extraneous Regulator and Target ID (and Target type since all target types are gene)
    cols_to_keep = ["Regulator Symbol", "Target Symbol", "Regulator Type"]
    filtered_df = df.loc[
        df["Regulator Type"].isin(selected_regulators), cols_to_keep
    ].copy()
    print(filtered_df.describe())
    """
           Regulator Symbol Target Symbol Regulator Type
    count           1881859       1881859        1881859
    unique            15603         20842              2
    top             Gm25239         Celf1          miRNA
    freq              16410          2549        1203978
    """
    print(filtered_df["Regulator Type"].value_counts())
    """
    Regulator Type
    miRNA     1203978
    lncRNA     677881
    """

    # Write filtered data to new CSV file
    filtered_filepath = f"../filtered_{file_name}.csv"
    filtered_df.to_csv(filtered_filepath, header=cols_to_keep, index=False)


if __name__ == "__main__":
    main()
