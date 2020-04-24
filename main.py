import numpy as np
import pandas as pd
import csv

if __name__ == '__main__':
    pd.options.display.max_rows
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

    df = pd.read_csv("sequences.csv")

    file = open("sequences.fasta")
    file = file.read().split(">")
    result = []
    for line in file:
        if "complete genome" in line:
            line = line.replace("complete genome","")
            tmp = line.split("|")
            tmp[2] = tmp[2].replace("\n","")
            result.append(tmp)

    df2 = pd.DataFrame(result, columns=["Accession", "extra_info", "DNA"])
    df2.to_csv("reee.csv")
    convert_dict = {'Accession': object}
    df['Accession'] = df['Accession'].astype(str)
    df2['Accession'] = df2['Accession'].astype(str)

    # df = pd.merge(df, df2,  on="Accession")
    print(df.)