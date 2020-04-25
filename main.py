import numpy as np
import pandas as pd
from efficient_apriori import apriori
import csv

def create_dataframe():
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
            line = line.replace("complete genome", "").replace(" |", "|")
            tmp = line.split("|")
            tmp[2] = tmp[2].replace("\n", "")
            result.append(tmp)

    df2 = pd.DataFrame(result, columns=["Accession", "extra_info", "DNA"])
    df2.to_csv("reee.csv")
    convert_dict = {'Accession': object}
    df['Accession'] = df['Accession'].astype(str)
    df2['Accession'] = df2['Accession'].astype(str)

    df = pd.merge(df, df2, how="right", on="Accession")
    del df["extra_info"]
    del df["Publications"]
    del df["Authors"]
    del df["Genotype"]
    del df["Segment"]
    del df["Species"]
    del df["Genus"]
    del df["Family"]
    del df["BioSample"]
    del df["GenBank_Title"]
    del df["Isolation_Source"]
    del df["Nuc_Completeness"]
    return df

class DNA():
    def __init__(self, dna):
        self.dna = dna

    def __str__(self):
        return self.dna

if __name__ == '__main__':
    df = create_dataframe()
    print(df)
    dna = {}
    for i in df['Geo_Location'].unique():
        dna[i] = [{df['DNA'][j]: df['Collection_Date'][j]} for j in df[df['Geo_Location'] == i].index]
     print(dna)
    # transactions = []
    # for item in dna:
    #     tmp_trans = []
    #     for n in range(3,8):
    #         chunks = [item[i:i+n] for i in range(0, len(item), n)]
    #         tmp_trans += chunks
    #     transactions.append(tmp_trans)

    # result = apriori(transactions)