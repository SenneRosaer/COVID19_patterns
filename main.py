import numpy as np
import pandas as pd
from efficient_apriori import apriori
from own_apriori import apriori2
import os
import pickle
from sklearn import tree
import graphviz


def create_dataframe():
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

    df = pd.read_csv("sequences.csv")

    file = open("MT372482.1")
    file = file.read().replace("\n", "")
    file = file.replace("gb|", "")
    file = file.split(">")
    result = []
    for line in file:
        tmp = line.split(".")
        if len(tmp) > 1:
            tmp[1] = tmp[1][1:]
            result.append(tmp)
        else:
            pass
    df2 = pd.DataFrame(result, columns=["Accession", "DNA"])
    convert_dict = {'Accession': object}
    df['Accession'] = df['Accession'].astype(str)
    df2['Accession'] = df2['Accession'].astype(str)

    df = pd.merge(df, df2, how="right", on="Accession")
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


def create_y():
    df = pd.read_csv("owid-covid-data.csv")
    df = df[df.location == "China"]
    del df["new_tests_per_thousand"]
    del df["tests_units"]
    del df["iso_code"]
    del df["new_tests"]
    del df["total_tests"]
    del df["new_deaths_per_million"]
    del df["total_deaths_per_million"]
    del df["new_cases_per_million"]
    del df["total_cases_per_million"]
    del df["new_deaths"]
    del df["new_cases"]
    print(df)
    return df


def checkSame(first, second):
    first_bases = list(first)
    second_bases = list(second)
    if len(first_bases) != len(second_bases):
        return False
    for i in range(len(first_bases)):
        if first_bases[i] != second_bases[i]:
            if first_bases[i] == "-" or first_bases[i] == "N":
                continue
            elif second_bases[i] == "-" or second_bases[i] == "N":
                continue
            else:
                return False
    return True


def create_chunks(dna_list, chunk_min=3, chunk_max=8):
    transactions = list()
    for dna_sequence in dna_list:
        sequence_transactions = list()
        for n in range(chunk_min, chunk_max):
            chunks = [dna_sequence[i:(i + n)] for i in range(0, len(dna_sequence), n)]
            sequence_transactions += chunks
        transactions.append(tuple(sequence_transactions))
    return transactions


def filter_transactions(transactions):
    new_trans = [list() for _ in range(len(transactions))]  # != new_trans = [list()] * len(transactions) MEM BULLSHIT
    for index1 in range(len(transactions[0])):
        is_same = True
        same_val = None
        frequency_dict = {}
        for index2 in range(len(transactions)):
            current = transactions[index2][index1]
            if current not in frequency_dict.keys():
                frequency_dict[current] = 1
            else:
                frequency_dict[current] += 1

            if same_val is None:
                same_val = transactions[index2][index1]
            if not checkSame(same_val, transactions[index2][index1]) and transactions[index2][index1]:
                is_same = False

        if not is_same:
            for index2 in range(len(transactions)):
                current_item = transactions[index2][index1]

                if frequency_dict[current_item] / len(transactions) < 0.9:
                    new_trans[index2].append((str(index1) + ":" + transactions[index2][index1]))
    return new_trans


def cache_transactions(transactions, file_name='cache.txt'):
    with open(file_name, 'wb') as fp:
        pickle.dump(transactions, fp)


def write_apriori_results(transactions, file_name='test.txt'):
    result = apriori2(transactions, min_support=0.7, min_confidence=0.85, max_length=5)
    with open(file_name, 'w') as fp:
        for freq_dict in result:
            for key in result[freq_dict]:
                string = ""
                string += str(key)
                string += ": " + str(result[freq_dict][key]) + "\n"
                fp.write(string)
            fp.write("\n\n\n=======================\n\n\n\n")


def frequent_itemsets_apriori(df, cache_results=True):
    if not os.path.isfile("cache.txt"):
        transactions = create_chunks(df['DNA'])
        transactions = filter_transactions(transactions)
        if cache_results:
            cache_transactions(transactions)
    else:
        with open('cache.txt', 'rb') as fp:
            transactions = pickle.load(fp)
    write_apriori_results(transactions)


if __name__ == '__main__':
    df = create_dataframe()
    frequent_itemsets_apriori(df)
    print("?")
