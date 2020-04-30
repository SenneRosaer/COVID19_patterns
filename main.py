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


if __name__ == '__main__':
    y = create_y()

    if not os.path.isfile("cache.txt"):
        df = create_dataframe()
        print(df)

        tmp = list(df['DNA'])
        transactions = []
        for item in tmp:
            #item = item[:5835]
            tmp_trans = []
            for n in range(3,8):
                chunks = [item[i:i + n] for i in range(0, len(item), n)]
                tmp_trans += chunks
            transactions.append(tuple(tmp_trans))

        new_trans = []
        for i in range(len(transactions)):
            new_trans.append([])
        for index in range(len(transactions[0])):
            same = True
            same_val = None

            frequency_dict = {

            }

            for index2 in range(len(transactions)):
                current = transactions[index2][index]
                if current not in frequency_dict.keys():
                    frequency_dict[current] = 1
                else:
                    frequency_dict[current] += 1

                if same_val is None:
                    same_val = transactions[index2][index]
                if not checkSame(same_val, transactions[index2][index]) and transactions[index2][index]:
                    same = False

            if not same:
                for index2 in range(len(transactions)):
                    current_item = transactions[index2][index]

                    if frequency_dict[current_item] / len(transactions) < 0.9:
                        new_trans[index2].append((str(index) + ":" + transactions[index2][index]))

        with open('cache.txt', 'wb') as fp:
            pickle.dump(new_trans, fp)
    else:
        with open('cache.txt', 'rb') as fp:
            new_trans = pickle.load(fp)

    for index in range(len(new_trans)):
        new_trans[index] = tuple(new_trans[index])
    print("done")

    print(new_trans)

    clf = tree.DecisionTreeClassifier()

    ###Write to file###
    # file = open("transactions.txt", "w")
    # for trans in new_trans:
    #     tmp = ""
    #     for item in trans:
    #         tmp += item +", "
    #
    #     tmp = tmp[:-2]
    #     tmp += "\n"
    #     file.write(tmp)
    # file.close()

    result = apriori2(new_trans,min_support=0.7, min_confidence=0.85,max_length=5)
    with open('test.txt', 'w') as fp:
        for freq_dict in result:
            for key in freq_dict:
                string = ""
                string += str(key)
                string += ": " + freq_dict[key]
                fp.write(string)
            fp.write("=======================")

    print("?")