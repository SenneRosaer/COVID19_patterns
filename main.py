import numpy as np
import pandas as pd
from efficient_apriori import apriori
from  prefixspan import PrefixSpan
import csv
from Bio import AlignIO

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


def transactionsToTidlist(transactions):
    """ Converts transactions matrix to tidlist.
        Return: List of the form [(item1, {tids1}), (item2, {tids2})]
        (Hint: Store them in a dict d (item -> set) and return list(d.items())
    """
    # TODO: Implement
    count = 0
    returndict = dict()
    for trans in transactions:
        for item in trans:
            if item not in returndict:
                tempset = set()
                tempset.add(count)
                returndict[item] = tempset
            else:
                returndict[item].add(count)
        count += 1
    print("yet")
    return list(returndict.items())


def eclat(df, minsup):
    tidlist = transactionsToTidlist(transactions)
    return _eclat([], tidlist, minsup)


def _eclat(prefix, tidlist, minsup):
    """ Implement the Eclat algorithm recursively.
        prefix: items in this depth first branch (the set alpha).
        tidlist: tids of alpha-conditional db.
        minsup: minimum support.
        return: list of itemsets with support > minsup. Format: [({item1, item2}, supp1), ({item1}, supp2)]
    """
    returnval = []
    if prefix == []:
        newtidlist = []
        for item, tidl in tidlist:
            if len(tidl) >= minsup:
                returnval.append(({item}, len(tidl)))
                newtidlist.append(({item}, tidl))
        if returnval == []:
            return []
        tidlist = newtidlist

    for i in range(0, len(tidlist)):
        n_tidlist = []
        item = tidlist[i][0]
        tidl = tidlist[i][1]
        for j in range(i + 1, len(tidlist)):
            item2 = tidlist[j][0]
            tidl2 = tidlist[j][1]
            tempitems = item.union(item2)
            temptidl = tidl.intersection(tidl2)
            if len(tempitems) == len(item) + 1:
                if len(temptidl) >= minsup:
                    n_tidlist.append((tempitems, temptidl))
                    returnval.append((tempitems, len(temptidl)))
        result_list = []
        if len(n_tidlist) > 1:
            result_list = _eclat([item], n_tidlist, minsup)
        returnval = returnval + result_list

    return returnval


if __name__ == '__main__':
    df = create_dataframe()
    print(df)

    geo = df.loc[df.Geo_Location == 'China']

    tmp = list(df['DNA'])
    transactions = []
    for item in tmp:
        item = item[:100]
        tmp_trans = []
        for n in range (3,8):
            chunks = [item[i:i+n] for i in range(0, len(item), n)]
            tmp_trans += chunks
        transactions.append(tuple(tmp_trans))

    new_trans = []
    for i in range(len(transactions)):
        new_trans.append([])
    for index in range(len(transactions[0])):
        same = True
        same_val = None
        for index2 in range(len(transactions)):
            if same_val == None:
                same_val = transactions[index2][index]
            if same_val != transactions[index2][index]:
                same = False
                break

        if not same:
            for index2 in range(len(transactions)):
                new_trans[index2].append(transactions[index2][index])


    print("done")
    result = apriori(transactions, min_support=0.55)
    print("?")

