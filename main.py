import numpy as np
import pandas as pd
from efficient_apriori import apriori


def create_dataframe():
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

    df = pd.read_csv("sequences.csv")

    file = open("MT372482.1")
    file = file.read().replace("\n","")
    file = file.replace("gb|","")
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
    if not prefix:
        newtidlist = []
        for item, tidl in tidlist:
            if len(tidl) >= minsup:
                returnval.append(({item}, len(tidl)))
                newtidlist.append(({item}, tidl))
        if not returnval:
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


def checkSame(first,second):
    tmp1 = list(first)
    tmp2 = list(second)
    if len(tmp1) != len(tmp2):
        return False
    else:
        ret = False
        for index in range(len(tmp1)):
            if tmp1[index] != tmp2[index]:
                if tmp1[index] == "-" or tmp1[index] == "N":
                    ret = True
                    continue
                else:
                    ret = False

                if tmp2[index] == "-" or tmp2[index] == "N":
                    ret = True
                    continue
                else:
                    ret = False

                if not ret:
                    break
            else:
                ret = True
        return ret
if __name__ == '__main__':
    df = create_dataframe()
    print(df)


    tmp = list(df['DNA'])
    transactions = []
    for item in tmp:
        # item = item[:100]
        tmp_trans = []
        for n in range(3, 8):
            chunks = [item[i:i + n] for i in range(0, len(item), n)]
            tmp_trans += chunks
        transactions.append(tuple(tmp_trans))

    new_trans = []
    for i in range(len(transactions)):
        new_trans.append([])
    for index in range(len(transactions[0])):
        same = True
        same_val = None
        for index2 in range(len(transactions)):
            if same_val is None:
                same_val = transactions[index2][index]
            if not checkSame(same_val,transactions[index2][index]) and transactions[index2][index] :
                same = False
                break

        if not same:
            for index2 in range(len(transactions)):
                new_trans[index2].append(transactions[index2][index])

    print("done")
    result = apriori(transactions, min_support=0.55)
    print("?", result)
