import numpy as np
import pandas as pd
from efficient_apriori import apriori
import os
import datetime
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
    tidlist = transactionsToTidlist(df)
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
    mortality = create_y()

    mortality['date'] = mortality['date'].apply(lambda x: datetime.datetime.strptime(x, "%Y-%m-%d"))
    if not os.path.isfile("cache.txt") and not os.path.isfile("final_list_cache.txt"):
        df = create_dataframe()
        print(df)

        del df["Accession"]
        del df["Length"]
        del df["Sequence_Type"]
        del df["Host"]
        del df["Release_Date"]


        date_dna_list = list(df.itertuples(index=False))
        transactions = []
        for i in date_dna_list:
            item = i[2]
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
                        new_trans[index2].append((transactions[index2][index]))

        trans_tuples = []
        for index, item in enumerate(date_dna_list):
            if str(item[1]) != "nan":
                if len(str(item[1])) == 7:
                    d = datetime.datetime.strptime(item[1], "%Y-%m")
                else:
                    d = datetime.datetime.strptime(item[1], "%Y-%m-%d")
                trans_tuples.append((new_trans[index], d ))

        final_list = []
        for t in trans_tuples:
            y = mortality[mortality["date"] == (t[1]+datetime.timedelta(days=7))]
            a =  list(y["total_deaths"])
            b = list(y["total_cases"])
            if len(a) and len(b):
                final_list.append((t[0], a[0] / b[0]))
        print(final_list)

        with open('cache.txt', 'wb') as fp:
            pickle.dump(new_trans, fp)

        with open('final_list_cache.txt', 'wb') as fp:
            pickle.dump(final_list, fp)
    else:
        with open('cache.txt', 'rb') as fp:
            new_trans = pickle.load(fp)

        with open('final_list_cache.txt', 'rb') as fp:
            final_list = pickle.load(fp)

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

    result = apriori(new_trans,min_support=0.7, min_confidence=0.85,max_length=5)
    # print("?", result)