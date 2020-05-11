import numpy as np
import pandas as pd
from efficient_apriori import apriori
from own_apriori import apriori2
import os
import datetime
import pickle
from sklearn.tree import export_graphviz, DecisionTreeRegressor
from sklearn.feature_extraction.text import CountVectorizer
import graphviz
import copy

def create_dataframe():
    """
    Creates dataframe from sequences.csv which contains data about location and time and MT.. file which contains
    DNA strings alligned
    Also removes unnecesary columns
    :return: Pandas dataframe with Dna strings, location, time, ...
    """
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

    df = pd.read_csv("input/sequences.csv")

    file = open("input/allignment.aln")
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
    df['Accession'] = df['Accession'].astype(str)
    df2['Accession'] = df2['Accession'].astype(str)

    df = pd.merge(df, df2, how="right", on="Accession")
    return df


def create_y():
    """
    This function is used to create a dataframe containing the mortality data
    :return: Dataframe with mortality data
    (columns: country, total deaths, amount of cases, date, percentage of population that is tested)
    """
    df = pd.read_csv("./input/owid-covid-data.csv")
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

    df['date'] = df['date'].apply(lambda x: datetime.datetime.strptime(x, "%Y-%m-%d"))

    return df


def checkSame(first, second):
    """
    Controls if two parts of a DNA string are the same, if one of them has - or N we consider this as a symbol that can
    be anything
    :param first: DNA part as string
    :param second: DNA part as string
    :return: Boolean for equal or not
    """
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


def create_chunks(date_dna_list, chunk_min=3, chunk_max=8):
    """
    Splits a list of strings of DNA in chunks of certain lenghts
    :param dna_list: list of strings of DNA
    :param chunk_min: minimum length of chunk
    :param chunk_max: maximum length of chunk
    :return: List of lists that contain the chunks
    """
    transactions = list()
    for item in date_dna_list:
        sequence_transactions = list()
        for n in range(chunk_min, chunk_max):
            chunks = []
            ATG_pos = 0
            for i in range(0, len(item[2]), n):
                if i < len(item[2])-2:
                    if item[2][i:i+3] == "ATG":
                        ATG_pos=i
                chunks.append((item[2][i:(i + n)], str(i) + "-ATG@" + str(i - ATG_pos - item[2].count("-", ATG_pos, i))))
            sequence_transactions += chunks
        transactions.append(tuple(sequence_transactions))
    return transactions


def create_chunks2(date_dna_list, chunk_min=3, chunk_max=8):
    """
    Splits a list of strings of DNA in chunks of certain lenghts
    :param dna_list: list of strings of DNA
    :param chunk_min: minimum length of chunk
    :param chunk_max: maximum length of chunk
    :return: List of lists that contain the chunks
    """
    transactions = list()
    for item in date_dna_list:
        sequence_transactions = list()
        for n in range(chunk_min, chunk_max):
            chunks = []
            ATG_pos = 0
            for i in range(0, len(item), n):
                if i < len(item) - 2:
                    if item[i:i + 3] == "ATG":
                        ATG_pos = i
                chunks.append((item[i:(i + n)], str(i) + "-ATG@" + str(i - ATG_pos - item.count("-", ATG_pos, i))))
            sequence_transactions += chunks
        transactions.append(tuple(sequence_transactions))
    return transactions


def filter_transactions(transactions):
    """
    Filters the transactions since they contain way to much data
    1) If every transaction (of the transaction list) has the same item at a certain index we remove this since it is
    duplicate ==> Done by taking every possible index in a transaction and check if it is the same for every transaction

    2) If at a certain index in the transaction we have a pattern with a very high frequency we remove this since we have
    more use for the less frequent ones ==> Done by using a dict to count the occurrence

    :param transactions: List of transactions that we need to filter
    :return: List of transactions without unnecessary data
    """
    new_trans = [list() for _ in range(len(transactions))]  # != new_trans = [list()] * len(transactions) MEM BULLSHIT
    for index1 in range(len(transactions[0])):
        is_same = True
        same_val = None
        frequency_dict = {}
        for index2 in range(len(transactions)):
            current = transactions[index2][index1][0]
            if current not in frequency_dict.keys():
                frequency_dict[current] = 1
            else:
                frequency_dict[current] += 1

            if same_val is None:
                same_val = transactions[index2][index1][0]
            if not checkSame(same_val, current) and current:
                is_same = False

        if not is_same:
            for index2 in range(len(transactions)):
                current_item, position = transactions[index2][index1]

                if frequency_dict[current_item] / len(transactions) < 0.9:
                    new_trans[index2].append((current_item + "@" + str(position)))
    return new_trans


def cache(transactions, file_name):
    """
    Cache certain transactions to use later since filtering can take some time
    :param transactions: what we will cache
    :param file_name: file to cache to
    """
    with open("cache/" + file_name, 'wb') as fp:
        pickle.dump(transactions, fp)


def uncache(file_name):
    """
    Load a cache file
    :param file_name: cache file to load
    :return: list of transactions that were in the cache file
    """
    with open("cache/" + file_name, 'rb') as fp:
        transactions = pickle.load(fp)
    return transactions


def is_cached(file_name):
    """
    Check whether a file is in the cache folder
    :param file_name: file to check for
    :return: True if file exists in cache folder, else False
    """
    return os.path.isfile("cache/" + file_name)


def write_apriori_results(results, file_name='test'):
    """
    Write the results of executing apriori to a file for easier usage
    :param results: Results of the apriori run
    :param file_name: File name we want to output (time will be added to name so we always have a new file)
    :return: results
    """
    apr = apriori2(results, min_support=0.7, min_confidence=0.85, max_length=5)
    result = apr[0]
    with open('output/' + file_name + datetime.datetime.now().strftime('%Y-%m-%d_%Hh%M') + '.txt', 'w') as fp:
        for freq_dict in result:
            for key in result[freq_dict]:
                string = ""
                string += str(key)
                string += ": " + str(result[freq_dict][key]) + "\n"
                fp.write(string)
            fp.write("\n\n\n=======================\n\n\n\n")

    with open('output/' + file_name + datetime.datetime.now().strftime('%Y-%m-%d_%Hh%M') + '-rules.txt', 'w') as fp:
        for rule in apr[1]:
            boolean = False
            for i in rule.rhs:
                if '[' in i:
                    boolean = True
            for i in rule.lhs:
                if '[' in i:
                    boolean = True
            if boolean:
                fp.write(str(rule) + "\n")
                boolean = False
    return apr


def create_final_list(mortality, transactions, date_dna_list, useTrans=False):
    """
    This function combines mortality with the transactions.
    The date is used to combine the transactions and the mortality.

    Only returns data that has mortality.

    :param mortality: mortality data
    :param transactions: list of transactions
    :param date_dna_list: list of dates
    :param useTrans: use the transactions or not
    :return: list that links DNA strings with mortality
    """
    trans_tuples = []
    for index, item in enumerate(date_dna_list):
        if str(item[2]) != "nan":
            if len(str(item[2])) == 7:
                d = datetime.datetime.strptime(item[2], "%Y-%m")
            else:
                d = datetime.datetime.strptime(item[2], "%Y-%m-%d")
            trans_tuples.append((item[3], d, item[1]))  # DNA, Date, Location

    final_list = []
    for index, t in enumerate(trans_tuples):
        y = mortality[mortality["date"] == (t[1] + datetime.timedelta(days=7))]
        y = y[y["location"] == t[2]]
        a = list(y["total_deaths"])
        b = list(y["total_cases"])
        if len(a) and len(b):
            if useTrans:
                m = 0
                if b[0] != 0:
                    m = a[0] / b[0]
                m = m * 100
                if m < 3:
                    m = "[0 - 3]"
                else:
                    m = "[3 - ..]"
                transactions[index].append(m)
            else:
                m = 0
                if b[0] != 0:
                    m = a[0] / b[0]
                m = m * 100
                final_list.append((t[0], m))

    return final_list


def make_date_dna_list(df):
    """
    Creates a list of the dates.
    :param df: Dataframe to create list from
    :return: list of dates
    """
    print(df)
    return list(df.itertuples(index=False))


class Tokenizer(object):
    """Custom tokenizer used to create tokens"""
    def __call__(self, doc):
        return set(doc.split(","))


def make_tree(list):
    """
    Creates a tree using DNA data in a list
    :param list: list containing DNA
    """
    string_X = []
    Y = []
    for i in list:
        string_X.append(i[0])
        Y.append(i[1])

    if not is_cached("string_x"):
        string_X = create_chunks2(string_X)
        string_X = filter_transactions(string_X)
        tmp = []
        for item in string_X:
            string = str(item) + ","
            tmp.append(string)
        string_X = tmp
        cache(string_X, "string_x")
    else:
        string_X = uncache("string_x")

    if not is_cached("vect"):
        vect = CountVectorizer(tokenizer=Tokenizer())
        vect.fit(string_X)
        X = vect.transform(string_X)
        cache(X, "vect_X")
        feature_names = vect.get_feature_names()
        cache(feature_names, "ft_names")
    else:
        X = uncache("vect_X")
        feature_names = uncache("ft_names")

    print("is fitting")
    if not is_cached("tree"):
        cls = DecisionTreeRegressor(min_samples_leaf=2)
        cls.fit(X, Y)
        cache(cls, "tree")
    else:
        cls = uncache("tree")

    # export to dot
    dot_data = export_graphviz(cls, out_file=None)
    for val in range(0, len(feature_names) - 1):
        dot_data = dot_data.replace("X[" + str(val) + "] <= 0.5", str(feature_names[val]))
    dot_data = dot_data.replace("True", "Does not contain")
    dot_data = dot_data.replace("False", "Does contain")
    dot_data = dot_data.replace("\'", "")
    graph = graphviz.Source(dot_data)
    graph.render("tree")


def frequent_itemsets_apriori(df, cache_results=True):
    """
    Runs all the apriori algorithm edited to only compute frequent sets
    :param df: Dataframe we want to use for calculation
    :param cache_results: Boolean to see if we want to use cached transactions
    :return: None
    """
    if not is_cached("cache.txt") and not is_cached("final_list_cache.txt"):
        mortality = create_y()
        date_dna_list = make_date_dna_list(df)
        del df["Geo_Location"]
        transactions = create_chunks(list(df.itertuples(index=False)))
        transactions = filter_transactions(transactions)
        create_final_list(mortality, transactions, date_dna_list, True)
        final_list = create_final_list(mortality, transactions, date_dna_list)
        if cache_results:
            cache(transactions, 'cache.txt')
            cache(final_list, 'final_list_cache.txt')
    else:
        transactions = uncache('cache.txt')
        final_list = uncache('final_list_cache.txt')
    make_tree(final_list)
    # write_apriori_results(transactions)


def frequent_itemsets_apriori_by_month(df, cache_results=True):
    """
    Get the frequent itemsets by month
    :param df: dataframe to get frequent itemsets from
    :param cache_results: if true, results get cached
    :return: list of frequent itemsets
    """
    df.info()
    df['Collection_Date'] = pd.to_datetime(df['Collection_Date'])
    df.info()
    list_of_dataframes = [df.loc[(df.Collection_Date.dt.month == _)] for _ in range(1, 10)]

    result_list = []
    for index, new_df in enumerate(list_of_dataframes):
        if len(list(new_df['DNA'])):
            transactions = create_chunks(make_date_dna_list(new_df))
            transactions = filter_transactions(transactions)
            results = write_apriori_results(transactions)
            result_list.append(results)

    final = []
    fin_dups = []
    for item in result_list:
        set1 = set()
        item = item[0]
        for tmp in item:
            t = set(item[tmp].items())
            set1 = set1.union(t)

        dups_tmp = set()
        for item2 in result_list:
            item2 = item2[0]
            if item != item2:
                set2 = set()
                for tmp in item2:
                    t = set(item2[tmp].items())
                    set2 = set2.union(t)
                duplicates = set1.intersection(set2)
                set1 = set1 - set2
                dups_tmp = dups_tmp.union(duplicates)
        fin_dups.append(dups_tmp)
        final.append(set1)

    with open("output/reee.txt", "w") as f:
        for index, item in enumerate(fin_dups):
            f.write("Dups in month " + str(index+1) + ": \n\n")
            for item2 in item:
                f.write(str(item2) + "\n")
            f.write("\n\n\n\n\n")

    with open("output/testresult.txt", "w") as f:
        for index, set_ in enumerate(final):
            f.write("Month " + str(index + 1) + ":\n")
            for item in set_:
                f.write(str(item) + "\n")
            f.write("\n\n--------\n\n")


if __name__ == '__main__':
    df = create_dataframe()
    print(df)
    frequent_itemsets_apriori(df)
