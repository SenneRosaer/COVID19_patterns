"""
Small adaptation of the efficient apriori implementation.
This version does not generate any rules.
"""

import typing
from efficient_apriori.itemsets import itemsets_from_transactions, ItemsetCount


def apriori2(
        transactions: typing.Union[typing.List[tuple], typing.Callable],
        min_support: float = 0.5,
        min_confidence: float = 0.5,
        max_length: int = 8,
        verbosity: int = 0,
        output_transaction_ids: bool = False,
):
    itemsets, num_trans = itemsets_from_transactions(
        transactions,
        min_support,
        max_length,
        verbosity,
        output_transaction_ids,
    )

    if itemsets and isinstance(next(iter(itemsets[1].values())), ItemsetCount):
        itemsets_for_rules = _convert_to_counts(itemsets)
    else:
        itemsets_for_rules = itemsets

    for itemset in itemsets.items():
        for itemset2 in itemset[1]:
            tmp = itemset[1][itemset2] / num_trans
            if tmp > 0.7 and tmp < 0.75:
                itemset[1][itemset2] = "[0.7 - 0.75]"
            elif tmp < 0.8:
                itemset[1][itemset2] = "[0.75 - 0.8]"
            elif tmp < 0.85:
                itemset[1][itemset2] = "[0.8 - 0.85]"
            elif tmp < 0.90:
                itemset[1][itemset2] = "[0.85 - 0.9]"
            else:
                itemset[1][itemset2] = "[0.9 - 1]"
    return itemsets, []


def _convert_to_counts(itemsets):
    itemsets_counts = {}
    for size, sets in itemsets.items():
        itemsets_counts[size] = {i: c.itemset_count for i, c in sets.items()}
    return itemsets_counts
