# COVID19_patterns
Project for the course 'Introduction to Computational Biology'.<br>
By `Jarne Colman`, `Senne Rosaer`, `Toon Meynen` and `Freek De Sagher`

The two most important functions are the frequent_itemsets_apriori and the frequent_itemsets_apriori_by_month.
The first one runs apriori on the filtered data and also computes a regression tree which gives more clear data.
The second one is also apriori but instead of running it on all the data is runs seperate for every month for better comparisons.

Input files will always be in the input directory
Selecting the data to be used is done in create_dataframe() and create_y(). In create_dataframe() there needs to be a csv with the id, collection date, geo-location etc. (taken from NCBI database) and also an alignment file. In create_y() there needs to be a csv with the all the information about counted cases.
