# COVID19_patterns
Project for the course 'Introduction to Computational Biology'.<br>
By `Jarne Colman`, `Senne Rosaer`, `Toon Meynen` and `Freek De Sagher`

## Installation
In order to install this project on your local machine, follow the next steps:
1. Clone or download this repository
2. Navigate to the root directory of the project
2. Install the requirements with `pip install -r requirements.txt`
3. The program can be executed using `python main.py`

## Data Used in the Report
This repository contains the datasets that were used for the report.
The datasets are located in the [`input`](https://github.com/SenneRosaer/COVID19_patterns/tree/master/input) directory.
This folder contains the mortality data and the aligned sequences of both China and the United States.

## Brief Explanation of the Code
The two most important functions are the `frequent_itemsets_apriori` and the `frequent_itemsets_apriori_by_month`.
The first one runs the apriori algorithm on the filtered data and also computes a regression tree which gives a clear representation of the data.
The second function also performs apriori, but instead of running it on all the data is runs seperate for every month for better comparisons.

## Run this Code with Newer Data
Input files will always be in the input directory.
In order to run this code with other/newer data, the following functions need to be adapted.

#### `create_dataframe()`
The `create_dataframe()` function requires a `.csv` file that requires the following columns: `Accession`, `Geolocation` and `Collection_Date` (line 24). This data is taken from the NCBI database. Further, this function also needs an `.aln` alignment file which contains the aligned DNA sequences of the country/countries you want to investigate (line 26).

#### `create_y()`
This function uses the mortality data from `owid-covid-data.csv`. This data comes from the [Github Repository](https://github.com/owid/covid-19-data/tree/master/public/data) of 'Our World In Data'. This file contains all information about the counted COVID-19 cases. 
