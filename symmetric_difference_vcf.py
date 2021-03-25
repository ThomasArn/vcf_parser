#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
import io

filename = sys.argv[1]

def read_vcf(filename):
    with open(filename, 'r') as f:
        lines = [line for line in f if not line.startswith('##')]
         dataframe = pd.read_csv(
            io.StringIO(''.join(lines)),sep='\t').drop(
            ['#CHROM', 'POS','ID', 'REF', 'ALT', 'QUAL', 'INFO', 'FORMAT'], axis=1)
        df_filtered  = dataframe[dataframe["FILTER"] == "PASS"].drop('FILTER',axis=1)
    return df_filtered

def symmetric_difference_pd(x, column1, column2):
    # compare the two columns and make an array with boolean values )
    comparison = np.where(x[column1] == x[column2], True, False)
    # calculate the numbers of different cases between the two columns
    number_diff_snps = comparison.size - np.count_nonzero(comparison)
    return number_diff_snps

vcf_dataframe = read_vcf(filename)

columns = list(vcf_dataframe)
dist_matrix = pd.DataFrame(columns=columns, index=columns)

for num_1, column_ref in enumerate(columns, start=1):
    for num_2, column_compare in enumerate(columns, start=1):

        diff_snps = symmetric_difference_pd(vcf_dataframe, column_ref, column_compare)
        #print(column_ref + "\t" + column_compare + "\t" + str(diff_snps)
        dist_matrix.at[[column_ref], [column_compare]] = diff_snps

dist_matrix.to_csv("distance_matrix.csv" )
