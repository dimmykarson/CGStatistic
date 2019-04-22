import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from util import kruskal_test
from rpy2.robjects.packages import importr


class CGStatistic():
    def __init__(self, dataframe_path, data_frame_columns,  figure_dir="", R_lib="PMCMR"):
        """
            Parameters
            ----------
            dataframe_path : str
                Path for dataframe results
            data_frame_columns : list
                List of dataframe column name
            figure_dir : str
                Path for output graphs
                

        """
        self.dataframe_path = dataframe_path
        self.data_frame_columns = data_frame_columns
        self.figure_dir = figure_dir
        self.dataframe = pd.read_csv(self.dataframe_path)
        self.dataframe.columns = self.data_frame_columns
        self.R_lib = R_lib
        self.post_test = self.R_lib.posthoc_kruskal_nemenyi_test


    def statistical_test_kruskal(self, unique_column, filter_columns, final_columns, result_metric_column, 
        filter_clause=None, group_by=None, aggregation=None, sort_by=True, verbose=False):
        """
            Run statistical test kruskal
            Parameters
            ----------
            unique_column : str
                column to filter into unique results
            filter_columns : list
                list of columns to filter in dataframe to plot
            filter_clause : clause
                clause for filter dataframe
            group_by : list
                List of columns in which dataframe that should be aggregated
            
        """
        df = self.dataframe.copy()
        if unique_column not in self.data_frame_columns:
            raise f"{unique_column} not in dataframe columns"

        unique_filter = df[unique_column].unique()
        for dt in unique_filter:
            x = None
            if filter_clause:
                x = df[[filter_columns]][filter_clause]
            else:
                x = df[[filter_columns]]

            if final_columns:
                x = x[final_columns]
            if not aggregation:
                aggregation = {final_columns:['mean', 'std']}
            if not group_by:
                raise "When need group_by column and aggregation model"
            mean = x.groupby([group_by], as_index=False).agg(aggregation)
            mean.columns = [result_metric_column, 'mean', 'std']
            if sort_by:
                mean = mean.sort_values(by='mean', ascending=False)
            x = x.sort_values(by=group_by)
            x.to_csv(f'{self.figure_dir}/stats_{dt}.csv', sep=";", index=False)

            if verbose:
                print("Means")
                print(mean)
                print(f"\nBest {result_metric_column}")
                bestp = mean.loc[mean['mean'].idxmax()]
                print(f"{bestp[result_metric_column]}: {bestp['mean']:.3f}")
            k = kruskal_test(x, result_metric_column, group_by)
            kruskal, posthoc = k.apply(filename=f'{self.figure_dir}/kruskal_{dt}')

            if verbose:
                if posthoc:
                    print(posthoc[0])
                    print(posthoc[1])
                    print("\nEffect-Size")
                    print(posthoc[2])












            






