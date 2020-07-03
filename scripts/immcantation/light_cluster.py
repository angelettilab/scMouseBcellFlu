#!/usr/bin/env python3
"""
Corrects IGH only cloning with IGK/L annotations (VERSION 1)
"""

# Imports
import os
import pandas as pd
import sys
from argparse import ArgumentParser

# Presto and changeo imports
from changeo.Gene import getGene


def clusterLinkage(cell_series, group_series):
    """
    Returns a dictionary of {cell_id : cluster_id} that identifies clusters of cells by analyzing their shared
    features (group_series) using single linkage. 

    Arguments:
      cell_series (iter): iter of cell ids.
      group_series (iter): iter of group ids.

    Returns:
      dict:  dictionary of {cell_id : cluster_id}.
    """

    # assign initial clusters
    # initial_dict = {cluster1: [cell1], cluster2: [cell1]}
    initial_dict = {}
    for cell, group in zip(cell_series, group_series):
        try:    
            initial_dict[group].append(cell)
        except KeyError:
            initial_dict[group] = [cell]
               
    # naive single linkage clustering (ON^2 best case, ON^3 worst case) ...ie for cells with multiple light chains
    # cluster_dict = {cluster1: [cell1, cell2]}, 2 cells belong in same group if they share 1 light chain 
    while True:
        cluster_dict = {}
        for i, group in enumerate(initial_dict.keys()):
            cluster_dict[i] = initial_dict[group]
            for cluster in cluster_dict:
                # if initial_dict[group] and cluster_dict[cluster] share common cells, add initial_dict[group] to cluster
                if cluster != i and any(cell in initial_dict[group] for cell in cluster_dict[cluster]):
                    cluster_dict[cluster] = cluster_dict[cluster] + initial_dict[group]
                    del cluster_dict[i]
                    break
        # break if clusters stop changing, otherwise restart 
        if len(cluster_dict.keys()) == len(initial_dict.keys()):
            break
        else:
            initial_dict = cluster_dict.copy()
    
    # invert cluster_dict for return
    assign_dict = {cell:k for k,v in cluster_dict.items() for cell in set(v)}
    
    return assign_dict


def lightCluster(heavy_file, light_file, out_file, doublets='drop', format='airr'):
    """
    Split heavy chain clones based on light chains

    Arguments:
      heavy_file (str): heavy chain input file.
      light_file (str): light chain input file.
      out_file (str): heavy chain output file.
      doublets (str): method for handling multiple heavy chains per cell. one of 'drop' or 'count'.
      format (str): file format. one of 'changeo' or 'airr'.
    """
    # Set column names
    if format == 'changeo':
        cell_id = 'CELL'
        clone_id = 'CLONE'
        v_call = 'V_CALL'
        j_call = 'J_CALL'
        junction_length = 'JUNCTION_LENGTH'
        umi_count = 'UMICOUNT'
    elif format == 'airr':
        cell_id = 'cell_id'
        clone_id = 'clone_id'
        v_call = 'v_call'
        j_call = 'j_call'
        junction_length = 'junction_length'
        umi_count = 'umi_count'
    else:
        sys.exit("Invalid format %s" % format)

    # read in heavy and light DFs
    heavy_df = pd.read_csv(heavy_file, dtype='object', na_values=['', 'None', 'NA'], sep='\t')
    light_df = pd.read_csv(light_file, dtype='object', na_values=['', 'None', 'NA'], sep='\t')

    # column checking
    expected_heavy_columns = [cell_id, clone_id, v_call, j_call, junction_length, umi_count]
    if set(expected_heavy_columns).issubset(heavy_df.columns) is False:
        raise ValueError("Missing one or more columns in heavy chain file: " + ", ".join(expected_heavy_columns))
    expected_light_columns = [cell_id, v_call, j_call, junction_length, umi_count]
    if set(expected_light_columns).issubset(light_df.columns) is False:
        raise ValueError("Missing one or more columns in light chain file: " + ", ".join(expected_light_columns))

    # Fix types
    heavy_df[junction_length] = heavy_df[junction_length].astype('int')
    light_df[junction_length] = light_df[junction_length].astype('int')

    # filter multiple heavy chains
    if doublets == 'drop':
        heavy_df = heavy_df.drop_duplicates(cell_id, keep=False)
        if heavy_df.empty is True:
            raise ValueError("Empty heavy chain data, after doublets drop. Are you combining experiments in a single file? If so, split your data into multiple files.")
    elif doublets == 'count':
        heavy_df[umi_count] = heavy_df[umi_count].astype('int')
        heavy_df = heavy_df.groupby(cell_id, sort=False).apply(lambda x: x.nlargest(1, umi_count))

    # transfer clone IDs from heavy chain df to light chain df
    clone_dict = {v[cell_id]:v[clone_id] for k, v in heavy_df[[clone_id, cell_id]].T.to_dict().items()}
    light_df = light_df.loc[light_df[cell_id].apply(lambda x: x in clone_dict.keys()), ]
    light_df[clone_id] = light_df.apply(lambda row: clone_dict[row[cell_id]], axis = 1)

    # generate a "cluster_dict" of CELL:CLONE dictionary from light df  (TODO: use receptor object V/J gene names)
    cluster_dict = clusterLinkage(light_df[cell_id],
                                  light_df.apply(lambda row:
                                                 getGene(row[v_call]) + ',' + \
                                                 getGene(row[j_call]) + ',' + \
                                                 str(row[junction_length]) + ',' + row[clone_id], axis=1))

    # add assignments to heavy_df
    heavy_df = heavy_df.loc[heavy_df[cell_id].apply(lambda x: x in cluster_dict.keys()), :]
    heavy_df[clone_id] = heavy_df[clone_id] + '_' + heavy_df.apply(lambda row: str(cluster_dict[row[cell_id]]), axis=1)

    # write heavy chains
    heavy_df.to_csv(out_file, sep='\t', index=False)


if __name__ == "__main__":
    """
    Parses command line arguments and calls main
    """
    # Define arguments
    parser = ArgumentParser()
    parser.add_argument('-d', dest='heavy_file', required=True,
                        help='Cloned heavy chain Change-O or AIRR TSV file.')
    parser.add_argument('-e', dest='light_file', required=True,
                        help='Corresponding Light chain Change-O or AIRR TSV file.')
    parser.add_argument('-o', dest='out_file', required=True,
                        help='Output file name.')
    parser.add_argument('--doublets', dest='doublets', default='drop', choices=('drop', 'count'),
                        help='Either drop cells with multiple heavy chains or keep the best one my UMI count.')
    parser.add_argument('--format', dest='format', default='airr', choices=('airr','changeo'),
                        help='File format.')

    # Parse arguments and call main
    args = parser.parse_args()

    # Check that files exist
    for f in [args.heavy_file, args.light_file]:
        if not os.path.isfile(f):  sys.exit('File %s does not exist.' % f)

    lightCluster(**args.__dict__)
