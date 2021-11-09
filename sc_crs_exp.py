'''
This script only deals with input quint containing preassigned clusters based on 
the transcriptome data.
'''
import os,sys
from numpy.core.defchararray import count
import pandas as pd
import numpy as np
import gzip
from scipy import stats
import pickle
import argparse

def return_abundant_cells(df, min_umi = 100):
    umi_dict = {}
    for _, row in df.iterrows():
        cell_bc = row['cellBC']
        if cell_bc not in umi_dict:
            umi_dict[cell_bc] = 1
        else:
            umi_dict[cell_bc] +=1
    abundant_list = []
    for key in umi_dict:
        if umi_dict[key] > min_umi:
            abundant_list.append(key)
    # Slice the dataset
    pop_df = df[df['cellBC'].isin(abundant_list)]
    return pop_df

def attach_promoter_to_quad(quad_df, info_df):
    pop_list = []
    for _, line in quad_df.iterrows():
        pBC = line['pBC']
        if len(info_df[info_df['pBC'] == pBC]['prom_id']) != 0:
            prom_id = str(info_df[info_df['pBC'] == pBC]['prom_id'].values[0])
            pop_list.append([line['cellBC'], line['umi'], prom_id, line['pBC'], line['rBC'], line['counts']])
    sexa = pd.DataFrame(pop_list,columns= ['cellBC', 'umi', 'prom_id', 'pBC', 'rBC', 'counts'])
    return sexa

def extract_sc_pBC_exp(sexa_df):
    '''
    Input: a haxa-column def contains scMRPA information
    Output: df with [cellBC, prom_id, pBC, norm_exp]
    '''
    pop_dict = {} 
    # Get the clusters
    clusters = set(sexa_df['cluster'].values)
    # Iterate through different clusters -- aka different cell type/states
    for cluster in clusters:
        slice_df = sexa_df[sexa_df['cluster'] == cluster]
        for _,line in slice_df.iterrows():
            temp_key = line[['cellBC', 'prom_id','pBC']]
            key = tuple(temp_key.values)
            rBC = str(line[['rBC']].values)
            if key not in pop_dict:
                pop_dict[key] = {}
                if rBC not in pop_dict[key]:
                    pop_dict[key][rBC] = 1
                else:
                    pop_dict[key][rBC] += 1
            else:
                if rBC not in pop_dict[key]:
                    pop_dict[key][rBC] = 1
                else:
                    pop_dict[key][rBC] += 1
        pop_list = []
        for key in pop_dict:
            info = list(key)
            num_plasmid = len(pop_dict[key])
            exp = sum(pop_dict[key].values())
            norm_exp = exp/num_plasmid
            info.append(norm_exp)
            info.append(exp)
            pop_list.append(info)
    sc_measurement_df = pd.DataFrame(pop_list,  columns= ['cellBC', 'prom_id', 'pBC', 'norm_exp' , 'direct_exp'])
    return sc_measurement_df

def extract_sc_pBC_normed_exp(sexa_df):
    '''
    Input: a haxa-column def contains scMRPA information
    Output: df with [cellBC, prom_id, pBC, norm_exp]
    '''
    pop_dict = {} 
    # Get the clusters
    #clusters = set(sexa_df['cluster'].values)
    cell_umi_dict = {}
    slide_df_cells = list(set(sexa_df['cellBC'].to_list()))
    for cell in slide_df_cells:
        total_umis = len(sexa_df[sexa_df['cellBC'] == cell])
        total_plasmid = len(set(sexa_df[sexa_df['cellBC'] == cell]['rBC'].to_list()))
        cell_umi_dict[cell] = total_umis/total_plasmid
    # Iterate through different clusters -- aka different cell type/states
    slice_df = sexa_df
    # A helper dict the records the total number of UMI for each cell
    for _,line in slice_df.iterrows():
        temp_key = line[['cellBC', 'prom_id','pBC']]
        key = tuple(temp_key.values)
        rBC = str(line[['rBC']].values)
        if key not in pop_dict:
            pop_dict[key] = {}
            if rBC not in pop_dict[key]:
                pop_dict[key][rBC] = 1
            else:
                pop_dict[key][rBC] += 1
        else:
            if rBC not in pop_dict[key]:
                pop_dict[key][rBC] = 1
            else:
                pop_dict[key][rBC] += 1
    pop_list = []
    for key in pop_dict:
        info = list(key)
        current_cell = info[0]
        normalization_factor = cell_umi_dict[current_cell]
        num_plasmid = len(pop_dict[key])
        exp = sum(pop_dict[key].values())
        # The proper normalization for the expression in a cell is the total
        # number of UMI divided by the number of plasmids it has 
        norm_exp = ((exp/num_plasmid)/normalization_factor)*1000
        info.append(norm_exp)
        info.append(exp)
        pop_list.append(info)
    sc_measurement_df = pd.DataFrame(pop_list,  columns= ['cellBC', 'prom_id', 'pBC', 'norm_exp' , 'direct_exp'])
    return sc_measurement_df

def filter_based_on_umi(sext_df, min_count = 2):
    '''
    Input：sext_df, a six-column table of cellBC, prom_id, pBC, rBC, umi, counts, and clusters the cell come from.
    Output: the unique sext that are supported by more than 1 read. 
    '''
    cell_list = list(set(sext_df['cellBC'].values))
    pop_list = []
    for cells in cell_list:
        cell_df = sext_df[sext_df['cellBC'] == cells]
        cell_list = list(set(sext_df['cellBC'].values))
        high_count_slice = cell_df[cell_df['counts'] > min_count]
        high_count_rBC = set(high_count_slice['rBC'].values)
        for _, row in cell_df.iterrows():
            if row['rBC'] in high_count_rBC:
                pop_list.append([row['cellBC'], row['prom_id'],row['pBC'],row['rBC'], row['umi'], row['counts'] ])
            else:
                possible_rBC = calculate_hamming(row['rBC'], high_count_rBC)
                if possible_rBC != 0:
                    pop_list.append([row['cellBC'], row['prom_id'],row['pBC'],row['rBC'], row['umi'], row['counts']])
    return pd.DataFrame(pop_list, columns = ['cellBC', 'prom_id', 'pBC', 'rBC', 'umi', 'counts'])


def calculate_hamming(rBC1, good_list):
    for rBC in good_list:
        hamming = hammingDist(rBC1, rBC)
        if hamming < 2:
            return rBC
        else:
            return 0


def hammingDist(str1, str2):
    '''
    Calculating hamming distance of two strings
    https://www.geeksforgeeks.org/hamming-distance-two-strings/
    '''
    i = 0
    count = 0
    while(i < len(str1)):
        if(str1[i] != str2[i]):
            count += 1
        i += 1
    return count

def extract_exp_distribution(sc_df):
    '''
    Input: dataframe with columns: cellBC, prom_id, pBC, norm_exp
    Output: dataframe with dictionary: prom_id: {expression per cell}
    '''
    pop_dict = {}
    prom_list = list(set(sc_df['prom_id'].values))
    for prom in prom_list:
        prom_slice = sc_df[sc_df['prom_id'] == prom]
        cellBC_list = list(set(prom_slice['cellBC'].values))
        exp_list = []
        for cellBC in cellBC_list:
            cell_slice = prom_slice[prom_slice['cellBC'] == cellBC]
            pBC_list = list(cell_slice['norm_exp'].values)
            exp_list.append(np.mean(pBC_list))
        pop_dict[prom] = exp_list
    return pop_dict

def return_pBC_mean_exp_percell(cell_pBC_exp):
    prom_list = list(set(cell_pBC_exp['prom_id'].to_list()))
    pop_list = []
    for inx, prom in enumerate(prom_list):
        if inx%50 == 0: 
            print(f'We are dealing with {inx}. ')
        prom_df = cell_pBC_exp[cell_pBC_exp['prom_id'] == prom]
        cell_list = list(set(prom_df['cellBC'].to_list()))
        cell_exp_list = []
        for cell in cell_list:
            cell_slice = prom_df[prom_df['cellBC'] == cell]
            cell_mean = np.sum(cell_slice['norm_exp'])
            cell_exp_list.append(cell_mean)
        prom_pBC_auc = np.sum(cell_exp_list)
        prom_pBC_mean = np.mean(cell_exp_list)
        prom_pBC_mad = stats.median_abs_deviation(cell_exp_list)
        prom_pBC_var = np.var(cell_exp_list)
        pop_list.append([prom, prom_pBC_mean, prom_pBC_auc, prom_pBC_mad, prom_pBC_var])
    pop_df = pd.DataFrame(pop_list, columns = ['prom_id', 'mean', 'auc', 'mad', 'var'])
    return pop_df

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    # Add the preprocess fileds calld 
    parser.add_argument('--quint', help='path to quint file with clusters', required=True)
    # Add the information about which promoter BC correlates with which promoter (genomic coordinates)
    parser.add_argument('--promlib', help='path to file named: prom_lib_master_info', required=True)
    # Add the name for the output file
    parser.add_argument('--out', help='output file name', required=True)
    parser.add_argument('--exp', help = 'name of the experiment', required = True)
    # Grab input arguments
    args= parser.parse_args()
    # Read in the quint file with columns as cellBC, umi, pBC,rBC, counts, cluster
    prom_quint = pd.read_csv(args.quint, sep= '\t')
    # Read in the promoter lib infor
    prom_lib_info = pd.read_csv(args.promlib, sep= '\t')
    # Add promoter id to the quint file
    prom_sext = attach_promoter_to_quad(prom_quint, prom_lib_info)
    # Return abundant cells more than x read per umi
    prom_filtered = filter_based_on_umi(prom_sext, min_count =0)
    # Return cells with more than 10 UMIs 
    prom_filtered_abundant = return_abundant_cells(prom_filtered, min_umi = 100)
    # Return the file containing the single-cell expression for each cBC
    sc_exp = extract_sc_pBC_normed_exp(prom_filtered_abundant)
    sc_exp.to_csv('./' + args.exp + '_cBC_exp_' + args.out + '.csv',index = False)
    # Return the expression mean, auc, median absolution deviation, variance for each CRS
    sc_stats = return_pBC_mean_exp_percell(sc_exp)
    sc_exp.to_csv('./'+ args.exp + '_cBC_stats_' + args.out + '.csv',index = False)
if __name__ == "__main__":
    main()


