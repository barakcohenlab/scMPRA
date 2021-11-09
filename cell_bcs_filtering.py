import os, sys
import pandas as pd
import numpy as np
import pickle
import argparse

'''
It feels like the progress in this analysis pipeline is a bit 
too slow, but it is always slow.
'''

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

def ec_cellBC(cellBC, cr_cellBC_list):
    pop_list = []
    for cr_cellBC in cr_cellBC_list:
        hamming = hammingDist(cellBC, cr_cellBC)
        if hamming <= 3:
            pop_list.append([cr_cellBC, hamming])
    if len(pop_list) == 1:
        return pop_list[0][0], pop_list[0][1]
    else:
        #print(f'{cellBC} is a bad cell!')
        return 0, 16

def cell_bc_ec(quad, cr_cellBC_list):
    '''
    Input: quad: tsv file with columns: cellBC, UMI, pBC, rBC, and read count
        cellBC: a tsv file with cellranger corrected cell barcodes.
    Output: cellranger error corrected quad file with the same structure.
    '''
    barcode_mapping_dict = {}
    # First extract the quad_cellBC set from quad file
    quad_cellBC = list(set(quad['cellBC'].values))
    # Then error correct the quad_cellBCs if it is in the cellBC list or within
    # If the cellBC is the exact match to the cr_cellBC, then the cellBC is kept
    for cellBC in quad_cellBC:
        if cellBC in cr_cellBC_list:
            barcode_mapping_dict[cellBC] = cellBC
        else:
            close_mapping, hamming = ec_cellBC(cellBC, cr_cellBC_list)
            if hamming <= 3:
                if close_mapping != 0:
                    barcode_mapping_dict[cellBC] = close_mapping
    # Then we go through all the list of cellBCs in the quad file and only record 
    # the ones whose cellBC is in the barcode_mapping_dict!
    pop_list = []
    for _, line in quad.iterrows(): 
        cellBC = line['cellBC']
        if cellBC in barcode_mapping_dict:
            pop_list.append([barcode_mapping_dict[cellBC], line['umi'], line['pBC'], line['rBC'], line['counts']])
    pop_quad = pd.DataFrame(pop_list,columns= ['cellBC', 'umi', 'pBC', 'rBC', 'counts'])
    return pop_quad

def main():
    # Input arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quad', help='path to quad file')
    parser.add_argument('--cellBC', help='path to cellranger tsv.gz cell barcode')
    parser.add_argument('--doc_name', help='name of document')
    args = parser.parse_args()
    # Main function
    # Message cellranger cellBC file
    cellBC = pd.read_csv(args.cellBC, sep = '\t' , header = None)
    # Remove the -1 thing in cellBC from cellranger
    cellBC = cellBC.replace('-1', '', regex=True).astype(str)
    cellBC = list(cellBC[0].values)
    # Read in the quad file 
    quad = pd.read_csv(args.quad, sep = '\t')
    corrected_quad = cell_bc_ec(quad, cellBC)
    # Save correct quad as a tsv file
    corrected_quad.to_csv('../cellranger_ec_' + args.doc_name + '.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()
