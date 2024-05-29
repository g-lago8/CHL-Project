import sys
import pandas as pd
import numpy as np

def ss_to_array(ss):
     # count number of each type of secondary structure
    vocab = list(set(df['secondary_structure'][0]))
    dict_vector = {v:0 for v in vocab}
    for i in ss:
        dict_vector[i] += 1   
    return dict_vector

def count_structures(ss):
    # count number of each type of secondary structure e.g. HHHHEEEEHHHH -> H:2, E:1
    vocab = list(set(df['secondary_structure'][0]))
    dict_vector = {v:0 for v in vocab}
    last_ss = ss[0]
    dict_vector[last_ss] = 1
    for new_ss in ss[1:]:
        if new_ss != last_ss:   
            last_ss = new_ss
            dict_vector[new_ss] += 1
    return dict_vector

if __name__ == '__main__':
    df=pd.read_csv('../datasets/pdb_files_ss.csv')
    df['secondary_structure_count'] = df['secondary_structure'].apply(count_structures)
    for key in df['secondary_structure_count'][0]:
        df[key] = df['secondary_structure_count'].apply(lambda x: x[key])
    df.drop(columns=['secondary_structure_count'], inplace=True)
    df.to_csv('../datasets/pdb_files_ss_count.csv', index=False)