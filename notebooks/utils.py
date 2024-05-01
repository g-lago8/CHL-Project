
import pandas as pd
from Bio.Seq import MutableSeq
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
from Bio.PDB import PDBParser 


def replace(sequence: MutableSeq, mutation:str) -> MutableSeq:
    """
    takes a MutableSequence and a string  of type "A000B" and performs
    """

    original = mutation[0]
    mut = mutation[-1]

    try:
        position = int(mutation[1:-1])
    except:
        print('Invalid mutation format: should be AnB where n is an integer, representing the position, A and B two letters representing an aminoacid')
        return
    
    if len(sequence)< position-1:
        print('Sequence too short')
        return
    print
    if sequence[position-1] != original.upper():
        
        print(f'Original aminoacid in the sequence and in the mutation do not correspond: original aminoacid: {sequence[position-1]}')
        return
    
    sequence[position-1] = mut.upper()
    return sequence


def create_graph_df():
    config = ProteinGraphConfig()
    df = pd.read_csv("../datasets/pdb_files.csv")
    graphs = {}
    structures = {}
    parser = PDBParser()
    for i, row in df.iterrows():
        structures[row['mutation']] = parser.get_structure(row['mutation'], row['pdb_file'])
    for i, row in df.iterrows():
        #print(row['mutation'])
        graphs[row['mutation'] ] = construct_graph(path = row['pdb_file'], config= config, verbose=False)
    df_patients =pd.read_excel('../datasets/aku_prin_v2.0.xlsx')
    df_patients = df_patients[['Protein change allele 1 ', 'Protein change allele 2']]

    df_patients['graph_allele1'] = [graphs[mut] if mut in graphs else None for mut in df_patients['Protein change allele 1 '] ]
    df_patients['graph_allele2'] = [graphs[mut] if mut in graphs else None for mut in df_patients['Protein change allele 2'] ]
    df_patients['structure_allele1'] = [structures[mut] if mut in structures else None for mut in df_patients['Protein change allele 1 '] ] 
    df_patients['structure_allele2'] = [structures[mut] if mut in structures else None for mut in df_patients['Protein change allele 2'] ]
    return df_patients