
import pandas as pd
from Bio.Seq import MutableSeq
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
from Bio.PDB import PDBParser 


def replace(sequence: MutableSeq, mutation:str, verbose = False) -> MutableSeq:
    """
    takes a MutableSequence and a string  of type "A000B" and performs
    """

    original = mutation[0]
    mut = mutation[-1]
    sequence_copy = sequence[:]
    try:
        position = int(mutation[1:-1])
    except:
        if verbose:
            print('Invalid mutation format: should be AnB where n is an integer, representing the position, A and B two letters representing an aminoacid')
        return
    
    if len(sequence_copy)< position-1:
        if verbose:
            print('Sequence too short')
        return
    print
    if sequence_copy[position-1] != original.upper():
        if verbose:
            print(f' aminoacid in the sequence and in the mutation do not correspond: original aminoacid: {sequence_copy[position-1]}')
        return
    
    sequence_copy[position-1] = mut.upper()
    return sequence_copy


def create_graph_df(pdb_path ="../datasets/pdb_files.csv", akussy_path ='../datasets/aku_prin_v2.0.xlsx', config=ProteinGraphConfig() ):

    df = pd.read_csv(pdb_path)
    graphs = {}
    structures = {}
    parser = PDBParser()
    for i, row in df.iterrows():
        structures[row['mutation']] = parser.get_structure(row['mutation'], row['pdb_file'])
    for i, row in df.iterrows():
        #print(row['mutation'])
        graphs[row['mutation'] ] = construct_graph(path = row['pdb_file'], config= config, )
    df_patients =pd.read_excel(akussy_path)
    df_patients = df_patients[['Protein change allele 1 ', 'Protein change allele 2']]

    df_patients['graph_allele1'] = [graphs[mut] if mut in graphs else None for mut in df_patients['Protein change allele 1 '] ]
    df_patients['graph_allele2'] = [graphs[mut] if mut in graphs else None for mut in df_patients['Protein change allele 2'] ]
    df_patients['structure_allele1'] = [structures[mut] if mut in structures else None for mut in df_patients['Protein change allele 1 '] ] 
    df_patients['structure_allele2'] = [structures[mut] if mut in structures else None for mut in df_patients['Protein change allele 2'] ]
    return df_patients


def get_subgraph(g, nodes):
    SG = g.__class__()
    SG.add_nodes_from((n, g.nodes[n]) for n in nodes)
    if SG.is_multigraph():
        SG.add_edges_from(
            (n, nbr, key, d)
            for n, nbrs in g.adj.items()
            if n in nodes
            for nbr, keydict in nbrs.items()
            if nbr in nodes
            for key, d in keydict.items()
        )
    else:
        SG.add_edges_from(
            (n, nbr, d)
            for n, nbrs in g.adj.items()
            if n in nodes
            for nbr, d in nbrs.items()
            if nbr in nodes
        )
    SG.graph.update(g.graph)
    return SG

def get_edges_from_nodes(g, nodes):
    active_edges = {}
    nodes_name = []
    for i in nodes:
        active_edges[i] = []

    for node in g.nodes:
        #get the int of the node at the end of the string
        nnode = int(node.split(':')[2])
        if nnode in nodes:
            nodes_name.append(node)
            # get the edges of the node
            for edge in g.edges(node, data=True):
                active_edges[nnode].append(edge)

    return nodes_name, active_edges

def get_edge_list_from_edgeData(node_list ,edgeData):
    active_edges_num = {}
    for i in node_list:
        active_edges_num[i] = []

    for key, value in edgeData.items():
        for edge in value:
            active_edges_num[key].append(edge[1].split(':')[2])
