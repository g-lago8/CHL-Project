from utils import replace

#H292R (no), R321P (no), R330S (ok), P332R (ok), R336K (no), R336T (no), N337D (no), R347P (no),
# Y350C (no), K353Q (ok), P359L (no), H371R (no), G372R (no), P373L (no), D374H (no), E401Q (ok)
# mutations with no pdb files: H292R, R321P, R336K, R336T, N337D, R347P, Y350C, P359L, H371R, G372R, P373L, D374H
import Bio
from Bio import SeqIO

from Bio.Seq import Seq, MutableSeq

def read_faa_file(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(record)
    return sequences

# Specify the path to your .faa file
file_path = "..\\datasets\\HGD_datasets\\ncbi_dataset\\data\\protein.faa"


# Read the .faa file
sequences = read_faa_file(file_path)

hgd = sequences[0].seq
# convert the sequence to a mutable sequence
hgd = MutableSeq(str(hgd))
mutations_to_get = ['H292R', 'R321P', 'R336K', 'R336T', 'N337D', 'R347P', 'Y350C', 'P359L', 'H371R', 'G372R', 'P373L', 'D374H']

mutations=[]
for mutation in mutations_to_get:
    mutations.append(replace(hgd, mutation))

# save the mutations to a file
with open('mutations.txt', 'w') as f:
    for mutation in mutations:
        f.write("%s\n" % mutation)
