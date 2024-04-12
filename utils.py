

from Bio.Seq import MutableSeq

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