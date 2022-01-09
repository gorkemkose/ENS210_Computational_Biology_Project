import pandas as pd     

def fastareader(filename):
    seqDict = {}
    infile = open(filename, 'r')
    for line in infile:
        if line[0] == ">":
            header = line.strip()[1:]
            seqDict[header] = ""
        else:
            seqDict[header] += line.strip()
    return seqDict

def map_mutation_position(aligned_filename):
    filename = "../mnx1_proteinSequence.fas"
    original_sequence = list(fastareader(filename).values())[0]
    sequence_dict = fastareader(aligned_filename)
    aligned_sequence = sequence_dict['NP 005506.3 motor neuron and pancreas homeobox protein 1 isoform 1 Homo sapiens']
    originalPos_alignedPos_dict = {}
    modified_aligned = aligned_sequence
    for i in range(len(original_sequence)):
        aa = original_sequence[i]
        position = modified_aligned.index(aa) 
        originalPos_alignedPos_dict[i] = position
        modified_aligned = modified_aligned[:position] + "-" + modified_aligned[position+1:]
    return originalPos_alignedPos_dict

