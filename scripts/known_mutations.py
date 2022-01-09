import map_mutation_positions
import pandas as pd     
import numpy as np

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

def mutation_pos_cons_score(mutation_position_list, pos_alignedPos_dict, conservation_table):
    cs_list = []
    for position in mutation_position_list:
        position -= 1 
        aligned_position = pos_alignedPos_dict[position]
        row = conservation_table.loc[aligned_position]  
        conservation_score = row[2]
        cs_list.append(conservation_score)
    return cs_list



def main():
    filename = "../paper_mutations.txt"
    mutations = open(filename, 'r').read()
    mutations = mutations.split(",")
    position_list = []
    #
    filename = "../mnx1_proteinSequence.fas"
    original_sequence = list(fastareader(filename).values())[0]
    #
    for mutation in mutations:
        substr = mutation[3:]
        if substr[1].isalpha(): #Then the mutation position is a 1-digit number
            position = substr[0]
        elif substr[2].isalpha(): #Then the mutation position is a 2-digit number
            position = substr[:2]
        else: #Otherwise
            position = substr[:3]
        if int(position) <= len(original_sequence):
            position_list.append(int(position))
    aligned_filename = "../aligned/fas/aligned-500.fas"
    position_dict = map_mutation_positions.map_mutation_position(aligned_filename)
    conservation_scores = pd.read_csv("../conservation/conservation_scores_s=500.tsv", sep = "\t")
    cons_score_list = mutation_pos_cons_score(position_list, position_dict, conservation_scores)
    #
    aligned_filename = "../aligned/fas/pruned-500.fas"
    position_dict = map_mutation_positions.map_mutation_position(aligned_filename)
    conservation_scores = pd.read_csv("../conservation/conservation_scores_pruned_s=500.tsv", sep = "\t")
    cons_score_list2 = mutation_pos_cons_score(position_list, position_dict, conservation_scores)
    #
    aligned_filename = "../aligned/fas/realigned-pruned-500.fas"
    position_dict = map_mutation_positions.map_mutation_position(aligned_filename)
    conservation_scores = pd.read_csv("../conservation/conservation_scores_pruned_realigned_s=500.tsv", sep = "\t")
    cons_score_list3 = mutation_pos_cons_score(position_list, position_dict, conservation_scores)
    #
    thresholds = open("../conservation_thresholds.txt", "w")
    thresholds.write(str(np.mean(cons_score_list))+" " + str(np.mean(cons_score_list2)) + " " +str(np.mean(cons_score_list3)))
    thresholds.close()

if __name__ == "__main__":
	main()