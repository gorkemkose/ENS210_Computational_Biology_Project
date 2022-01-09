import pandas as pd  
import numpy as np
import known_mutations
import map_mutation_positions


variants = pd.read_csv("../gnomad_variants.csv")
cols_to_keep = ['HGVS Consequence', 'VEP Annotation', 'ClinVar Clinical Significance', 'Allele Count', 'Allele Number']
variants = variants[cols_to_keep]

missense_variants = variants.loc[variants['VEP Annotation'] == 'missense_variant']
missense_variants = missense_variants.loc[missense_variants['ClinVar Clinical Significance'].isnull()]

position_list = []
for i in list(missense_variants.index):
    row = missense_variants.loc[i]  
    mutation = row[0]
    substr = mutation[5:]
    if substr[1].isalpha(): #Then the mutation position is a 1-digit number
        position = substr[0]
    elif substr[2].isalpha(): #Then the mutation position is a 2-digit number
        position = substr[:2]
    else: #Otherwise
        position = substr[:3]
    position_list.append(int(position))

conservation_thresholds = open("../conservation_thresholds.txt", 'r').read().strip().split(" ")
conservation_thresholds = [float(score) for score in conservation_thresholds]

aligned_filename = "../aligned/fas/aligned-500.fas"
position_dict = map_mutation_positions.map_mutation_position(aligned_filename)
conservation_scores = pd.read_csv("../conservation/conservation_scores_s=500.tsv", sep = "\t")
cons_score_list = known_mutations.mutation_pos_cons_score(position_list, position_dict, conservation_scores)
threshold1 = conservation_thresholds[0]
isPathogenic = []
for score in cons_score_list: 
    if score >= threshold1:
        isPathogenic.append(True)
    else:
        isPathogenic.append(False)

missense_variants1 = missense_variants
missense_variants1['cons_score'] = cons_score_list
missense_variants1['isPathogenic'] = isPathogenic
missense_variants1['Allele_frequency'] = missense_variants1['Allele Count']/missense_variants1['Allele Number']
missense_variants1 = missense_variants1[['HGVS Consequence', 'cons_score', 'isPathogenic', 'Allele_frequency']]
missense_variants1.to_csv("../classification/classification_s=500.csv")

aligned_filename = "../aligned/fas/aligned-pruned-500.fas"
position_dict = map_mutation_positions.map_mutation_position(aligned_filename)
conservation_scores = pd.read_csv("../conservation/conservation_scores_pruned_s=500.tsv", sep = "\t")
cons_score_list = known_mutations.mutation_pos_cons_score(position_list, position_dict, conservation_scores)
threshold2 = conservation_thresholds[1]
isPathogenic = []
for score in cons_score_list: 
    if score >= threshold2:
        isPathogenic.append(True)
    else:
        isPathogenic.append(False)

missense_variants2 = missense_variants
missense_variants2['cons_score'] = cons_score_list
missense_variants2['isPathogenic'] = isPathogenic
missense_variants2['Allele_frequency'] = missense_variants2['Allele Count']/missense_variants2['Allele Number']
missense_variants2 = missense_variants2[['HGVS Consequence', 'cons_score', 'isPathogenic', 'Allele_frequency']]
missense_variants2.to_csv("../classification/classification_pruned_s=500.csv")


aligned_filename = "../aligned/fas/realigned-pruned-500.fas"
position_dict = map_mutation_positions.map_mutation_position(aligned_filename)
conservation_scores = pd.read_csv("../conservation/conservation_scores_pruned_realigned_s=500.tsv", sep = "\t")
cons_score_list = known_mutations.mutation_pos_cons_score(position_list, position_dict, conservation_scores)
threshold3 = conservation_thresholds[2]

isPathogenic = []
for score in cons_score_list: 
    if score >= threshold3:
        isPathogenic.append(True)
    else:
        isPathogenic.append(False)

missense_variants3 = missense_variants
missense_variants3['cons_score'] = cons_score_list
missense_variants3['isPathogenic'] = isPathogenic
missense_variants3['Allele_frequency'] = missense_variants3['Allele Count']/missense_variants3['Allele Number']
missense_variants3 = missense_variants3[['HGVS Consequence', 'cons_score', 'isPathogenic', 'Allele_frequency']]
missense_variants3.to_csv("../classification/classification_pruned_realigned_s=500.csv")
