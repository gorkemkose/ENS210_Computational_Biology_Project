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

filename = "../aligned/fas/aligned-500.fas"
#filename = "./aligned/fas/aligned-pruned-500.fas"
#filename = "./aligned/fas/realigned-pruned-500.fas"
seqDict = fastareader(filename)
sequences = list(seqDict.values())
aaList = ["A", "R", "N" , "D" , "C" , "Q" , "E" , "G" , "H" , "I" , "L" , "K" , "M" , "F" , "P" , "S" , "T" , "W" ,"Y" ,"V"]

consensus = ""
consensus_aminoacid_score = {}
for pos in range(len(sequences[0])):
    aa_percent_dict = dict.fromkeys(aaList, 0)
    aminoacids_onsamepos = ""
    for seq in sequences:
        aminoacids_onsamepos += seq[pos]
    for aa in aaList:
        count_aa = aminoacids_onsamepos.count(aa)
        aa_percent_dict[aa] = count_aa / len(sequences)
    aa_percent_dict = dict(sorted(aa_percent_dict.items(), key=lambda x: x[1], reverse=True))
    consensus_aa = list(aa_percent_dict.keys())[0]
    max_score = list(aa_percent_dict.values())[0]
    consensus += consensus_aa
    consensus_aminoacid_score[pos] = {consensus_aa:max_score}

conservation_scores_file = open("../conservation/conservation_scores_pruned_s=500.tsv", "w")
conservation_scores_file.write("Position"+"\t"+"Consensus Aminoacid"+"\t" +"Conservation Score" + "\n")

for k in consensus_aminoacid_score:
    inner_dict = consensus_aminoacid_score[k]
    conservation_scores_file.write(str(k) + "\t" + list(inner_dict.keys())[0] + "\t" + str(list(inner_dict.values())[0])+ "\n")

conservation_scores_file.close()