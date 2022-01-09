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
#filename = "./hits/hits-500.fas"
seqDict = fastareader(filename)

inFile = open("../trees/clade_sequences.txt")
clade_sequences = [sequence.strip() for sequence in inFile.readlines()]
clade_sequences = [sequence[: sequence.find("_", 3)] for sequence in clade_sequences]

for i  in range(len(clade_sequences)):
    sequence = clade_sequences[i]
    if "LOW QUALITY" in sequence or "PREDICTED" in sequence:
        word_list = sequence.split(" ") 
        identifier = word_list[0] + "_" + word_list[1] 
        clade_sequences[i] = identifier 

clade_sequences = [sequence.replace("_", " ") for sequence in clade_sequences]

outFile = open("../hits/hits-pruned-500.fas", "w")
for header in seqDict:
    word_list = header.split(" ")
    #match = word_list[0]
    match = word_list[0] + " " + word_list[1]
    if match in clade_sequences:
        outFile.write(">" + header + "\n" + seqDict[header] + "\n")

outFile.close()