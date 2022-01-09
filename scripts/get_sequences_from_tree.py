from ete3 import Tree

t = Tree("../trees/500_midpoint_rooted.nwk", format=1, quoted_node_names=True)

homo_sapiens = ['NP_001158727.1_motor_neuron_and_pancreas_homeobox_protein_1_isoform_2_Homo_sapiens', \
                 'NP_005506.3_motor_neuron_and_pancreas_homeobox_protein_1_isoform_1_Homo_sapiens']

nodes_to_keep = t.get_common_ancestor(homo_sapiens).get_leaves()
t.prune(nodes_to_keep)
single_hs = t.get_children()[1].get_leaves()
other_hs = t.get_children()[0].get_leaves()
if len(other_hs) > len(single_hs):
    single_hs = other_hs
t.prune(single_hs)

t.write(format=1, outfile="pruned_500.nwk")

inFile = open("../trees/clade_sequences.txt", "w")
for sequence in list(single_hs):
    inFile.write(sequence.name + "\n")

inFile.close()