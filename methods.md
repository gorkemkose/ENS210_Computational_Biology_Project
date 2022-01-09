# Methods 
All Python and R scripts are located inside the folder scripts. 

After retrieving homologous sequences from BlastP, you can use MEGA11 to align the protein sequences. Then, again using MEGA11, you can construct the phylogenetic tree and apply midpoint rooting in FigTree. 

By running the get_sequences_from_tree.py, you can subset the clade of interest and by supplying the pruned tree to subset_clade_hits.py, you can construct fasta files which only contains the sequences in the pruned tree. 

You can use the conservation_score.py script to compute conservation scores. You should use the fasta files located at ./aligned/fas directory. All the files in this directory are constructed using the steps above. The conservation_score.py script generates conservation tables which are stored inside ./conservation directory. 

You can use conservation_histogram.py script in order to observe the distribution of conservation scores, which the results are also located at ./conservation directory. You can use conservation_plot.R script to see the conservation scores by position, which the results are also located at ./conservation directory. 

Using map_mutation_positions.py script in known_mutations.py, you can get the conservation scores for the mutations that are known to be pathogenic. The thresholds for the classification purpose are stored inside conservation_thresholds.txt inside ./classification directory.

Using this conservation_thresholds.txt file on gnomad_variation_analysis.py, you can classify the mutations of unknown pathogenic significance as either pathogenic or neutral and these classifications are stored inside ./classification directory. 

In order to assess if there is any statistical difference between the allele frequencies of pathogenic and neutral mutations, you can conduct a t-test using t-test.py script. Also using the allele_frequency_correlation.R, you can plot a scatterplot of conservation scores and allele frequencies in order to observe any correlation between those variables.  