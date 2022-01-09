import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

d = pd.read_csv("../conservation/conservation_scores_pruned_realigned_s=500.tsv", sep="\t")
scores = d['Conservation Score']

q25, q75 = np.percentile(scores, [25, 75])
bin_width = 2 * (q75 - q25) * len(scores) ** (-1/3)
bins = round((scores.max() - scores.min()) / bin_width)

plt.hist(scores, bins=bins)
plt.draw()
plt.xlabel("Conservation score")
plt.ylabel("Number of positions")
plt.savefig("../conservation/conservation_histogram_pruned_realigned_s=500.png")
plt.close()

#Reference for the correct number of bins
#https://stackoverflow.com/questions/33203645/how-to-plot-a-histogram-using-matplotlib-in-python-with-a-list-of-data