from scipy import stats as st
import pandas as pd 

df = pd.read_csv("../classification/classification_pruned_realigned_s=500.csv")

a = df.loc[df['isPathogenic'] == True, 'Allele_frequency'].to_numpy()
b = df.loc[df['isPathogenic'] == False, 'Allele_frequency'].to_numpy()

pvalue = st.ttest_ind(a=a, b=b, equal_var = True).pvalue

