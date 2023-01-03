from email.headerregistry import UniqueUnstructuredHeader
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter
from scipy import stats
#opens the excel sheet as a dataframe so we can work with it in Pandas
NTonly = pd.read_excel("new_nprotein_seq-2.xlsx")
NTseq = NTonly.iloc[:, [0]]
NTseq

#loops through each and makes it into a matrix so we can check by index if nucleotides differ/ each row appended to one list
NTseqM = []  # list of sequences
for i in range(len(NTseq)):
    NTseqM.append(NTseq.iat[i,0])


count = 0
for i in range(len(NTseqM)):
    if(NTseqM[0] == NTseqM[i]):
        count += 1
print(count)

#loop through seq, counting unique sequences in n_terminal and putting into dict
n_terminal_count = dict()
for seq in NTseqM:
    n_terminal = seq[45:167]
    n_terminal_count[n_terminal] = n_terminal_count.get(n_terminal, 0) + 1
print(n_terminal_count.values())

#loop through seq, counting unique sequences in c_terminal and putting into dict
c_terminal_count = dict()
for seq in NTseqM:
    c_terminal = seq[246:354]
    c_terminal_count[c_terminal] = c_terminal_count.get(c_terminal, 0) + 1
print(c_terminal_count.values())

#convert dict values to list
n_terminal = list(n_terminal_count.values())
c_terminal = list(c_terminal_count.values())

results = stats.levene(n_terminal, c_terminal, center = 'trimmed')
print("Levene Test: " + str(results))

#equalize the lengths of list
for x in range(len(n_terminal) - len(c_terminal)):
    c_terminal.append(0)

print(n_terminal)
print(c_terminal)

#perform ranksum
results2 = stats.ranksums(c_terminal, n_terminal, 'less')
print("RankSums Test: " + str(results2))