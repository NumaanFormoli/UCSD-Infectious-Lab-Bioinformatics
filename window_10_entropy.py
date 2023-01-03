from email.headerregistry import UniqueUnstructuredHeader
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter
from scipy import stats
#opens the excel sheet as a dataframe so we can work with it in Pandas
NTonly = pd.read_excel("new_nprotein_seq_2.xlsx")
NTseq = NTonly.iloc[:, [0]]
NTseq
#normalize lengths by adding 'N' to sequences shorter than the longest sequence
# lengths = []
# for i in range(len(NTseq)):
#     lengths.append(len(NTseq.iat[i,0]))

# for i in range(len(NTseq)):
#     if len(NTseq.iat[i,0]) < max(lengths):
#         toadd = max(lengths) - len(NTseq.iat[i,0])
#         for k in range(toadd):
#             NTseq.iat[i,0] += "Z"
#             k +=1
#loops through each and makes it into a matrix so we can check by index if nucleotides differ/ each row appended to one list
NTseqM = []  # list of sequences
for i in range(len(NTseq)):
    NTseqM.append(NTseq.iat[i,0])

# start used for iteration
start = 0 
current = []
unique_windows = []
value = 0
entropy_values = []
sum = 0
seqCount = len(NTseqM)


# loops through column
while start < (len(NTseqM[0])-9): 
    end = start + 10

    # appends each 100 amino acid window in all the rows
    for i in range(len(NTseqM)):
        current.append(NTseqM[i][start:end])    
        
    # finds number of occurrence of all unique 100 amino acid sequences in that specfic window for all the samples
    unique_windows = Counter(current).values()  

    #calculation for entropy values
    for x in unique_windows:
            #Value (A) = (# of this particular unique sequence A present / total # samples ) log base 2 (# of this particular unique sequence A present  / total # samples)
        value = (x/seqCount) * np.log2(x/seqCount)
        sum = sum + value
        print(sum)
    entropy_values.append(-1 * sum)

    #reset to start on the next window
    start = start + 1
    current = []
    unique_windows = []
    index = []                          # each position in one of the sequences
    sum = 0

# puts data into csv file
entropy_table = {"Entropy Values": entropy_values}
df = pd.DataFrame(entropy_table)
df.to_csv("/Users/numaanformoli/Documents/EntropyAnalysis/20221128_10_window_entropy_values.csv")

# puts starting point at midpoint 5
index.extend([-4, -3, -2, -1])

for i in range(len(NTseqM[0])):
    index.append(i)

# creates graph
plt.figure(figsize = (20, 7))
plt.bar(index[0:(len(NTseqM[0])-9)], entropy_values);
plt.plot()
plt.title("Shannon's Entropy")
plt.ylabel('Entropy Values')
plt.xlabel('Amino Acid Window 10')
plt.show() 


#calculates rank sum and compares N_terminal and C_terminal
def rank_sum_test_window10(entropy_values):
    values = entropy_values
    print(len(values[246:354]))
    n_terminal_values = values[45:167]
    print("N_terminal: " + str(len(n_terminal_values)))
    c_terminal_values = values[246:355]
    print("C_terminal: " + str(len(c_terminal_values)))
    results = stats.ranksums(n_terminal_values, c_terminal_values, 'greater')
    #Z-score indicates how much a given value differs from the standard deviation.
    print(results)

# n terminal (45 to 167)
# c terminal (246 to 355)

rank_sum_test_window10(entropy_values)