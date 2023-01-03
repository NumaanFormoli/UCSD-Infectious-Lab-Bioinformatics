
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
m = ['G', "A", "V", "C", "P", "L", "I", "M", "W", "F", "K", "R", "H", "S", "T", "Y", "N", "Q", "D", "E"]
#opens the excel sheet as a dataframe so we can work with it in Pandas
NTonly = pd.read_excel("amino_sequences.xlsx")
NTseq = NTonly.iloc[:, [0]]
NTseq
#normalize lengths by adding 'N' to sequences shorter than the longest sequence
lengths = []
for i in range(len(NTseq)):
    lengths.append(len(NTseq.iat[i,0]))
#print(max(lengths))

for i in range(len(NTseq)):
    if len(NTseq.iat[i,0]) < max(lengths):
        toadd = max(lengths) - len(NTseq.iat[i,0])
        for k in range(toadd):
            NTseq.iat[i,0] += "Z"
            k +=1
#loops through each and makes it into a matrix so we can check by index if nucleotides differ
NTseqM = []  # list of sequences
for i in range(len(NTseq)):
    NTseqM.append(NTseq.iat[i,0])
    
print(len(NTseqM))

p1 = []

for j in range(len(NTseqM[1])):    # loop through first sequence length in p1 list
    cc = []
    st = []
    for k in range(len(NTseqM)):   # loop through length of sequences (348) add each amino acid in the different rows same coulm
        st.append(NTseqM[k][j]) # adds all amino acids in the first position to one list
    p1.append(st)    # adds list to another list of each position
    print(p1[0])
            
#len(NTseqM[0])
#len(NTseqM)
#cc
print(len(p1[0]))
p2 = []
for j in p1:        # each list of one position
    counter = []
    for i in range(len(m)):     # length of amino acids 
        counter.append( (1+ j.count(m[i]))/(20^420)) #add to counter list number of certain amino acid in first list
    p2.append(counter)  # add to another list for each position for all amino acids
u2 = np.log2(p2)    # calculation

arr2 = []
#print(sum(p[2]))
for i in range(len(p2)):
    Sum = 0
    for j in range(len(p2[i])):
        Sum +=p2[i][j] * u2[i][j]
    arr2.append(-1*Sum)
plt.figure(figsize = (20, 7))
plt.bar(range(len(NTseqM[0])),arr2)
plt.xlabel("Amino acid positions")
plt.show()

###########


from collections import Counter
WuKabat = [] #list of values
seqCount = len(NTseqM) #length 

def mostCommon(lst):
    return Counter(lst).most_common(1)[0][0]

start = 0
current = []
unique = []
count = 0

while start < (len(NTseqM[0])-99): #anywhere you see 49 and 50 can be changed to adjust the window
    end = start + 100
    for i in range(len(NTseqM)):
        current.append(NTseqM[i][start:end])    # adds 100 windows to list
    unique = Counter(current).keys()        # finds unique keys and puts them in list unique
    for x in current:
        if (x == mostCommon(current)):
            count = count + 1
    WuKabat.append(seqCount*len(unique)/count)
    start = start + 1
    current = []
    unique = []
    count = 0
#print(WuKabat)
from matplotlib import pyplot as plt
xs = [1, 200]
plt.figure(figsize = (20, 7))
plt.bar(index[0:(len(NTseqM[0])-99)],WuKabat) 
plt.title("Wu-Kabat Variability Coefficient")
plt.ylabel('Variability ')
plt.xlabel('Amino Acids')
#plt.vlines(x = [46, 127], ymin = 0, ymax = max(xs),
           #colors = 'purple')
#plt.vlines(x = [247, 315], ymin = 0, ymax = max(xs),
           #colors = 'red')
plt.show() 