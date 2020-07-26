#Plotting of the dense problem (squared)

from pyplot.matplotlib import *
import numpy as np

 # To display the densified CSR matrix

from numpy import *
import matplotlib.pyplot as plt

ia_f = loadtxt("data/matrices/ia.txt")
ja_f = loadtxt("data/matrices/ja.txt")
a_f = loadtxt("data/matrices/a.txt")

ia = ia_f[:]
ja = ja_f[:]
a = a_f[:]

dense = []

#creating a zero-filled dense matrix LoL
for i in range(len(ia) + 1):
    line = []
    for j in range(len(ia) + 1):
        line.append(0.0)
    dense.append(line)

for i in range(len(ia)):
    if(i < len(ia) - 1) :
        for k in range(ia[i].astype(int64), ia[i+1].astype(int64)):
            dense[i][ja[k].astype(int64)] = 10000.0

#displaying results
plt.matshow(dense)
plt.show()
