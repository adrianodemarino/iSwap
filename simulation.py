
import pandas as pd 
import numpy as np 
import random as rd
import matplotlib.pyplot as plt
import scipy.special as sps

hg = pd.read_csv("hg38.size",sep="\t")
sample = pd.read_csv("sample_name.tsv",sep="\t")

#Sequence of Interest Sofi Creator
def Sofi():
	selection = hg.sample(1)
	chrom = selection.CHR.values
	position = rd.randrange(selection.Size.values)
	Construct = chrom+"_"+str(position)
	return Construct[0]

def GammaDistribution():	
	shape, scale = 2., 2. 
	s = np.random.gamma(shape, scale, 1000)
	count, bins, ignored = plt.hist(s, 50, density=True)
	y = bins**(shape-1)*(np.exp(-bins/scale)/(sps.gamma(shape)*scale**shape))
	plt.plot(bins, y, linewidth=2, color='r')
	plt.show()

def ShearsiteDistribution(len=1000000):
	shape, scale = 2., 2. 
	s = np.random.gamma(shape, scale, len)
	s = (s*100).astype(int)
	s = s[s<=1000]
	return s

order_i = {1:10,2:100,3:1000,4:10000,5:100000,6:1000000,7:10000000,8:100000000,9:1000000000,10:10000000000}
order_f = {1:99,2:999,3:9999,4:99999,5:999999,6:9999999,7:99999999,8:999999999,9:999999999,10:9999999999} 

def getSimuNumber(number): 
    number = magnitude(number) 
    return np.random.randint(1*order[number-2],1*orderf[number-2]) 
