
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