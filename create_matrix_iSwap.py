import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import seaborn as sns
import glob


filename = glob.glob("*multiplexing")
for file in filename:
	df = pd.read_csv(filei,sep="\t",header=None,usecols=[1,2])
	df[[2,3]] = df[2].str.split(".",expand=True)[[2,3]] 

	df = df.groupby([2,3]).sum().reset_index() 

	pivot = pd.pivot_table(df, columns=2,index=3,values=1)

	matrix = pivot.values
	rows = np.shape(matrix)[0] 
	columns = np.shape(matrix)[1]

	threshold = df[1].quantile(0.984)

	sumcol = matrix.sum(1)
	sumrow = matrix.sum(0)

	sumcol_noTP = pivot[pivot <= threshold].fillna(0).values.sum(1)
	sumrow_noTP = pivot[pivot <= threshold].fillna(0).values.sum(0)

	final = []

	for i in range(rows):
		for j in range(columns):
			if matrix[i][j] >= threshold:
				final.append([matrix[i][j],sumcol_noTP[j]+sumrow_noTP[i]])


	pd.DataFrame(final)