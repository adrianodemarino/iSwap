import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import seaborn as sns
import glob


filename = glob.glob("*multiplexing")
for num,file in enumerate(filename):
	try:
		df = pd.read_csv(file,sep="\t",header=None,usecols=[1,2])
	except Exception as e:
		continue
	df[[2,3]] = df[2].str.split(".",expand=True)[[2,3]] 
	threshold = df[1].quantile(0.984)

	df = df.groupby([2,3]).sum().reset_index() 
	prob = list(dict(df[2]+"."+df[3]).values())
	LTR = df[df[1] >= threshold][2].value_counts().to_dict() 
	LC = df[df[1] >= threshold][3].value_counts().to_dict() 
	
	pivot = pd.pivot_table(df, columns=2,index=3,values=1)

	matrix = pivot.values
	rows = np.shape(matrix)[0] 
	columns = np.shape(matrix)[1]


	sumcol = matrix.sum(1)
	sumrow = matrix.sum(0)

	sumrows_noTP = pivot[pivot <= threshold].fillna(0).values.sum(1)
	sumcols_noTP = pivot[pivot <= threshold].fillna(0).values.sum(0)

	final = []

	for i in range(rows):
		for j in range(columns):
			if matrix[i][j] >= threshold:
				final.append([matrix[i][j],sumcols_noTP[j]+sumrows_noTP[i],pivot.columns[j],pivot.index[i]])


	final = pd.DataFrame(final)
	final[4] = final[2].map(LTR)
	final[5] = final[3].map(LC)

	final.to_csv("dataset.{x}.tsv".format(x=num),sep="\t",index=False,header=None)