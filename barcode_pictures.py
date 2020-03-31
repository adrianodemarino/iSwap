#!/usr/bin/python3
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as colors
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument('-f', dest="file", help="", action="store", required=True)

args = parser.parse_args()
file = args.file

import re
def humanSorted(l):
    def tryint(s):
        try:
            return int(s)
        except:
            return s
    def alphanum_key(s):
        return [ tryint(c) for c in re.split('([0-9]+)', s) ]
    def sort_nicely(l):
        return sorted(l, key=alphanum_key)
    return sort_nicely(l)

df = pd.read_csv(file,sep='\t').drop_duplicates()

#First Picture#################################################################
a = pd.pivot_table(df,values="used",columns="bc2",index="bc1").astype(int)
#a = a.reindex(index=humanSorted(df.bc1.unique()), columns=humanSorted(df.bc2.unique()))
f, ax = plt.subplots(figsize=(8, 8))

colors = ["gray", "cyan", "orange"]
cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

# Draw the heatmap with the mask and correct aspect ratio
_ = sns.heatmap(a, cmap=cmap,square=True,  linewidths=.5, cbar_kws={"shrink": .5})

# Set the colorbar labels
colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0.35,0.95,1.65])
colorbar.set_ticklabels(['Unused combination', 'Expected combination',"Swapped fraction"])
#ax.invert_xaxis()

ax.set_xlabel("Barcode 2",fontsize=20)
ax.set_ylabel("Barcode 1",fontsize=20)
ax.tick_params(labelsize=11)
_.set_yticklabels(_.get_yticklabels(), rotation = 0, fontsize = 10)
# fix for mpl bug that cuts off top/bottom of seaborn viz
b, t = plt.ylim() # discover the values for bottom and top
b += 0.5 # Add 0.5 to the bottom
t -= 0.5 # Subtract 0.5 from the top
plt.ylim(b, t) # update the ylim(bottom, top) values


ax.get_figure().savefig(file.split('.')[0]+'_barcodeStrategy.png',format='png', dpi=500)
##########################################################################################


#Second Picture#################################################################
b = pd.pivot_table(df,values="reads",columns="bc2",index="bc1").fillna(0)
b = b.reindex(index=humanSorted(df.bc1.unique()), columns=humanSorted(df.bc2.unique()))
fig, bx = plt.subplots(figsize=(10, 10))

log_norm = LogNorm()
# Draw the heatmap with the mask and correct aspect ratio
# color bar keyword arguments
cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(df.reads.min())), 1+math.ceil(math.log10(df.reads.max())))]
cbar_kws = {"shrink":.5,'label':'Log10 read counts',"ticks": cbar_ticks}

pcm = sns.heatmap(b, cmap="viridis",square=True,  linewidths=.5, cbar_kws=cbar_kws,norm=log_norm)

bx.invert_xaxis()
bx.set_xlabel("Barcode 2",fontsize=20)
bx.set_ylabel("Barcode 1",fontsize=20)
bx.tick_params(labelsize=11)

bx.get_figure().savefig(file.split('.')[0]+'_LogReadsCount.png',format='png', dpi=500)

