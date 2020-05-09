#!/usr/bin/python3
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as colors
import seaborn as sns
from matplotlib.colors import LogNorm
import math


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
'''
mapper = {'N701':"A",
 'N702':"B",
 'N703':"C",
 'N704':"D",
 'N705':"E",
 'N706':"F",
 'N707':"G",
 'N710':"H",
 'N711':"I",
 'N712':"K",
 'N714':"L",
 'N715':"M",
 'N716':"N",
 'N718':"O",
 'N719':"P",
 'N720':"Q",
 'N721':"R",
 'N722':"S",
 'N723':"T",
 'N724':"U",
 'N726':"V",
 'N727':"X",
 'N728':"Y",
 'N729':"Z"}

mapper2 = {'S502':1,
 'S503':2,
 'S505':3,
 'S506':4,
 'S507':5,
 'S508':6,
 'S510':7,
 'S511':8,
 'S513':9,
 'S515':10,
 'S516':11,
 'S517':12,
 'S518':13,
 'S520':14,
 'S521':15,
 'S522':16}
df.bc1 = df.bc1.map(mapper)     
df.bc2 = df.bc2.map(mapper2)
''' 
#First Picture#################################################################
a = pd.pivot_table(df,values="expected",columns="bc2",index="bc1").astype(int)
#a = a.reindex(index=humanSorted(df.bc1.unique()), columns=humanSorted(df.bc2.unique()))
f, ax = plt.subplots(figsize=(12, 10))

colors = ["gray", "cyan"]
cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

# Draw the heatmap with the mask and correct aspect ratio
_ = sns.heatmap(a, cmap=cmap,square=True,  linewidths=.5, cbar_kws={"shrink": .5})

# Set the colorbar labels
colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0.25,0.75])
colorbar.set_ticklabels(['Unused combination', 'Expected combination'])
#ax.invert_xaxis()

ax.set_xlabel("Barcode 2",fontsize=20)
ax.set_ylabel("Barcode 1",fontsize=20)
ax.tick_params(labelsize=11)
_.set_yticklabels(_.get_yticklabels(), rotation = 0, fontsize = 13)
_.set_xticklabels(_.get_xticklabels(), rotation = 0, fontsize = 13)
# fix for mpl bug that cuts off top/bottom of seaborn viz
b, t = plt.ylim() # discover the values for bottom and top
b += 0.5 # Add 0.5 to the bottom
t -= 0.5 # Subtract 0.5 from the top
plt.ylim(b, t) # update the ylim(bottom, top) values


ax.get_figure().savefig(file.split('.')[0]+'_barcodeStrategy.png',format='png', dpi=500)
##########################################################################################


#Second Picture#################################################################
b = pd.pivot_table(df,values="reads",columns="bc2",index="bc1")
#b = b.reindex(index=humanSorted(df.bc1.unique()), columns=humanSorted(df.bc2.unique()))
fig, bx = plt.subplots(figsize=(12, 10))

log_norm = LogNorm()
# Draw the heatmap with the mask and correct aspect ratio
# color bar keyword arguments
cbar_ticks = [math.pow(10, i) for i in range(math.floor(math.log10(df.reads.min())), 1+math.ceil(math.log10(df.reads.max())))]
cbar_kws = {"shrink":.5,'label':'Log10 read counts',"ticks": cbar_ticks}

pcm = sns.heatmap(b, cmap="viridis",square=True,  linewidths=.5, cbar_kws=cbar_kws,norm=log_norm)

bx.set_xlabel("Barcode 2",fontsize=20)
bx.set_ylabel("Barcode 1",fontsize=20)
bx.tick_params(labelsize=11)
pcm.set_yticklabels(pcm.get_yticklabels(), rotation = 0, fontsize = 13)
pcm.set_xticklabels(pcm.get_xticklabels(), rotation = 0, fontsize = 13)
# fix for mpl bug that cuts off top/bottom of seaborn viz
b, t = plt.ylim() # discover the values for bottom and top
b += 0.5 # Add 0.5 to the bottom
t -= 0.5 # Subtract 0.5 from the top
plt.ylim(b, t) # update the ylim(bottom, top) values


bx.get_figure().savefig(file.split('.')[0]+'_LogReadsCount.png',format='png', dpi=500)

#terza Picture#################################################################
medfrac_mapped = df.mapped[~df.expected].median()/df.mapped[df.expected].median()
totfrac_mapped = df.mapped[~df.expected].sum()/df.mapped[df.expected].sum()
medfrac_all = df.reads[~df.expected].median()/df.reads[df.expected].median()
totfrac_all = df.reads[~df.expected].sum()/df.reads[df.expected].sum()

'''
The distribution of total number of mapped reads (i.e., library sizes) for all 
combinations are shown in Figure 4 and Figure 5. The impossible combinations have a 
median mapped-read library size that is 1.5% of the median size of the expected combinations. 
The total number of mapped reads assigned to impossible combinations is 1.1% of that assigned
to expected combinations. Note that the empty well is still considered as an expected barcode 
combination due to the presence of ERCC spike-in transcripts.
'''


nex = df[df.Expected == True][["Reads","Mapped","Expected"]]
nex.Expected = "Expected barcode combinations"
Mapped = nex.Mapped.copy()
nex.Mapped = "Total Reads"
nex = nex.rename({"Mapped":"Read Type"},axis=1)
Mapped = Mapped.to_frame()
Mapped = Mapped.rename({"Mapped":"Reads"},axis=1)
Mapped["Read Type"] = "Mapped Reads"
Mapped["Expected"] = "Expected barcode combinations"
new = pd.concat([nex,Mapped])

nex = df[df.Expected == False][["Reads","Mapped","Expected"]]
nex.Expected = "Unused Combinations"
Mapped = nex.Mapped.copy()
nex.Mapped = "Total Reads"
nex = nex.rename({"Mapped":"Read Type"},axis=1)
Mapped = Mapped.to_frame()
Mapped = Mapped.rename({"Mapped":"Reads"},axis=1)
Mapped["Read Type"] = "Mapped Reads"
Mapped["Expected"] = "Unused Combinations"
new2 = pd.concat([nex,Mapped])

new = pd.concat([new,new2])

fig, cx = plt.subplots(figsize=(12, 10))
my_pal = {"Mapped Reads": "blue", "Total Reads": "cyan"}
x = sns.boxplot(x="Expected",y="Reads",hue="Read Type",data=new, palette=my_pal)
plt.yscale("log") 
cx.set_xlabel("Combinations",fontsize=20)
cx.set_ylabel("Reads",fontsize=20)
cx.tick_params(labelsize=11)
x.set_xticklabels(x.get_xticklabels(), rotation = 0, fontsize = 13)

cx.get_figure().savefig('3_log_expe_unused.png',format='png', dpi=500)

#Quarta Picture#################################################################
#################################################################


fig, cx = plt.subplots(figsize=(12, 10))
#df.rename({"bc1":"Barcode 1","bc2":"Barcode 2","reads":"Reads","expected":"Expected"},axis=1,inplace=True) 
plt.yscale("log")
ax = sns.swarmplot(x="Expected",y="Reads",hue="Barcode 1",data=df,size=5,dodge=False,order=[True,False],palette="Greens")
#ax = sns.swarmplot(x="Expected",y="Reads",hue="Barcode 2",data=df,size=5,dodge=False,order=[True,False],palette="Reds")
cx.set_xlabel("Combinations",fontsize=20)
cx.set_ylabel("Reads",fontsize=20)
ax.tick_params(labelsize=14)
cx.get_figure().savefig('4_barcode1_pos.png',format='png', dpi=500)


#Quarta Picture#################################################################
#################################################################


pvt_unused = pd.pivot_table(index="Barcode 2",columns="Barcode 1",values="Reads",data=df[df.Expected == False]) 
pvt_used = pd.pivot_table(index="Barcode 2",columns="Barcode 1",values="Reads",data=df[df.Expected == True]) 

x_swapped_reads = {}

for i in pvt_unused.columns:
    for j in pvt_unused.index:
        if np.isnan(pvt_unused[i][j]):
            continue
        else:
            avai_fraction = pvt_used[i].sum()+pvt_used.loc[j].sum()
            x_swapped_reads[i+str(j)] = (pvt_unused[i][j],avai_fraction)

mydata = pd.DataFrame.from_dict(x_swapped_reads).T 
mydata.rename({1:"Avaiable Swapping Reads", 0:"Observed Swapped Reads"},axis=1,inplace=True)

ax = sns.regplot(y=mydata["Observed Swapped Reads"], x=mydata["Avaiable Swapping Reads"], marker="+")

ax.get_figure().savefig('5_swapped_fraction.png',format='png', dpi=500)
