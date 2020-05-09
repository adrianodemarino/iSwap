import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as colors
import seaborn as sns


LC = {'LC10': np.nan,'LC12': np.nan,'LC14': np.nan,'LC16': np.nan,'LC18': np.nan,'LC2': np.nan,'LC20': np.nan,'LC22': np.nan,'LC24': np.nan,'LC26': np.nan,'LC28': np.nan,'LC30': np.nan,'LC32': np.nan,'LC34': np.nan,'LC36': np.nan,'LC38': np.nan,'LC4': np.nan,'LC40': np.nan,'LC42': np.nan,'LC44': np.nan,'LC46': np.nan,'LC48': np.nan,'LC50': np.nan,'LC52': np.nan,}
LTR = ['LTR1','LTR3','LTR5','LTR7','LTR9','LTR11','LTR13','LTR15','LTR17','LTR19','LTR21','LTR23']
df = pd.DataFrame(LC, index=LTR).fillna(0).astype(int)
df = df.reindex(index=humanSorted(df.index), columns=humanSorted(df.columns))


f, ax = plt.subplots(figsize=(10, 10))
colors = ["gray", "cyan"]
cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

# Draw the heatmap with the mask and correct aspect ratio
_ = sns.heatmap(df, cmap=cmap,square=True,  linewidths=.3, cbar_kws={"shrink": .5})

# Set the colorbar labels
colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0.25,0.75])
colorbar.set_ticklabels(['Unused combination', 'Expected combination'])


ax.set_xlabel("Barcode 2",fontsize=20)
ax.set_ylabel("Barcode 1",fontsize=20)
ax.tick_params(labelsize=11)

ax.get_figure().savefig('iSwap_logo.png',format='png', dpi=500)


df['LC44']['LTR7'] =1
df['LC44']['LTR9'] =1
df['LC44']['LTR11'] =1
df['LC44']['LTR13'] =1
df['LC44']['LTR15'] =1
df['LC44']['LTR17'] =1
df['LC46']['LTR7'] =0
df['LC46']['LTR9'] =1
df['LC48']['LTR9'] =1
df['LC48']['LTR11'] =0
df['LC50']['LTR9'] =1
df['LC50']['LTR11'] =1
df['LC50']['LTR13'] =1
df['LC48']['LTR13'] =1
df['LC46']['LTR13'] =1



#a
df['LC32']['LTR5'] =1
df['LC34']['LTR5'] =1
df['LC36']['LTR5'] =1
df['LC38']['LTR5'] =1
df['LC38']['LTR7'] =1
df['LC38']['LTR9'] =1
df['LC38']['LTR11'] =1
df['LC38']['LTR13'] =1
df['LC40']['LTR13'] =1
df['LC36']['LTR13'] =1
df['LC34']['LTR13'] =1
df['LC32']['LTR13'] =1
df['LC32']['LTR11'] =1
df['LC32']['LTR9'] =1
df['LC34']['LTR9'] =1
df['LC36']['LTR9'] =1



#w
df['LC20']['LTR13'] =0
df['LC20']['LTR11'] =1
df['LC20']['LTR9'] =1
df['LC20']['LTR7'] =1
df['LC20']['LTR5'] =1
df['LC20']['LTR3'] =0
df['LC22']['LTR11'] =0
df['LC22']['LTR9'] =0
df['LC22']['LTR13'] =1
df['LC24']['LTR13'] =0
df['LC24']['LTR11'] =1
df['LC24']['LTR9'] =1
df['LC24']['LTR7'] =0
df['LC24']['LTR5'] =0
df['LC24']['LTR3'] =0
df['LC26']['LTR13'] =1
df['LC28']['LTR11'] =1
df['LC28']['LTR9'] =1
df['LC28']['LTR7'] =1
df['LC28']['LTR5'] =1
df['LC28']['LTR3'] =0


#### A
df['LC20']['LTR13'] =0
df['LC20']['LTR11'] =0
df['LC20']['LTR9'] =0
df['LC20']['LTR7'] =0
df['LC20']['LTR5'] =0
df['LC22']['LTR3'] =0
df['LC24']['LTR5'] =0
df['LC24']['LTR7'] =0
df['LC24']['LTR9'] =0
df['LC24']['LTR11'] =0
df['LC24']['LTR13'] =0
df['LC22']['LTR9'] =0