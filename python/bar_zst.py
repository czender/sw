# Purpose: Plot effects of lossy codecs on final compression ratio
# History: Created for E3SM Newsletter 20220428
# Usage:
# cd ~/anaconda/bin;conda activate
# cd ~/sw/python;python bar_zst.py
# conda deactivate

import numpy as np
import matplotlib.pyplot as plt

# Bar width
barWidth=0.25
fgr=plt.subplots(figsize=(12,8))

# Bar bottom
bar_btm=0.0

# Bar height
dat_zst=[1.61,1.61,1.61]
dat_shf=[1.83,1.83,1.83]
dat_bgr=[2.68,3.21,4.22]
dat_gbr=[3.18,4.01,6.51]
dat_btr=[3.23,3.94,5.87]

# Set bar positions on X-axis
br0=0.25
br1=1+np.arange(len(dat_bgr))
br2=[x+barWidth for x in br1]
br3=[x+barWidth for x in br2]

# Make plot
plt.bar(br0,dat_zst,color='c',width=0.75,label='Zstandard (level = 3)',bottom=bar_btm)
plt.bar(br0,dat_shf,color='y',width=0.75,label='+ Shuffle',bottom=bar_btm)
plt.bar(br0,dat_zst,color='c',width=0.75,bottom=bar_btm)
plt.bar(br1,dat_bgr,color ='r',width=barWidth,edgecolor='grey',label='+ BitGroom (NSD)',bottom=bar_btm)
plt.bar(br2,dat_gbr,color ='g',width=barWidth,edgecolor='grey',label='+ GranularBR (NSD)',bottom=bar_btm)
plt.bar(br3,dat_btr,color ='b',width=barWidth,edgecolor='grey',label='+ BitRound (NSB)',bottom=bar_btm)
plt.bar(br1+barWidth,dat_shf,color='y',width=0.75,edgecolor='grey',bottom=bar_btm)
plt.bar(br1+barWidth,dat_zst,color='c',width=0.75,edgecolor='grey',bottom=bar_btm)

# Decorate
plt.title('Lossless & Lossy Contributions to Compression Ratio',fontweight='bold',fontsize=20)
plt.xlabel('Number of Significant Digits or Bits Retained',fontweight='bold',fontsize=20)
plt.ylabel('Final Compression Ratio',fontweight='bold',fontsize=20)
plt.xticks([r+barWidth for r in range(1+len(dat_bgr))],['NSD=7, NSB=23','NSD=4, NSB=12','NSD=3, NSB=9','NSD=2, NSB=6'],fontsize=15)
plt.yticks(fontsize=15)

plt.legend(fontsize=20)
plt.axis([-.5, 4, 1, 7]) # works
plt.show()
