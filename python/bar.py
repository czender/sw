# Purpose: Plot effects of lossy codecs on final compression ratio
# History: Created for E3SM Newsletter 20220428

import numpy as np
import matplotlib.pyplot as plt

# Bar width
barWidth=0.25
fig=plt.subplots(figsize=(12,8))

# Bar bottom
bar_btm=0.0

# Bar height
dat_dfl=[1.78,1.78,1.78]
dat_bgr=[2.55,3.02,3.79]
dat_gbr=[2.96,3.69,5.73]
dat_btr=[2.99,3.64,5.29]

# Set bar positions on X-axis
br0=0.25
br1=1+np.arange(len(dat_bgr))
br2=[x+barWidth for x in br1]
br3=[x+barWidth for x in br2]

# Make plot
plt.bar(br0,dat_dfl,color='c',width=0.75,label='DEFLATE (level = 1)',bottom=bar_btm)
plt.bar(br1,dat_bgr,color ='r',width=barWidth,edgecolor='grey',label='+ BitGroom (NSD)',bottom=bar_btm)
plt.bar(br2,dat_gbr,color ='g',width=barWidth,edgecolor='grey',label='+ GranularBR (NSD)',bottom=bar_btm)
plt.bar(br3,dat_btr,color ='b',width=barWidth,edgecolor='grey',label='+ BitRound (NSB)',bottom=bar_btm)
plt.bar(br1+barWidth,dat_dfl,color='c',width=0.75,edgecolor='grey',bottom=bar_btm)

# Decorate
plt.title('Lossless & Lossy Contributions to Compression Ratio',fontweight='bold',fontsize=20)
plt.xlabel('Number of Significant Digits or Bits Retained',fontweight='bold',fontsize=20)
plt.ylabel('Final Compression Ratio',fontweight='bold',fontsize=20)
plt.xticks([r+barWidth for r in range(1+len(dat_bgr))],['NSD=7, NSB=23','NSD=4, NSB=12','NSD=3, NSB=9','NSD=2, NSB=6'],fontsize=15)
plt.yticks(fontsize=15)

plt.legend(fontsize=20)
plt.axis([-.5, 4, 1, 6]) # works
plt.show()
