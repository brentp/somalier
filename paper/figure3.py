import sys

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.font_manager

import matplotlib as mpl
mpl.rcParams['font.family'] = 'arial'

colors = sns.color_palette("Set2")

df = pd.read_csv(sys.argv[1], sep="\t") # 1kg.pairs.tsv

s = 25
if df.shape[0] > 1000:
    s = 10

print(s)


for i, ex_rel in enumerate(df.expected_relatedness.unique()):
    sub = df.loc[df.expected_relatedness == ex_rel, :]
    plt.scatter(sub.ibs0, sub.ibs2, label="Unrelated" if ex_rel == -1 else
            "Identical", ec="none", c=colors[i], s=s)

sub = df.loc[df.ibs0 < 450]
plt.scatter(sub.ibs0, sub.ibs2, label="Unexpected relatedness",
        ec="k", c=colors[0], s=s+4)

plt.xlabel("Number of IBS0 sites")
plt.ylabel("Number of IBS2 sites")

sns.despine()
plt.legend(frameon=False)
#plt.savefig("somalier-figure3.eps")
plt.savefig("somalier-figure3.png")
#plt.show()
