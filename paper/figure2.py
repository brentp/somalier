import sys

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

colors = sns.color_palette("Set2")

fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)

dfa = pd.read_csv(sys.argv[1], sep="\t")
dfb = pd.read_csv(sys.argv[2], sep="\t")

s = 30

titles = ["Original Samples", "After Correction"]

for i, df in enumerate((dfa, dfb)):
    ax = axes[i]
    ax.set_title(titles[i])

    for j, ex_rel in enumerate(df.expected_relatedness.unique()):
        sub = df.loc[df.expected_relatedness == ex_rel, :]
        ax.scatter(sub.ibs0, sub.ibs2, label="unrelated" if ex_rel == -1 else "identical",
                ec="none", c=colors[j], s=s)

    ax.set_xlabel("Number of IBS0 sites", fontsize=13)
    if i == 0:
        ax.set_ylabel("Number of IBS2 sites", fontsize=13)
    print(i)

    if i == 0:
        ax.legend(title="expected relatedness")

    sns.despine()
plt.savefig("somalier-figure2.eps")
plt.show()
