import sys

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

colors = sns.color_palette("Set2")

fig, ax = plt.subplots(1, 2, figsize=(10, 5))

df = pd.read_csv(sys.argv[1], sep="\t") # 1kg.samples.tsv
#sample	pedigree_sex	gt_depth_mean	gt_depth_sd	depth_mean	depth_sd	ab_mean	ab_std	n_hom_ref	n_het	n_hom_alt	n_unknown	p_middling_ab	X_depth_mean	X_n	X_hom_ref	X_het	X_hom_alt	Y_depth_mean	Y_n

for i, ex_sex in enumerate(df.pedigree_sex.unique()):
    sub = df.loc[df.pedigree_sex == ex_sex, :]
    ax[1].scatter(sub.X_hom_alt, sub.Y_depth_mean, label=ex_sex, ec="none", c=colors[i])
    ax[0].scatter(sub.X_hom_alt, sub.X_het, label=ex_sex, ec="none", c=colors[i])

ax[0].set_xlabel("Number of 1/1 sites on X chromosome")
ax[0].set_ylabel("Number of 0/1 sites on X chromosome")
ax[0].legend(title="reported sex")

ax[1].set_xlabel("Number of 1/1 sites on X chromosome")
ax[1].set_ylabel("Mean depth on Y chromosome")

sns.despine()
plt.tight_layout()
plt.savefig("somalier-figure4.eps")
plt.show()
