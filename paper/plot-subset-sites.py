import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import sys
sns.set_palette("Set2")

df = pd.read_csv(sys.argv[1], sep="\t")

fig, axes = plt.subplots(1, 1, figsize=(8, 4))

try:
    axes[0]
except (IndexError, TypeError):
    axes = (axes,)


sns.barplot(x="n", y="fp", hue="strict", data=df, ax=axes[0])

axes[0].set_xlabel("Number of sites")
axes[0].set_ylabel("False-positive rate")

#sns.barplot(x="n", y="tp", hue="strict", data=df, ax=axes[1])
#axes[1].set_xlabel("Number of sites")
#axes[1].set_ylabel("True-positive rate")

plt.savefig("subset-sites.png")
plt.show()

"""
n	tp	fp	fn	strict
10	0.816905	0.183080	0.000014	false
20	0.859777	0.140218	0.000005	false
40	0.925964	0.074012	0.000024	false
100	0.985616	0.014384	0.000000	false
200	0.997638	0.002362	0.000000	false
400	0.999724	0.000276	0.000000	false
1000	0.999986	0.000014	0.000000	false
2000	1.000000	0.000000	0.000000	false
4000	1.000000	0.000000	0.000000	false


"""
