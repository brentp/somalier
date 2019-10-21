from pathlib import Path
import sys
import time


def read_somalier(path):
    """
    take a path to a single .somalier file and return a simple datastructure
    containing the sample information
    """
    data = Path(path).read_bytes()
    version = int.from_bytes(data[:1], byteorder="little")
    assert version == 2, ("bad version for:", path)
    data = data[1:]
    sample_L = int.from_bytes(data[:1], byteorder="little")
    data = data[1:]
    sample = data[:sample_L].decode()
    data = data[sample_L:]

    nsites = int.from_bytes(data[:2], byteorder="little")
    data = data[2:]
    nxsites = int.from_bytes(data[:2], byteorder="little")
    data = data[2:]
    nysites = int.from_bytes(data[:2], byteorder="little")
    data = data[2:]

    sites = np.frombuffer(data[:nsites * 3 * 4], dtype=np.uint32).reshape((nsites, 3))
    data = data[nsites * 3 * 4:]
    x_sites = np.frombuffer(data[:nxsites * 3 * 4], dtype=np.uint32).reshape((nxsites, 3))
    data = data[nxsites * 3 * 4:]
    y_sites = np.frombuffer(data[:nysites * 3 * 4], dtype=np.uint32).reshape((nysites, 3))

    return dict(sample=sample, sites=sites, x_sites=x_sites)

def to_gt(ab):
    result = np.zeros_like(ab) - 1
    result[(ab >= 0) & (ab < 0.015)] = 0
    result[(ab >= 0.2) & (ab <= 0.8)] = 1
    result[(ab > 0.985)] = 2
    return result

if __name__ == "__main__":

    from sklearn.decomposition import PCA
    import seaborn as sns
    from sklearn import svm
    from sklearn.metrics import confusion_matrix, accuracy_score
    from sklearn.pipeline import make_pipeline
    import numpy as np
    import pandas as pd
    from matplotlib import pyplot as plt
    import argparse

    sns.set_palette(sns.color_palette("Set1", 12))
    colors = sns.color_palette()

    p = argparse.ArgumentParser("predict ancestry given a labelled background set")
    p.add_argument("--labels", required=True, help="tsv file of sample => ancestry for background set first column must be sample-id")
    p.add_argument("--label-column", help="column name with population label from --labels", default="superpop")
    p.add_argument("--backgrounds", nargs="+", help="path to background *.somalier files matching those specified in labels")
    p.add_argument("--samples", nargs="+", help="path to sample *.somalier for ancestry prediction")
    p.add_argument("--plot", help="path to save figure. if not specified, the plot is `show`n")

    args = p.parse_args()

    label = args.label_column

    bg_samples = []
    bg_ABs = []
    test_samples = []
    test_ABs = []
    for f in args.backgrounds:
        s = read_somalier(f)
        depth = s["sites"].sum(axis=1)
        ab = s["sites"][:, 0] / np.maximum(depth, 1).astype(float)
        ab[depth < 5] = -1
        bg_ABs.append(ab)
        bg_samples.append(s["sample"])

    if args.samples is None: args.samples = []
    for f in args.samples:
        s = read_somalier(f)
        depth = s["sites"].sum(axis=1)
        ab = s["sites"][:, 0] / np.maximum(depth, 1).astype(float)
        ab[depth < 5] = -1
        test_ABs.append(ab)
        test_samples.append(s["sample"])

    bg_sample_df = pd.read_csv(args.labels, sep="\t", escapechar='#', index_col=0)

    n_components = 5
    clf = make_pipeline(PCA(n_components=n_components, whiten=True, copy=True, svd_solver="randomized"),
                svm.SVC(C=3, probability=True, gamma="auto"))


    # convert labels to integers
    bg_samples = np.array(bg_samples)
    test_samples = np.array(test_samples)
    bg_ABs = np.array(bg_ABs, dtype=float)
    test_ABs = np.array(test_ABs, dtype=float)

    unk = (bg_ABs == -1).sum(axis=0)
    rm = unk / float(len(bg_samples)) > 0.3
    if len(test_ABs) > 0:
        unk = (test_ABs == -1).sum(axis=0)
        rm |=  (unk / float(len(test_samples)) > 0.4)

    bg_ABs = bg_ABs[:, ~rm]

    np.save("thousandG.npy", bg_ABs)
    if len(test_ABs) > 0:
        test_ABs = test_ABs[:, ~rm]

    bg_sample_df = bg_sample_df.loc[bg_samples, :]
    targetL = list(bg_sample_df[label].unique())


    target = np.array([targetL.index(p) for p in bg_sample_df[label]])

    clf.fit(bg_ABs, target)

    bg_reduced = clf.named_steps["pca"].transform(bg_ABs)

    # generate PC labels for n components
    labels_pc = [("PC%d" % i) for i in range(1,(n_components+1))]

    # export csv
    bg_reduced_df = pd.DataFrame(data=bg_reduced, index=bg_sample_df[label].values, columns=labels_pc)
    bg_reduced_df.index.name = "ancestry"
    bg_reduced_df.to_csv(path_or_buf="somalier.background_pca.csv")
    
    # export json
    bg_reduced_df.reset_index(level=0, inplace=True)
    bg_reduced_df.to_json(path_or_buf="somalier.background_pca.json",)

    if len(test_ABs) > 0:
        test_reduced = clf.named_steps["pca"].transform(test_ABs)
        test_pred = clf.predict(test_ABs)
        test_prob = clf.predict_proba(test_ABs)

        # export prediction and PC's as csv
        sample_df = pd.DataFrame(data=np.hstack((test_prob.round(2), test_reduced)), index=test_samples, columns=[*targetL, *labels_pc])
        sample_df.index.name = "#sample"
        sample_df.to_csv(path_or_buf="somalier.sample_prediction_pca.csv")


    fig, axes = plt.subplots(1) #, len(targetL) + 1, figsize=(22, 12))
    axes = (axes,)
    for i, l in enumerate(sorted(set(bg_sample_df[label]))):
        sel = bg_sample_df[label] == l
        ibg_sample_df = bg_sample_df.loc[sel, :]
        axes[0].scatter(bg_reduced[sel, 0], bg_reduced[sel, 1], label=l, s=8,
                alpha=0.15, c=[colors[i % len(colors)]], ec='none')

        if len(test_ABs) > 0:
            test_sel = (test_pred == targetL.index(l))
            axes[0].scatter(test_reduced[test_sel, 0], test_reduced[test_sel, 1], s=4,
                alpha=0.9, c=[colors[i % len(colors)]])

    #for ax in axes: ax.legend()
    plt.legend()
    plt.tight_layout()
    if args.plot in ["", None]:
        plt.show()
    else:
        plt.savefig(args.plot)
