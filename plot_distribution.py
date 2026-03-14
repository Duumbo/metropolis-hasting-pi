import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import collections
import os

#cmap = sns.color_palette("mako_r", as_cmap=True)
#cmap = sns.color_palette("rocket_r", as_cmap=True)
#cmap = sns.color_palette("rocket", as_cmap=True)
#cmap = sns.light_palette("seagreen", as_cmap=True)
#cmap = sns.color_palette("light:salmon", as_cmap=True)
cmap = sns.color_palette("YlOrBr", as_cmap=True)

fps = os.listdir("data")
with open(f"data/{fps[0]}", 'r') as fp:
    for ln in fp:
        if ln.startswith("#"):
            metadata = ln
parameters = metadata.split(" ")
size = int(parameters[1])

freq = np.zeros((6, size, size))
for (ii, simulation) in enumerate(fps):

    with open(f"data/{simulation}", 'r') as fp:
        for ln in fp:
            if ln.startswith("#"):
                metadata = ln
    parameters = metadata.split(" ")
    size = int(parameters[1])

    data = np.genfromtxt(f"data/{simulation}")
    xi = range(0, size)
    yi = range(0, size)
    x, y = np.meshgrid(xi, yi)

    c = collections.Counter(data)
    c = sorted(c.items())
    for i in c:
        posx = int(i[0]) // size
        posy = int(i[0]) % size
        freq[ii, posx, posy] = int(i[1])

freq_max = np.max(freq)

for (ii, simulation) in enumerate(fps):
    with open(f"data/{simulation}", 'r') as fp:
        for ln in fp:
            if ln.startswith("#"):
                metadata = ln
    parameters = metadata.split(" ")
    size = int(parameters[1])
    fig, ax = plt.subplots()

    sns.heatmap(freq[ii, :, :], square=True, vmax=freq_max, cmap=cmap)
    ax.set_title(fr"$N_{{WARMUP}}={parameters[3]}, \pi\approx{np.round(float(parameters[5]),4)}$", fontsize=22)
    ax.set_xlabel(r"$x$", fontsize=20)
    ax.set_ylabel(r"$y$", fontsize=20)
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.tick_params(axis="both", which="minor", labelsize=14)
    xticks = ax.get_xticks()
    xticks_labels = ax.get_xticklabels()
    new_ticks = []
    new_labels = []
    for (ll, tick) in enumerate(xticks):
        if ll % 2 == 0:
            new_ticks.append(xticks[ll])
            new_labels.append(xticks_labels[ll])
        if ll == len(xticks) - 1 & ll % 2 != 0:
            new_ticks.append(xticks[ll])
            new_labels.append(xticks_labels[ll])
    ax.set_xticks(new_ticks)
    ax.set_xticklabels(new_labels)
    xticks = ax.get_yticks()
    xticks_labels = ax.get_yticklabels()
    new_ticks = []
    new_labels = []
    for (ll, tick) in enumerate(xticks):
        if ll % 2 == 0:
            new_ticks.append(xticks[ll])
            new_labels.append(xticks_labels[ll])
        if ll == len(xticks) - 1 & ll % 2 != 0:
            new_ticks.append(xticks[ll])
            new_labels.append(xticks_labels[ll])
    ax.set_yticks(new_ticks)
    ax.set_yticklabels(new_labels)

    fig.tight_layout()
    fig.savefig(f"{simulation}.png")

plt.show()
