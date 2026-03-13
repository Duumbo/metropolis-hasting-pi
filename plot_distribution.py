import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import collections

with open('distribution', 'r') as fp:
    metadata = fp.readline()
parameters = metadata.split(" ")
size = int(parameters[1])

data = np.genfromtxt("distribution")
xi = range(0, size)
yi = range(0, size)
x, y = np.meshgrid(xi, yi)

freq = np.zeros((size, size))
c = collections.Counter(data)
c = sorted(c.items())
for i in c:
    posx = int(i[0]) // size
    posy = int(i[0]) % size
    freq[posx, posy] = int(i[1])


print(freq)
mean = np.mean(freq)
print(parameters)

fig, ax = plt.subplots()

sns.heatmap(freq)
#ax.hlines(mean, np.min(pos), np.max(pos), color = "black", linestyle = "dashed")
ax.set_title(fr"$N_{{WARMUP}}={parameters[3]}$", fontsize=22)
ax.set_xlabel(r"Position (indice)", fontsize=20)
ax.set_ylabel(r"Nombre de visite", fontsize=20)
ax.tick_params(axis="both", which="major", labelsize=14)
ax.tick_params(axis="both", which="minor", labelsize=14)
fig.tight_layout()
plt.show()
fig.savefig(f"{parameters[3]}.png")
