import matplotlib.pyplot as plt
import numpy as np
import collections

SIZE = 50

with open('distribution', 'r') as fp:
    metadata = fp.readline()
parameters = metadata.split(" ")

data = np.genfromtxt("distribution")
indices = range(0, SIZE * SIZE)

c = collections.Counter(data)
c = sorted(c.items())
pos = [int(i[0]) for i in c]
freq = [int(i[1]) for i in c]

print(pos)
print(freq)
mean = np.mean(freq)
print(parameters)

fig, ax = plt.subplots()

ax.bar(pos, freq)
ax.hlines(mean, np.min(pos), np.max(pos), color = "black", linestyle = "dashed")
ax.set_title(fr"$N_{{WARMUP}}={parameters[3]}$", fontsize=22)
ax.set_xlabel(r"Position (indice)", fontsize=20)
ax.set_ylabel(r"Nombre de visite", fontsize=20)
ax.tick_params(axis="both", which="major", labelsize=14)
ax.tick_params(axis="both", which="minor", labelsize=14)
fig.tight_layout()
plt.show()
fig.savefig(f"{parameters[3]}.png")
