import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

color_setup = ['#5a6855','red','color']
setup = color_setup
cpus_list = [1, 8, 16, 32, 64]
split_list = [x/60 for x in [2571, 457, 309, 230, 216]]
nosplit_list= [x/60 for x in [3980, 638, 398, 295, 279]]
values = ['1', '8', '16', '32','64']

fig, ax = plt.subplots()

plt.ylim(3,70)
plt.xticks(cpus_list, values)

ax.scatter(cpus_list, split_list, label="GTDB-Tk v2", marker="s", s=30, color=setup[0])
ax.plot(cpus_list, split_list,linestyle='dashed', color=setup[0])
ax.scatter(cpus_list, nosplit_list, label="GTDB-Tk v1",color=setup[1])
ax.plot(cpus_list, nosplit_list,linestyle=':', color=setup[1])

ax.set_yscale('log')
ax.set_yticks([3,5,7,10,20,30,50,70])
ax.yaxis.set_major_formatter(ScalarFormatter())

plt.ylabel('Runtime (h)')
plt.xlabel('Number of CPUs')
plt.title('Fig.1: GTDB-Tk runtime for 1000 genomes (log scale)')

# show a legend on the plot
plt.legend(loc=1, prop={'size': 12},frameon=False)
plt.show()