.. _performance/Performance:

Performance
===========


| GTDB-Tk v2 runs 22% to 35% faster when processing 1000 genomes with 1 to 64 CPUs (**Fig.1**) and is >40% faster when processing 5,000 genomes using 32 CPUs (**Fig.2**).
| The tests below were run on a machine with 4 AMD EPYC 7402 24-Core Processor and 512 GB of RAM.

.. plot::

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

**Fig. 1**: Processing time (log scale) for 1,000 randomly selected GEM MAGs for increasing numbers of CPUs.

.. plot::

      color_setup = ['#5a6855','red','color']
      setup = color_setup
      pool_size_list=[10, 50, 100, 200, 500, 1000, 2000, 5000]
      split_list=[x/60 for x in [88, 150, 160, 169, 195, 235, 312, 558]]
      nosplit_list=[x/60 for x in [137, 151, 163, 180, 219, 280, 416, 934]]
      values = [10,500,1000,1500,2000,2500,3000,3500,4000,4500,5000]

      plt.scatter(pool_size_list, split_list, label="GTDB-Tk v2" , marker="s", s=30, color=setup[0])
      plt.plot(pool_size_list, split_list,linestyle='dashed', color=setup[0])
      plt.scatter(pool_size_list, nosplit_list, label="GTDB-Tk v1",color=setup[1])
      plt.plot(pool_size_list, nosplit_list,linestyle=':', color=setup[1])

      # naming the x axis

      plt.xticks(values)
      plt.yticks([1,3,5,7,9,11,13,15])

      # naming the axis
      plt.ylabel('Runtime (h)')
      plt.xlabel('Number of genomes')
      # giving a title to my graph
      plt.title('Fig.2: GTDB-Tk runtime with 32CPUs')

      # show a legend on the plot
      plt.legend(loc=2, prop={'size': 12},frameon=False)

      # function to show the plot
      plt.show()
**Fig. 2**: Processing time with 32 CPUs on increasing numbers of randomly selected GEM MAGs.
