import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np

num_procs = np.array([256, 512, 1024, 2048])
compression = np.array([374, 174, 101, 64])
ulv = np.array([21, 13, 8.6, 6.0])
sol = np.array([5.3, 5.7, 8, 9.6])
matvec = np.array([13.7, 7.0, 4.3, 1.8])
ideal = np.array([np.max(compression) / 2**i for i in range(4)])
ideal_ulv = np.array([np.max(ulv) / 2**i for i in range(4)])

datasets = [compression, ulv, sol, matvec, ideal, ideal_ulv]
plotstyle = ['o-', 's-', '^-', 'v-', 'k--', 'k--']
plots=[]
plt.clf()
rcParams.update({'font.size': 10, 'legend.fontsize': 8})
for (i, d) in enumerate(datasets[:-2]):
    p, = plt.loglog(num_procs, d, plotstyle[i])
    plots.append(p)
for (i, d) in enumerate(datasets[-2:]):
    plt.loglog(num_procs, d, 'k--')
plt.legend(plots, ['compr.', 'ULV', 'solution', 'matvec'])
plt.xlim([250, 2200])
plt.ylim([1, 1000])
plt.ylabel('t/s')
plt.xlabel('# MPI procs.')
plt.gcf().subplots_adjust(left=0.15,bottom=0.15,right=0.98,top=0.97)
plt.gcf().set_size_inches(4, 3)
plt.xticks(num_procs, num_procs)
plt.savefig('strong_scaling.pdf')

