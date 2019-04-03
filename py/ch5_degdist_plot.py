from graph_tools import *
from nem_graph import *
import matplotlib.pyplot as plt
import numpy as np

#Set parameters
p = 100
q = 1000
mu = 0.20

#Number of graphs to be simulated
n = 1000

#Get theoretical distribution
th_deg_dist = nem_th_deg_dist(p, q, mu)
while (th_deg_dist[-1] == 0):
    th_deg_dist.pop()

#Get the degree distribution of n Nem. graphs
graphs = []
for i in range(n):
    graphs.append( degree_distribution ( nem_random_graph (p, q, mu) ) )

#Fill with zero to have the same lenght
l = max ([len(item) for item in graphs])
graphs = np.array([ np.array(item + ([0] * (l - len(item)))) for item in graphs ])
#Get the mean degree distruibution
mean = np.mean(graphs, axis=0)
#Fill with zeros
mean = np.pad(mean, (0, len(th_deg_dist) - l), 'constant', constant_values=(0,0))
#std = np.std(graphs, axis=0)

#Plot both theoretical and mean distribution
plt.figure(1, figsize=(6,4))
plt.plot(mean, 'bo')
plt.plot(th_deg_dist, 'ro')
#Set axis
plt.xlabel(r'$k$')
plt.ylabel(r'$\mathbb{P}(deg(v) = k)$')

#Save plot to file
#svg
filesvg = 'fig/ch5fig2.svg'
plt.savefig(filesvg, format='svg', bbox_inches='tight')
#pdf
filepdf = 'fig/ch5fig2.pdf'
plt.savefig(filepdf, format='pdf', bbox_inches='tight')

#Show plot
#plt.show()
