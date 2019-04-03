import matplotlib.pyplot as plt
import numpy as np
from nem_graph import *
import networkx as nx
import math

#mu for graphs
mu = [0.05, 0.15, 0.30]

#order of the graphs
p = 250
halfp = math.floor(p/2)

#size of the graphs
q = 500

#Start plot
plt.figure(1, figsize=(12,4))
#list of subplots
subplots = [131, 132, 133]

#Generating random positions for vertexes
my_pos_A = nx.random_layout(range(halfp))
my_pos_B = nx.random_layout(range(halfp, p))

#Find the needed shift to separete community
x_shift = (max(item[0] for item in my_pos_A.values()) - min(item[0] for item in my_pos_B.values()))*0.25
y_shift = (max(item[1] for item in my_pos_B.values()) - min(item[1] for item in my_pos_A.values()))*1.5

#stretch graph
for i in my_pos_A:
    my_pos_A[i][0] = my_pos_A[i][0] * 2
    my_pos_A[i][1] = (my_pos_A[i][1] + y_shift) * 0.25
for i in my_pos_B:
    my_pos_B[i][0] = (my_pos_B[i][0] + x_shift) * 2
    my_pos_B[i][1] = my_pos_B[i][1] * 0.25

#get final positions of nodes
my_pos = {}
my_pos.update(my_pos_A)
my_pos.update(my_pos_B)

#Create the graphs
for i in range(len(mu)):
    #Create the graph
    G = nem_random_graph(p, q, mu[i], 1)
    plt.subplot(subplots[i])
    #Set nodes color
    color = ['r']*halfp + ['b']*(p-halfp)
    #Draw graph
    nx.draw_networkx(G, pos = my_pos, node_size=45, with_labels=False, node_color=color)
    #Caption
    plt.xlabel(r'$\mu = $' + str(mu[i]))

#Save plot to file
#svg
filesvg = 'fig/ch5fig1.svg'
plt.savefig(filesvg, format='svg', bbox_inches='tight')
#pdf
filepdf = 'fig/ch5fig1.pdf'
plt.savefig(filepdf, format='pdf', bbox_inches='tight')

#Show plot
#plt.show()
