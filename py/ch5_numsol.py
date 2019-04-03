import matplotlib.pyplot as plot
import numpy as np
from nem_graph import *
import math
from graph_tools import *
import networkx as nx

#Set parameters
number_of_mus = 61
mus = np.linspace(0.0,1.0, num=number_of_mus)

rhoA0 = 0.4
theta = 0.4

p = 100
half_p = math.floor(p/2)
q = 1500

rhos = []

#For each value of mu
for mu in mus:
    degdist = nem_th_deg_dist(p, q, mu)
    rhoB0 = 0.0
    rhoA = 0.0
    rhoB = 0.0
    check = True
    #Iterate untill...
    while check:
        rhoA_old = rhoA
        rhoB_old = rhoB
        qA = (1 - mu)*rhoA + mu * rhoB
        qB = mu * rhoA + (1 - mu) * rhoB
        rhoA = rhoA0 + (1 - rhoA0) * sum([degdist[k] * sum( [binomial(k, m, qA) for m in range(math.ceil(theta*k), k + 1)] ) for k in range(1, p)])
        rhoB = rhoB0 + (1 - rhoB0) * sum([degdist[k] * sum( [binomial(k, m, qB) for m in range(math.ceil(theta*k), k + 1)] ) for k in range(1, p)])
        #...difference between two step is less than 0.000001 in both A e B
        if abs(rhoA - rhoA_old) < 1e-6 and abs(rhoB - rhoB_old) < 1e-6:
            check = False
    rho = (rhoA + rhoB)/2
    #save the results
    rhos.append([rhoA, rhoB, rho])

#Prepare plot
plt.plot(mus, rhos)
#Restrict plotted mu values
plt.xlim(right = 0.5, left=0.0)
plt.xlabel(r'$\mu$')
plt.legend([r'$\rho_\infty^A$', r'$\rho_\infty^B$', r'$\rho_\infty$'], loc='upper right')
plt.title(r'$\rho_0^A = 0.4$ $\theta = 0.4$')

#Save plot to file
#svg
filesvg = 'fig/ch5fig2b.svg'
plt.savefig(filesvg, format='svg', bbox_inches='tight')
#pdf
filepdf = 'fig/ch5fig2b.pdf'
plt.savefig(filepdf, format='pdf', bbox_inches='tight')

#Show plot
#plt.show()
