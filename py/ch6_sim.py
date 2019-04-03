import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import numpy as np
import math
from graph_tools import *
import csv
import networkx as nx
import pandas as pd
import scipy.stats

#Set parameters
number_of_pis = 26
number_of_thetas = 41
number_of_graphs = 50

p = 100
half_p = math.floor(p/2)

#Set the values of variables
pis = np.concatenate((np.linspace(0.0,0.1, num=(number_of_pis - 1), endpoint=False), np.linspace(0.1,0.5, num=number_of_pis, endpoint=True)), axis=None)
thetas = np.linspace(0.1,0.5, num=number_of_thetas)

#Select operations to be performed
simulate = False
num_res = False
parse = False
plot = True

#Simulation
if simulate:
    #Select file in which save data
    data = open('data/watts_sim.csv', 'w')
    #Write columns' name in file
    fields = ['p', 'pi', 'theta', 'graph', 'n', 'rho', 'G2_z']
    writer = csv.writer(data, delimiter=",")
    writer.writerow(fields)
    for pi in pis:
        #print('pi=' + str(pi))
        for i in range(number_of_graphs):
            #Create Erdos Renyi graph
            G = nx.fast_gnp_random_graph(p, pi)
            for theta in thetas:
                #find max degree of vulnerable vertices
                vuln_degree = math.floor(1.0 / theta)
                #Find vulnerable vertices
                vulns = [v for (v,k) in G.degree() if k <= vuln_degree]
                if vulns != []:
                    #Find connected componets of the subgraph induced by vulnerable vertices
                    conn_comp = [list(item) for item in list(nx.connected_components(G.subgraph(vulns)))]
                    #and their lenght
                    conn_comp_size = [len(item) for item in conn_comp]
                    n = max(conn_comp_size)
                    #choose as initial seed a vertex in the biggest component
                    seed = [conn_comp[conn_comp_size.index(n)][0]]
                else:
                    #If there aren't vulnerable vertices initial seed is arbitrary
                    n = 0
                    seed = [0]
                #Find the fraction of active nodes at the equilibrium
                rho = len(evolve(G, seed, [theta]*p))/p
                #Get the degree degree_distribution of the graph
                deg_dist = degree_distribution (G)
                #Calculate G''(1)/z
                if G.number_of_edges() != 0:
                    G2_z = sum([(k*(k-1)*deg_dist[k]) for k in range(2, min(vuln_degree + 1, len(deg_dist)))]) / (2*G.number_of_edges() / float(G.number_of_nodes()))
                else:
                    G2_z = None
                #print (n, len(rho), G2_z)
                #Save results to file
                writer.writerow([p, pi, theta, i, n, rho, G2_z])
    data.close()

#To get numerical approximation
if num_res:
    #Load file to save data
    data = open('data/watts_num.csv', 'w')
    fields = ['p', 'pi', 'theta', 'n', 'G2_z']
    writer = csv.writer(data, delimiter=",")
    #write columns' name
    writer.writerow(fields)
    for pi in pis:
        print('pi=' + str(pi))
        #Get the theoretical degree distribution of an Erdos Reyi graph
        degdist = [scipy.stats.binom.pmf(k, p-1, pi) for k in range (p)]
        for theta in thetas:
            #Find max degree of vulnerable nodes
            vuln_degree = math.floor(1.0 / theta)
            #Compute G''(1)/z
            G2_z = sum([k*(k-1)*degdist[k]/(pi*(p-1)) for k in range(2, vuln_degree + 1)])
            #compute mean n
            n = sum([degdist[k] for k in range (vuln_degree + 1)]) + sum([k*degdist[k] for k in range (vuln_degree + 1)])**2 / ((1 - G2_z)*(pi*(p-1)))
            #save data to file
            writer.writerow([p, pi, theta, n, G2_z])
    data.close()

#Prepare data to plot
if parse:
    #load numerical data
    num_data = pd.read_csv('data/watts_num.csv', sep=',', header=0)
    #drop unutilized data
    num_data = num_data.drop(columns=['p'])
    num_data = num_data[num_data.pi != 0]
    #Transform to matrix
    num_data = num_data.pivot(index='pi', columns='theta')
    #Save data
    num_data.to_pickle(('data/watts_num.pkl'))
    del num_data

    #Load data
    sim_data = pd.read_csv('data/watts_sim.csv', sep=',', header=0)
    #drop unutilized data
    sim_data = sim_data.drop(columns=['p', 'graph'])
    sim_data = sim_data[sim_data.pi != 0]
    #Aggregate data and transform to matrix
    sim_data = sim_data.groupby(['pi', 'theta']).mean()
    sim_data = sim_data.unstack(level=-1)
    #Save to file
    sim_data.to_pickle(('data/watts_sim.pkl'))
    del sim_data

if plot:
    #Load data
    data_sim = pd.read_pickle(('data/watts_sim.pkl'))
    data_num = pd.read_pickle(('data/watts_num.pkl'))
    #Prepare plot
    plt.figure(1, figsize=(20,10))
    plot, subplot = plt.subplots(nrows = 2, ncols = 3, num=1)
    plot.subplots_adjust(hspace=0.3)
    plt.rcParams.update({'font.size': 12})

    #phase diagrams plots
    tmp = subplot[0][0].pcolormesh(thetas, pis[1:],  data_sim['G2_z'], cmap='Blues')
    plt.colorbar(tmp, ax=subplot[0][0])
    subplot[0][0].set_xlabel(r'$\theta$')
    subplot[0][0].set_ylabel(r'$\pi$')
    subplot[0][0].set_title(r"$G''(1) / z$ - Simulation")
    subplot[0][0].set_xlim(right = 0.4)
    subplot[0][0].set_ylim(top = 0.3)

    tmp = subplot[0][1].pcolormesh(thetas, pis[1:],  data_sim['n'], cmap='Blues')
    plt.colorbar(tmp, ax=subplot[0][1])
    subplot[0][1].set_xlabel(r'$\theta$')
    subplot[0][1].set_ylabel(r'$\pi$')
    subplot[0][1].set_title(r'$n_{max}$ - Simulation')
    subplot[0][1].set_xlim(right = 0.4)
    subplot[0][1].set_ylim(top = 0.3)

    tmp = subplot[0][2].pcolormesh(thetas, pis[1:],  data_sim['rho'], cmap='Blues')
    plt.colorbar(tmp, ax=subplot[0][2])
    subplot[0][2].set_xlabel(r'$\theta$')
    subplot[0][2].set_ylabel(r'$\pi$')
    subplot[0][2].set_title(r'$|s_{\infty}|$ - Simulation')
    subplot[0][2].set_xlim(right = 0.4)
    subplot[0][2].set_ylim(top = 0.3)

    tmp = subplot[1][0].pcolormesh(thetas, pis[1:],  data_num['G2_z'], cmap='Blues')
    plt.colorbar(tmp, ax=subplot[1][0])
    subplot[1][0].set_xlabel(r'$\theta$')
    subplot[1][0].set_ylabel(r'$\pi$')
    subplot[1][0].set_title(r"$G''(1) / z$ - Model approximation")
    subplot[1][0].set_xlim(right = 0.4)
    subplot[1][0].set_ylim(top = 0.3)

    tmp = subplot[1][1].pcolormesh(thetas, pis[1:],  data_num['n'], cmap='RdBu', norm = mplcolors.Normalize(vmin=-50.,vmax=50.))
    plt.colorbar(tmp, ax=subplot[1][1])
    subplot[1][1].set_xlabel(r'$\theta$')
    subplot[1][1].set_ylabel(r'$\pi$')
    subplot[1][1].set_title(r'$\langle n \rangle$ - Model approximation')
    subplot[1][1].set_xlim(right = 0.4)
    subplot[1][1].set_ylim(top = 0.3)

    #Delete unused plot space
    plot.delaxes(subplot[1][2])

    #Save plot to file
    #svg
    filesvg = 'fig/ch6fig1.svg'
    plt.savefig(filesvg, format='svg', bbox_inches='tight')
    #pdf
    filepdf = 'fig/ch6fig1.pdf'
    plt.savefig(filepdf, format='pdf', bbox_inches='tight')

    #Show plot
    #plt.show()
