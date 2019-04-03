import matplotlib.pyplot as plot
import numpy as np
from nem_graph import *
import math
from graph_tools import *
import csv
import networkx as nx
import pandas as pd
import random

#Set parameters
number_of_mus = 41
number_of_rho_zero = 41
number_of_graphs = 10
number_of_repetition = 10
number_of_thetas = 5

p = 100
half_p = math.floor(p/2)
q = 1500

#Set the values of variables
mus = np.linspace(0.0,1.0, num=number_of_mus)
rho_A_zero = np.linspace(0.0, 1.0, num=number_of_rho_zero)
thetas = np.linspace(0.1,0.5, num=number_of_thetas)

#Select operations to be performed
simulate = False
num_res = False
parse = False
plot = False
slide_plot = True

#Simulation
if simulate:
    #Select file in which save data
    data = open('data/nem_sim.csv', 'w')
    #Write columns' name in file
    fields = ['p', 'q', 'theta', 'mu', 'graph', 'rhoA0', 'rep', 'rhoA', 'rhoB', 'rho']
    writer = csv.writer(data, delimiter=",")
    writer.writerow(fields)
    for mu in mus:
        print('mu=' + str(mu))
        for i in range(number_of_graphs):
            #Create connected Nem. graph
            G = nem_random_graph(p, q, mu)
            while not (nx.is_connected(G) or mu == 0.0):
                G = nem_random_graph(p, q, mu)
            for theta in thetas:
                for seed_size in rho_A_zero:
                    for j in range(number_of_repetition):
                        #Select a random seed with the given dimension
                        seed = random.sample(range(half_p), math.floor(half_p*seed_size))
                        #Find the equilibrium
                        active = evolve(G, seed, [theta]*p)
                        #Measure values of interest
                        rho_A = len([v for v in active if v < half_p])/half_p
                        rho_B = len([v for v in active if v >= half_p])/(p-half_p)
                        rho = (rho_A + rho_B) / 2
                        #Save data to file
                        writer.writerow([p, q, theta, mu, i, seed_size, j, rho_A, rho_B, rho])
    data.close()

#To get numerical approximation
if num_res:
    #Load file to save data
    data = open('data/nem_num.csv', 'w')
    fields = ['p', 'q', 'theta', 'mu', 'rhoA0', 'rhoA', 'rhoB', 'rho']
    writer = csv.writer(data, delimiter=",")
    #write columns' name
    writer.writerow(fields)
    for mu in mus:
        print('mu=' + str(mu))
        #Get theoretical degree distribution
        degdist = nem_th_deg_dist(p, q, mu)
        for theta in thetas:
            for seed_size in rho_A_zero:
                rhoA0 = seed_size
                rhoB0 = 0.0
                rhoA = 0.0
                rhoB = 0.0
                #Iterate untill...
                check = True
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
                #save the results to file
                writer.writerow([p, q, theta, mu, seed_size, rhoA, rhoB, rho])
    data.close()

#Prepare data to plot
if parse:
    #load numerical data
    num_data = pd.read_csv('data/nem_num.csv', sep=',', header=0)
    #drop unutilized data
    num_data = num_data.drop(columns=['p','q'])
    #Divide data for theta's value
    theta_data = [num_data[num_data.theta == a_theta] for a_theta in thetas]
    #drop unutilized data
    theta_data = [df.drop(columns=['theta']) for df in theta_data]
    #Transform to matrix
    theta_data = [df.pivot(index='rhoA0', columns='mu') for df in theta_data]
    for i in range(len(thetas)):
        #Save data
        theta_data[i].to_pickle(('data/nem_num_th' + str(thetas[i]) + '.pkl'))
    del num_data
    del theta_data

    #Load data
    sim_data = pd.read_csv('data/nem_sim.csv', sep=',', header=0)
    #drop unutilized data
    sim_data = sim_data.drop(columns=['p','q', 'graph', 'rep'])
    #Divide data for theta's value
    theta_data = [sim_data[sim_data.theta == a_theta] for a_theta in thetas]
    #drop unutilized data
    theta_data = [df.drop(columns=['theta']) for df in theta_data]
    #mean the data
    theta_data = [df.groupby(['rhoA0', 'mu']).mean() for df in theta_data]
    #Transform to matrix
    theta_data = [df.unstack(level=-1) for df in theta_data]
    for i in range(len(thetas)):
        #Save data
        theta_data[i].to_pickle(('data/nem_sim_th' + str(thetas[i]) + '.pkl'))
    del sim_data
    del theta_data


if plot:
    #Select values of rhoA_0 for the plots in fourth columns
    cuts = [(1,3), (3,6), (10,12), (15,19), (24,28)]
    for theta in thetas:
        #Laod data
        data_sim = pd.read_pickle(('data/nem_sim_th' + str(theta) + '.pkl'))
        data_num = pd.read_pickle(('data/nem_num_th' + str(theta) + '.pkl'))
        #Prepare plot
        plt.figure(1, figsize=(20,10))
        plot, subplot = plt.subplots(nrows = 2, ncols = 4, num=1)
        plot.subplots_adjust(hspace=0.3, wspace=0.2)
        plt.rcParams.update({'font.size': 12})
        #plt.suptitle(r'$\theta = $' + str(round(theta, 1)), fontsize = 20)

        #Define rows and columns contents
        cols = ['rhoA', 'rhoB', 'rho']
        rows = [data_num, data_sim]

        for i in range(2):
            #phase diagrams plots
            for j in range(4):
                if j != 3:
                    tmp = subplot[i][j].pcolormesh(mus, rho_A_zero, rows[i][cols[j]], cmap='Reds')
                    plt.colorbar(tmp, ax=subplot[i][j])
                    subplot[i][j].set_xlabel(r'$\mu$')
                    subplot[i][j].set_ylabel(r'$\rho_0^A$')
            #fourth columns, section plot
            for data in np.split(rows[i].iloc[cuts[list(thetas).index(theta)][i], :],3):
                subplot[i][3].plot(mus, data)
            subplot[i][3].set_xlim(right = 0.5)
            subplot[i][3].set_xlabel(r'$\mu$')
            subplot[i][3].set_title(r'Section at $\rho_0^A =$' + str(round(rho_A_zero[cuts[list(thetas).index(theta)][i]],3)))
            subplot[i][3].legend([r'$\rho_\infty^A$', r'$\rho_\infty^B$', r'$\rho_\infty$'], loc='upper right')

        #Plots' titles
        subplot[0][0].set_title(r'$\rho_{\infty}^A$ - Numeric approximation')
        subplot[0][1].set_title(r'$\rho_{\infty}^B$ - Numeric approximation')
        subplot[0][2].set_title(r'$\rho_{\infty}$ - Numeric approximation')
        subplot[1][0].set_title(r'$\rho_{\infty}^A$ - Simulation')
        subplot[1][1].set_title(r'$\rho_{\infty}^B$ - Simulation')
        subplot[1][2].set_title(r'$\rho_{\infty}$ - Simulation')

        #Save plot to file
        #svg
        filesvg = 'fig/ch5fig3th'+ str(round(theta, 1)).replace('.', '') +'.svg'
        plt.savefig(filesvg, format='svg', bbox_inches='tight')
        #pdf
        filepdf = 'fig/ch5fig3th'+ str(round(theta, 1)).replace('.', '') +'.pdf'
        plt.savefig(filepdf, format='pdf', bbox_inches='tight')

        #Show plot
        #plt.show()

        #Clear plot
        plt.clf()

if slide_plot:
    #fix theta
    theta = thetas[2]
    #Load data
    data_sim = pd.read_pickle(('data/nem_sim_th' + str(theta) + '.pkl'))
    data_num = pd.read_pickle(('data/nem_num_th' + str(theta) + '.pkl'))
    #Prepare plot
    plt.figure(1, figsize=(20,10))
    plot, subplot = plt.subplots(nrows = 2, ncols = 3, num=1)
    plot.subplots_adjust(hspace=0.3, wspace=0.2)
    plt.rcParams.update({'font.size': 12})
    #plt.suptitle(r'$\theta = $' + str(round(theta, 1)), fontsize = 20)

    #Define rows and columns contents
    cols = ['rhoA', 'rhoB', 'rho']
    rows = [data_num, data_sim]

    #phase diagrams plots
    for i in range(2):
        for j in range(4):
            if j != 3:
                tmp = subplot[i][j].pcolormesh(mus, rho_A_zero, rows[i][cols[j]], cmap='Reds')
                plt.colorbar(tmp, ax=subplot[i][j])
                subplot[i][j].set_xlabel(r'$\mu$')
                subplot[i][j].set_ylabel(r'$\rho_0^A$')

    #Plots' titles
    subplot[0][0].set_title(r'$\rho_{\infty}^A$ - Numeric approximation')
    subplot[0][1].set_title(r'$\rho_{\infty}^B$ - Numeric approximation')
    subplot[0][2].set_title(r'$\rho_{\infty}$ - Numeric approximation')
    subplot[1][0].set_title(r'$\rho_{\infty}^A$ - Simulation')
    subplot[1][1].set_title(r'$\rho_{\infty}^B$ - Simulation')
    subplot[1][2].set_title(r'$\rho_{\infty}$ - Simulation')

    #Save plot to file
    #svg
    filesvg = 'fig/ch5slide.svg'
    plt.savefig(filesvg, format='svg', bbox_inches='tight')
    #pdf
    filepdf = 'fig/ch5slide.pdf'
    plt.savefig(filepdf, format='pdf', bbox_inches='tight')

    #plt.show()
