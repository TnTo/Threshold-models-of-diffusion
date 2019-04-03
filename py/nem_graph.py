import networkx as nx
import numpy.random as rnd
import math
import matplotlib.pyplot as plt
from graph_tools import *

def nem_random_graph(n, m, mu, seed = None):
    #Initialization
    rnd.seed(seed)
    G = nx.Graph()

    #Create n vertexes, half in A and half in B
    vertexes = range(n)
    half_n = math.floor(n/2)
    vertexes_A = range(half_n)
    vertexes_B = range(half_n, n)
    G.add_nodes_from(vertexes)

    #Create m edges, mu*M between A and B, the other in A or in B
    inter_edges_num = math.floor(mu * m)
    intra_edges_num = m - inter_edges_num

    #create the inter-community edges
    while (G.number_of_edges() < inter_edges_num):
        v = rnd.choice(vertexes)
        if v in vertexes_A:
            u = rnd.choice(vertexes_B)
        else:
            u = rnd.choice(vertexes_A)
        if v != u:
            G.add_edge(v,u)
            #Already existing edges are managed by NetworkX

    #create the remaining intra-community edges
    while (G.number_of_edges() < m):
        v = rnd.choice(vertexes)
        if v in vertexes_A:
            u = rnd.choice(vertexes_A)
        else:
            u = rnd.choice(vertexes_B)
        if v != u:
            G.add_edge(v,u)
            #Already existing edges are managed by NetworkX

    #return the graph
    return G

#Theoretical degree distribution of Nem. graphs
def nem_th_deg_dist(p, q, mu):
    s = 4.0 * mu * q / (p ** 2.0)
    r = (4.0 * (1.0 - mu) * q) / (p * (p - 2.0))
    deg_dist = []
    for k in range(q+1):
        sum = 0.0
        for n in range(k+1):
            b = binomial(math.floor(p/2), n, s)
            c = binomial(math.floor(p/2)-1, k-n, r)
            sum = sum + b*c
        deg_dist.append(sum)
    return deg_dist
    #Returned as list
