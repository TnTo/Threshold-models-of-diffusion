from scipy.special import binom

#Compute (n k) p^k (1-p)^(n-k)
def binomial (n, k, p):
    if n < k:
        return 0
    else:
        return (binom(n, k) * (p ** k) * ((1 - p) ** (n - k)))

#Return the degree distribution of a graph as list
def degree_distribution (G):
    p = number_of_nodes(G)
    degree_sequence = [d for n, d in G.degree()]
    degree_dist = []
    for i in range(max(degree_sequence)+1):
        degree_dist.append(degree_sequence.count(i)/p)
    return degree_dist

#This function simulate the dynamics of the diffusion
def evolve (G, seed, theta):
    #The seed must be a list of verteices, represented by thei indices
    s_0 = seed
    #Iterate until equilibrium is reached
    check = True
    while (check):
        #List of active nodes at t+1
        s_1 = []
        for v in G.nodes:
            if v in s_0:
                #Active nodes at t will are active at t+1
                s_1.append(v)
            else:
                #Set of neighbours
                act_neigh = [s for s in G.neighbors(v) if s in s_0]
                if len(act_neigh) >= (theta[v] * G.degree[v]):
                    #Active neighbours are more than threshold?
                    s_1.append(v)
        if set(s_1) == set(s_0):
            #Equilibrium reached (the process is markovian)
            check = False
        else:
            #t+1 becomes new t
            s_0 = s_1
    return s_0
