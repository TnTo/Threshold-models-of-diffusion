import matplotlib.pyplot as plt
import numpy as np

#order of graph
p = 10

#list of vertices
vertexes = np.arange(0, p)

#thresholds for 0, 1/p-1 ... p-2/p-1, 1 distribution (case A)
threshold_A = vertexes / (p-1)

#thresholds for 1/p-1, 1/p-1 ... p-2/p-1, 1 distribution (case B)
threshold_B = threshold_A.copy()
threshold_B[0] = threshold_B [1]

#define T
def T (n, threshold):
    return sum(1 for x in threshold if x*(p-1) <= n)

#We know from the model that T(n) channge in this case only for integer numbers
#Compute T for case A
T_values_A = []
for x in range(0, p+1):
    T_values_A.append(T(x, threshold_A))
T_values_A = np.array(T_values_A)
#Compute T for case B
T_values_B = []
for x in range(0, p+1):
    T_values_B.append(T(x, threshold_B))
T_values_B = np.array(T_values_B)

#Start plot
plt.figure(1, figsize=(12,5))
plt.subplots_adjust(hspace=0.5)

#Case A
plt.subplot(121)
#Set axis labels
plt.xlabel('x')
plt.ylabel('T(x)')
#Plot y=x
plt.plot(np.array([0, p]), np.array([0, p]),'b:')
#Plot T
plt.plot(np.arange(0, p+1), T_values_A, 'ro')
for i in range(0,p):
    plt.plot(np.array([i, i+1]), np.array([T(i, threshold_A), T(i, threshold_A)]), 'r-')

#case B
plt.subplot(122)
#Set axis labels
plt.xlabel('x')
plt.ylabel('T(x)')
#Plot y=x
plt.plot(np.array([0, p]), np.array([0, p]),'b:')
#Plot T
plt.plot(np.arange(0, p+1), T_values_B, 'ro')
for i in range(0,p):
    plt.plot(np.array([i, i+1]), np.array([T(i, threshold_B), T(i, threshold_B)]), 'r-')

#Save plot to file
#svg
filesvg = 'fig/ch4fig1.svg'
plt.savefig(filesvg, format='svg', bbox_inches='tight')
#pdf
filepdf = 'fig/ch4fig1.pdf'
plt.savefig(filepdf, format='pdf', bbox_inches='tight')

#Show plot
#plt.show()
