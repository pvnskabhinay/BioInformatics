# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 18:37:40 2019

@author: pvnsk
"""
###############################################################################
                                #1.1
###############################################################################

import pyedflib
import numpy as np
import connectivipy as cp
import matplotlib.pyplot as plt

####### Eyes open ########
filename1 = "C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\S064\\S064R01.edf"
eyes_open = pyedflib.EdfReader(filename1)

k = eyes_open.signals_in_file
N = eyes_open.getNSamples()[0]

signalOpen = np.zeros((k, N, 1))
for i in np.arange(k):
    signalOpen[i, :, 0] = eyes_open.readSignal(i)

labelsOpen = eyes_open.getSignalLabels()

labelsOpen  = list(map(lambda x: x.replace('.',''),  labelsOpen ))
labelsD = dict(enumerate(labelsOpen))
numericD = {key: idx for idx, key in enumerate(labelsOpen)}
####### Eyes closed ########

filename2 = "C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\S064\\S064R02.edf"
eyes_closed = pyedflib.EdfReader(filename2)

k = eyes_closed.signals_in_file
N = eyes_closed.getNSamples()[0]

signalClosed = np.zeros((k, N, 1))
for i in np.arange(k):
    signalClosed[i, :, 0] = eyes_closed.readSignal(i)

labelsClosed = eyes_closed.getSignalLabels()

labelsClosed  = list(map(lambda x: x.replace('.',''),  labelsClosed))

########### DTF #############
# Eyes open
mv = cp.Mvar # static class mvar
# Transforming matrix to connectivipy 
dataOpen = cp.Data(data=signalOpen, fs=eyes_open.getSampleFrequency(0),chan_names=labelsOpen)

#https://connectivipy.readthedocs.io/en/latest/tutorial.html
# using this link try to find the best fit model using Vieira-Morf algorithm
bestOpen, crit = mv.order_akaike(signalOpen, 15, 'vm')
plt.plot(1+np.arange(len(crit)), crit, 'g')
plt.show()
plt.title('Best fit for eyes-open case')
print(bestOpen)
# use Yule Walker algorithm to fit the mvar model
dataOpen.fit_mvar(bestOpen, 'yw')

#DTF
dtf = cp.conn.DTF()
avOpen, vfOpen = dataOpen.mvar_coefficients
dtfval_open = dtf.calculate(avOpen, vfOpen, eyes_open.getSampleFrequency(0))


# Eyes closed
mv = cp.Mvar # static class mvar
# Transforming matrix to connectivipy 
dataClosed = cp.Data(data=signalClosed, fs=eyes_closed.getSampleFrequency(0),chan_names=labelsClosed)

#https://connectivipy.readthedocs.io/en/latest/tutorial.html
# using this link try to find the best fit model using Vieira-Morf algorithm
bestClosed, crit = mv.order_akaike(signalClosed, 15, 'vm')
plt.plot(1+np.arange(len(crit)), crit, 'g')
plt.show()
plt.title("Best fit for eyes-open case")
print(bestClosed)
# use Yule Walker algorithm to fit the mvar model
dataClosed.fit_mvar(bestClosed, 'yw')

#DTF
dtf = cp.conn.DTF()
avClosed, vfClosed = dataClosed.mvar_coefficients
dtfval_closed = dtf.calculate(avClosed, vfClosed, eyes_closed.getSampleFrequency(0))



# Desired density 20%
# eyes open
tmpmatrix20 = dtfval_open[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix20, 0)
desired_density=0.2
tol=1e-3
binSpace = [0.0, 1.0]
finalVal20, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal20 != tmp:
    finalVal20 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix20 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

# eyes closed
tmpmatrix20_1 = dtfval_closed[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix20_1, 0)
desired_density=0.2
tol=1e-3
binSpace = [0.0, 1.0]
finalVal20_1, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal20_1 != tmp:
    finalVal20_1 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix20_1 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

# Binary adjacency matrix
# eyes open

binaryMatrixOpen = np.array(tmpmatrix20 > finalVal20)
plt.imshow(binaryMatrixOpen)
plt.title('Graphical representation of the binary adjacent matrix - eyes open case')
plt.show()

binaryMatrixClosed = np.array(tmpmatrix20_1 > finalVal20_1)
plt.imshow(binaryMatrixClosed)
plt.title('Graphical representation of the binary adjacent matrix - eyes closed case')
plt.show()


###############################################################################
                                #1.3
###############################################################################

# Desired density 1%
# eyes open
tmpmatrix1 = dtfval_open[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix1, 0)
desired_density=0.01
tol=1e-3
binSpace = [0.0, 1.0]
finalVal1, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal1 != tmp:
    finalVal1 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix1 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

# eyes closed
tmpmatrix1_1 = dtfval_closed[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix1_1, 0)
desired_density=0.01
tol=1e-3
binSpace = [0.0, 1.0]
finalVal1_1, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal1_1 != tmp:
    finalVal1_1 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix1_1 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

# Binary adjacency matrix
# eyes open

binaryMatrixOpen = np.array(tmpmatrix1 > finalVal1)
plt.imshow(binaryMatrixOpen)
plt.title('Graphical representation of the binary adjacent matrix - eyes open case - 1%')
plt.show()

binaryMatrixClosed = np.array(tmpmatrix1_1 > finalVal1_1)
plt.imshow(binaryMatrixClosed)
plt.title('Graphical representation of the binary adjacent matrix - eyes closed case - 1%')
plt.show()

# Desired density 5%
# eyes open
tmpmatrix5 = dtfval_open[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix5, 0)
desired_density=0.05
tol=1e-3
binSpace = [0.0, 1.0]
finalVal5, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal5 != tmp:
    finalVal5 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix5 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

# eyes closed
tmpmatrix5_1 = dtfval_closed[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix5_1, 0)
desired_density=0.05
tol=1e-3
binSpace = [0.0, 1.0]
finalVal5_1, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal5_1 != tmp:
    finalVal5_1 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix5_1 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

binaryMatrixOpen = np.array(tmpmatrix5 > finalVal5)
plt.imshow(binaryMatrixOpen)
plt.title('Graphical representation of the binary adjacent matrix - eyes open case - 5%')
plt.show()

binaryMatrixClosed = np.array(tmpmatrix5_1 > finalVal5_1)
plt.imshow(binaryMatrixClosed)
plt.title('Graphical representation of the binary adjacent matrix - eyes closed case - 5%')
plt.show()

# Desired density 10%
# eyes open
tmpmatrix10 = dtfval_open[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix10, 0)
desired_density=0.1
tol=1e-3
binSpace = [0.0, 1.0]
finalVal10, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal10 != tmp:
    finalVal10 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix10 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

# eyes closed
tmpmatrix10_1 = dtfval_closed[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix10_1, 0)
desired_density=0.1
tol=1e-3
binSpace = [0.0, 1.0]
finalVal10_1, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal10_1 != tmp:
    finalVal10_1 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix10_1 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

binaryMatrixOpen = np.array(tmpmatrix10 > finalVal10)
plt.imshow(binaryMatrixOpen)
plt.title('Graphical representation of the binary adjacent matrix - eyes open case - 10%')
plt.show()

binaryMatrixClosed = np.array(tmpmatrix10_1 > finalVal10_1)
plt.imshow(binaryMatrixClosed)
plt.title('Graphical representation of the binary adjacent matrix - eyes closed case - 10%')
plt.show()

# Desired density 20%
# eyes open
tmpmatrix20 = dtfval_open[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix20, 0)
desired_density=0.2
tol=1e-3
binSpace = [0.0, 1.0]
finalVal20, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal20 != tmp:
    finalVal20 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix20 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

# eyes closed
tmpmatrix20_1 = dtfval_closed[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix20_1, 0)
desired_density=0.2
tol=1e-3
binSpace = [0.0, 1.0]
finalVal20_1, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal20_1 != tmp:
    finalVal20_1 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix20_1 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

binaryMatrixOpen = np.array(tmpmatrix20 > finalVal20)
plt.imshow(binaryMatrixOpen)
plt.title('Graphical representation of the binary adjacent matrix - eyes open case - 20%')
plt.show()

binaryMatrixClosed = np.array(tmpmatrix20_1 > finalVal20_1)
plt.imshow(binaryMatrixClosed)
plt.title('Graphical representation of the binary adjacent matrix - eyes closed case - 20%')
plt.show()

# Desired density 30%
# eyes open
tmpmatrix30 = dtfval_open[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix30, 0)
desired_density=0.3
tol=1e-3
binSpace = [0.0, 1.0]
finalVal30, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal30 != tmp:
    finalVal30 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix30 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

# eyes closed
tmpmatrix30_1 = dtfval_closed[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix30_1, 0)
desired_density=0.3
tol=1e-3
binSpace = [0.0, 1.0]
finalVal30_1, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal30_1 != tmp:
    finalVal30_1 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix30_1 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

binaryMatrixOpen = np.array(tmpmatrix30 > finalVal30)
plt.imshow(binaryMatrixOpen)
plt.title('Graphical representation of the binary adjacent matrix - eyes open case - 30%')
plt.show()

binaryMatrixClosed = np.array(tmpmatrix30_1 > finalVal30_1)
plt.imshow(binaryMatrixClosed)
plt.title('Graphical representation of the binary adjacent matrix - eyes closed case - 30%')
plt.show()

# Desired density 50%
# eyes open
tmpmatrix50 = dtfval_open[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix50, 0)
desired_density=0.5
tol=1e-3
binSpace = [0.0, 1.0]
finalVal50, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal50 != tmp:
    finalVal50 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix50 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

# eyes closed
tmpmatrix50_1 = dtfval_closed[11, :, :].reshape((64, 64))
np.fill_diagonal(tmpmatrix50_1, 0)
desired_density=0.5
tol=1e-3
binSpace = [0.0, 1.0]
finalVal50_1, tmp, best_n = -1, 0, 4096
target = int(k*(k-1)*desired_density)

while finalVal50_1 != tmp:
    finalVal50_1 = tmp
    for alpha in np.arange(binSpace[0], binSpace[1], (binSpace[1] - binSpace[0])/100.0):
        n = len(np.where(tmpmatrix50_1 > alpha)[0])
        if abs(n-target) < abs(best_n-target) and n-target>0:
            tmp, best_n = alpha, n
    binSpace[0], binSpace[1] = tmp-tol, tmp+tol

best_n
best_n/(k*(k-1))
print("Network density for eyes open graph after the threshold is:"+str(best_n/(k*(k-1))))

binaryMatrixOpen = np.array(tmpmatrix50 > finalVal50)
plt.imshow(binaryMatrixOpen)
plt.title('Graphical representation of the binary adjacent matrix - eyes open case - 50%')
plt.show()

binaryMatrixClosed = np.array(tmpmatrix50_1 > finalVal50_1)
plt.imshow(binaryMatrixClosed)
plt.title('Graphical representation of the binary adjacent matrix - eyes closed case - 50%')
plt.show()

###############################################################################
                                #1.5
###############################################################################
import pandas as pd
import networkx as nx


# we created this text file using the following process:
# 1st we used the library eegkit, which has the coordinates for 87 channels
# then we write these coordinates into a csv file. 
# The R script for this can be found in the file: eegCoordfromR.R
# Then using basic excel functions index and match, we get the coordinates for our 64 channels
# and save it in a text file
# https://cran.r-project.org/web/packages/eegkit/eegkit.pdf
position_data = pd.read_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\channel_locations.txt', sep="\t", index_col=0)
position_data = {label.replace('.', ''): (x, y) for label, x, y in zip(position_data['label'], position_data['x'], position_data['y'])}

#position_data = create_position_graph('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\channel_locations.txt')

def create_directed_graph(input_matrix, alpha, position, weighted=False):
    G = nx.DiGraph()
    x, y = np.where(input_matrix > alpha)
    G.add_nodes_from(position_data.keys(), pos=position)
    G.add_edges_from(list(zip(map(lambda idx: labelsD[idx], x), map(lambda idy: labelsD[idy], y))))
    return G



directed_graph_open = create_directed_graph(input_matrix=tmpmatrix20,
                                            alpha=finalVal20,
                                            position=position_data)

#directed_graph_open = create_directed_graph(input_matrix=tmpmatrix1,
#                                            alpha=finalVal1,
#                                            position=position_data)

directed_graph_closed = create_directed_graph(input_matrix=tmpmatrix20_1,
                                              alpha=finalVal20_1,
                                              position=position_data)

nx.draw_networkx(directed_graph_open, position_data)
plt.title('Topological network representation - eyes open case')
plt.show()

nx.draw_networkx(directed_graph_closed, position_data)
plt.title('Topological network representation - eyes closed case')
plt.show()

###############################################################################
                                #2.1
###############################################################################
import igraph as ig

def global_indices(graph):
    coeffs = nx.average_clustering(graph.to_undirected())
    pathLen = nx.average_shortest_path_length(graph.to_undirected())
    return coeffs, pathLen

def local_indices(graph):
    g = ig.Graph.TupleList(graph.edges(), directed=True)
    Vertices = ig.VertexSeq(g)
    node, degree, indegree, outdegree = [], [], [], []        
    for i in Vertices:
        node.append(i['name'])
        indegree.append(g.degree(i['name'], type="in"))
        outdegree.append(g.degree(i['name'], type="out"))
        degree.append(g.degree(i['name']))
    Datall = {'node': node, 'in-degree': indegree, 'out-degree': outdegree,'degree': degree}
    deg = pd.DataFrame(data=Datall)
    deg = deg[['node', 'in-degree', 'out-degree', 'degree']]
    deg=deg.sort_values(by=['degree'], ascending=False)
    indeg=deg.sort_values(by=['in-degree'], ascending=False)
    outdeg=deg.sort_values(by=['out-degree'], ascending=False)
    return(deg,indeg,outdeg)

'''For Eyes Open'''
coeff,path_len = global_indices(directed_graph_open)
degreeOpen,indegreeOpen,outdegreeOpen = local_indices(directed_graph_open)
print ('Average Clustering Coefficient for Eyes Open Graph is :'+str(coeff))
print ('Average Shortest Path Length for Eyes Open Graph is :'+str(path_len))
print ('Top Ten Channels for Local indices Are :')
degreeOpen.head(10)
print ('Top Ten Channels for in-degree Are :')
indegreeOpen.head(10)
print ('Top Ten Channels for out-degree Are :')
outdegreeOpen.head(10)

'''For Eyes Closed'''
coeff,path_len = global_indices(directed_graph_closed)
degreeClosed,indegreeClosed,outdegreeClosed = local_indices(directed_graph_closed)
print ('Average Clustering Coefficient for Eyes Closed Graph is :'+str(coeff))
print ('Average Shortest Path Length for Eyes Closed Graph is :'+str(path_len))
print ('Top Ten Channels for Local indices Are :')
degreeClosed.head(10)
print ('Top Ten Channels for in-degree Are :')
indegreeClosed.head(10)
print ('Top Ten Channels for out-degree Are :')
outdegreeClosed.head(10)


###############################################################################
                                #2.2
###############################################################################
coeffOpen,path_lenOpen = global_indices(directed_graph_open)
coeffClose,path_lenClose = global_indices(directed_graph_closed)
nodesNumOpen = len(directed_graph_open.nodes)
edgesNumOpen = len(directed_graph_open.edges)
nodesNumClose = len(directed_graph_closed.nodes)
edgesNumClose = len(directed_graph_closed.edges)

# https://stats.stackexchange.com/a/211047
Grandom = nx.gnm_random_graph(nodesNumOpen, edgesNumOpen)
coeffGr, path_lenGr = global_indices(Grandom)
# coeffGr = 0.40088298988557414
# path_lenGr = 1.599702380952381
lambda1 = path_lenOpen/path_lenGr
gamma1 = coeffOpen/coeffGr
# there are 2 conditions: 
# 1. lambda≈1 and 2. gamma>1
# both are satisfied, so this network can be called a small-world network


# trying r library
import rpy2
print(rpy2.__version__)
import os
os.environ['R_HOME'] = 'C:\Program Files\R\R-3.5.1' #path to your R installation
os.environ['R_USER'] = 'c:\program files (x86)\microsoft visual studio\shared\anaconda3_64\lib\site-packages (2.9.5)' #path depends on where you installed Python. Mine is the Anaconda distribution

from rpy2.robjects.packages import importr
base = importr('base')

import rpy2.robjects.packages as rpackages
utils = rpackages.importr('utils')
utils.install_packages("qgraph")
qgraph = importr('qgraph')

utils.install_packages("igraph")
igraph = importr('igraph')

import rpy2.robjects.numpy2ri
from rpy2.robjects.numpy2ri import numpy2ri
ro.conversion.py2ri = numpy2ri
rpy2.robjects.numpy2ri.activate()


from itertools import chain
nodeList = list(directed_graph_open.nodes)
edgeList = list(directed_graph_open.edges)
edgeList = list(chain.from_iterable(edgeList))
np.savetxt("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\edgeListSmallWorld.csv", edgeList, delimiter=",", fmt='%s')

edgeListnp = np.asarray(edgeList)

shape1 = (int(len(edgeListnp)/2),2)
edgeListnp = edgeListnp.reshape(shape1)

g2 = igraph.graph_from_edgelist(edgeListnp)
dta = qgraph.smallworldIndex(g2)

'''
the output is as follows:
R object with classes: ('list',) mapped to:
<ListVector - Python:0x000002B061E59D48 / R:0x000002B05F7B7F90>
[FloatVe..., FloatVe..., FloatVe..., FloatVe..., FloatVe...]
  transitivity: <class 'rpy2.robjects.vectors.FloatVector'>
  R object with classes: ('numeric',) mapped to:
<FloatVector - Python:0x000002B061E59808 / R:0x000002B05FB34B10>
[0.498335]
  transitivity_random: <class 'rpy2.robjects.vectors.FloatVector'>
  R object with classes: ('numeric',) mapped to:
<FloatVector - Python:0x000002B061E59108 / R:0x000002B05FFA0FE0>
[0.197021]
  APL: <class 'rpy2.robjects.vectors.FloatVector'>
  R object with classes: ('numeric',) mapped to:
<FloatVector - Python:0x000002B061E59548 / R:0x000002B05FFA0DA0>
[1.616667]
  APL_random: <class 'rpy2.robjects.vectors.FloatVector'>
  R object with classes: ('numeric',) mapped to:
<FloatVector - Python:0x000002B061E59188 / R:0x000002B06024D428>
[1.913198]
  index: <class 'rpy2.robjects.vectors.FloatVector'>
  R object with classes: ('numeric',) mapped to:
<FloatVector - Python:0x000002B061E59A08 / R:0x000002B061EF9C68>
[2.993280]
'''

transitivity = dta[0][0] # 0.498335
transitivity_random = dta[1][0] # 0.197021
APL = dta[2][0] # 1.616667
APL_random = dta[3][0] # 1.913198
SWI = dta[4][0] # 2.993280

# we double-checked these values in R also, they were similar
print("The Small worldness Index of this network is: ", SWI)

'''
# latticize the graph to get C_l and L_l
latG = nx.lattice.Graph(incoming_graph_data=nx.Graph.(nodesNumOpen, edgesNumOpen))
C_l, L_l = global_indices(latG) # 0.40088298988557414, 1.599702380952381
C_r = coeffGr # 0.40088298988557414
L_r = path_lenGr # 1.599702380952381
C = coeffOpen #  0.7370509093949917
L = path_lenOpen # 1.6364087301587302
# Now according to the formula given here: https://en.wikipedia.org/wiki/Small-world_network
# The small world index SWI is given as 
SWI = ((L-L_l)/(L_r-L_l))*((C-C_r)/C_l-C_r)'''


###############################################################################
                                #2.4
###############################################################################
                        
directed_graph_open1 = create_directed_graph(input_matrix=tmpmatrix1,
                                            alpha=finalVal1,
                                            position=position_data)

directed_graph_closed1 = create_directed_graph(input_matrix=tmpmatrix1_1,
                                              alpha=finalVal1_1,
                                              position=position_data)

directed_graph_open5 = create_directed_graph(input_matrix=tmpmatrix5,
                                            alpha=finalVal5,
                                            position=position_data)

directed_graph_closed5 = create_directed_graph(input_matrix=tmpmatrix5_1,
                                              alpha=finalVal5_1,
                                              position=position_data)

directed_graph_open10 = create_directed_graph(input_matrix=tmpmatrix10,
                                            alpha=finalVal10,
                                            position=position_data)

directed_graph_closed10 = create_directed_graph(input_matrix=tmpmatrix10_1,
                                              alpha=finalVal10_1,
                                              position=position_data)

directed_graph_open20 = create_directed_graph(input_matrix=tmpmatrix20,
                                            alpha=finalVal20,
                                            position=position_data)

directed_graph_closed20 = create_directed_graph(input_matrix=tmpmatrix20_1,
                                              alpha=finalVal20_1,
                                              position=position_data)

directed_graph_open30 = create_directed_graph(input_matrix=tmpmatrix30,
                                            alpha=finalVal30,
                                            position=position_data)

directed_graph_closed30 = create_directed_graph(input_matrix=tmpmatrix30_1,
                                              alpha=finalVal30_1,
                                              position=position_data)

directed_graph_open50 = create_directed_graph(input_matrix=tmpmatrix50,
                                            alpha=finalVal50,
                                            position=position_data)

directed_graph_closed50 = create_directed_graph(input_matrix=tmpmatrix50_1,
                                              alpha=finalVal50_1,
                                              position=position_data)



'''Network Indices as a Measure of Density for Eyes Open'''
ClusteringCoeff_Open, PathLength_Open, ClusteringCoeff_close, PathLength_close = [], [], [], []
#ClusteringCoeff_Open.append(global_indices(directed_graph_open1)[0])
ClusteringCoeff_Open.append(global_indices(directed_graph_open5)[0])
ClusteringCoeff_Open.append(global_indices(directed_graph_open10)[0])
ClusteringCoeff_Open.append(global_indices(directed_graph_open20)[0])
ClusteringCoeff_Open.append(global_indices(directed_graph_open30)[0])
ClusteringCoeff_Open.append(global_indices(directed_graph_open50)[0])
#PathLength_Open.append(global_indices(directed_graph_open1)[1])
PathLength_Open.append(global_indices(directed_graph_open5)[1])
PathLength_Open.append(global_indices(directed_graph_open10)[1])
PathLength_Open.append(global_indices(directed_graph_open20)[1])
PathLength_Open.append(global_indices(directed_graph_open30)[1])
PathLength_Open.append(global_indices(directed_graph_open50)[1])
'''Network Indices as a Measure of Density for Eyes Close'''
#ClusteringCoeff_close.append(global_indices(directed_graph_closed1)[0])
#ClusteringCoeff_close.append(global_indices(directed_graph_closed5)[0])
ClusteringCoeff_close.append(global_indices(directed_graph_closed10)[0])
ClusteringCoeff_close.append(global_indices(directed_graph_closed20)[0])
ClusteringCoeff_close.append(global_indices(directed_graph_closed30)[0])
ClusteringCoeff_close.append(global_indices(directed_graph_closed50)[0])
#PathLength_close.append(global_indices(directed_graph_closed1)[1])
#PathLength_close.append(global_indices(directed_graph_closed5)[1])
PathLength_close.append(global_indices(directed_graph_closed10)[1])
PathLength_close.append(global_indices(directed_graph_closed20)[1])
PathLength_close.append(global_indices(directed_graph_closed30)[1])
PathLength_close.append(global_indices(directed_graph_closed50)[1])

# Clustering Coefficients vs Densities
cco = [.05,.10,.20,.30,.50] 
ccc = [.10,.20,.30,.50] # for 5% density the graph was not connected
plt.plot(cco,ClusteringCoeff_Open,c='red')
plt.xlabel('Density')
plt.ylabel('Clustering Coefficient')
plt.plot(ccc,ClusteringCoeff_close,c='blue')
plt.xlabel('Density')
plt.ylabel('Clustering Coefficient')
plt.title('Density Vs Clustering Coefficient')
plt.legend(['Eyes Open Clustering Coeff','Eyes Close Clustering Coeff'])
plt.show()

# Path Lengths vs Densities
plo = [.05,.10,.20,.30,.50]
plc = [.10,.20,.30,.50] # for 5% density the graph was not connected
plt.plot(plo,PathLength_Open,c='red')
plt.xlabel('Density')
plt.ylabel('Avg Path Length')
plt.plot(plc,PathLength_close,c='blue')
plt.xlabel('Density')
plt.ylabel('Avg Path Length')
plt.title('Density Vs Avg Path Length')
plt.legend(['Eyes Open Avg Path','Eyes CloseAvg Path'])
plt.show()

###############################################################################
                                #3.1
###############################################################################

########### here we use the software mfinder

with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_open.txt', 'w') as f:
    for i, j in directed_graph_open20.edges():
        f.write("\t".join(map(str, [numericD[i]+1, numericD[j]+1, 1]))+"\n")

with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_closed.txt', 'w') as f:
    for i, j in directed_graph_closed20.edges():
        f.write("\t".join(map(str, [numericD[i]+1, numericD[j]+1, 1]))+"\n")

import subprocess

# we need to pass parameters to mfinder to open it (as in cmd). The sytax is as follows
#[mfinder.exe,input txt file, -s (motif size), size, -r, number of random networks to generate, number, -f, output filename]
mfinderLink = 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\mfinder1.2\\mfinder1.2.exe'
parameters = [mfinderLink, 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_open.txt', '-omem','-s', '3', '-r', '200', '-f', 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\outputsMfinderOpen']
p = subprocess.Popen(parameters)
parameters = [mfinderLink, 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_closed.txt', '-omem', '-s', '3', '-r', '200', '-f', 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\outputsMfinderClosed']
p = subprocess.Popen(parameters)

# from this we get the results in txt files, we have copied those results into csv files, which we will load to display the results

motifstatsDictOpen = pd.read_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\motifstats-open.csv', ).to_dict()
motifstatsDictOpen.keys()
#['MOTIF ID', 'Frequency', 'Interval', 'z score', 'p - value']
pd.DataFrame(motifstatsDictOpen)

motifstatsDictClosed = pd.read_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\motifstats-closed.csv', ).to_dict()
motifstatsDictClosed.keys()
pd.DataFrame(motifstatsDictClosed)

def calculateInterval(dct):
    for i in range(len(dct['Interval'])):
        dct['lowerLimit'][i] = float(dct['Interval'][i].split('+-')[0]) - float(dct['Interval'][i].split('+-')[1])
        dct['upperLimit'][i] = float(dct['Interval'][i].split('+-')[0]) + float(dct['Interval'][i].split('+-')[1])
    return dct
motifstatsDictOpen['lowerLimit']=dict()
motifstatsDictOpen['upperLimit']=dict()
motifstatsDictClosed['lowerLimit']=dict()
motifstatsDictClosed['upperLimit']=dict()
motifstatsDictOpen = calculateInterval(motifstatsDictOpen)
motifstatsDictClosed = calculateInterval(motifstatsDictClosed)
pd.DataFrame(motifstatsDictOpen)
pd.DataFrame(motifstatsDictClosed)

# To figure out the statistical significance of each row, we apply the following:
# if frequency of a row > upper limit = motif
# if frequency of a row < lower limit = anti-motif

def statSignificance(dct):
    for i in range(len(dct['MOTIF ID'])):
        if dct['Frequency'][i] > dct['upperLimit'][i]:
            dct['Statistical Significance'][i] = ['MOTIF']
        elif dct['Frequency'][i] < dct['lowerLimit'][i]:
            dct['Statistical Significance'][i] = ['ANTI-MOTIF']
        else:
            dct['Statistical Significance'][i] = ['NA']
    return dct

motifstatsDictOpen['Statistical Significance'] = dict()
motifstatsDictClosed['Statistical Significance'] = dict()

motifstatsDictOpen = statSignificance(motifstatsDictOpen)
motifstatsDictClosed = statSignificance(motifstatsDictClosed)
pd.DataFrame(motifstatsDictOpen)
pd.DataFrame(motifstatsDictClosed)
pd.DataFrame(motifstatsDictOpen).to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\motifstats-open1.csv')
pd.DataFrame(motifstatsDictClosed).to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\motifstats-closed1.csv')

###############################################################################
                                #3.2
###############################################################################
# For this question we are supposed to consider only the configuration: A -> B <- C
# This configuration is described in the mFinder manual (motif dictionary) as id36 under 3 node subgraphs
# for this task we have to use an extra parameter as mentioned in the manual: -ospmem (which outputs members of specific subgraph only)
                
mfinderLink = 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\mfinder1.2\\mfinder1.2.exe'
parameters = [mfinderLink, 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_open.txt', '-s', '3', '-r', '200', '-f', 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\id36open', '36', '-ospmem', '36']
p = subprocess.Popen(parameters)   
parameters = [mfinderLink, 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_closed.txt', '-s', '3', '-r', '200', '-f', 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\id36closed','36', '-ospmem', '36',]
p = subprocess.Popen(parameters)   

# eyes open
file36open = open("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\id36open_MEMBERS.txt",'r')
edges36open = []
rowsAll = ''
for rows in file36open:
    rowsAll += rows
    
for line in rowsAll.split('\n'):
    if '\t' in line:
        edges36open += [(int(line.split('\t')[0]), int(line.split('\t')[2]))]
        edges36open += [(int(line.split('\t')[1]), int(line.split('\t')[2]))]
        
# to remove duplicates
edges36open = list(set(edges36open))
len(edges36open)
edges36opena = [["",""]]*len(edges36open)
for i in range(len(edges36open)):
    print(i,labelsD[edges36open[i][0]-1],labelsD[edges36open[i][1]-1])
    edges36opena[i] = [labelsD[edges36open[i][0]-1],labelsD[edges36open[i][1]-1]]

nodesOpen = set()
for i,j in edges36opena:
    nodesOpen.update([i,j])

Gopen36 = nx.DiGraph()
Gopen36 .add_edges_from(edges36opena)
positionOpen = {(x, y) for (x, y) in position_data.items() if x in list(nodesOpen)}
nx.draw_networkx(Gopen36 , dict(positionOpen ))
plt.title('Topological representation for eyes open case of networks with this configuration')
plt.show()

# eyes closed
file36closed = open("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\id36closed36_MEMBERS.txt",'r')
edges36closed = []
rowsAll = ''
for rows in file36closed:
    rowsAll += rows
    
for line in rowsAll.split('\n'):
    if '\t' in line:
        edges36closed += [(int(line.split('\t')[0]), int(line.split('\t')[2]))]
        edges36closed += [(int(line.split('\t')[1]), int(line.split('\t')[2]))]
        
edges36closed
len(edges36closed)
# to remove duplicates
edges36closed = list(set(edges36closed))
len(edges36closed)
edges36closeda = [["",""]]*len(edges36closed)
for i in range(len(edges36closed)):
    print(i,labelsD[edges36closed[i][0]-1],labelsD[edges36closed[i][1]-1])
    edges36closeda[i] = [labelsD[edges36closed[i][0]-1],labelsD[edges36closed[i][1]-1]]

nodesClosed = set()
for i,j in edges36closeda:
    nodesClosed.update([i,j])
Gclosed36 = nx.DiGraph()
Gclosed36 .add_edges_from(edges36opena)
positionClosed = {(x, y) for (x, y) in position_data.items() if x in list(nodesClosed)}
nx.draw_networkx(Gclosed36 , dict(positionClosed ))
plt.title('Topological representation for eyes closed case of networks with this configuration')
plt.show()

###############################################################################
                                #3.3
###############################################################################
# Determining the motifs involved with a selected channel from parieto-occipital scalp region
# we found a research paper mentioning chnnaels in parietooccipital region as (PO; PO7, PO3, O1)
# So, we considered the channel Po3
numericD['Po3'] # the id is 56
'''
with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_openPO.txt', 'w') as f:
    for i, j in directed_graph_open.edges():
        if i == 'Po3' or j == 'Po3':
            f.write("\t".join(map(str, [numericD[i]+1, numericD[j]+1, 1]))+"\n")

with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_closedPO.txt', 'w') as f:
    for i, j in directed_graph_closed.edges():
        if i == 'Po3' or j == 'Po3':
            f.write("\t".join(map(str, [numericD[i]+1, numericD[j]+1, 1]))+"\n")

import subprocess

# we need to pass parameters to mfinder to open it (as in cmd). The sytax is as follows
#[mfinder.exe,input txt file, -s (motif size), size, -r, number of random networks to generate, number, -f, output filename]
mfinderLink = 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\mfinder1.2\\mfinder1.2.exe'
parameters = [mfinderLink, 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_openPO.txt', '-omem','-s', '3', '-r', '200', '-f', 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\outputsMfinderOpenPO']
p = subprocess.Popen(parameters)
parameters = [mfinderLink, 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_closedPO.txt', '-omem', '-s', '3', '-r', '200', '-f', 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\outputsMfinderClosedPO']
p = subprocess.Popen(parameters)                

# motif analysis
motifstatsDictOpenPO = pd.read_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\motifstats-openPO.csv', ).to_dict()
motifstatsDictOpenPO.keys()
#['MOTIF ID', 'Frequency', 'Interval', 'z score', 'p - value']
pd.DataFrame(motifstatsDictOpenPO)

motifstatsDictClosedPO = pd.read_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\motifstats-closedPO.csv', ).to_dict()
motifstatsDictClosedPO.keys()
pd.DataFrame(motifstatsDictClosedPO)

def calculateInterval(dct):
    for i in range(len(dct['Interval'])):
        dct['lowerLimit'][i] = float(dct['Interval'][i].split('+-')[0]) - float(dct['Interval'][i].split('+-')[1])
        dct['upperLimit'][i] = float(dct['Interval'][i].split('+-')[0]) + float(dct['Interval'][i].split('+-')[1])
    return dct
motifstatsDictOpenPO['lowerLimit']=dict()
motifstatsDictOpenPO['upperLimit']=dict()
motifstatsDictClosedPO['lowerLimit']=dict()
motifstatsDictClosedPO['upperLimit']=dict()
motifstatsDictOpenPO = calculateInterval(motifstatsDictOpenPO)
motifstatsDictClosedPO = calculateInterval(motifstatsDictClosedPO)
pd.DataFrame(motifstatsDictOpenPO)
pd.DataFrame(motifstatsDictClosedPO)

# To figure out the statistical significance of each row, we apply the following:
# if frequency of a row > upper limit = motif
# if frequency of a row < lower limit = anti-motif

def statSignificance(dct):
    for i in range(len(dct['MOTIF ID'])):
        if dct['Frequency'][i] > dct['upperLimit'][i]:
            dct['Statistical Significance'][i] = ['MOTIF']
        elif dct['Frequency'][i] < dct['lowerLimit'][i]:
            dct['Statistical Significance'][i] = ['ANTI-MOTIF']
        else:
            dct['Statistical Significance'][i] = ['NA']
    return dct

motifstatsDictOpenPO['Statistical Significance'] = dict()
motifstatsDictClosedPO['Statistical Significance'] = dict()

motifstatsDictOpenPO = statSignificance(motifstatsDictOpenPO)
motifstatsDictClosedPO = statSignificance(motifstatsDictClosedPO)
pd.DataFrame(motifstatsDictOpenPO)
pd.DataFrame(motifstatsDictClosedPO)
'''
############################

# eyes open
filePOopen = open("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\outputsMfinderOpenPO_MEMBERS.txt",'r')
edgesPOopen = []
rowsAll = ''
for rows in filePOopen:
    rowsAll += rows
    
for line in rowsAll.split('\n'):
    if '\t' in line:
        edgesPOopen += [(int(line.split('\t')[0]), int(line.split('\t')[2]))]
        edgesPOopen += [(int(line.split('\t')[1]), int(line.split('\t')[2]))]
        
# to remove duplicates
edgesPOopen = list(set(edgesPOopen))
len(edgesPOopen)
edgesPOopena = [["",""]]*len(edgesPOopen)
for i in range(len(edgesPOopen)):
    print(i,labelsD[edgesPOopen[i][0]-1],labelsD[edgesPOopen[i][1]-1])
    edgesPOopena[i] = [labelsD[edgesPOopen[i][0]-1],labelsD[edgesPOopen[i][1]-1]]

nodesOpen = set()
for i,j in edgesPOopena:
    nodesOpen.update([i,j])

GopenPO = nx.DiGraph()
GopenPO.add_edges_from(edgesPOopena)
positionOpen = {(x, y) for (x, y) in position_data.items() if x in list(nodesOpen)}
nx.draw_networkx(GopenPO , dict(positionOpen ))
plt.title('All the motifs involved with the channel - Po3 - eyes open')
plt.show()
print("The total number of nodes in this network are: ",len(GopenPO.nodes)) # 49
print("The total number of edges in this network are: ",len(GopenPO.edges)) # 473

####
# eyes closed
filePOclosed = open("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\outputsMfinderClosedPO_MEMBERS.txt",'r')
edgesPOclosed = []
rowsAll = ''
for rows in filePOclosed:
    rowsAll += rows
    
for line in rowsAll.split('\n'):
    if '\t' in line:
        edgesPOclosed += [(int(line.split('\t')[0]), int(line.split('\t')[2]))]
        edgesPOclosed += [(int(line.split('\t')[1]), int(line.split('\t')[2]))]
        
# to remove duplicates
edgesPOclosed = list(set(edgesPOclosed))
len(edgesPOclosed)
edgesPOcloseda = [["",""]]*len(edgesPOclosed)
for i in range(len(edgesPOclosed)):
    print(i,labelsD[edgesPOclosed[i][0]-1],labelsD[edgesPOclosed[i][1]-1])
    edgesPOcloseda[i] = [labelsD[edgesPOclosed[i][0]-1],labelsD[edgesPOclosed[i][1]-1]]

nodesClosed = set()
for i,j in edgesPOcloseda:
    nodesClosed.update([i,j])

GclosedPO = nx.DiGraph()
GclosedPO .add_edges_from(edgesPOcloseda)
positionClosed = {(x, y) for (x, y) in position_data.items() if x in list(nodesClosed)}
nx.draw_networkx(GclosedPO , dict(positionClosed ))
plt.title('All the motifs involved with the channel - Po3 - eyes closed')
plt.show()
print("The total number of nodes in this network are: ",len(GclosedPO.nodes)) # 42
print("The total number of edges in this network are: ",len(GclosedPO.edges)) # 333


###############################################################################
                                #3.4
###############################################################################
########### here we use the software mfinder
# motifs analysis for 4-node motifs
with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_open-4node.txt', 'w') as f:
    for i, j in directed_graph_open20.edges():
        f.write("\t".join(map(str, [numericD[i]+1, numericD[j]+1, 1]))+"\n")

with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_closed-4node.txt', 'w') as f:
    for i, j in directed_graph_closed20.edges():
        f.write("\t".join(map(str, [numericD[i]+1, numericD[j]+1, 1]))+"\n")

import subprocess

# we need to pass parameters to mfinder to open it (as in cmd). The sytax is as follows
#[mfinder.exe,input txt file, -s (motif size), size, -r, number of random networks to generate, number, -f, output filename]
mfinderLink = 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\mfinder1.2\\mfinder1.2.exe'
# processing time: 1.06 hours
parameters = [mfinderLink, 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_open-4node.txt', '-omem','-s', '4', '-r', '200', '-f', 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\outputsMfinderOpen4nodes']
p = subprocess.Popen(parameters)
# processing time: 1.08 hours
parameters = [mfinderLink, 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\directed_closed-4node.txt', '-omem', '-s', '4', '-r', '200', '-f', 'C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\outputsMfinderClosed4nodes']
p = subprocess.Popen(parameters)

# from this we get the results in txt files, we have copied those results into csv files, which we will load to display the results

motifstatsDictOpen4nodes = pd.read_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\motifstats-open4nodes.csv', ).to_dict()
motifstatsDictOpen4nodes.keys()
#['MOTIF ID', 'Frequency', 'Interval', 'z score', 'p - value']
pd.DataFrame(motifstatsDictOpen4nodes)

motifstatsDictClosed4nodes = pd.read_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\motifstats-closed4nodes.csv', ).to_dict()
motifstatsDictClosed4nodes.keys()
pd.DataFrame(motifstatsDictClosed4nodes)

def calculateInterval(dct):
    for i in range(len(dct['Interval'])):
        dct['lowerLimit'][i] = float(dct['Interval'][i].split('+-')[0]) - float(dct['Interval'][i].split('+-')[1])
        dct['upperLimit'][i] = float(dct['Interval'][i].split('+-')[0]) + float(dct['Interval'][i].split('+-')[1])
    return dct
motifstatsDictOpen4nodes['lowerLimit']=dict()
motifstatsDictOpen4nodes['upperLimit']=dict()
motifstatsDictClosed4nodes['lowerLimit']=dict()
motifstatsDictClosed4nodes['upperLimit']=dict()
motifstatsDictOpen4nodes = calculateInterval(motifstatsDictOpen4nodes)
motifstatsDictClosed4nodes = calculateInterval(motifstatsDictClosed4nodes)
pd.DataFrame(motifstatsDictOpen4nodes)
pd.DataFrame(motifstatsDictClosed4nodes)

# To figure out the statistical significance of each row, we apply the following:
# if frequency of a row > upper limit = motif
# if frequency of a row < lower limit = anti-motif

def statSignificance(dct):
    for i in range(len(dct['MOTIF ID'])):
        if dct['Frequency'][i] > dct['upperLimit'][i]:
            dct['Statistical Significance'][i] = ['MOTIF']
        elif dct['Frequency'][i] < dct['lowerLimit'][i]:
            dct['Statistical Significance'][i] = ['ANTI-MOTIF']
        else:
            dct['Statistical Significance'][i] = ['NA']
    return dct

motifstatsDictOpen4nodes['Statistical Significance'] = dict()
motifstatsDictClosed4nodes['Statistical Significance'] = dict()

motifstatsDictOpen4nodes = statSignificance(motifstatsDictOpen4nodes)
motifstatsDictClosed4nodes = statSignificance(motifstatsDictClosed4nodes)
pd.DataFrame(motifstatsDictOpen4nodes)
pd.DataFrame(motifstatsDictClosed4nodes)
pd.DataFrame(motifstatsDictOpen4nodes).to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\motifstats-open1-4nodes.csv')
pd.DataFrame(motifstatsDictClosed4nodes).to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\motifstats-closed1-4nodes.csv')




###############################################################################
                                #4.1
###############################################################################

import community

###Louvain Clustering 
# Eyes Open 
com = community.best_partition(directed_graph_open.to_undirected())
#community.generate_dendrogram(directed_graph_open.to_undirected())
clustersOpen = []
tempdf=pd.DataFrame()
tempdf['Channel'] = pd.Series(list(com.keys()))
tempdf['Community'] = pd.Series(list(com.values()))
clustersOpen.append(tempdf[tempdf['Community']==0])
clustersOpen.append(tempdf[tempdf['Community']==1])
clustersOpen.append(tempdf[tempdf['Community']==2])
clustersOpen.append(tempdf[tempdf['Community']==3])

print ('Number of Communities Eyes Open case:'+str(len(set(com.values())))) #4
len(clustersOpen[0]),len(clustersOpen[1]),len(clustersOpen[2]),len(clustersOpen[3])

# Eyes Closed 
com = community.best_partition(directed_graph_closed.to_undirected())
clustersClosed = []
tempdf=pd.DataFrame()
tempdf['Channel'] = pd.Series(list(com.keys()))
tempdf['Community'] = pd.Series(list(com.values()))
clustersClosed.append(tempdf[tempdf['Community']==0])
clustersClosed.append(tempdf[tempdf['Community']==1])
clustersClosed.append(tempdf[tempdf['Community']==2])
print ('Number of Communities Eyes Open case:'+str(len(set(com.values())))) #3
len(clustersClosed[0]),len(clustersClosed[1]),len(clustersClosed[2])

clustersOpen[0].to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\1ClusterOpen.csv')
clustersOpen[1].to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\2ClusterOpen.csv')
clustersOpen[2].to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\3ClusterOpen.csv')
clustersOpen[3].to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\4ClusterOpen.csv')

clustersClosed[0].to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\1ClusterClose.csv')
clustersClosed[1].to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\2ClusterClose.csv')
clustersClosed[2].to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\3ClusterClose.csv')

###############################################################################
                                #4.3
###############################################################################

import igraph

g1 = ig.Graph.TupleList(directed_graph_open.edges(), directed=True)

# Open
#directed
ceb = g1.community_infomap()
ceb.summary()
ceb = g1.community_label_propagation()
ceb.summary()
ceb = g1.community_leading_eigenvector()
ceb.summary()
ceb = g1.community_spinglass()
ceb.summary()
# undirected
g1 = ig.Graph.TupleList(directed_graph_open.edges(), directed=False)
ceb = g1.community_multilevel()
ceb.summary()

# Close
g1 = ig.Graph.TupleList(directed_graph_closed.edges(), directed=True)

#directed
ceb = g1.community_infomap()
ceb.summary()
ceb = g1.community_label_propagation()
ceb.summary()
ceb = g1.community_leading_eigenvector()
ceb.summary()
ceb = g1.community_spinglass()
ceb.summary()
# undirected
g1 = ig.Graph.TupleList(directed_graph_closed.edges(), directed=False)
ceb = g1.community_multilevel()
ceb.summary()






