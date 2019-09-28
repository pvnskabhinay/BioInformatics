# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 11:51:49 2019

@author: pvnsk
"""
# Bioinformatics NetBio-2
import pandas as pd
import networkx as nx

SGI = pd.read_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\4-1-SeedGeneInteractome.csv")
I = pd.read_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\4-3-intersectionInteractome.csv")
U = pd.read_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\4-2-UnionInteractome.csv")
Genes = ["ADCY3", "CXCR2", "DNMT3B", "FAP", "FOS", "FUT2", "GPR35", "IFIH1", "IL23R", "IL27", "IL2RA", "KEAP1", "LCE3B", "NFKB1", "NXPE1", "OR5B21", "OSMR", "PTPN2", "RAVER1", "RNF186", "RPS6KB1", "SH2B3", "TNFRSF6B", "TYK2"]

SGI.head()
I.head()
U.head()

# Seed Gene Interactome
nodesSGI = set()
tmp = set(SGI['Interactor A Gene Symbol'])
tmp1 = set(SGI['Interactor B Gene Symbol'])
nodesSGI = tmp.union(tmp1)
edgesSGI = []
for i in range(len(SGI)):
    edgesSGI.append([SGI['Interactor A Gene Symbol'][i],SGI['Interactor B Gene Symbol'][i]])
G_SGI = nx.Graph()
G_SGI.add_edges_from(edgesSGI)
nx.draw_random(G_SGI,with_labels = True)
nx.draw_kamada_kawai(G_SGI,with_labels = True)

# Intersection Interactome
nodesI = set()
tmp = set(I['Interactor A Gene Symbol'])
tmp1 = set(I['Interactor B Gene Symbol'])
nodesI = tmp.union(tmp1)
edgesI = []
for i in range(len(I)):
    edgesI.append([I['Interactor A Gene Symbol'][i],I['Interactor B Gene Symbol'][i]])
G_I = nx.Graph()
G_I.add_edges_from(edgesI)
nx.draw_random(G_I,with_labels = True)
nx.draw_kamada_kawai(G_I,with_labels = True)

# Union Interactome
nodesU = set()
tmp = set(U['Interactor A Gene Symbol'])
tmp1 = set(U['Interactor B Gene Symbol'])
nodesU = tmp.union(tmp1)
edgesU = []
for i in range(len(U)):
    edgesU.append([U['Interactor A Gene Symbol'][i],U['Interactor B Gene Symbol'][i]])
G_U = nx.Graph()
G_U.add_edges_from(edgesU)
#nx.draw_kamada_kawai(G_U,node_size = 5)
nx.draw_random(G_U, node_size = 5)

'''
1. Calculate the main network measures for SGI, I and U
'''

#1.1 Global measures of SGI, U and I
#A. Number of nodes and links

# SGI
print("The number of nodes for SGI graph is: ", len(G_SGI.nodes))
print("The number of links for SGI graph is: ", len(G_SGI.edges))
# Intersection interactome
print("The number of nodes for I graph is: ", len(G_I.nodes))
print("The number of links for I graph is: ", len(G_I.edges))
# Union interactome
print("The number of nodes for U graph is: ", len(G_U.nodes))
print("The number of links for U graph is: ", len(G_U.edges))

# B. Number of Connected components
print("The number of connected components for SGI graph is: ", nx.number_connected_components(G_SGI))
print("The number of connected components for I graph is: ", nx.number_connected_components(G_I))
print("The number of connected components for U graph is: ", nx.number_connected_components(G_U))

# C. Number of isolated nodes
print("The number of isolated nodes for SGI graph is: ", nx.number_of_isolates(G_SGI))
print("The number of isolated nodes for I graph is: ", nx.number_of_isolates(G_I))
print("The number of isolated nodes for U graph is: ", nx.number_of_isolates(G_U))

# D. Average Path length

# SGI
PathLenSGI=[]
for g in nx.connected_component_subgraphs(G_SGI):
    PathLenSGI.append(nx.average_shortest_path_length(g)) 
print("The average path length for SGI graph is: ", PathLenSGI)
print(sum(PathLenSGI)/len(PathLenSGI))
# I
PathLenI=[]
for g in nx.connected_component_subgraphs(G_I):
    PathLenI.append(nx.average_shortest_path_length(g)) 
print("The average path length for I graph is: ", PathLenI)
print(sum(PathLenI)/len(PathLenI))
# U
PathLenU=[]
for g in nx.connected_component_subgraphs(G_U):
    PathLenU.append(nx.average_shortest_path_length(g)) 
print("The average path length for U graph is: ", PathLenU)
print(sum(PathLenU)/len(PathLenU))

# E. Average Degree
nx.average_degree_connectivity(G_SGI)
degreeSGI = list(G_SGI.degree())
tmp = 0
for i in range(len(degreeSGI)):
    tmp = tmp+degreeSGI[i][1]
print("Average degree from SGI graph is: ", tmp/len(degreeSGI))

nx.average_degree_connectivity(G_I)
degreeI = list(G_I.degree())
tmp = 0
for i in range(len(degreeI)):
    tmp = tmp+degreeI[i][1]
print("Average degree from SGI graph is: ", tmp/len(degreeI))

nx.average_degree_connectivity(G_U)
degreeU = list(G_U.degree())
tmp = 0
for i in range(len(degreeU)):
    tmp = tmp+degreeU[i][1]
print("Average degree from SGI graph is: ", tmp/len(degreeU))

# F. Average Clustering coefficient

print("The Average Clustering Coefficient of SGI Graph is: ", nx.average_clustering(G_SGI))
print("The Average Clustering Coefficient of I Graph is: ", nx.average_clustering(G_I))
print("The Average Clustering Coefficient of U Graph is: ", nx.average_clustering(G_U))

# G. Network diameter and radius
# As the graphs are not connected, we consider all connected subgraphs

DiaSGI, RadSGI = [], []
g1 = [] # all the subgraphs
for g in nx.connected_component_subgraphs(G_SGI):
    d=nx.diameter(g)
    r=nx.radius(g)
    DiaSGI.append(d)
    RadSGI.append(r)
    g1.append(g)
    if d != 0:
        print("diameter: ", d)
    if r != 0:
        print("radius: ", r)

DiaI, RadI = [], []
g2 = [] # all the subgraphs
for g in nx.connected_component_subgraphs(G_I):
    d=nx.diameter(g)
    r=nx.radius(g)
    DiaI.append(d)
    RadI.append(r)
    g2.append(g)
    if d != 0:
        print("diameter: ", d)
    if r != 0:
        print("radius: ", r)

DiaU, RadU = [], []
g3 = [] # all the subgraphs
for g in nx.connected_component_subgraphs(G_U):
    d=nx.diameter(g)
    r=nx.radius(g)
    DiaU.append(d)
    RadU.append(r)
    g3.append(g)
    if d != 0:
        print("diameter: ", d)
    if r != 0:
        print("radius: ", r)

# H. Centralization
# https://stackoverflow.com/a/35248202

N_SGI = G_SGI.order()
centralitySGI = nx.degree_centrality(G_SGI).values()
maxCSGI = max(centralitySGI)
CentralizationSGI = ((N_SGI*maxCSGI)-sum(centralitySGI))/(N_SGI - 1)**2

N_I = G_I.order()
centralityI = nx.degree_centrality(G_I).values()
maxCI = max(centralityI)
CentralizationI = ((N_I*maxCI)-sum(centralityI))/(N_I - 1)**2

N_U = G_U.order()
centralityU = nx.degree_centrality(G_U).values()
maxCU = max(centralityU)
CentralizationU = ((N_U*maxCU)-sum(centralityU))/(N_U - 1)**2

print("The centralization value for SGI graph is: ", CentralizationSGI)
print("The centralization value for I graph is: ", CentralizationI)
print("The centralization value for U graph is: ", CentralizationU)

#1.2 Global and local measures for Largest Connected Component (LCC) 

# Isolating the LCC of I and U, (we also did SGI for comparison)
LCC_I = max(nx.connected_component_subgraphs(G_I), key=len)
LCC_U = max(nx.connected_component_subgraphs(G_U), key=len)
LCC_SGI = max(nx.connected_component_subgraphs(G_SGI), key=len)

#A. Global measures

def global_measures(G):
    nodesG=len(G.nodes())
    edgesG=len(G.edges())
    PathLen=nx.average_shortest_path_length(G)
    AvgDeg=nx.average_degree_connectivity(G)
    AvgCC=nx.average_clustering(G)
    Dia=nx.diameter(G)
    Rad=nx.radius(G)
    N = G.order()
    centrality = nx.degree_centrality(G).values()
    maxC = max(centrality)
    Centralization = ((N*maxC)-sum(centrality))/(N - 1)**2
    deg1 = list(G.degree())
    tmp = 0
    for i in range(len(deg1)):
        tmp = tmp+deg1[i][1]
    deg = tmp/len(deg1)
    colnames = ['Nodes','Links','Avg Path Length','Avg Degree','Avg Degree Val','Avg Clustering Coeff','Network Diameter','Network Radius','Centralization']
    Gmeasures = pd.DataFrame(columns = colnames)
    row1 = [nodesG, edgesG, PathLen, AvgDeg, deg, AvgCC, Dia, Rad, Centralization]    
    Gmeasures.loc[len(Gmeasures)] = row1
    return Gmeasures

global_I = global_measures(LCC_I)
global_U = global_measures(LCC_U)
global_SGI = global_measures(LCC_SGI)

global_I.transpose().to_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\1\\global_I.csv")
global_U.transpose().to_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\1\\global_U.csv")
global_SGI.transpose().to_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\1\\global_SGI.csv")

print("The global indices for I graph are: ", global_I.transpose())
print("The global indices for U graph are: ", global_U.transpose())
print("The global indices for SGI graph are: ", global_SGI.transpose())
    
# B. Local measures
def local_measures(G):
    degC = nx.degree_centrality(G)
    betweenC = nx.betweenness_centrality(G)
    eigenVectC= nx.eigenvector_centrality_numpy(G)
    closenessC = nx.closeness_centrality(G)
    ratio = {k: betweenC[k]/degC[k] for k in degC}
    colnames = ['Node Degree Centrality','Betweenness Centrality','EigenVector Centrality','Closeness Centrality','Ratio']
    Lmeasures = pd.DataFrame(columns = colnames)
    row1 = [degC, betweenC, eigenVectC, closenessC, ratio]
    Lmeasures.loc[len(Lmeasures)] = row1
    return Lmeasures    

local_I = local_measures(LCC_I)
local_U = local_measures(LCC_U)
local_SGI = local_measures(LCC_SGI)

print("The local indices for I graph are: ", local_I.transpose())
print("The local indices for U graph are: ", local_U.transpose())
print("The local indices for SGI graph are: ", local_SGI.transpose())

local_I.transpose().to_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\1\\local_I.csv")
local_U.transpose().to_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\1\\local_U.csv")
local_SGI.transpose().to_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\1\\local_SGI.csv")

# Plotting the graphs of LCC_I, LCC_U and LCC_SGI
nx.draw_kamada_kawai(LCC_I,with_labels = True)
nx.draw_random(LCC_I, with_labels = True)

nx.draw_kamada_kawai(LCC_U,node_size = 5)
nx.draw_random(LCC_U, node_size = 5)

nx.draw_kamada_kawai(LCC_SGI,with_labels = True)
nx.draw_random(LCC_SGI, with_labels = True)
'''
2. Apply clustering methods for disease modules discovery
'''
# MCL
# we use the library https://github.com/guyallard/markov_clustering
import markov_clustering as mc

def markov_clustering(G):
    matrix = nx.to_scipy_sparse_matrix(G)
    result = mc.run_mcl(matrix)
    clusters = mc.get_clusters(result)  
    mc.draw_graph(matrix, clusters,  node_size=50, with_labels=False, edge_color="silver")
    return clusters

mcI = markov_clustering(LCC_I)
mcU = markov_clustering(LCC_U)
numericI = {key: idx for idx, key in enumerate(LCC_I.nodes)}
labelsI = dict(enumerate(LCC_I.nodes))
valI = {}
numericU = {key: idx for idx, key in enumerate(LCC_U.nodes)}
labelsU = dict(enumerate(LCC_U.nodes))

markovDict = {"mcI": mcI, "mcU": mcU}

mcI_1 = [[] for i in range(len(mcI))]
for i in range(len(mcI)):
    tmp = []
    for j in range(len(mcI[i])):
        tmp.append(labelsI[mcI[i][j]])
    mcI_1[i] = tmp


mcU_1 = [[] for i in range(len(mcU))]
for i in range(len(mcU)):
    tmp = []
    for j in range(len(mcU[i])):
        tmp.append(labelsU[mcU[i][j]])
    mcU_1[i] = tmp

markovDict_1 = {"mcI": mcI_1, "mcU": mcU_1}
        
# Louvain 
# we first convert this nx graph to igraph, as the inbuilt function community Louvain
# works only on igraph
import igraph
import community
#pip install cairocffi
#import cairocffi
#import cairo
import louvain
import igraph as ig
import matplotlib.pyplot as plt
igLCC_I = igraph.Graph.TupleList(LCC_I.edges(), directed=True)
igLCC_U = igraph.Graph.TupleList(LCC_U.edges(), directed=True)
partitionI = louvain.find_partition(igLCC_I, louvain.CPMVertexPartition,resolution_parameter = 0.05)
ig.plot(partitionI, with_labels=True)
partitionU = louvain.find_partition(igLCC_U, louvain.CPMVertexPartition,resolution_parameter = 0.05)
ig.plot(partitionU)
partitionIa = [() for i in range(len(partitionI))]
for i in range(len(partitionI)):
    partitionIa[i] = tuple(partitionI[i])

partitionUa = [() for i in range(len(partitionU))]
for i in range(len(partitionU)):
    partitionUa[i] = tuple(partitionU[i])

louvainDict = {"lcI": list(partitionIa), "lcU": list(partitionUa)}
#louvainDict.keys()
#louvainDict['lcI']

#https://python-louvain.readthedocs.io/en/latest/
# I
partition = community.best_partition(LCC_I)

size = float(len(set(partition.values())))
pos = nx.spring_layout(LCC_I)
count = 0.
for com in set(partition.values()) :
    count = count + 1.
    list_nodes = [nodes for nodes in partition.keys()
                                if partition[nodes] == com]
    nx.draw_networkx_nodes(LCC_I, pos, list_nodes, node_size = 20)


nx.draw_networkx_edges(LCC_I, pos, alpha=0.5)
plt.show()


# U
partition = community.best_partition(LCC_U)

size = float(len(set(partition.values())))
pos = nx.spring_layout(LCC_U)
count = 0.
for com in set(partition.values()) :
    count = count + 1.
    list_nodes = [nodes for nodes in partition.keys()
                                if partition[nodes] == com]
    nx.draw_networkx_nodes(LCC_U, pos, list_nodes, node_size = 20)


nx.draw_networkx_edges(LCC_U, pos, alpha=0.5)
plt.show()

partitionI = community.best_partition(LCC_I)
partitionU = community.best_partition(LCC_U)

commsI = set(partitionI.values())
tmp0, tmp1, tmp2 = [], [], []
for key in partitionI.keys():    
    if partitionI[key] == 0:
        tmp0.append(key)
    elif partitionI[key] == 1:
        tmp1.append(key)
    else:
        tmp2.append(key)

lcI = [tmp0, tmp1, tmp2]

commsU = set(partitionU.values())
tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12 = [], [], [], [], [], [], [], [], [], [], [], [], []
for key in partitionU.keys():    
    if partitionU[key] == 0:
        tmp0.append(key)
    elif partitionU[key] == 1:
        tmp1.append(key)
    elif partitionU[key] == 2:
        tmp2.append(key)
    elif partitionU[key] == 3:
        tmp3.append(key)
    elif partitionU[key] == 4:
        tmp4.append(key)
    elif partitionU[key] == 5:
        tmp5.append(key)
    elif partitionU[key] == 6:
        tmp6.append(key)
    elif partitionU[key] == 7:
        tmp7.append(key)
    elif partitionU[key] == 8:
        tmp8.append(key)
    elif partitionU[key] == 9:
        tmp9.append(key)
    elif partitionU[key] == 10:
        tmp10.append(key)
    elif partitionU[key] == 11:
        tmp11.append(key)
    elif partitionU[key] == 12:
        tmp12.append(key)

lcU = [tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12 ]

louvainDict = {"lcI": lcI, "lcU": lcU}
markovDict_1.keys()
louvainDict.keys()

####################
from scipy.stats import hypergeom

M, N = 0, 0
d = {}
mid = 0
for key in markovDict_1.keys():
    for cluster in markovDict_1[key]:
        mid = mid+1
        M = M + len(cluster)
#        print(cluster)
        for i in range(len(cluster)):
            c = cluster[i]
            if c in Genes:
                N = N + 1
        
        for i in range(len(cluster)):
            n = len(cluster[i])
            k = 0
            sg = []
            if cluster[i] in Genes:
#                print(i, cluster[i])
                k = k+1
                sg.append(cluster[i])
                pval = 1-hypergeom.cdf(k,M,n,N)
                d[i] = [sg,cluster,cluster[i],pval,key,mid,k,n] # SeedGene List (if >1), Module, SeedGene, pVal, algorithm, module ID, Num of seed genes, num of genes in the module
        

M, N = 0, 0
d1 = {}
mid = 0
for key in louvainDict.keys():
    for cluster in louvainDict[key]:
        mid = mid+1
        M = M + len(cluster)
        for i in range(len(cluster)):
            c = cluster[i]
            if c in Genes:
                N = N + 1
        for i in range(len(cluster)):
            n = len(cluster[i])
            k = 0
            sg = []
            if cluster[i] in Genes:
                k = k+1
                sg.append(cluster[i])
                pval = 1-hypergeom.cdf(k,M,n,N)
                d1[i] = [sg,cluster,cluster[i],pval,key,mid,k,n] # SeedGene List (if >1), Module, SeedGene, pVal, algorithm, module ID, Num of seed genes, num of genes in the module
                
M, N # 3312, 26

p = 0.05
pdmMarkov = {}
pdmLouvain = {}
for i in d:
    if d[i][3] < p:
       pdmMarkov[i] = d[i] 

for i in d1:
    if d1[i][3] < p:
       pdmLouvain[i] = d1[i] 

pdmMarkovDF = pd.DataFrame(columns = ["SeedGeneList","ModuleGenes","SeedGene","p-value","AlgorithmUsed","module ID","Num of SeedGenes","Num of Genes in Module"])
pdmLouvainDF = pd.DataFrame(columns = ["SeedGeneList","ModuleGenes","SeedGene","p-value","AlgorithmUsed","module ID","Num of SeedGenes","Num of Genes in Module"])

len(d1)
for i in pdmMarkov.keys():
    pdmMarkovDF.loc[i] = pdmMarkov[i]

for i in pdmLouvain.keys():
    pdmLouvainDF.loc[i] = pdmLouvain[i]

pdmLouvainDF.head()
pdmLouvainDF.to_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\2\\pdmLouvain.csv")
pdmMarkovDF.to_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\2\\pdmMarkov.csv")

pdmU = pd.DataFrame(columns = ["SeedGeneList","ModuleGenes","SeedGene","p-value","AlgorithmUsed","module ID","Num of SeedGenes","Num of Genes in Module"])
pdmI = pd.DataFrame(columns = ["SeedGeneList","ModuleGenes","SeedGene","p-value","AlgorithmUsed","module ID","Num of SeedGenes","Num of Genes in Module"])

for i in pdmMarkov.keys():
    if pdmMarkov[i][4] == 'mcU':    
        pdmU.loc[len(pdmU)] = pdmMarkov[i]
    else:
        pdmI.loc[len(pdmI)] = pdmMarkov[i] 


for i in pdmLouvain.keys():
    if pdmLouvain[i][4] == 'lcU':    
        pdmU.loc[len(pdmU)] = pdmLouvain[i]
    else:
        pdmI.loc[len(pdmI)] = pdmLouvain[i] 

pdmU.to_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\2\\pdmU.csv")
pdmI.to_csv("C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\2\\pdmI.csv")

'''
3. Enrichment Analysis on the disease modules
'''

GeneListI = set()
for i in range(len(pdmI)):
    for j in range(len(pdmI['ModuleGenes'][i])):
        GeneListI.add(pdmI['ModuleGenes'][i][j])

GeneListU = set()
for i in range(len(pdmU)):
    for j in range(len(pdmU['ModuleGenes'][i])):
        GeneListU.add(pdmU['ModuleGenes'][i][j])

print(len(GeneListI),len(GeneListU))
import csv
with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\3\\GeneListI.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(list(GeneListI))

with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 3\\3\\GeneListU.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(list(GeneListU))

# Save these 2 csv files as xls in the same folder manually, and using the file uniprot-list.xlsx from the 1st homework
# we use basic excel functions index and match to get the corresponding UNIPROT IDs of these genes and use the links: 
# (http://www.innatedb.com/redirect.do?go=batchPw and http://www.innatedb.com/redirect.do?go=batchGo)
# The Pathway and GO analysis for Intersection and Union interactomes are saved as PDM-I-GO.xls, PDM-I-Pathways.xls, PDM-U-GO.xls, and PDM-U-Pathways.xls

'''
4. Finding Putative Disease protiens using the DIAMOnD tool
'''

#python .\DIAMOnD.py .\biogridPPI.txt .\seed_genes.txt 200 .\out.txt
# we run this line in the command prompt to use the DIAMOnD tool. The output is list of 200 DIAMOnD nodes.
# then we save these nodes in an excel file to find out the Gene names using basic excel functions using
# the list of genes and their corresponding numbers that we allocated. And then again using excel functions
# we find the UNIPROT ids of those genes and save them in an xls file and then use this file to get the 
# pathways and GO analysis.
# The final GO analysis output is saved in files: PDM-DIAMOnD-Pathways.xls and PDM-DIAMOnD-GO.xls







