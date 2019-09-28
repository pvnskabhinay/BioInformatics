# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 09:39:54 2018

@author: pvnsk
"""

import csv
from bioservices import PSICQUIC
import pandas as pd

# intializations
gene_species = []
c=0
j=0
df = pd.DataFrame(columns=('ID','Seed_Gene','GeneA','GeneB','Column1','Column2','Column3'))


# 1st type
# no data:  "LCE3B", "NXPE1",  "OR5B21"
#Genes = ["ADCY3", "CXCR2", "DNMT3B", "FAP", "FOS", "FUT2", "GPR35", "IFIH1", "IL23R", "IL27", "IL2RA", "KEAP1", "LCE3B", "NFKB1", "NXPE1", "OR5B21", "OSMR", "PTPN2", "RAVER1", "RNF186", "RPS6KB1", "SH2B3", "TNFRSF6B", "TYK2"]
Genes = ["ADCY3", "CXCR2", "DNMT3B", "FAP", "FOS", "FUT2", "GPR35", "IFIH1", "IL23R", "IL27", "IL2RA", "KEAP1", "NFKB1","OSMR", "PTPN2", "RAVER1", "RNF186", "RPS6KB1", "SH2B3", "TNFRSF6B", "TYK2"]

s = PSICQUIC(verbose=False)

for gene in Genes:
    gene_species.append((gene+" AND species:9606"))


for gene in gene_species:
    print(gene)
    data = s.query("biogrid", gene)
    for row in data:
        df.loc[c] = [j,Genes[j],row[0],row[1],row[2],row[3],row[4]]
        #print(Genes[j])
        c+=1
        #print(c)
    j+=1
    #print(j)
    df.to_csv('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\BioGrid.csv',encoding='utf8')


# 2nd type
from bioservices import BioGRID
#b = BioGRID(query=["ADCY3", "CXCR2", "DNMT3B", "FAP", "FOS", "FUT2","LCE3B","NXPE1","OR5B21", "GPR35", "IFIH1", "IL23R", "IL27", "IL2RA", "KEAP1", "NFKB1", "OSMR", "PTPN2", "RAVER1", "RNF186", "RPS6KB1", "SH2B3", "TNFRSF6B", "TYK2"
 #                 ],taxId="9606")
# no data: "LCE3B", "NXPE1", "OR5B21"
#Genes = ["ADCY3", "CXCR2", "DNMT3B", "FAP", "FOS", "FUT2", "GPR35", "IFIH1", "IL23R", "IL27", "IL2RA", "KEAP1", "LCE3B", "NFKB1", "NXPE1", "OR5B21", "OSMR", "PTPN2", "RAVER1", "RNF186", "RPS6KB1", "SH2B3", "TNFRSF6B", "TYK2"]
Genes = ["ADCY3", "CXCR2", "DNMT3B", "FAP", "FOS", "FUT2", "GPR35", "IFIH1", "IL23R", "IL27", "IL2RA", "KEAP1", "NFKB1",  "OSMR", "PTPN2", "RAVER1", "RNF186", "RPS6KB1", "SH2B3", "TNFRSF6B", "TYK2"]
len(Genes)
totalgene =[]
for i in Genes :
    print(i)
    b = BioGRID(query= i,taxId="9606")
    totalgene.append(b.biogrid.interactors)

len(b.biogrid.interactors)
len(totalgene)

sum1 = 0
for i in range(0,len(totalgene)):
    sum1 = sum1 + len(totalgene[i])
    
sum1


with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\BioGrid1.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(totalgene)

##### Plotting network graph
import networkx as nx
import matplotlib.pyplot as plt

Genes = ["ADCY3", "CXCR2", "DNMT3B", "FAP", "FOS", "FUT2", "GPR35", "IFIH1", "IL23R", "IL27", "IL2RA", "KEAP1", "LCE3B", "NFKB1", "NXPE1", "OR5B21", "OSMR", "PTPN2", "RAVER1", "RNF186", "RPS6KB1", "SH2B3", "TNFRSF6B", "TYK2"]

#nx.test()
G = nx.Graph()
#pos=nx.spring_layout(G) # positions for all nodes
G.add_nodes_from(Genes)
nx.draw(G, node_color = 'b',with_labels = True)

for i in range(0,len(totalgene)):
    G.add_edges_from(totalgene[i])

color_map = []
for node in G:
    if node in Genes:
        color_map.append('blue')
    else: color_map.append('yellow')

size_map = []
for node in G:
    if node in Genes:
        size_map.append(800)
    else: size_map.append(100)

nx.draw_kamada_kawai(G, node_color = color_map)

nx.draw_random(G,node_color = color_map, node_size = size_map)
nx.draw_kamada_kawai(G, node_color = color_map,  node_size = size_map)

plt.draw() 
G.number_of_edges()
G.number_of_nodes()
list(G.nodes)
list(G.edges)
list(G.adj["ADCY3"])

G1 = nx.Graph()

# To check if non-seed genes are interacting
for gene in list(G.adj):
    if gene not in Genes:
        for i in range(0,len(list(G.adj[str(gene)]))):
            if list(G.adj[str(gene)])[i] not in Genes:
                G1.add_edge(str(gene),list(G.adj[str(gene)])[i])
                 
nx.draw_kamada_kawai(G1)
nx.draw_random(G1,node_color = 'pink')

len(list(G1.nodes))
len(list(G1.edges))

list(G1.nodes)
list(G1.edges)
            
G2 = nx.Graph()
for gene in list(G.adj):
    if gene in Genes:
        for i in range(0,len(list(G.adj[str(gene)]))):
            if list(G.adj[str(gene)])[i] in Genes:
                G2.add_edge(str(gene),list(G.adj[str(gene)])[i])            
            
len(list(G2.nodes))
len(list(G2.edges))

list(G2.nodes)
list(G2.edges)    
            
            
##################### IID
iidData = pd.read_excel('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\IID.xlsx')
# 1st method (considering all the interactions downloaded from the iid website)
G3 = nx.Graph()
edgeList = []
for i in range(0,len(iidData)):
    #tmpList = []
    edgeList.append([iidData['Query Symbol'][i],iidData['Partner Symbol'][i]])
    #edgeList.append(tmpList)
G3.add_edges_from(edgeList)

len(list(G3.nodes))
len(list(G3.edges))

list(G3.nodes)
list(G3.edges)    

color_map = []
for node in G3:
    if node in Genes:
        color_map.append('blue')
    else: color_map.append('yellow')

size_map = []
for node in G3:
    if node in Genes:
        size_map.append(800)
    else: size_map.append(100)


#IID
nx.draw_random(G3,node_color = color_map, node_size = size_map)

# 2nd method: considering only the interactions which have iid (in the sources column)
iidData = pd.read_excel('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\IID.xlsx')
import numpy as np
num1 = 0
sLength = len(iidData['Sources'])
iidData['iid'] = pd.Series(np.random.randn(sLength), index=iidData.index)
for i in range(len(iidData)):
    if "iid" in iidData['Sources'][i]:
#        print(iidData['Sources'][i])
        num1 = num1 + 1
        iidData['iid'][i] = 1
    else:
        iidData['iid'][i] = 0
        
len(iidData[iidData['iid'] == 1])
# making a subset from this data which contain iid in the source
iidOnly = iidData[iidData['iid'] == 1]
iidOnly = iidOnly.reset_index(drop=True)
iidOnly.head()

G4 = nx.Graph()
edgeList = []
for i in range(0,len(iidOnly)):
    #tmpList = []
#    print([iidOnly['Query Symbol'][i],iidOnly['Partner Symbol'][i]])
    print(i)
    edgeList.append([iidOnly['Query Symbol'][i],iidOnly['Partner Symbol'][i]])
    #edgeList.append(tmpList)
G4.add_edges_from(edgeList)

len(list(G4.nodes))
len(list(G4.edges))

list(G4.nodes)
list(G4.edges)    

color_map = []
for node in G4:
    if node in Genes:
        color_map.append('blue')
    else: color_map.append('yellow')

size_map = []
for node in G4:
    if node in Genes:
        size_map.append(800)
    else: size_map.append(100)

nx.draw_random(G3,node_color = color_map, node_size = size_map)

# Writing this into a file
with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\iidOnly.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(list(G4.edges))



F = nx.compose(G,G4)

color_map = []
for node in F:
    if node in Genes:
#        print(node)
        color_map.append('blue')
    else: color_map.append('yellow')

size_map = []
for node in F:
    if node in Genes:
        size_map.append(800)
    else: size_map.append(25)

nx.draw_random(F,node_color = color_map, node_size = size_map)

# To check if non-seed genes are interacting
G1 = nx.Graph()
for gene in list(F.adj):
    if gene not in Genes:
        #print(gene)
        for i in range(0,len(list(F.adj[gene]))):
            if list(F.adj[gene])[i] not in Genes:
                G1.add_edge(str(gene),list(F.adj[str(gene)])[i])

len(G1.nodes)
nx.draw_random(G1,with_labels = True)

with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\nonSeedGeneInteractions.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(G1.edges)
