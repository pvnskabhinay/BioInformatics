# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 13:33:57 2018

@author: pvnsk
"""
import pandas as pd
import csv

combinedData = pd.read_excel('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\BioGRIDandIIDcombined.xlsx')
combinedData1 = pd.read_excel('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\BioGRIDandIIDcombined1.xlsx')

Genes = ["ADCY3", "CXCR2", "DNMT3B", "FAP", "FOS", "FUT2", "GPR35", "IFIH1", "IL23R", "IL27", "IL2RA", "KEAP1", "LCE3B", "NFKB1", "NXPE1", "OR5B21", "OSMR", "PTPN2", "RAVER1", "RNF186", "RPS6KB1", "SH2B3", "TNFRSF6B", "TYK2"]

#4.1 Seed Genes Interactome
tmpList = []
for i in range(0, len(combinedData)):
    if combinedData['Gene A'][i] in Genes:
        if combinedData['Gene B'][i] in Genes:
#            print(combinedData['Database'][i])
            tmpList.append([combinedData['Gene A'][i],combinedData['Gene B'][i],combinedData['Database'][i]])
            
len(tmpList)
with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\4-1-SeedGeneInteractome.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(tmpList)

# IID only

tmpList = []
for i in range(0, len(combinedData1)):
    if combinedData1['Gene A'][i] in Genes:
        if combinedData1['Gene B'][i] in Genes:
#            print(combinedData['Database'][i])
            tmpList.append([combinedData1['Gene A'][i],combinedData1['Gene B'][i],combinedData1['Database'][i]])
            
len(tmpList)
with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\4-1-SeedGeneInteractome1.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(tmpList)
     
# 4.2 Union Interactome
tmpList = []
for i in range(0, len(combinedData)):
    if (combinedData['Gene A'][i] in Genes) or (combinedData['Gene B'][i] in Genes):
#        if combinedData['Gene B'][i] in Genes:
#            print(combinedData['Database'][i])
            tmpList.append([combinedData['Gene A'][i],combinedData['Gene B'][i],combinedData['Database'][i]])

len(tmpList)
with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\4-2-UnionInteractome.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(tmpList)

# IID only

tmpList = []
for i in range(0, len(combinedData1)):
    if (combinedData1['Gene A'][i] in Genes) or (combinedData1['Gene B'][i] in Genes):
#        if combinedData['Gene B'][i] in Genes:
#            print(combinedData['Database'][i])
            tmpList.append([combinedData1['Gene A'][i],combinedData1['Gene B'][i],combinedData1['Database'][i]])

len(tmpList)
with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\4-2-UnionInteractome1.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(tmpList)




# 4.3 Intersection Interactome
# Using the data from the file with only IID
tmpList1 = []
tmpList2 = []
for i in range(0,len(tmpList)):
    for j in range(0,len(tmpList)):
        if [tmpList[i][0],tmpList[i][1]] == [tmpList[j][0],tmpList[j][1]] or [tmpList[i][0],tmpList[i][1]] == [tmpList[j][1],tmpList[j][0]]:
            if tmpList[i][2] != tmpList[j][2]:
                tmpList1.append([tmpList[i][0],tmpList[i][1]])
                tmpList2.append([tmpList[j][0],tmpList[j][1]])
    
len(tmpList1)
len(tmpList2)

len(list(G.edges))
len(list(G4.edges))
lst1 = list(G.edges)
lst2 = list(G4.edges)

 
tmpList1 = []
for i in range(0,len(list(G.edges))):
    for j in range(0,len(list(G4.edges))):
        
        if [list(G.edges)[i][0],list(G.edges)[i][1]] == [list(G4.edges)[j][0],list(G4.edges)[j][1]]:
            print(i,j)
            tmpList1.append([list(G.edges)[i][0],list(G.edges)[i][1]])

'''tmpList1 = []
for i in range(0,len(lst1)):
    for j in range(0,len(lst2)):
        
        if [lst1[i][0],lst1[i][1]] == [lst2[j][0],lst2[j][1]]:
            print(i,j)
            tmpList1.append([lst1[i][0],lst1[i][1]])
'''

len(tmpList1)

with open('C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\4-3-intersectionInteractome.csv', 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(tmpList1)

    