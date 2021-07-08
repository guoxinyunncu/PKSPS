# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 10:12:23 2021

@author: 20170426-2
"""
import numpy as np
import pandas as pd
import networkx as nx
import random

def MakeGraph(S,subname_list,kinname_list,sweight,kweight):#kinase-substrate network 
    G1 = nx.Graph()
    G1.add_nodes_from(subname_list)
    G1.add_nodes_from(kinname_list)
    G1.add_weighted_edges_from(sweight)
    G1.add_weighted_edges_from(kweight)
    esmall=[(u,v) for (u,v,d) in G1.edges(data=True) if d['weight'] <=0.7]
    G1.remove_edges_from(esmall)
    k_sinter=[]
    for i in range(len(kinname_list)):
        for j in range(len(subname_list)):
            if S[i][j]!=0:
                a=[]
                a.append(kinname_list[i])
                a.append(subname_list[j])
                a.append(S[i][j])
                k_sinter.append(a)
    G1.add_weighted_edges_from(k_sinter)    
    return G1

def requry_graph(S,substrate,subname_list,kinname_list,sweight,kweight):
    A= S.copy()
    b = subname_list.index(substrate)
    c = np.zeros(shape=(len(kinname_list)))
    A[:, b]=c  #The Kinase information of the requry substrate protein was removed
    G1= MakeGraph(A,subname_list,kinname_list,sweight,kweight)
    
    return G1,A
    



def MWBM(sub,kin,G1):# matching algorithm
    Ns = G1[sub]
    Nk = G1[kin]
    Nks=[]
    for i in Ns:
        Nks.append(i)
    for j in Nk:
        Nks.append(j)
    sub_graph = G1.subgraph(Nks)#Neighbors subnet
    maxbig=nx.max_weight_matching(sub_graph, maxcardinality=False, weight='weight')#maximum matching
    wij=0
    for line in  maxbig:
        wij += G1.get_edge_data(line[0],line[1])['weight']
    return wij


def findkin(kinname,data3):#Look for  the substrates of the kinase
    kin_gene =[line[2] for line in data3]
    kin_gene=np.array(kin_gene)
    index = np.argwhere(kin_gene==kinname)
    kinase_set =[]
    for i in range(len(index)):
        kinase_set.append(data3[index[i][0]])
    return kinase_set


def get_match_score(s1,s2,score_matrix):#The values of amino acids S1 and S2 in the substitution matrix
    score = score_matrix[s1][s2]
    return score


def symmetric_normalization(kinname_list,subname_list,S):#normalization Incidence matrix S 
    sums=S.sum(axis=0)
    sumk=S.sum(axis=1)
    for i in range(len(kinname_list)):
        for j in range(len(subname_list)):
            if (sumk[i]==0 and sums[j]!=0):
                S[i][j]=((S[i][j]/sums[j]))/2
            if (sums[j]==0 and sumk[i]!=0):
                S[i][j]=((S[i][j]/sumk[i]))/2
            if (sums[j]==0 and sumk[i]==0):
                S[i][j]=0
            if (sums[j]!=0 and sumk[i]!=0):
                S[i][j]=((S[i][j]/sums[j])+(S[i][j]/sumk[i]))/2                    
    return S

def location_sort(list_in): #Sort by value
    max_location =sorted(enumerate(list_in), key=lambda y:y[1],reverse=True)
    list_out = []
    for i in range(len(max_location)):
        list_out.append(max_location[i][0])#output location list
 
    return (list_out)

def PKSPS_NET(substrate,A,subname_list,kinname_list,G1,sweight,kweight):
    for kin in kinname_list:#   
        b = subname_list.index(substrate)
        a = kinname_list.index(kin)
        A[a][b]=MWBM(substrate,kin,G1)
    sub_pre = A[:,b]
    sort = location_sort(sub_pre)    
    sort_kin = [kinname_list[line] for line in sort]
    sort10 =sort_kin[0:10]
    P_net10 = []
    for i in range(len(sort10)):
        value = []
        value.append(sort10[i])
        value.append((i+1)/317)
        P_net10.append(value)
        
    return P_net10, sort10



def PKSPS_seq(substrate,query,sort10,data3,data5,score_matrix):
    P = []
    for kin in sort10:
        sp = [ line[4] for line in findkin(kin,data3) if line[1]!= substrate ]
        sb1 = [seq[4] for seq in data3 if seq not in sp and seq[1]!= substrate]
        sub = [line[0] for line in findkin(kin,data3) if line[1]!= substrate]
        sb2 = [seq[4] for seq in data5 if seq[0] in sub]
        if len(sb2)>=len(sb1):
            sb = sb1 + random.sample(sb2, len(sb1))
        else:
            sb = random.sample(sb1, len(sb2)) + sb2
        rsp=len(sp)
        rsn=len(sb)
        comparison = sb.copy()
        dist_list=[]       
        for q in range(len(comparison)):#计算查询蛋白和背景集的相似度
            sequence1 = query
            sequence2 = comparison[q]
            dist=0
            for i in range(len(sequence2)):
                match = get_match_score(sequence2[i].upper(),sequence1[i].upper(),score_matrix)
                sim = match-(-4)/11-(-4)
                dist = dist + sim
            dist = dist/len(sequence2)
            dist_list.append(dist)
        sort_dist = location_sort(dist_list)
        D = sum(dist_list[:rsp])
        runsum = []
        sumnum=0
        for line in sort_dist:
            if line<=rsp:
                sumnum=sumnum+(dist_list[line])/D
            else:
                sumnum= sumnum+(-1/rsn)
            runsum.append(sumnum)    
        esp = max(runsum)    
        ES=[]
        ES.append(esp)
        for i in range(999):
            R=random.sample(sort_dist, rsp)
            sumnum=0
            runsum = []
            D=0
            for line in R:
                D=D+dist_list[line]
            for line in sort_dist:
                if line in R:
                    sumnum=sumnum+(dist_list[line])/D
                else:
                    sumnum= sumnum+(-1/rsn)
                runsum.append(sumnum)    
            es = max(runsum)
            ES.append(es)
        sort_ES = location_sort(ES)
        p =(sort_ES.index(0)+1)/1000    
        P.append(p)
    return P

def PKSPS(l,P_net10,P_seq):
    value = P_net10.copy()
    for i in range(len(P_net10)):
        value[i].append(P_seq[i])
        p_combine = l*P_net10[i][1]+(1-l)*P_seq[i]
        value[i].append(p_combine)
    return value
        
    
