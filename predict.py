# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 10:22:59 2021

@author: 20170426-2
"""

import Data_read 
import PKSPS
import numpy as np

def main():
    print('PKSPS: a novel method for predicting kinase-specific phosphorylation sites combined with protein-protein interaction information and local sequence information.')
    substrate =input("Please enter the substrate protein for inquiry:")
    query = input("Please enter a sequence of queries:")
    threshold = input("Please select a threshold for the output:")
    threshold = float(threshold)
    data = Data_read.data_read('data_all.csv')
    data_neg = Data_read.data_read('data_neg.csv')
    kinname_list = Data_read.name_read('kinname_list.txt')
    subname_list = Data_read.name_read('subname_list.txt')
    sweight = Data_read.simmat_read('sweight.csv')
    kweight = Data_read.simmat_read('kweight.csv')
    score_matrix = Data_read.read_score_matrix('blosum62.xlsx')
    
    S= np.zeros(shape=(len(kinname_list),len(subname_list)))#连接矩阵
    for line in data:
        if line[1] in kinname_list:
            line[1] = line[1]+'S'
        i = kinname_list.index(line[2])
        j = subname_list.index(line[1])
        S[i][j] = S[i][j]+1
    S=PKSPS.symmetric_normalization(kinname_list,subname_list,S)
    
    G1,A = PKSPS.requry_graph(S,substrate,subname_list,kinname_list,sweight,kweight)
    P_net10, sort10 = PKSPS.PKSPS_NET(substrate,A,subname_list,kinname_list,G1,sweight,kweight)
    P_seq = PKSPS.PKSPS_seq(substrate,query,sort10,data,data_neg,score_matrix)
    P_combine= PKSPS.PKSPS(0.3,P_net10,P_seq)
    result=[]
    for i in range(len(P_combine)):
        re = []
        re.append(P_combine[i][0])
        re.append(P_combine[i][3])
        result.append(re)   
    output=[line for line in result if line[1] < threshold]
    print(output)
    
        
    
if __name__ == "__main__" :
    main()    