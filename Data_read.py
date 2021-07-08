# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 10:28:46 2021

@author: 20170426-2
"""
import numpy as np
import pandas as pd

def data_read(filepath):#positive nagetive
    data = pd.read_csv(filepath)
    data = np.array(data).tolist()
    
    return data


def simmat_read(mat):#substrate-substrate similarity matrix;kinase-kinase similarity matrix
    weight = pd.read_csv(mat,header = None,)
    weight = np.array(weight).tolist()
    return weight

def name_read(namelist):
    name = pd.read_csv(namelist,header = None,)
    name = np.array(name).tolist()
    name_list=[]
    for line in name:
        name_list.append(line[0])
    
    return name_list

def read_score_matrix(file):
    score_matrix = pd.read_excel(file,index_col=0)    
    return score_matrix