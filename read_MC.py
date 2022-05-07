# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 23:45:30 2019

@author: waterhorse
"""
import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt

def main():
    E = read('E_MC')
    g = readG2('RDF_MC')
    E = np.array(E)
    E_per_part = E
    plt.rc('text', usetex=True)
    PlotEnergy(E)
    mean_g = []
    err_g = []
    for val in g:
        mean_g.append(np.mean(val))
        err_g.append(np.sqrt(np.var(val)/4000))
    PlotRDF(mean_g,np.array(err_g))
#    
    t_b = np.arange(1,800,1)
    s = 20000
    s = int(s)
    mean,var,err = blockAverage(t_b,E[5000:s])
    mean2,var2,err2 = blockAverage(t_b,E**2)
    T = 1.4
    C_V = var[199]/(500*T**2) + 3/2
    print(C_V,mean[500],np.sqrt(err[500]),2*np.sqrt(var[500])*err[500]/(500*T**2))
    plt.figure(figsize =[1.25*6.4,1.25*4.8])
    plt.plot(t_b,np.sqrt(err))
    plt.axhline(y = 1.6,color='g',linestyle = '--')
    plt.axvline(x = 500,color='g',linestyle = '--')
    plt.ylabel('Err$(<E_P > )[\epsilon]$')
    plt.xlabel('Box length')
    
    
    
def blockAverage(t_b,A):
    Number_blocks = len(t_b)
    blockMean = [0]*Number_blocks
    blockVar = [0]*Number_blocks
    blockErr =[0]*Number_blocks
    cnt = 0
    for block_len in t_b:
        Num_of_block = int(np.floor(len(A)/block_len))
        temp = [0]*Num_of_block
        for i in range(Num_of_block):
            block_begin = i*block_len
            block_end = block_begin + block_len
            temp[i] = np.mean(A[block_begin:block_end])
        
        blockMean[cnt] = np.mean(temp)
        blockErr[cnt] = np.var(temp)/(Num_of_block)
        blockVar[cnt] = np.var(temp)
        cnt += 1
        
    return blockMean,blockVar,blockErr
    
def PlotRDF(g,err):
    fig = plt.figure(figsize =[1.25*6.4,1.25*4.8])
    x_g = np.linspace(0,len(g)*0.0075,len(g))
    plt.errorbar(x_g,g,err*20,fmt = 'b--',ecolor = 'r')
    plt.axhline(y=1,color ='g',linestyle = '--')
    plt.ylabel('$g(r)$')
    plt.xlabel('$r/\sigma$')
    plt.legend(['reference: $1$','g(r)'])
    plt.show()
    
def PlotEnergy(E):
    fig = plt.figure(figsize =[1.25*6.4,1.25*4.8])
    x = np.linspace(0,len(E)/(1e4),len(E));    
    plt.plot(x,E)
    plt.text((max(x)-min(x))/2,(max(E)+min(E))/2,'$\Delta = 0.4\sigma$'+ ' ,acceptance ratio '+ '$32.97 \%$')
    plt.ylabel('$E_P[\epsilon]$')
    plt.xlabel('Number of steps $(10^{6})$')
    plt.savefig('text.png')
  
    
def readG2(file_name):
    file_name = file_name +'.txt'
    g  = []
    for i in range(400):
        g.append([])
        
    cnt =0
    c2=1
    frames =1
    f = open(file_name,"r")
    for line in f:
        if line.strip():
            res = line.split(" ")
            g[cnt].append(float(res[0]))
            cnt += 1
        else:
            cnt = 0
            
#    for i in range(400):
#       g[i] /= frames
    return g
    
    
def read(file_name):
    file_name += '.txt'
    f= open(file_name,'r')
    E = []
    for line in f:
        s = line.split(" ")
        E.append(float(s[0]))
    return E

if __name__ =='__main__':
    main()