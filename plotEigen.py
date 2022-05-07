# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

def main():
    EigenVectors = readEigenVectors("EigenVector")
    EigenValues = readEigenValues("EigenValues.txt")
    EigenValuesI = readEigenValues("EigenValuesI.txt")
    EigenVectorsI = readEigenVectorsI()
    EigenVectorsHO2 = readEigenVectors("EigenVector_ho2")
    
    
    
    
    plt.rc('text', usetex=True)
    plotOverLap(EigenVectors,EigenVectorsHO2)
    for i in range(4):
        print(np.linalg.norm(EigenVectorsI[i]),np.dot(EigenVectorsI[i],EigenVectorsI[i+1]))
    plotEigenVectors(EigenVectors,EigenVectorsI,3)
    ##plt.figure(figsize = (8,6))
    ##xx = np.linspace(-10,10,len(EigenVectorsI[0]))
    ##plt.plot(xx,EigenVectorsI[0])
    plt.figure(figsize = (8,6))
    x = np.linspace(-5 ,5,1000)
    y = HO(x)
    plt.plot(x,y)
    for i in range(9):
        plt.axhline(EigenValuesI[i],color='g',linestyle ='-.')
        plt.axhline(EigenValues[i],color = 'r',linestyle = '--')
    plt.ylim([0,12])
    plt.xlabel('$z[\sqrt{m\omega/\hbar}]$',fontsize =14)
    plt.ylabel('$V(z)$ and $E$ $[\hbar \omega]$')
    plt.legend(['$V(z) = 0.5z^2$','$E$-Power iteration','$E$-EigenSolver'])
    
def plotOverLap(EigenVectors,EigenVectorsHO2):
    fig = plt.figure(figsize=(6,8))
    x = np.linspace(-10,10,len(EigenVectors[0]))
    labels=[]
    for i in range(3):
        plt.plot(x,EigenVectors[i])
        plt.plot(x,EigenVectorsHO2[i],'--')
        labels.append('$\psi_{'+ str(i) +'HO}$')
        labels.append('$\psi_{'+ str(i) +'Bump}$')
    plt.legend(labels)
    plt.ylabel('$\psi(z)_i$',fontsize =14)
    plt.xlabel('$z[\sqrt{m\omega/\hbar}]$',fontsize =14)
        
    
def HO(x):
    return 0.5*x**2 + 4*np.exp(-4*x**2)

def plotEigenVectors(EigenVectors,EigenVectorsI,numOfEigenvectors):
    fig,(ax_1,ax_2) = plt.subplots(1,2,sharex = True,sharey = True)
    
    fig.set_figheight(6)
    fig.set_figwidth(10)
    x = np.linspace(-10,10,len(EigenVectors[0]))
    labels = []
    for i in range(numOfEigenvectors):
        ax_1.plot(x,np.array(EigenVectors[i]))
        ax_2.plot(x,(np.array(EigenVectorsI[i])/np.linalg.norm(EigenVectorsI[i])))
        labels.append('$\psi_' +str(i)+ '$')
    ax_1.set_title('EigenSolver-Built in',fontsize =18)
    ax_2.set_title('Iterative Power Method',fontsize =18)
    ax_1.set_xlabel('$z[\sqrt{m\omega/\hbar}]$',fontsize = 14)
    ax_1.set_ylabel('$\psi(z)_i$',fontsize = 14)
    ax_2.set_xlabel('$z[\sqrt{m\omega/\hbar}]$',fontsize = 14)
    fig.legend(labels)
    
    #ax_1.set_xlabel('$z[\sqrt{m\omega/\hbar}]$',fontsize = 14)
    #ax_1.set_ylabel('$\psi(z)_i$',fontsize = 14)
    #plt.text(-7.5,0.2,'Using EigenSolver')
    

def readEigenValues(file_name):

    EigenValues = []
    f = open(file_name)
    for line in f:
        split = line.split(",")
        EigenValues.append(float(split[0]))
    return EigenValues
        
def readEigenVectorsI():
    EigenVectors = []
    EigenVector = []
    f = open("EigenVectorsI.txt","r")
    for line in f:
        if line !="\n":
            EigenVector.append(float(line))
        if not line.strip():
            EigenVectors.append(EigenVector)
            EigenVector = []
    return EigenVectors

def readEigenVectors(file_name):
    file_name = file_name+'.txt'
    EigenVectors = []
    EigenVector = []
    f = open(file_name,"r")
    for line in f:
        split = line.split(",")
        if len(split) == 2: 
            EigenVector.append(float(split[0]))
        if not line.strip():
            EigenVectors.append(EigenVector)
            EigenVector = []
    return EigenVectors
    
if __name__ =="__main__":
    main()