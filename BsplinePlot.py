# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 20:21:06 2019

@author: waterhorse
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

def main():
    bspline = readSplines()
    phi = readPhi("Phi.txt")
    phi2 = readPhi("Phi2.txt")
    eigenPhi = readEigen2()
    IonEnergy = readIonEnergy()
    #Phi = readPhi("Phi1.txt")
    
#    for val in IonEnergy:
#        print("test! ",val)
    

    #plt.ylim([0,1])
    
    plotIonEnergy(IonEnergy)
    xx = np.linspace(0.01,10.5,len(bspline[0]))
    plt.figure(figsize = (10,5))
    test = [0,0.875,1.01,1.75,2.625,3.5,4.375,5.25,6.125,7,7.875,8.75,9.625,10.5]
#    for i in range(100  -4):
#        plt.plot(xx,bspline[i])
    #for i in range(len(test)):
    #    plt.axvline(x=test[i],linestyle='--',color = 'g')
   # plt.title("B-splines for $N = 22$ knot points  with $k=4$",fontsize = 18)
   # plt.ylabel('$B^{k}_n(r)$',fontsize = 14)
   # plt.xlabel('$r$[arb. units]',fontsize = 14)
   # plt.figure(figsize = (10,5))
    #x = np.logspace(-2,3,len(phi))
    x = np.logspace(-7,3,len(eigenPhi[0]))
    h = (x[2]-x[1])/len(x)
    #N = h*sum(Phi) 
    # - min(np.array(phi)/x)
    #plt.axvline(x=1)
    
   # plt.plot(x,(2/(4*np.pi))*(np.array(phi)/x)**2)
    labels = []
    cnt = 0
    for phi_e in eigenPhi: 
        plt.plot(x,4*np.pi*(x**2)*(np.array(phi_e)))
        labels.append("Iteration number : " + str(cnt))
        cnt +=1
    #plt.plot(x,4*np.pi*(x**2)*(np.array(eigenPhi[1])))
    plt.xlim([0,5])
   # plt.plot(x,np.array(Phi)/x)
    
    r = np.linspace(0.0,10.5,1000)
    plt.legend(labels)
    plt.title("Convergence of charge density[Neon]",fontsize = 22)
    plt.xlabel('$r[a_0]$',fontsize = 14)
    plt.ylabel('$4\pi r^2 rho (r)$',fontsize = 14)
    
    
    [Eval,Evec] = readEigen("EigenVals0.txt","eigenvector.txt")
    #[Eval2,Evec2] =readEigen("EigenVals.txt","eigenvector2.txt")
    Eval = np.array(Eval)
    Evec = np.array(Evec)
   # Eval2 = np.array(Eval2)
   # Evec2 = np.array(Evec2)
   
    #x = np.linspace(min(Eval),max(Eval),len(Eval))
    converter = 2*13.60562953
    #plotEvals(Eval)
    plotEvecs(Evec)
    #Eval = np.array(sorted(Eval))*converter
    #plt.scatter(x[0:8],Eval[0:8])
    
    #idx = Eval.argsort()[::1]   
    #Eval = Eval[idx]
    test = Eval
    arr = []
    for i in range(len(Eval)):
        print(abs(test[i] -(- 1/(2*(i+1)**2))),i,test[i])
        #arr.append(abs(test[i] -(- 13.60562953/(i+2)**2)))
        if test[i] < 0:
            arr.append(abs(test[i] -(- 1/(2*(i+1)**2))))
            #arr.append(test[i])
        
    plt.figure(figsize = (10,5))
    xx = np.arange(len(arr))
    plt.scatter(xx,arr)
    plt.title("Energy difference for All bounded states",fontsize = 22)
    plt.ylabel('$|\Delta E|$',fontsize = 14)
    plt.xlabel('Energy-level',fontsize = 14)
    print(len(Evec[0]))
    
    
def plotIonEnergy(IonEnergy):
    converter = 2*13.60562953
    plt.figure(figsize = (10,5))
    xIon = np.arange(2,19,1)
    plt.plot(xIon,IonEnergy*converter,'r^')
    plt.xticks(np.arange(min(xIon), max(xIon)+1, 1.0))
    plt.grid(linewidth = 1)
    plt.title("Ionization energy for different elements",fontsize = 22)
    plt.ylabel("$E_I[eV]$",fontsize = 14)
    plt.xlabel("Atomic number $Z$",fontsize = 14)
    plt.text(2.5,23,"He",fontsize = 12)
    plt.text(10.5,17,"Ne",fontsize = 12)
    plt.text(17.5,11,"Ar",fontsize = 12)
    
    
    
def readIonEnergy():
    f = open("IonEnergy.txt")
    diff = []
    EIon = []
    for line in f:
        if not line.strip():
            EIon.append(abs(diff[1]-diff[0]))
            print(diff)
            diff =[]
        else:
            diff.append(float(line))
    return np.array(EIon)
            

def plotEvecs(Evec):
    plt.figure(figsize = (10,5))
    xx = np.logspace(-4,2,len(Evec[0]))
    labels = []
    for i in range(6):
        ##plt.plot(xx,(Evec2[i+1]**2))
        plt.plot(xx,(Evec[i]**2))
        labels.append('$|P_{' + str(i+1) + 's}(r)|^2$')
    plt.ylabel('Probability density',fontsize = 14)
    plt.xlabel('$r[a_0]$',fontsize = 14)
    plt.legend(labels,prop = {'size': 14})
    #plt.legend(['$|P_{2s}(r)|^2$','$|P_{2p}(r)|^2$'],prop = {'size': 14})
    plt.title("Probability distribution for hydrogen atom with $\ell = 0 $ ",fontsize =22)
    plt.grid()
    plt.xlim([0,15])
    
    
    
def plotEvals(Eval):
    plt.figure(figsize = (10,5))
    plt.axhline(y = 0 , color = 'g')
    for i in range(6):
        plt.axhline(y = -1/(2*(i+1)**2),linestyle ='-.',color = 'r')
        plt.axhline(y = Eval[i],linestyle ='--')

    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    plt.ylabel("$E[eV]$")
    #plt.ylim([-14,5])
    #plt.text(0.35,2.5,"Non-bound states",fontsize = 14)
    plt.legend(['0 reference','$E_0/ i^2$','Solved eigen-values'],prop = {'size': 14})
    plt.title("Energy levels for  hydrogen atom with $\ell = 1$",fontsize =22)


def Phi(r):
    phi = []
    Q = 1.0
    R = 5.0
    R2 = 6.0
    eps = 1/(4*np.pi)
    for val in r:
        if val <R:
            phi.append(Q/R)  
        elif R<= val <=R2:
            phi.append()
        #else:
        else:
            phi.append(0)
    return phi

def readEigen2():
    f = open("phiEigen.txt")
    Evec = []
    Evecs = []
    for line in f:
        if not line.strip():
            Evecs.append(Evec)
            Evec = []
        else:
            Evec.append(float(line))
    return Evecs

def readEigen(val,vec):
    f = open(val)
    F = open(vec)
    Evals = []
    Evecs = []
    Evec = []
    for line in f:
        Evals.append(float(line))
    for line in F:
        if not line.strip():
            Evecs.append(Evec)
            Evec = []
        else:
            Evec.append(float(line))
    return [Evals,Evecs]
        
    
def readPhi(file_name):
    f = open(file_name)
    phi = []
    for line in f:
        phi.append(float(line))
    return phi 

def readSplines():
    f = open("Bsplines.txt")
    cnt = 0
    bspline = []
    for r in range(1200):
        bspline.append([])
    for line in f:
            split = line.split(" ")
            if not line.strip():
                cnt = 0
            else:
                bspline[cnt].append(float(line))
                cnt += 1
    return bspline
if __name__=='__main__':
    main()