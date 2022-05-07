# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 22:27:48 2019

@author: waterhorse
"""
import numpy as np
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import Symbol


def main():
    h = 0.1
    g = -0.707976139657665
    k = 0.1
    s = 10
    x = np.arange(s)
    h_l = [h]
    g_l = [g]
    k_l = [k]
    
    for r in range(s-1):
       h = rec_h(h,0.1)
       k = rec_k(k,0.1)
       g = rec_g(g,0.1,0.1)
       h_l.append(h)
       k_l.append(k)
       g_l.append(g)
       
       
    plt.rc('text', usetex=True)
      
    plt.plot(x,h_l)
    plt.plot(x,k_l)
    #plt.plot(x,g_l)
    #plt.ylim([0,0.5])
#    x = np.linspace(-5,5,1000)
#    y = 2*x + (1/4)*np.log(16*(np.cosh(0.1)**2)*np.cosh(0.3)*np.cosh(0.1))
#    lin = x
#    plt.plot(x,lin)
#    plt.plot(x,y)
    plt.xlabel('iterations')
    plt.ylabel('parameters')
    plt.legend(['$h\'$','$k\'$','$k_0$'])
    X = Symbol('x')
    print(solve((1/2)*np.log(np.cosh(0.2+X)/np.cosh(0.2-X))))
    
def rec_h(h_old,k):
    return h_old + (1/2)*np.log(np.cosh(2*k + h_old)/np.cosh(2*k - h_old))

def rec_g(g_old,k,h):
    return 2*g_old + (1/4)*np.log(16*np.cosh(2*k + h)*np.cosh(2*k -h)*np.cosh(h)**2)

def rec_k(k_old,h):
    return (1/4)*np.log((np.cosh(2*k_old-h)*np.cosh(2*k_old+h)/np.cosh(h)**2))
    
    
if __name__ == '__main__':
    main()