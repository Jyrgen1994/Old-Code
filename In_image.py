# -*- coding: utf-8 -*-
"""
Created on Mon May 27 22:15:51 2019

@author: waterhorse
"""
import numpy as np
import matplotlib.pyplot as plt
import ImageLib as im

def main():
    N = 64
    m = N/2
    Im = im.ImageLib(N)
    Im.im_mat[30:34,30:34] = 1
    Im.im_mat[31:33,2:4] = 1
    Im.im_mat_noisy = Im.im_mat;
    #Im.defaultTransform(N,m)
    #Im.NosiyDefault(N,m,False,0.9)
    Im.pictureFilteringBW('mya_china.jpg',kernel_type ='log',sigma =2.5)
    
    
    
    
    

        

    
if __name__ == '__main__':
    main()