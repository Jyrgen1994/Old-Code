# -*- coding: utf-8 -*-
"""
Created on Mon May 27 22:19:21 2019

@author: waterhorse
"""
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import cv2
from mpl_toolkits.mplot3d import Axes3D
import scipy.ndimage as nd
import math

class ImageLib:
    im_mat = np.zeros((1,1))
    im_mat_noisy = np.zeros((1,1))
    laplace_filter = np.zeros((1,1))
    gaussian_blur = np.zeros((1,1))
    
    def __init__(self,im_size):
        self.im_mat = np.zeros((im_size,im_size))
        self.im_mat_noisy = np.zeros((im_size,im_size))
        self.laplace_filter = np.zeros((im_size,im_size))
        self.gaussian_blur = (1/256)*np.array([[1,4,6,4,1],[4,16,24,16,4],[6,24,36,24,6],[4,16,24,16,4],[1,4,6,4,1]])

    def defaultTransform(self,N,m):
        pic = self.im_mat;
        transform = np.fft.fft2(pic)
        transform_shift = np.fft.fftshift(transform)
        
        k = np.arange(1,m+1)
        X,Y = np.meshgrid(k,k)
        done = False
        plt.figure()
        while(not done):
            percentage = float(input("Percentace to remove: "))
          
            if(percentage <= 0 or percentage >=1):#break if percentage quote 0> p or p>1
                done =True
                break
            
            r = percentage*m
            Hq  = np.power(X,2) +np.power(Y,2) < r**2; H = np.hstack((Hq,np.fliplr(Hq))); H= np.vstack((H,np.flipud(H)));
            freq_transform = (H*transform); zapp = np.fft.ifft2(freq_transform);
            self.plotImages(pic,transform,transform_shift,zapp,percentage)
            
    def NosiyDefault(self,N,m,get_high_freq,p):
        pic = self.im_mat_noisy
        pic +=  0.2*np.random.rand(N,N);#Add noise to picture
        transform = np.fft.fft2(pic);transform_shift = np.fft.fftshift(transform);
        max_freq = abs(transform).max();
        print(max_freq)
        if(get_high_freq):
            H = abs(transform/max_freq)>p
        else:
            H = abs(transform/max_freq)<p
        
        freq_transform = H*transform;freq_transform_shift = np.fft.fftshift(freq_transform);zapp = np.fft.ifft2(freq_transform)
        
        
        k = np.arange(1,m+1)
        X,Y = np.meshgrid(k,k)
        r = p*m;
        
        plt.gray()
        plt.subplot(1,4,1)
        plt.imshow(pic)
        plt.subplot(1,4,2)
        plt.imshow(abs(transform_shift))
        plt.subplot(1,4,3)
        plt.imshow(abs(freq_transform_shift))
        plt.subplot(1,4,4)
        plt.imshow(abs(zapp))
        
    def pictureFilteringBW(self,picture_name,kernel_type,sigma = 1.4):
        pic = cv2.imread(picture_name,0)
        rows,cols = pic.shape
        #pic += 0.2*np.random.rand(rows,cols)
        transform = np.fft.fft2(pic)
        shift = np.fft.fftshift(transform)
        mag_freq = 20*np.log(np.abs(shift))
        print(len(pic))
        
        crow,ccol = int(rows/2),int(cols/2)
        
        filter_const = 20
        shift[crow-filter_const:crow+filter_const,ccol-filter_const:ccol+filter_const]= 0
        i_shift = np.fft.ifftshift(shift)
        pic_back = np.fft.ifft2(i_shift)
        pic_back = abs(pic_back)
        title = ''
        if kernel_type == 'log':
            kernel =self.createLoGKernel2(sigma,9)
            title = 'Filter: Laplace of Gaussian with pixel deviation σ = ' + str(sigma)
        elif kernel_type == 'Gaussian':
            kernel = self.gaussian_blur
        else:
            kernel = self.laplacianDefault()
            title = 'Filter: Laplacian'
        
        pic_filter = self.Convolution(kernel,pic)
        pic_zero = self.zeroCrossing(pic_filter)
        pic_enhanced = self.edgeEnhancement(pic_zero,pic)
        plt.figure(figsize =(12,15))
        plt.subplot(2,2,1)
        plt.title('Original Image')
        plt.imshow(pic,cmap = 'gray')
        plt.xticks([]);plt.yticks([]);
        plt.subplot(2,2,2)
        plt.title(title)
        plt.imshow(pic_enhanced,cmap='gray')
        plt.xticks([]);plt.yticks([]);
        plt.subplot(2,2,3)
        plt.title('Absolute value of filter after zero cross')
        plt.imshow(pic_zero,cmap ='gray')
        plt.xticks([]);plt.yticks([]);        
        plt.subplot(2,2,4)
        plt.title('Absolute value of filter before zero cross')    
        plt.imshow(np.abs(pic_filter),cmap = 'gray')
        plt.xticks([]);plt.yticks([]);        
        plt.show()
        
        
        plt.figure(figsize = (12,15))
        
        plt.subplot(2,2,1)
        plt.title('Original Image')
        plt.imshow(pic,cmap='gray')
        plt.xticks([]);plt.yticks([]);
        plt.subplot(2,2,2)
        plt.title('log of shifted FFT')
        plt.imshow(mag_freq,cmap='gray')
        plt.xticks([]);plt.yticks([]);
        plt.subplot(2,2,3)
        plt.title('Edge detection filter')
        plt.imshow(pic_back,cmap ='gray')
        plt.xticks([]);plt.yticks([]);
        plt.subplot(2,2,4)
        plt.title('log of shifted FFT with middle removed')
        plt.imshow(20*np.log(abs(shift)),cmap='gray')
        plt.xticks([]);plt.yticks([]);
        plt.show()
        
    def zeroCrossing(self,LoG):
        threshold = np.abs(LoG).mean()*0.75
        result = np.zeros(LoG.shape)
        rows,cols = LoG.shape
        for r in range(1,rows-1):
            for c in range(1,cols-1):
                diff_patch = LoG[r-1:r+2,c-1:c+2]
                peak = LoG[r,c]
                max_peak = diff_patch.max()
                min_peak = diff_patch.min()
                if(peak>0):
                    zero_cross = True if min_peak < 0 else False
                else:
                    zero_cross = True if max_peak > 0 else False
                
                if ((max_peak-min_peak)>threshold) and zero_cross:
                    #result[r,c] = abs(LoG[r,c])
                    result[r,c] = 1
        print(result.max())
        return result
        
    def edgeEnhancement(self,filtered_image,image):
        enhanced_image = np.zeros(image.shape)
        rows,cols = image.shape
        filtered_image = np.abs(filtered_image)
        filtered_image *= image.max()/filtered_image.max()
        for r in range(rows):
            for c in range(cols):
                enhanced_image[r,c] += (-filtered_image[r,c] + image[r,c])

        print('hej ',image.max(),np.abs(filtered_image).max(),np.abs(filtered_image).min())
        return np.abs(enhanced_image)
        
    def Convolution(self,kernel,image):
        rows,cols = image.shape
        krows,kcols = kernel.shape
        padded_image = np.zeros(shape = (rows+krows,cols+kcols))
        padded_image[krows//2:-krows//2,kcols//2:-kcols//2] = image.copy()
        filtered_image = np.zeros(shape = image.shape)
        for r in range(rows):
            for c in range(cols):
                for kr in range(krows):
                    for kc in range(kcols):
                        filtered_image[r,c] += padded_image[r + kr,c + kc]*kernel[kr,kc]
        #filtered_image -= kcols/2
        return filtered_image
        
        
    def laplacianDefault(self):
        #return np.array([[-1,-1,-1],[-1,8,-1],[-1,-1,-1]])
        return np.array([[1,1,1,1,1],[1,1,1,1,1],[1,1,-24,1,1],[1,1,1,1,1],[1,1,1,1,1]])
    
    def gaussLaplaceFilter(self,x,y,sigma):
        return -1.0/(np.pi*sigma**4)*(1.0-(x**2 + y**2)/(2.0*sigma**2))*np.exp(-(x**2+y**2)/(2.0*sigma**2))
    
    def createLoGKernel2(self,sigma,N):
        _min = int(math.ceil(-N/2))
        _max = int(math.floor(N/2))
        print(_min,_max)
        FoG_kernel = np.zeros((N,N))
        cnt_x = 0
        sum_all = 0
        sum_pos = 0
        correct_factor = 487
        #correct_factor = 400
        for x in range(_min,_max+1):
            for y in range(_min,_max+1):
                deriv_gauss = self.gaussLaplaceFilter(x,y,sigma)
                if deriv_gauss >0:
                    sum_pos += deriv_gauss
                sum_all += deriv_gauss
                
        
        if sum_pos == 0:
            sum_pos = 1
        print('hoo ',sum_all,sum_pos,abs(sum_all/sum_pos))
        for x in range(_min,_max+1):
            cnt_y = 0
            for y in range(_min,_max+1):
                deriv_gauss = self.gaussLaplaceFilter(x,y,sigma)
                if deriv_gauss >0 :
                    FoG_kernel[cnt_x,cnt_y] = ((deriv_gauss+ deriv_gauss*abs(sum_all/sum_pos))*correct_factor )
                else:
                    FoG_kernel[cnt_x,cnt_y] = (deriv_gauss*correct_factor)
                
                
                cnt_y += 1
            cnt_x += 1
        #print(sum_all,sum_pos,N*abs(sum_all/sum_pos))
        print(FoG_kernel,FoG_kernel.sum())
        return FoG_kernel
        
    
    def plotImages(self,pic,transform,transform_shift,zapp,percentage):
        plt.gray()
        plt.subplot(2,2,1) 
        plt.title("Original Image")
        plt.imshow(pic)
        plt.subplot(2,2,2)
        plt.title("FFT")
        plt.imshow(abs(transform))
        plt.subplot(2,2,3)
        plt.title("shifted FFr")
        plt.imshow(abs(transform_shift))
        plt.subplot(2,2,4)
        plt.title(str(int(100*(1-percentage))) + "removed")
        plt.imshow(abs(zapp))
        plt.show()
        
    
    