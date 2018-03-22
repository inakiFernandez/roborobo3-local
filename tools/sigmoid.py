# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 14:29:59 2017

@author: fernandi
"""

import matplotlib 
matplotlib.use('svg')
import matplotlib.pyplot as plt
import os, time
import numpy as np
import math

if __name__ == '__main__':
    maxV = 2
    minV = -maxV
    xdata = np.arange(minV,maxV,0.01)

    ydata = [np.tanh(x) for x in xdata]
    dpi = 300
    labelFontSize = 30
    bgcolor="white"
    figAll = plt.figure(1,figsize=[8,6])
    axisAll = plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)
    
    plt.plot(xdata,ydata)
    #plotRuns.plot_one_curve(fitnessAllR, "orange", axisAll, "Random", True)
    axisAll.yaxis.grid(color="#BBBBBB", linestyle='-', linewidth=1)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.yticks(np.arange(-1, 1.1, 0.25))
    axisAll.tick_params(axis='both', which='major', labelsize=labelFontSize-10)
    
    bgcolor="white"
    axisAll.tick_params(axis='both', which='major', labelsize=labelFontSize-10)    

    axisAll.set_xlim((minV,maxV))
    axisAll.set_ylim((-1.1,1.1))
    
    #axisFitCumulAll.set_xticklabels(["B", "R"],fontsize=12)    
            
          
        
    axisAll.set_xlabel("Weighted sum",fontsize=labelFontSize)
    axisAll.set_ylabel("Effector activation",fontsize=labelFontSize)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)    

    colorValues = np.arange(-0.999,0.999,0.2499)    
    
    inverseTanh = []    
    for v in colorValues:
        inverseTanh.append(math.atanh(v))
    nVal = 100
    for i,invV in enumerate(inverseTanh):
        axisAll.plot([invV]*nVal,np.linspace(start=0.0,stop=colorValues[i],num=nVal), 
                     color="#999999", linestyle='--')
    print("Sigmoid saved")
    time.sleep(2)
    figAll.savefig("/home/fernandi/svn/docs/ecal2017/presentationECAL/fig/sigmoid.svg"  , dpi=dpi*3)        
    plt.close(figAll)



