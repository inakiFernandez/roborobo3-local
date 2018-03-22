# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:13:52 2017

@author: fernandi
"""

#python3 tools/histogramsNbRobots.py robotsPerItem.log logs/robotsPerItemHist.png --png
import matplotlib 
matplotlib.use('svg')
import matplotlib.pyplot as plt
import os, time
import numpy as np
import brewer2mpl
import importlib
#import multirunFitness
import argparse
import glob
import multirunFitness as plotRuns
import scipy.stats as stats
from matplotlib import gridspec

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='To png or not')

    parser.add_argument('filenameHist', help='Filename (for histogram)')        
    parser.add_argument('outfile', help='output name for png')
    
      
    parser.add_argument('--png', action='store_true', help='output to png file?')
    args = parser.parse_args()

    isToPng = args.png
    #nbExp =  args.nbExp
        
    outfile= args.outfile
    filenameInput = args.filenameHist

    dpi = 100
    text_file = open(filenameInput, "r")
    lines = text_file.readlines()
    text_file.close()    
    data = []
    for line in lines:
        data.append(line.split(" ")[:-1])
    
    
    #print(data)
    #data = np.genfromtxt(filenameInput)
    valuesSet = np.unique(data)
    #print("Set of values: ", valuesSet)
    #histograms = {}
    bins = []
    #for val in valuesSet:
    for vec in data:
        binVector = np.bincount(vec)
        bins.append(binVector)
        #histograms.append({val:})
    #print("Bins:")    
    #print(bins)
    transposedBins = np.array(bins).T  
  
    ind = np.arange(len(bins))  # the x locations for the groups
    width = 0.6       # the width of the bars   
    fig, ax = plt.subplots()
    colors = ["red","green","blue","orange","pink"]
    for i, oneBin in enumerate(bins):
        for j,val in enumerate(oneBin[2:]):
            if j==0:
                ax.bar(i + width * j +0.2, val,width,color=colors[j],bottom=0)
            else:
                ax.bar(i + 0.2, val,width,color=colors[j],bottom=np.sum(oneBin[2:][:j]))
        #ax.bar(i+width+0.01, oneBin[1],width)
        #ax.bar(i+width*2+0.01, oneBin[2],width)
    #n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)
    fig.savefig(outfile,dpi=dpi)
    binsProportion = []
    for i, oneBin in enumerate(bins):
        totalGeneration = np.sum(oneBin)
        normalizedBin = np.array([float(val)/float(totalGeneration) for val in oneBin])
        binsProportion.append(normalizedBin)
    #print(np.array(binsProportion))
    fig, ax = plt.subplots()
    colors = ["red","green","blue","orange","pink"]
    
    for i, oneBin in enumerate(binsProportion):
        bars = []
        for j,val in enumerate(oneBin[2:]):
            if j==0:
                bars.append(ax.bar(i + width * j +0.2, val,width,color=colors[j],bottom=0))
            else:
                bars.append(ax.bar(i + 0.2, val,width,color=colors[j],bottom=np.sum(oneBin[2:][:j]))) #  + width * j 
    #legend
    import matplotlib.patches as mpatches

    patches = []
    for j, color in enumerate(colors[:3]):
        patches.append(mpatches.Patch(color=color, label=str(j+2)+" Robots"))
    
    plt.legend(handles=patches, loc=4)
    
    fig.savefig(outfile + "Proportion.png",dpi=dpi)    
    