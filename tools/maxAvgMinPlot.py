# -*- coding: utf-8 -*-
"""
script plot several runs: quartiles and median through time
"""
import matplotlib
matplotlib.use('svg')

import matplotlib.pyplot as plt
import numpy as np
import sys
from pylab import *
import brewer2mpl
from operator import itemgetter

def read_logfile(fname):
    d = []
    fh = open(fname, 'r')
    for line in fh:
        l = []
        data = line.split()
        for o in data:
            l.append(float(o))
        d.append(l)
    fh.close()
    return d


def perc(data_l):
    data = np.asarray(data_l)
    median = np.zeros(data.shape[0])
    perc_25 = np.zeros(data.shape[0])
    perc_75 = np.zeros(data.shape[0])
    
    for i in range(0, len(median)):
        median[i] = np.median(data[i, :])
        perc_25[i] = np.percentile(data[i, :], 75)
        perc_75[i] = np.percentile(data[i, :], 25)
        #perc_25[i] = np.percentile(data[i, :], 5)
        #perc_75[i] = np.percentile(data[i, :], 95)
    return median, perc_25, perc_75


def plot_mean_curve(data, color, axis, label):
    mean = np.mean(data, 1)
    axis.plot(mean, lw=8, label=label, color=color)
    
    axis.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
    axis.tick_params(axis='x', direction='out')
    axis.tick_params(axis='y', length=0)
    for spine in axis.spines.values():
        spine.set_position(('outward', 5))
    axis.set_axisbelow(True)


def plot_one_curve(data, color, axis, label, quartiles=False):
    med, perc_25, perc_75 = perc(data)
    #print(data)
    #print(med)
    #print(perc_25)
    #print(perc_75)    
    #if color[3] == 0.0:
    #    color = (color[0],color[1],color[2],0.1)
    if quartiles:
        axis.fill_between(np.arange(0, len(med)), perc_25, perc_75,
                          alpha=0.25, linewidth=0, color=color)
    lineWidth = 2
    handle = axis.plot(med, lw=lineWidth, label=label,               
              color=color,linestyle="-")
    gridcolor="#FFFFFF"    
    
    #axis.spines['top'].set_visible(False)
    #axis.spines['right'].set_visible(False)
    #axis.spines['left'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
    axis.tick_params(axis='x', direction='out')
    axis.tick_params(axis='y', length=0)
    #for spine in axis.spines.values():
    #    spine.set_position(('outward', 5))
    #axis.set_axisbelow(True)
    #axis.grid(color='red', linestyle='-', linewidth=1)  
    #plt.grid(color=gridcolor,linewidth=1,linestyle='-') 
    #legend()
    return handle

def plot_one_curve_precomp(data, color, axis, label):
    
    med = data[1]
    perc_25 = data[0]
    perc_75 = data[2]
    
    axis.fill_between(np.arange(0, len(med)), perc_25, perc_75, alpha=0.25, linewidth=0, color=color)
    lineWidth = 25
    axis.plot(med, lw=lineWidth, label=label,               
              color=color,linestyle="-")
    gridcolor="#FFFFFF"    
    
    #axis.spines['top'].set_visible(False)
    #axis.spines['right'].set_visible(False)
    #axis.spines['left'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
    axis.tick_params(axis='x', direction='out')
    axis.tick_params(axis='y', length=0)
    #for spine in axis.spines.values():
    #    spine.set_position(('outward', 5))
    #axis.set_axisbelow(True)
    #axis.grid(color='red', linestyle='-', linewidth=1)  
    #plt.grid(color=gridcolor,linewidth=1,linestyle='-') 


def taskIntervals(horSize,interv=25):
    #vertical coordinates for taskswitch
    isT1 = True
    xcoords =  np.arange(0,horSize,interv)
    plt.locator_params(axis='y',nbins=30)
    for xc in xcoords:
        plt.axvline(x=xc,color='gray')
        if(isT1):
            plt.axvspan(xc, xc+interv, facecolor='grey', alpha=0.07)
        else:
            plt.axvspan(xc, xc+interv, facecolor='yellow', alpha=0.07)
        isT1 = not isT1

def sortLegend(objs,labels):
    print(labels)
    typeComm = [x[1] for x in labels]
    sp = [float(x[20:-1]) for x in labels]
    print(sp)
    print(typeComm)
    result = ([x for _,_,x in sorted(zip(typeComm,sp,objs), key=itemgetter(0,1),reverse=True)],
               [x for _,_,x in sorted(zip(typeComm,sp,labels), key=itemgetter(0,1),reverse=True)])
    return result

if __name__ == "__main__":
    # args: file name
    if len(sys.argv) < 3:
        sys.exit("Error: wrong number of arguments\n" +
                 "Usage: python multirunFitness [-d] filename [filename2 [...]]")

    print(sys.argv)
    bgcolor = "white" 
    bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
    colors = bmap.mpl_colors
    
    #axis = subplot2grid((1, 1), (0, 0))
    fig = plt.figure("MaxAvgMin",figsize=(21, 12.05),dpi=100,facecolor=bgcolor)
    axis = plt.subplot2grid((1, 1), (0, 0)) #,facecolor=bgcolor) 
    yAxisName =  "Fitness" #"Global Diversity" #
    legendObjects = []
    for i in range(len(sys.argv) - 2):
        dat = []
        dat = read_logfile(sys.argv[2 + i])
        print(dat)
        print(i+2, ": ",sys.argv[i+2])
        datColumns = [list(x) for x in zip(*dat)]
        print(datColumns)
        #print(len(dat))
        #datRow = [list(x)  for x in zip(*dat)]
        color = colors[(i+1)%len(colors)]
        axis.fill_between(np.arange(0, len(datColumns[2])), datColumns[2], datColumns[0],  alpha=0.25, linewidth=0, color=color)
        lineWidth = 2
        axis.plot(datColumns[1], lw=lineWidth, label=r''+sys.argv[1], color=color,linestyle="-")
         
        #axis.spines['top'].set_visible(False)
        #axis.spines['right'].set_visible(False)
        #axis.spines['left'].set_visible(False)
        axis.get_xaxis().tick_bottom()
        axis.get_yaxis().tick_left()
        axis.tick_params(axis='x', direction='out')
        axis.tick_params(axis='y', length=0)
#    plt.ylabel(yAxisName,fontsize = 40)
#    plt.xlabel("Time",fontsize = 40)
    plt.xticks(fontsize = 35)
    plt.yticks(fontsize = 35)    
    #axis.set_ylim(ymin=-0.01,ymax=np.max(fitness))
    #axis.set_xlim(xmin=-0.01,xmax=200)
    #plt.title("Swarm Fitness over Generations",fontsize=50) 
    plt.tight_layout()

    savefig(sys.argv[1] + ".test.png", dpi=100 )
