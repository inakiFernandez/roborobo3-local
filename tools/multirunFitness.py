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

def list_length(L):
    if L:
        return 1 + list_length(L[1:])
    return 0
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

def read_offspringfile(fname):
    fh = open(fname, 'r')
    #print(fname)
    #print(fh)
    rawData = []
    for line in fh:
        rawData.append(line)
    fh.close()
    filteredData = []
    isInGen = False
    generation = -1
    l = []
    for rawDatum in rawData:
        if "Start" in rawDatum:
            generation = int(rawDatum.split()[1])
            isInGen = True         
            l = []  
        else:
            if "End" in rawDatum:
                filteredData.append(l) #[generation,l])
            else:                    
                if not isInGen:
                    print("Mal")
                    exit()
                l.append([float(x) for i,x in enumerate(rawDatum.split(',')) if i != 5])
    return filteredData


def perc(data_l):
    data = np.asarray(data_l)
    median = np.zeros(data.shape[0])
    perc_25 = np.zeros(data.shape[0])
    perc_75 = np.zeros(data.shape[0])
    
    for i in range(0, len(median)):
        median[i] = np.median(data[i])
        perc_25[i] = np.percentile(data[i], 75)
        perc_75[i] = np.percentile(data[i], 25)
        #perc_25[i] = np.percentile(data[i, :], 5)
        #perc_75[i] = np.percentile(data[i, :], 95)
    return median, perc_25, perc_75


def plot_mean_curve(data, color, axis, label):
    mean = np.mean(data, 1)
    axis.plot(mean, lw=4, label=label, color=color)
    
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
        shadedColor = (color[0],color[1],color[2],color[3] * 0.5)
        axis.fill_between(np.arange(0, len(med)), perc_25, perc_75,
                          #alpha=0.25, 
                          linewidth=0, color=shadedColor)
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
    fig = plt.figure("Cooperative Foraging",figsize=(18, 12),dpi=70,facecolor=bgcolor)
    axis = plt.subplot2grid((1, 1), (0, 0)) #,facecolor=bgcolor) 
    yAxisName =  "Swarm Fitness" #"Global Diversity" #
    #names = ["Distributed", "Centralized"]
    names = ["F", "L","M","S"]
    names = ["Full", "Large","Medium","Small"] # #,"Full,", "Rate=0.1,Large","Rate=0.1,Medium","Rate=0.1,Small"]

    #typeComm= ["D","C"]
    #selP = ["0.0","0.25","0.5","0.75","1.0"]
    excessParam = 2
    #if sys.argv[1]=="-d":
    #    excessParam+=1
    #    yAxisName = "Local Diversity"
    nbRuns = 20
    labelsPlots = []
    handlesPlots = []
    colorsPlots = []
    legendObjects = []
    for i in range(len(sys.argv) - excessParam +1):
        dat = []
        #colors = [(1.0, 0.50, 0.0, float(i % 4 )/4.0 - 0.10),(0.25, 0.0, 1.0, float(i % 4 )/4.0 - 0.10)]
        print(colors)
        #if i % 4 < 1:
        #    colors =[(1.0, 0.50, 0.0, 0.05), (0.25, 0.0, 1.0, 0.05)]
        #if sys.argv[1]=="-d":
        #    for j in range(30):    
        #        dat.append([np.average(x) for x in read_logfile(sys.argv[excessParam + i] + "run-"+str(j+1) + ".log.localDiv.log")])
        #    datRow = [list(x)  for x in zip(*dat)]
        #    dat = datRow
        #    colorsPlots.append(colors[(i)%len(colors)])
        #    p = Rectangle((0, 0), 1, 1, fc=colors[(i)%len(colors)])        
        #else:
        dat = read_logfile(sys.argv[1 + i])
        color = colors[i]
        print(color)
        colorsPlots.append(color)
        p = Rectangle((0, 0), 1, 1, fc=color)
        #print(i+2, ": ",sys.argv[i+2])
        #print(len(dat))
        #datRow = [list(x)  for x in zip(*dat)]
        name = names[i] #"$"+typeComm[(i) % len(typeComm)] + ", θ_{\mathrm{sp}}="+ selP[(i) //2]+"$"
        #print(dat)
        labelsPlots.append(name)
        legendObjects.append(p)
        handlesPlots.append(plot_one_curve(dat, colors[(i)%len(colors)], axis, r''+name ,True))#sys.argv[excessParam + i], True)  #names[i], True)
    
    plt.ylabel(yAxisName,fontsize = 40)
    #plt.ylabel("Global Diversity",fontsize = 40)
    #plt.ylabel("Avg. Local Diversity",fontsize = 40)
    plt.xlabel("Generations",fontsize = 40)
    plt.xticks(fontsize = 35)
    plt.yticks(fontsize = 35)
    print(colorsPlots)
    print(len(colorsPlots))
    
    #legSorted,lblSorted =sortLegend(legendObjects,labelsPlots)
    legSorted,lblSorted = legendObjects,labelsPlots
    legend = axis.legend(legSorted,lblSorted,loc='best', fontsize = 23,ncol=2)
    legend.get_frame().set_alpha(0.75)
    axis.set_ylim(ymin=-0.01,ymax=6)
    #axis.set_xlim(xmin=-0.01,xmax=200)
    title = "$Broadcast\ rate \propto rank$" #"$Broadcast\ every\ 10\ it.$" # "$Always\ Broadcast$" #
    plt.title(r''+title,fontsize=50) 
    print(sys.argv[excessParam - 1] + ".test.all.png")
    plt.tight_layout()
    savefig(sys.argv[excessParam - 1] + ".test.all.png", dpi=100 )
    #plt.show()
