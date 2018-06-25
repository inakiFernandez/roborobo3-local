# -*- coding: utf-8 -*-
"""
script plot several runs: quartiles and median through time
"""
import matplotlib
matplotlib.use('svg')
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import sys
from pylab import *
import brewer2mpl
from operator import itemgetter
import pandas

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



def plot_one_curve(x,data, color, axis, label, quartiles=False):
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
    lineWidth = 5
    handle = axis.plot_date(x,med, lw=lineWidth, label=label,               
              color=color,linestyle="-")
    gridcolor="#DDDDDD"    
    
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
    axis.grid(color=gridcolor, linestyle='-', linewidth=1)  
    #plt.grid(color=gridcolor,linewidth=1,linestyle='-') 
    
    return handle

if __name__ == "__main__":
    print(sys.argv)
    # args: file name
    if len(sys.argv) != 2:
        sys.exit("Error: wrong number of arguments\n" +
                 "Usage: python multirunFitness fileout")

    print(sys.argv)
    bgcolor = "white" 
    bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
    colors = bmap.mpl_colors
    
    #axis = subplot2grid((1, 1), (0, 0))
    fig = plt.figure("OptimCells",figsize=(10, 11),dpi=100,facecolor=bgcolor)
    axes = [plt.subplot2grid((3, 1), (0, 0)) , plt.subplot2grid((3, 1), (1, 0)),plt.subplot2grid((3, 1), (2, 0))]
    fnames = ["staff.log","incidents.log","clients.log"]
    fnamesS = ["staffS.log","incidentsS.log","clientsS.log"]
    expnames = ["#Deployed Agents", "#Failures", "#Clients without Power"]
    #years = matplotlib.dates.YearLocator()   # every year
    #months = matplotlib.dates.MonthLocator()  # every month
    #yearsFmt = matplotlib.dates.DateFormatter('%Y')
    dateFormat = '%H:%M %m/%d'
    
    datemin = np.datetime64(dt.datetime(2009, 1, 24, 7, 30))
    datemax = np.datetime64(dt.datetime(2009, 1, 27, 7, 0))
    print(datemin, datemax)        

    for i,fname in enumerate(fnames):
        dat = read_logfile(fname)
        datS = read_logfile(fnamesS[i])
        duration = datemax-datemin
        d = duration // len(dat)
        x = [datemin + d * i for i in range(len(dat))]
        plot_one_curve(x,dat, colors[1], axes[i], "GRN" ,False)
        plot_one_curve(x,datS, colors[2], axes[i], "Standard Policy",False)
        axes[i].set_title(expnames[i], fontsize = 22)
        
        axes[i].tick_params(labelsize = 18)
        #axes[i].set_yticklabels(labelsize = 18)
        #axes[i].xaxis.set_major_locator(years)
        #axes[i].xaxis.set_major_formatter(yearsFmt)
        #axes[i].xaxis.set_minor_locator(months)
        axes[i].set_xlim([datemin, datemax])
        
    #fig.autofmt_xdate()
    legends = [axes[0].legend(loc='best', fontsize = 22,ncol=2).get_frame().set_alpha(0.75),
               axes[1].legend(loc='best', fontsize = 22,ncol=2).get_frame().set_alpha(0.75),
               axes[2].legend(loc='best', fontsize = 22,ncol=2).get_frame().set_alpha(0.75)]
    plt.tight_layout()
    #axis.set_ylim(ymin=-0.01,ymax=np.max(fitness))
    #axis.set_xlim(xmin=-0.01,xmax=200)
    axes[0].format_xdata = matplotlib.dates.DateFormatter('%H:%M %m/%d')
    axes[1].format_xdata = matplotlib.dates.DateFormatter('%H:%M %m/%d')
    axes[2].format_xdata = matplotlib.dates.DateFormatter('%H:%M %m/%d')
            
    savefig(sys.argv[1] + ".png", dpi=100 )
    #plt.show()
