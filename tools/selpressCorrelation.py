# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 10:25:48 2017

@author: fernandi
For evostar19

"""
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

def plotMultipleNov18(exp,filetag,numberOfRobots,selectionPressures,controllerTypes,doVariance=True):
    parser = argparse.ArgumentParser(description='To png or not')
    parser.add_argument('outFolderName', help='Folder name out') 
    parser.add_argument('--png', action='store_true', help='output to png file?')
    args = parser.parse_args()
    isToPng = args.png
    outfolder = args.outFolderName  
    dpi = 100

    filesL1Subdiv = {25:{0.0:'R25.S0.5.C0.SP0.0.nbreg0.add0.del0.25', 
                         0.25:'R25.S0.5.C0.SP0.25.nbreg0.add0.del0.25',
                         0.5:'R25.S0.5.C0.SP0.5.nbreg0.add0.del0.25', 
                         0.75:'R25.S0.5.C0.SP0.75.nbreg0.add0.del0.25', 
                         1.0:'R25.S0.5.C0.SP1.0.nbreg0.add0.del0.25'},
                         
                   50:{0.0:'R50.S0.5.C0.SP0.0.nbreg0.add0.del0.25', 
                         0.25:'R50.S0.5.C0.SP0.25.nbreg0.add0.del0.25',
                         0.5:'R50.S0.5.C0.SP0.5.nbreg0.add0.del0.25', 
                         0.75:'R50.S0.5.C0.SP0.75.nbreg0.add0.del0.25', 
                         1.0:'R50.S0.5.C0.SP1.0.nbreg0.add0.del0.25'},      
               
                   100:{0.0:'R100.S0.5.C0.SP0.0.nbreg0.add0.del0.25', 
                        0.25:'R100.S0.5.C0.SP0.25.nbreg0.add0.del0.25', 
                        0.5:'R100.S0.5.C0.SP0.5.nbreg0.add0.del0.25',
                        0.75:'R100.S0.5.C0.SP0.75.nbreg0.add0.del0.25', 
                        1.0:'R100.S0.5.C0.SP1.0.nbreg0.add0.del0.25'}, 
               
                   200: {0.0:'R200.S0.5.C0.SP0.0.nbreg0.add0.del0.25', 
                         0.25:'R200.S0.5.C0.SP0.25.nbreg0.add0.del0.25', 
                         0.5:'R200.S0.5.C0.SP0.5.nbreg0.add0.del0.25', 
                         0.75:'R200.S0.5.C0.SP0.75.nbreg0.add0.del0.25', 
                         1.0:'R200.S0.5.C0.SP1.0.nbreg0.add0.del0.25'}
               }
    #print([str(i) + " " + str(j) for i in numberOfRobots for j in selectionPressures])
    filesSelected = [filesL1Subdiv[i][j] for i in numberOfRobots for j in selectionPressures]
    #print(filesSelected)
    
    dataFiles = [[[x.split('.')[0][1:], #x.split('.')[3][1],
                      (x.split('.')[4]+'.'+x.split('.')[5])[2:]],                     
                     x+'/all.log.'+filetag] for x in filesSelected if ".png" not in x]
    #print(dataFiles)
    labelFontSize = 14

    figAll = plt.figure(1,figsize=[8,6])
    gs = gridspec.GridSpec(1, 1, width_ratios=[1]) 
    axisAll = plt.subplot(gs[0]) 
    
    dataAll = [] 
    for i,fname in enumerate(dataFiles):
        dataAll.append([tuple(fname[0]),plotRuns.read_logfile(outfolder  + fname[1])])
 
    nColors = len(dataAll[0][1][0])
    print("NB colors: ",nColors)    

    colorArr = []
    colorArrShaded = []
    rawColors = {25:(1.0, 0.8, 0.4),50:(0.0,0.0,1.0),100:(0.1,0.55,0.18),200:(1.0,0.0,0.0)} # (1.0,0.55,0.18)
    #print(sorted(rawColors.keys()))
    usedRawColors = [rawColors[k] for k in sorted(rawColors.keys()) if k in numberOfRobots]        
        
    for i, nbR in enumerate(numberOfRobots):
        for j in range(len(selectionPressures)):
            colorArr += [usedRawColors[i] +tuple([1.0])]
            colorArrShaded += [usedRawColors[i]+tuple([0.7 * selectionPressures[j] +0.3])]
    
    count=0
    for i,v in enumerate(dataAll):
        plotRuns.plot_one_curve(v[1], colorArrShaded[count], axisAll,  r"$|\mathbf{S}|=$" + str(v[0][0]) + ", " + r"$\theta_{\mathrm{sp}=}$" + str(v[0][1]) , doVariance)
        count = (count+ 1) % len(colorArr)
    
    axisAll.tick_params(axis='both', which='major', labelsize=labelFontSize-3)
    legend = axisAll.legend(loc='lower right', #shadow=True, #bbox_to_anchor=[0, 1],
           ncol=4, fancybox=True) #, title="Legend")
    frame = legend.get_frame()
    frame.set_facecolor((0.8,0.8,0.8,0.025))
    frame.set_edgecolor('white')
    plt.xlabel("Generations",fontsize = 20) #,y=-0.06)
    plt.ylabel(exp,fontsize = 20) #,y=-0.06)
    for label in legend.get_texts():
        label.set_fontsize(6.5)
        label.set_color('black')
    for label in legend.get_lines():
        label.set_linewidth(4.6)  # the legend line width

    plt.tight_layout()
    plt.subplots_adjust(top=0.95)   
    strDoVariance = "All"
    if doVariance:
        strDoVariance = "MedAll"
    if isToPng:        
        print(exp+" saved to ", outfolder)
        time.sleep(0.5)
        figAll.savefig(outfolder+ "/"+ exp + strDoVariance #"MedAll" # "All" #
                       + "-R" + "-".join([str(x) for x in numberOfRobots]) 
                       + "-SP" + "-".join([str(x) for x in selectionPressures])
                       + "-C" + "-".join([str(x) for x in controllerTypes])
                       + ".png"  , dpi=dpi*3)        
        plt.close(figAll)
       
    traitIndexes = {"Fitness": 1, 
                    #"Genomes" : 3,  
                    #"Collisions": 2
                    } #, "Distance" : 'distance' }
    
    if exp not in traitIndexes.keys():
        return
    traitIndex = traitIndexes[exp]
    #print(outfolder + "/run*-logOffspring.txt",files)    
    allSpearData = {}
    for i, nbR in enumerate(numberOfRobots):
        v2={}
        for j,sp in enumerate(selectionPressures):
            fname = outfolder + filesL1Subdiv[nbR][sp]
            #print(fname)
            files = glob.glob(fname + "/run*-logOffspring.txt") #run*
            v1 = []
            for f in files:            
                v1.append(plotRuns.read_offspringfile(f))
            v2[sp] = v1

        allSpearData[nbR] = v2
    #print(allSpearData)
    #print(len(allSpearData[25][0.0][0]))
    #print(len(allSpearData[25]))
    #print(len(allSpearData[25][1.0]))
                
    fig = plt.figure(1,figsize=[8,6])
    gs = gridspec.GridSpec(1, 1, width_ratios=[1]) 
    ax = plt.subplot(gs[0]) 
    count = 0
    for i, nbR in enumerate(numberOfRobots):
        for j,sp in enumerate(selectionPressures):
            statAllVal = []             
            for run in allSpearData[nbR][sp]:
                #print(run)
                statVal = []
                offspring = [[x[0] for x in gen] for gen in run[:-1]] 
                fitness = [[x[int(traitIndex)] for x in gen] for gen in run[:-1]] 
                for k,gen in enumerate(offspring):
                    ranksO = stats.rankdata(offspring[k])
                    ranksF = stats.rankdata(fitness[k])
                    #kendall with number of discordances, spearman with deviations
                    tau = stats.kendalltau(ranksO,ranksF)[0]
                    statVal.append(tau) #spearmanr

                statAllVal.append(statVal)
            statistics = [list(x) for x in zip(*statAllVal)]
            plotRuns.plot_one_curve(statistics, 
                                    colorArrShaded[count], ax, 
                                    r"$|\mathbf{S}|=$" + str(nbR) + ", "  + r"$\theta_{\mathrm{sp}=}$" + str(sp) , 
                                    doVariance) #todo different traits: fit, collisions, genomes, items, distance #later with local diversity?
            count = (count+ 1) % len(colorArr)
            #print(statAllVal)            
    
    ax.tick_params(axis='both', which='major', labelsize=11)
    legend = ax.legend(loc='lower right', #shadow=True, #bbox_to_anchor=[0, 1],
                       ncol=4, fancybox=True) #, title="Legend")
    frame = legend.get_frame()
    frame.set_facecolor((0.8,0.8,0.8,0.025))
    frame.set_edgecolor('white')
    plt.xlabel("Generations",fontsize = 20) #,y=-0.06)
    plt.ylabel(r"$P_{\tau}$",fontsize = 20) #,y=-0.06)

    for label in legend.get_texts():
        label.set_fontsize(6.5)
        label.set_color('black')
    for label in legend.get_lines():
        label.set_linewidth(4.6)  # the legend line width

    plt.tight_layout()
    plt.subplots_adjust(top=0.95)    
    outfname = "/kendalltest"+ exp + strDoVariance \
               + "-R" + "-".join([str(x) for x in numberOfRobots]) \
               + "-SP" + "-".join([str(x) for x in selectionPressures]) \
               + "-C" + "-".join([str(x) for x in controllerTypes]) \
               + ".png"
    print("Spear saved to ", outfolder + outfname )
    #time.sleep(1)
    fig.savefig(outfolder+ outfname  , dpi=140)        
    plt.close(fig)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='To png or not')
    parser.add_argument('outFolderName', help='Folder name out') 
    parser.add_argument('--png', action='store_true', help='output to png file?')
    args = parser.parse_args()
    outfolder = args.outFolderName  


    doVariance = [False, 
                  True]
    #cheating here
    """
    """
    listExp = [
                [[25,50,100,200],[0,1.0]],     
                #[[25,50,100,200],[0]],
                #[[25,50,100,200],[1.0]],
                [[25,200],[0,1.0]],
                #[[25,200],[0,0.25,0.5]],
                #[[25,200],[0,0.5,1.0]],
                #[[25,200],[0.25,0.5,0.75]],
                #[[25,200],[0]],
                #[[25,200],[1.0]],
                #[[50,100,200],[0,0.25,0.5,0.75,1.0]],
                #[[200],[0,0.25,0.5,0.75,1.0]],
                #[[100],[0,0.25,0.5,0.75,1.0]],
                #[[50],[0,0.25,0.5,0.75,1.0]],
                #[[25],[0,0.25,0.5,0.75,1.0]],
                #[[25],[0,1.0]],
                #[[25],[0,0.25,0.5]],
                #[[25],[0,0.5,1.0]],
                #[[25],[0.25,0.5,0.75]],                
                #[[50],[0,1.0]],
                #[[50],[0,0.25,0.5]],
                #[[50],[0,0.5,1.0]],
                #[[50],[0.25,0.5,0.75]],
                #[[100],[0,1.0]],
                #[[100],[0,0.25,0.5]],
                #[[100],[0,0.5,1.0]],
                #[[100],[0.25,0.5,0.75]],
                #[[200],[0,1.0]],
                #[[200],[0,0.25,0.5]],
                #[[200],[0,0.5,1.0]],
                #[[200],[0.25,0.5,0.75]],              
                #[[100,200],[0,1.0]],
                #[[100,200],[0,0.25,0.5]],
                #[[100,200],[0,0.5,1.0]],
                #[[100,200],[0.25,0.5,0.75]],
                #[[50,200],[0,1.0]],
                #[[50,200],[0,0.25,0.5]],
                #[[50,200],[0,0.5,1.0]],
                #[[50,200],[0.25,0.5,0.75]],
                #[[50,200],[0,1.0]],
                #[[50,200],[0,0.25,0.5]],
                #[[25,50,100,200],[0,0.25,0.5,0.75,1.0]],                
                #[[50,100],[1.0]],
                #[[50,100],[0,0.25,0.5]],
                #[[50,100],[0,0.5,1.0]],
                #[[50,100],[0.25,0.5,0.75]],
                #[[25,100,200],[0,1.0]],
                #[[25,100,200],[0,0.25,0.5]],
                #[[25,100,200],[0,0.5,1.0]],
                #[[25,100,200],[0.25,0.5,0.75]],
                #[[25,200],[0,1.0]],
                #[[25,200],[0,0.25,0.5]],
                #[[25,200],[0,0.5,1.0]],
                #[[25,200],[0.25,0.5,0.75]],
                #[[25,100],[0,1.0]],
                #[[25,100],[0,0.25,0.5]],
                #[[25,100],[0,0.5,1.0]],
                #[[25,100],[0.25,0.5,0.75]]
               
              ]


    """
    

    """

    cTypes = [0,
              #1
              ]

    """
    nbOfRobots = [25,
                  50,
                  100,
                  200
                  ]
    selPressures = [0,
                    #0.25,
                    #0.5,
                    #0.75,
                    1.0
                    ]
    addRates = [0,
                #0.25
                ]


    
    filename = "R100.S0.5.C0.SP1.0.nbreg0.add0.del0.25/"
    filename = "run24-logOffspring.txt"  
    
    data = plotRuns.read_offspringfile(outfolder  + filename)

    print(len(data))    
    print(data)    
    
    offspring1 = [y[1] for y in data]    
    offspring = [[x[0] for x in y] for y in offspring1][:-1]

    fitness1 = [y[1] for y in data]    
    fitness = [[x[1] for x in y] for y in fitness1][:-1]
    
    #print(offspring)
    #print(len(offspring))
    #print(fitness)
    #print(len(fitness))
    
    #kendall with number of discordances, spearman with deviations
    spearmanVal = []
    for i,o in enumerate(offspring):
        f = fitness[i]
        ranksO = stats.rankdata(o)
        ranksF = stats.rankdata(f)
        spearmanVal.append(stats.kendalltau(ranksO,ranksF)[0]) #spearmanr
    
    
    fig = plt.figure(1,figsize=[8,6])
    gs = gridspec.GridSpec(1, 1, width_ratios=[1]) 
    ax = plt.subplot(gs[0]) 
    
    plotRuns.plot_one_curve(spearmanVal, (0.7,0.3,0.4,0.7), ax, "Spearman S.P coeff" , True)
    ax.tick_params(axis='both', which='major', labelsize=11)
    legend = ax.legend(loc='lower right', shadow=True, #bbox_to_anchor=[0, 1],
           ncol=5, fancybox=True) #, title="Legend")
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    
    
    for label in legend.get_texts():
        label.set_fontsize(5)
    for label in legend.get_lines():
        label.set_linewidth(0.6)  # the legend line width

    plt.tight_layout()
    plt.subplots_adjust(top=0.95)    
       
    print("saved to ", outfolder)
    time.sleep(2)
    fig.savefig(outfolder+ "/spearmantest.png"  , dpi=140)        
    plt.close(fig)
  
    print(spearmanVal)
    exit()
    """


    exp = "Fitness" # "Diversity" # "Genomes" # "Collisions" # "Distance" #
    filetags = {"Fitness": '4', "Diversity" : '3', "Genomes" : 'genomes',  "Collisions": 'collisions'} #, "Distance" : 'distance' }
    for b in doVariance:    
        for e in listExp:
            nbOfRobots = e[0]
            selPressures = e[1]
            for k in filetags.keys():      
                exp = k
                filetag = filetags[exp]
                plotMultipleNov18(exp,filetag,nbOfRobots,selPressures,cTypes,doVariance=b)
    
    
    exit()    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    parser = argparse.ArgumentParser(description='To png or not')
    
    #parser.add_argument('itemsPerColorFile', help='input filename items per color')
    #parser.add_argument('itemsGatheredAndPossibleFile', help='input filename items gathered and possible')
    #parser.add_argument('numberChangesFile', help='input filename number of changes')
    #parser.add_argument('givenAverageReward', help='input filename average individual reward')
    parser.add_argument('folderName', help='Folder name data') 
    parser.add_argument('outFolderName', help='Folder name out') 
    #parser.add_argument('nbExp', type=int, help='Number of variants of the experiments')
    
    parser.add_argument('--png', action='store_true', help='output to png file?')
    #parser.add_argument('itPerTask', type=int, help='Number of iterations per task')                
    #parser.add_argument('runId', type=int, help='ID of the current run')                
    args = parser.parse_args()

    isToPng = args.png
    #nbExp =  args.nbExp
        
    outfolder = args.outFolderName  
    folderData= args.folderName
    #cheating here
    numberOfRobots = [50,
                      #100,
                      #200
                      ]
    selectionPressures = [0,
                          #0.25,
                          #0.5,
                          #0.75,
                          1.0
                          ]
    controllerTypes = [0,
                       1
                       ]
    #runId = args.runId
    dpi = 100
    
    #os.path.dirname(os.path.abspath(args.folderNameBest))
    #dirFilesItemsPerColorB = glob.glob(datapathB + "/items-run*.log") # 
    filesL1 = os.listdir(outfolder)    
    filesL1 = [x for x in filesL1 if (".png" not in x) and (".log" not in x)]
    sortedFilesL1 = sorted(filesL1, key=str.lower)

    #cheating    
    
    filesL1 = ['R50.T1.S0.5.C0.SP0.0', 
               #'R50.T1.S0.5.C0.SP0.25',
               #'R50.T1.S0.5.C0.SP0.5', 
               #'R50.T1.S0.5.C0.SP0.75', 
               'R50.T1.S0.5.C0.SP1.0',
               
               'R100.T1.S0.5.C0.SP0.0', 
               #'R100.T1.S0.5.C0.SP0.25', 
               #'R100.T1.S0.5.C0.SP0.5',
               #'R100.T1.S0.5.C0.SP0.75', 
               'R100.T1.S0.5.C0.SP1.0', 
               'R200.T1.S0.5.C0.SP0.0', 
               #'R200.T1.S0.5.C0.SP0.25', 
               #'R200.T1.S0.5.C0.SP0.5', 
               #'R200.T1.S0.5.C0.SP0.75', 
               'R200.T1.S0.5.C0.SP1.0'
               ]
    filesL1Subdiv = {50:{0.0:'R50.T1.S0.5.C0.SP0.0', 
                         0.25:'R50.T1.S0.5.C0.SP0.25',
                         0.5:'R50.T1.S0.5.C0.SP0.5', 
                         0.75:'R50.T1.S0.5.C0.SP0.75', 
                         1.0:'R50.T1.S0.5.C0.SP1.0'},
               
                   100:{0.0:'R100.T1.S0.5.C0.SP0.0', 
                        0.25:'R100.T1.S0.5.C0.SP0.25', 
                        0.5:'R100.T1.S0.5.C0.SP0.5',
                        0.75:'R100.T1.S0.5.C0.SP0.75', 
                        1.0:'R100.T1.S0.5.C0.SP1.0'}, 
               
                   200: {0.0:'R200.T1.S0.5.C0.SP0.0', 
                         0.25:'R200.T1.S0.5.C0.SP0.25', 
                         0.5:'R200.T1.S0.5.C0.SP0.5', 
                         0.75:'R200.T1.S0.5.C0.SP0.75', 
                         1.0:'R200.T1.S0.5.C0.SP1.0'}
               }

    sortedFilesL1= filesL1

    filesL2  = [os.listdir(x) for x in sortedFilesL1 if os.path.isdir(x)]

    filesSelected = [filesL1Subdiv[i][j] for i in numberOfRobots for j in selectionPressures]
    print(filesSelected)
    
    fitnessFiles2 = [[[x.split('.')[0][1:],x.split('.')[4][1],
                      (x.split('.')[5]+'.'+x.split('.')[6])[2:]],                     
                     x+'/all.log.1'] for x in filesSelected if ".png" not in x]
    #print(x.split('.')[0][1:],x.split('.')[4][1],
    #                  (x.split('.')[5]+'.'+x.split('.')[6])[2:])
    fitnessFiles = [[[x.split('.')[0][1:],x.split('.')[4][1],
                      (x.split('.')[5]+'.'+x.split('.')[6])[2:]],                     
                     x+'/all.log.1'] for x in sortedFilesL1 if ".png" not in x]
    fitnessFiles = fitnessFiles2
    diversityFiles = [[[x.split('.')[0][1:],x.split('.')[4][1],
                      (x.split('.')[5]+'.'+x.split('.')[6])[2:]],                     
                     x+'/all.log.3'] for x in filesSelected if ".png" not in x]    
    #print(filesL1)
    #print(filesL2)
    #print(fitnessFiles)
    labelFontSize = 14

    #print(dirFilesColorChanges)
    figAll = plt.figure(1,figsize=[8,6])
    #gs = gridspec.GridSpec(1, 2, width_ratios=[8, 1]) 
    gs = gridspec.GridSpec(1, 1, width_ratios=[1]) 
    axisAll = plt.subplot(gs[0]) #,facecolor=bgcolor) #plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)
    #axisFitCumulAll= plt.subplot(gs[1], #facecolor=bgcolor,
    #                             sharey=axisAll) #plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)
    #exit()
    fitnessAll = [] #{}
    for i,fname in enumerate(fitnessFiles):
        #print(fname)
        fitnessAll.append([tuple(fname[0]),plotRuns.read_logfile(outfolder  + fname[1])])
    #print(fitnessAll)
    
    print(len(fitnessAll))
    print(sorted(fitnessAll)) #.keys())
    fitnessAllPerRun = fitnessAll
    print(len(fitnessAllPerRun[0]))    
    print(len(fitnessAllPerRun[0][0]))
    bgcolor="white"
    nColors = len(fitnessAllPerRun[0][1][0])
    print("NB colors: ",nColors)    
    colorArr = []
    colorArrShaded = []
    for j in range(nColors):
        absVal = (256/nColors) * j;
        colorArr += [(absVal/256.0, (255 - absVal)/256.0, 0)]
        colorArrShaded += [(absVal/256.0, (255 - absVal)/256.0, 0,0.6)]
    colorArr = []
    colorArrShaded = []
    rawColors = [(0.0,0.0,1.0),(0.1,0.55,0.18),(1.0,0.0,0.0)] # (1.0,0.55,0.18)
    rawColors = {50:(0.0,0.0,1.0),100:(0.1,0.55,0.18),200:(1.0,0.0,0.0)} # (1.0,0.55,0.18)
    print(sorted(rawColors.keys()))
    usedRawColors = [rawColors[k] for k in sorted(rawColors.keys()) if k in numberOfRobots]        
        
    for i, nbR in enumerate(numberOfRobots):
        for j in range(len(selectionPressures)):
            colorArr += [usedRawColors[i] +tuple([1.0])]
            colorArrShaded += [usedRawColors[i]+tuple([0.8 * selectionPressures[j] +0.1])]
    
    #change colors to rainbow [], .append
    #colorArr = [(0.0,0.0,1.0),(0.0,145.0/255.0,1.0),(0.0,1.0,218.0/255.0),(0.0,1.0,72.0/255.0),
    #            (72.0/255.0,1.0,0.0),(218.0/255.0,1.0,0.0),(1.0,145.0/255,0.0),(1.0,0.0,0.0)]
    #colorArrShaded = [(0.0,0.0,1.0,0.6),(0.0,145.0/255.0,1.0,0.6),(0.0,1.0,218.0/255.0,0.6),(0.0,1.0,72.0/255.0,0.6),
    #            (72.0/255.0,1.0,0.0,0.6),(218.0/255.0,1.0,0.0,0.6),(1.0,145.0/255,0.0,0.6),(1.0,0.0,0.0,0.6)]

    #exp,fitnessValues
    #for k,v in fitnessAllPerRun:
    #    print(k,v)
    #    fitnessAll[k] = [list(x) for x in zip(*v)]
    #print(fitnessAll)
    #print(len(fitnessAll))
    #print(fitnessAll.keys())
    #fitnessCumulated = [np.average(fitRun[int(0.8 * len(fitRun)):]) for fitRun in fitnessAllPerRun]    
    #exit()
    # to compare here: mannWhitCumulFit = stats.mannwhitneyu(fitnessCumulatedB,fitnessCumulatedR)
    #print("P-val mann whitney U Cumulated Fitness best vs. random" + str(mannWhitCumulFit[1]))
    count=0
    #for k in fitnessAll.keys():
    for i,v in enumerate(fitnessAll):
        #print(fitnessAll[tuple(k)])
        
        plotRuns.plot_one_curve(v[1], colorArrShaded[count], axisAll, v[0], True)
        count = (count+ 1) % len(colorArr)
        #plotRuns.plot_one_curve(fitnessAllR, "orange", axisAll, "Random", True)
    
        
    #axisAll.grid(color="#FFFFFF", linestyle='-', linewidth=1)
    
    axisAll.tick_params(axis='both', which='major', labelsize=labelFontSize-3)
    legend = axisAll.legend(loc='lower right', shadow=True, #bbox_to_anchor=[0, 1],
           ncol=5, fancybox=True) #, title="Legend")
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    
    
    for label in legend.get_texts():
        label.set_fontsize(5)
    for label in legend.get_lines():
        label.set_linewidth(0.6)  # the legend line width

    bgcolor="white"
    #figFitCumulAll = plt.figure(152)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)    
    if isToPng:        
        print("Fitness saved to ", outfolder)
        time.sleep(2)
        figAll.savefig(outfolder+ "/FitnessAll"
                       + "-R" + "-".join([str(x) for x in numberOfRobots]) 
                       + "-SP" + "-".join([str(x) for x in selectionPressures])
                       + "-C" + "-".join([str(x) for x in controllerTypes])
                       + ".png"  , dpi=dpi*3)        
        plt.close(figAll)
    
    figDiv = plt.figure(1,figsize=[8,6])
    #gs = gridspec.GridSpec(1, 2, width_ratios=[8, 1]) 
    gs = gridspec.GridSpec(1, 1, width_ratios=[1]) 
    axisDiv = plt.subplot(gs[0]) #,facecolor=bgcolor) #plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)
    divAll = [] #{}
    for i,fname in enumerate(diversityFiles):
        #print(fname)
        divAll.append([tuple(fname[0]),plotRuns.read_logfile(outfolder  + fname[1])])
    #print(fitnessAll)
    
    print(len(divAll))
    print(sorted(divAll)) #.keys())
    divAllPerRun = divAll
    print(len(divAllPerRun[0]))    
    print(len(divAllPerRun[0][0]))
    bgcolor="white"
    nColors = len(divAllPerRun[0][1][0])
    colorArr = []
    colorArrShaded = []
    rawColors = {50:(0.0,0.0,1.0),100:(0.1,0.55,0.18),200:(1.0,0.0,0.0)} # (1.0,0.55,0.18)
    usedRawColors = [rawColors[k] for k in sorted(rawColors.keys()) if k in numberOfRobots]        
        
    for i, nbR in enumerate(numberOfRobots):
        for j in range(len(selectionPressures)):
            colorArr += [usedRawColors[i] +tuple([1.0])]
            colorArrShaded += [usedRawColors[i]+tuple([0.8 * selectionPressures[j] +0.1])]
    
    count=0
    for i,v in enumerate(divAll):
        plotRuns.plot_one_curve(v[1], colorArrShaded[count], axisDiv, v[0], True)
        count = (count+ 1) % len(colorArr)
    
    axisDiv.tick_params(axis='both', which='major', labelsize=labelFontSize-3)
    legend = axisDiv.legend(loc='upper right', shadow=True, #bbox_to_anchor=[0, 1],
           ncol=5, fancybox=True) #, title="Legend")
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    axisDiv.set_xlabel("Diversity")    
    
    for label in legend.get_texts():
        label.set_fontsize(5)
    for label in legend.get_lines():
        label.set_linewidth(0.6)  # the legend line width

    bgcolor="white"
    #figFitCumulAll = plt.figure(152)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)    
    if isToPng:       
        outfileName = outfolder + "/diversityAll" + "-R" + "-".join([str(x) for x in numberOfRobots]) + "-SP" + "-".join([str(x) for x in selectionPressures]) + "-C" + "-".join([str(x) for x in controllerTypes]) + ".png"  
        print("Diversity saved to ", outfileName)
        time.sleep(2)
        figDiv.savefig(outfileName, dpi=dpi*3)        
        plt.close(figDiv)

    exit()    
    
        
    
    
    
    violin_parts = axisFitCumulAll.violinplot([fitnessCumulatedB,fitnessCumulatedR], positions=[0.5,1.0],
                                              points=len(fitnessCumulatedB), widths=0.25,
                                              showmeans=False, showextrema=True, showmedians=True)  
    
    violin_parts['bodies'][0].set_facecolor((0.0,0.0,1.0,0.5))
    violin_parts['bodies'][0].set_edgecolor("blue")
    violin_parts['bodies'][1].set_facecolor((1.0,0.75,0.0,0.5))
    violin_parts['bodies'][1].set_edgecolor("orange")
    
    #print([[1.0,v] for v in fitnessCumulatedB])
    axisFitCumulAll.plot([0.5 for v in fitnessCumulatedB], fitnessCumulatedB, 'x', color="blue",mew=3, ms=3)
    axisFitCumulAll.plot([1.0 for v in fitnessCumulatedR], fitnessCumulatedR, 'x', color="orange", mew=3, ms=3)
        
    axisFitCumulAll.set_xticks([0.5,1.0])
    axisAll.tick_params(axis='both', which='major', labelsize=labelFontSize-3)    

    axisFitCumulAll.set_xlim((0.2,1.3))
    #axisFitCumulAll.set_xticklabels(["Best", "Random"],rotation=90,fontsize=12)    
    axisFitCumulAll.set_xticklabels(["B", "R"],fontsize=12)    
    axisAll.set_xlabel("Generations",fontsize=labelFontSize)
    #axisAll.set_ylabel("Swarm Fitness (#Items)",fontsize=labelFontSize)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)    
    if isToPng:        
        print("Fitness saved to both", outfileBest)
        time.sleep(2)
        figAll.savefig(outfileBest + "FitnessAllBoth.png"  , dpi=dpi*3)        
        plt.close(figAll)
    
    nColors = 8
    colorArr = []
    for j in range(nColors):
        absVal = (256/nColors) * j;
        colorArr += [(absVal/256.0, (255 - absVal)/256.0, 0)]
    
    #change colors to rainbow [], .append
    colorArr = [(0.0,0.0,1.0),(0.0,145.0/255.0,1.0),(0.0,1.0,218.0/255.0),(0.0,1.0,72.0/255.0),
                (72.0/255.0,1.0,0.0),(218.0/255.0,1.0,0.0),(1.0,145.0/255,0.0),(1.0,0.0,0.0)]
    colorArrShaded = [(0.0,0.0,1.0,0.6),(0.0,145.0/255.0,1.0,0.6),(0.0,1.0,218.0/255.0,0.6),(0.0,1.0,72.0/255.0,0.6),
                (72.0/255.0,1.0,0.0,0.6),(218.0/255.0,1.0,0.0,0.6),(1.0,145.0/255,0.0,0.6),(1.0,0.0,0.0,0.6)]
    
    nbItems = 100 # 40 # 80 # 200 # 30
    nbSteps = 800 #500
    
    valuesItemsAggAllB = [[],[]]
    for fname in dirFilesGatheredB:
        with open(fname) as f:
            content = f.readlines()
    
        content = [x.strip() for x in content]
    
        contentSplit = [x.split(" ") for x in content]
                    
        #contentIncr= [list(x) for x in  zip(*contentSplit)]                         
    
        valuesItemsB = [ [(float(x[0]) / float(x[1])),  # / float(nbItems), 
                          float(x[1]) / float(nbItems)] if(x[1] != '0') else [0.0, 0.0] for x in contentSplit]
    
        valuesItemsB= [list(x) for x in zip(*valuesItemsB)]
    
        valuesItemsAggGenerationB = [[],[]]
        for i,x in enumerate(valuesItemsB):    
            nbGen = int(len(x) / nbSteps)
            aggPerIter = [np.average(x[nbSteps * i:nbSteps * (i + 1)]) 
                          for i in range(nbGen)]
            valuesItemsAggGenerationB[i] = aggPerIter
        valuesItemsAggAllB[0].append(valuesItemsAggGenerationB[0])
        valuesItemsAggAllB[1].append(valuesItemsAggGenerationB[1])    

    for i in range(len(valuesItemsAggAllB)):
        valuesItemsAggAllB[i] =   [list(x) for x in zip(*valuesItemsAggAllB[i])] 
        
    
    valuesItemsAggAllR = [[],[]]
    for fname in dirFilesGatheredR:
        with open(fname) as f:
            content = f.readlines()
    
        content = [x.strip() for x in content]
    
        contentSplit = [x.split(" ") for x in content]
                    
        #contentIncr= [list(x) for x in  zip(*contentSplit)]                         
    
        valuesItemsR = [ [(float(x[0]) / float(x[1])),  # / float(nbItems), 
                          float(x[1]) / float(nbItems)] if(x[1] != '0') else [0.0, 0.0] for x in contentSplit]
    
        valuesItemsR= [list(x) for x in zip(*valuesItemsR)]
    
        valuesItemsAggGenerationR = [[],[]]
        for i,x in enumerate(valuesItemsR):    
            nbGen = int(len(x) / nbSteps)
            aggPerIter = [np.average(x[nbSteps * i:nbSteps * (i + 1)]) 
                          for i in range(nbGen)]
            valuesItemsAggGenerationR[i] = aggPerIter
        valuesItemsAggAllR[0].append(valuesItemsAggGenerationR[0])
        valuesItemsAggAllR[1].append(valuesItemsAggGenerationR[1])    

    for i in range(len(valuesItemsAggAllR)):
        valuesItemsAggAllR[i] =   [list(x) for x in zip(*valuesItemsAggAllR[i])]
    #print(len(valuesItemsAggAll[0]))
    #print(len(valuesItemsAggAll[0][0]))
    possibleAllBPerRun = [list(x) for x in zip(*valuesItemsAggAllB[1])]
    gatheredAllBPerRun = [list(x) for x in zip(*valuesItemsAggAllB[0])]
    possibleAllRPerRun = [list(x) for x in zip(*valuesItemsAggAllR[1])]
    gatheredAllRPerRun = [list(x) for x in zip(*valuesItemsAggAllR[0])]
    possibleCumulatedB = [np.average(possibleRun[int(0.8 * len(possibleRun)):]) for possibleRun in possibleAllBPerRun]    
    possibleCumulatedR = [np.average(possibleRun[int(0.8 * len(possibleRun)):]) for possibleRun in possibleAllRPerRun]
    
    
    gatheredCumulatedB = [np.average(gatheredRun[int(0.8 * len(gatheredRun)):]) for gatheredRun in gatheredAllBPerRun]    
    gatheredCumulatedR = [np.average(gatheredRun[int(0.8 * len(gatheredRun)):]) for gatheredRun in gatheredAllRPerRun]
    #?Do per iteration?    
    mannWhitCumulPossible = stats.mannwhitneyu(possibleCumulatedB,possibleCumulatedR)
    mannWhitCumulGathered= stats.mannwhitneyu(gatheredCumulatedB,gatheredCumulatedR)       
    print("P-val mann whitney U Cumulated possible best vs. random" + str(mannWhitCumulPossible[1]))
    print("P-val mann whitney U Cumulated gathered best vs. random" + str(mannWhitCumulGathered[1]))

    bgcolor="white"
    fig2 = plt.figure(2,figsize=[8,6])
    gs = gridspec.GridSpec(1, 2, width_ratios=[8, 1])
    axis2 = plt.subplot(gs[0],facecolor=bgcolor) # plt.subplot2grid((1, 2), (0, 0), facecolor=bgcolor)
    axisBxpPossible = plt.subplot(gs[1],facecolor=bgcolor,sharey=axis2)
    
    #axis2 = plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)
    
    plotRuns.plot_one_curve(valuesItemsAggAllB[1], "blue", axis2, "Best", True)    
    plotRuns.plot_one_curve(valuesItemsAggAllR[1], "orange", axis2, "Random", True)    
    #plt.suptitle("Proportion of possible items to be collected per generation",fontsize=labelFontSize)
    axis2.set_xlabel("Generations",fontsize=labelFontSize)
    #axis2.set_ylabel("% of possible items",fontsize=labelFontSize)
    #axis2.grid(color="#FFFFFF", linestyle='-', linewidth=1)    
    axis2.tick_params(axis='both', which='major', labelsize=labelFontSize-3)    
    
    legend = axis2.legend(loc='upper left', shadow=True)
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(1.5)
    bgcolor="white"
    #figFitCumulAll = plt.figure(152)
    
    violin_parts = axisBxpPossible.violinplot([possibleCumulatedB,possibleCumulatedR], positions=[0.5,1.0],
                                              points=len(possibleCumulatedB), widths=0.25,
                                              showmeans=False, showextrema=True, showmedians=True)  
    
    violin_parts['bodies'][0].set_facecolor((0.0,0.0,1.0,0.5))
    violin_parts['bodies'][0].set_edgecolor("blue")
    violin_parts['bodies'][1].set_facecolor((1.0,0.75,0.0,0.5))
    violin_parts['bodies'][1].set_edgecolor("orange")
    
    #print([[1.0,v] for v in fitnessCumulatedB])
    axisBxpPossible.plot([0.5 for v in possibleCumulatedB], possibleCumulatedB, 'x', color="blue",mew=3, ms=3)
    axisBxpPossible.plot([1.0 for v in possibleCumulatedR], possibleCumulatedR, 'x', color="orange", mew=3, ms=3)
        
    axisBxpPossible.set_xticks([0.5,1.0])
    axisBxpPossible.set_xlim((0.2,1.3))
    axisAll.tick_params(axis='y', which='major', labelsize=labelFontSize-3)        
    #axisBxpPossible.set_xticklabels(["Best", "Random"],rotation=90,fontsize=12)    
    axisBxpPossible.set_xticklabels(["B", "R"],fontsize=12)    
    

    plt.tight_layout()    
    plt.subplots_adjust(top=0.93)  
    if isToPng:        
        print("Possible Best vs. Random saved to ", outfileBest)
        time.sleep(2)
        fig2.savefig(outfileBest + "PossibleProportionsBoth.png"  , dpi=dpi * 3)        
        plt.close(fig2)       
    
    
    bgcolor="white"
    fig2 = plt.figure(2,figsize=[8,6])
    gs = gridspec.GridSpec(1,2,width_ratios=[8,1])
    axis2 = plt.subplot(gs[0],facecolor=bgcolor)
    axisBxpGathered = plt.subplot(gs[1],facecolor=bgcolor,sharey=axis2)
    plotRuns.plot_one_curve(valuesItemsAggAllB[0], "blue", axis2, "Best", True)    
    plotRuns.plot_one_curve(valuesItemsAggAllR[0], "orange", axis2, "Random", True)    
    #plt.suptitle("Proportion of collected items over number of possible per generation",fontsize=labelFontSize)
    axis2.set_xlabel("Generations",fontsize=labelFontSize)
    #axis2.set_ylabel("% of collected items over possible",fontsize=labelFontSize)
    #axis2.grid(color="#FFFFFF", linestyle='-', linewidth=1)    
    axis2.tick_params(axis='both', which='major', labelsize=labelFontSize-3)    
    legend = axis2.legend(loc='upper left', shadow=True)
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width
    
    violin_parts = axisBxpGathered.violinplot([gatheredCumulatedB,gatheredCumulatedR], positions=[0.5,1.0],
                                              points=len(gatheredCumulatedB), widths=0.25,
                                              showmeans=False, showextrema=True, showmedians=True)  
    
    violin_parts['bodies'][0].set_facecolor((0.0,0.0,1.0,0.5))
    violin_parts['bodies'][0].set_edgecolor("blue")
    violin_parts['bodies'][1].set_facecolor((1.0,0.75,0.0,0.5))
    violin_parts['bodies'][1].set_edgecolor("orange")
    
    #print([[1.0,v] for v in fitnessCumulatedB])
    axisBxpGathered.plot([0.5 for v in gatheredCumulatedB], gatheredCumulatedB, 'x', color="blue",mew=3, ms=3)
    axisBxpGathered.plot([1.0 for v in gatheredCumulatedR], gatheredCumulatedR, 'x', color="orange", mew=3, ms=3)
        
    axisBxpGathered.set_xticks([0.5,1.0])
    axisBxpGathered.set_xlim((0.2,1.3))
    axis2.tick_params(axis='y', which='major', labelsize=labelFontSize-3)        
    #axisBxpGathered.set_xticklabels(["Best", "Random"],rotation=90,fontsize=12)    
    axisBxpGathered.set_xticklabels(["B", "R"],fontsize=12)    
    
    plt.tight_layout()    
    plt.subplots_adjust(top=0.93)  
               
    if isToPng:        
        print("Gathered Best vs. Random saved to ", outfileBest)
        time.sleep(2)
        fig2.savefig(outfileBest + "GatheredProportionsBoth.png"  , dpi=dpi * 3)        
        plt.close(fig2)       
    
    print(len(fitnessAllBPerRun))
    print(len(fitnessAllBPerRun[0]))
    exit
    #kendall with number of discordances, spearman with deviations
    spearmanBest = []
    for i in range(len(fitnessAllBPerRun)):    
        ranksChB = stats.rankdata(fitnessAllBPerRun[i])
        print("Changes size",len(ranksChB))
        #print(ranksChB)
        ranksFitB = stats.rankdata(colorChangesAvgAggGenerationB[i])
        print("Fit size",len(ranksFitB))
        #print(ranksFitB)
        spearmanBest.append([stats.spearmanr(ranksChB,ranksFitB[:-1])])
    
    spearmanRandom = []
    for i in range(len(fitnessAllRPerRun)):    
        ranksChR = stats.rankdata(fitnessAllRPerRun[i])
        ranksFitR = stats.rankdata(colorChangesAvgAggGenerationR[i])
        spearmanRandom.append([stats.spearmanr(ranksChR,ranksFitR[:-1])])
        
    bgcolor="white"
    figSpearman = plt.figure(3587)
    axisSpearman = plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)    
    
    #TODO 
    plotRuns.plot_one_curve([x[0] for x in spearmanBest], "blue", axisSpearman, "Best (ρ statistic)", True)
    #TODO  check following line
    plotRuns.plot_one_curve([x[0] for x in spearmanRandom], "orange", axisSpearman, "Random (ρ statistic)", True)            
    axisSpearman.set_title("Spearman-ρ coefficient between generations \n ranked by fitness and number of color changes")
    legend = axisSpearman.legend(loc='upper left', shadow=True)
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    
    for label in legend.get_texts():
        label.set_fontsize('large')
    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width
        
    #TODO axisSpearman.set_xlabel("Timesteps")
    #TODO axisSpearman.set_ylabel("% led color changes")
    if isToPng:        
        print("Spearman saved to ", outfileBest)
        #time.sleep(2)
        figSpearman.savefig(outfileBest + "SpearmanRho.png"  , dpi=dpi * 3)        
        plt.close(figSpearman)
    
    #TODO fitnessAllR
    #TODO colorChangesAvgAggB
    
    #TODO spearmanR.append(stats.spearmanr(x,y)) ranksG = stats.rankdata(data)
    #plt.plot([datum[1] for datum in spearmanR],label="Spearman-ρ p-val.")
    #plt.legend(loc='upper right', fontsize = 9)
    #plt.ylim(ymin=minPValRanks)
    #plt.grid(True,color=gridcolor)   
    #figureKendall.get_axes()[0].set_facecolor(bgcolor)

    
        
    normalizedItemsPerColorAllB = [np.sort(x/np.linalg.norm(x)) for x in itemsPerColorAllB]    
    
    correlationMatrixB = []
    for i,vals in enumerate(normalizedItemsPerColorAllB):
        line = []
        for j in range(i):
            line.append(np.dot(vals,normalizedItemsPerColorAllB[j]))
        correlationMatrixB.append(line)

    """    
    print("Correlation Matrix Best")
    for l in correlationMatrixB:
        for v in l:
            print("%.2f" % v, end='  ')
        print("")
    """
    flattenedCorrMatrixB = [item for sublist in correlationMatrixB for item in sublist]   

    itemsPerColorAllR = []
    for i,fname in enumerate(dirFilesItemsPerColorR):    
        with open(fname) as f:
            content = f.readlines()
        
        content = [x.strip() for x in content]
        
        contentSplitR = [x.split(" ") for x in content][1:]
        
        itemsPerColorRunR = [list(x) for x in zip(*contentSplitR)]

        itemsPerColorRunR = [np.sum([float(y) for y in x]) for x in itemsPerColorRunR]        
        itemsPerColorAllR.append(itemsPerColorRunR)        
        
        contentIncrR = []
        for it in contentSplitR:
            tmp = []
            aggreg = 0.0
            for val in it:
                aggreg += float(val)
                tmp += [aggreg]
            contentIncrR += [tmp]
        
        contentIncrR = [list(x) for x in zip(*contentIncrR)] # contentSplitR)]
        
        bgcolor="white"
        fig1 = plt.figure(1)
        axis1 = plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)
        for j,lColor in enumerate(contentIncrR):
            if j == 0:
                axis1.fill_between(np.arange(0, len(lColor)), np.zeros(len(lColor)) ,  lColor,
                                   alpha=0.55, linewidth=0, color=colorArr[j])
            else:
                axis1.fill_between(np.arange(0, len(lColor)), contentIncrR[j-1],  lColor,
                                   alpha=0.55, linewidth=0, color=colorArr[j])
        axis1.set_ylabel("#Items")  
        axis1.set_xlabel("Generations")  
        plt.title("Collected items per color over time of a Random run.")       
        if isToPng:        
            print("Items per color random saved to ", outfileRandom)
            #time.sleep(2)
            fig1.savefig(outfileRandom + "ItemsPerColorR" + str(i) + ".png"  , dpi=dpi * 3)        
            plt.close(fig1)      
             
    normalizedItemsPerColorAllR = [np.sort(x/np.linalg.norm(x)) for x in itemsPerColorAllR]    
    
    correlationMatrixR = []
    for i,vals in enumerate(normalizedItemsPerColorAllR):
        line = []
        for j in range(i):
            line.append(np.dot(vals,normalizedItemsPerColorAllR[j]))
        correlationMatrixR.append(line)
    """
    print("Correlation Matrix Random")
    for l in correlationMatrixR:
        for v in l:
            print("%.2f" % v, end='  ')
        print("")
    """
    flattenedCorrMatrixR = [item for sublist in correlationMatrixR for item in sublist]   
    print(flattenedCorrMatrixB)    
    print(np.median(flattenedCorrMatrixB))
    print(np.percentile(flattenedCorrMatrixB,25))
    print(np.percentile(flattenedCorrMatrixB,75))
    figBxpCorr = plt.figure(13322)
    axisBxpCorr = plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)

    violin_parts = axisBxpCorr.violinplot([flattenedCorrMatrixB, flattenedCorrMatrixR], 
                                          points=len(flattenedCorrMatrixR), widths=0.8,
                                          showmeans=False, showextrema=True, showmedians=True)  
    axisBxpCorr.set_ylim([0.0,1.0])
    violin_parts['bodies'][0].set_facecolor((0.0,0.0,1.0,0.5))
    violin_parts['bodies'][0].set_edgecolor("blue")
    violin_parts['bodies'][1].set_facecolor((1.0,0.75,0.0,0.5))
    violin_parts['bodies'][1].set_edgecolor("orange")
    
    axisBxpCorr.plot([1.0 for v in flattenedCorrMatrixB],flattenedCorrMatrixB, 'x', color="blue", mew=3, ms=3)
    axisBxpCorr.plot([2.0 for v in flattenedCorrMatrixR],flattenedCorrMatrixR, 'x', color="orange", mew=3, ms=3)
    
    plt.title("Collinearity measures of proportion of items per color",  fontsize=15)
    
    axisBxpCorr.set_xticks([1,2])
    axisBxpCorr.set_xticklabels(["Best", "Random"]) 
    for tick in axisBxpCorr.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    #axisBxpCorr.set_xlabel("")
    #axisBxpCorr.set_ylabel("")     
    
    plt.show()
    if isToPng:        
        print("Correlation values of gathered items between runs saved to ", outfileBest + "All")
        time.sleep(2)
        figBxpCorr.savefig(outfileBest + "CorrelationItemsPerColorBoth.png"  , dpi=dpi * 3)        
        plt.close(figBxpCorr)    
    
    itemsPerColorAllB = [list(x) for x in zip(*itemsPerColorAllB)]      
    itemsPerColorAllR = [list(x) for x in zip(*itemsPerColorAllR)]      
    
    figColorItems = plt.figure(122,figsize=(8, 6)) 
    gs = gridspec.GridSpec(1, 2, width_ratios=[8, 1]) 

    axisColorItems =  plt.subplot(gs[0],facecolor=bgcolor) # plt.subplot2grid((1, 2), (0, 0), facecolor=bgcolor,colspan=12)
    
    axisBxpEntropy =  plt.subplot(gs[1],facecolor=bgcolor) # plt.subplot2grid((1, 2), (0, 1), facecolor=bgcolor,colspan=1)    
    #print(itemsPerColorAll)
    pos = np.arange(1,len(itemsPerColorAllB), 1) #[1, 2, 4, 5, 7, 8]  #, pos
    #for colorNumbers in itemsPerColorAll:
    normalizeViolin = True
    if normalizeViolin:
        sumColorsB = []
        sumColorsR = []
        for x in [list(v) for v in zip(*itemsPerColorAllB)]:
            sumColorsB.append(np.sum(x))
        for x in [list(v) for v in zip(*itemsPerColorAllR)]:
            sumColorsR.append(np.sum(x))            
        
        itemsPerColorAllB = [[y/sumColorsB[i] for y in x] for i,x in enumerate([list(v) for v in zip(*itemsPerColorAllB)])]
        itemsPerColorAllB = [list(v) for v in zip(*itemsPerColorAllB)]
        itemsPerColorAllR = [[y/sumColorsR[i] for y in x] for i,x in enumerate([list(v) for v in zip(*itemsPerColorAllR)])]
        itemsPerColorAllR = [list(v) for v in zip(*itemsPerColorAllR)]
    
    itemsPerColorRunsB = [list(v) for v in zip(*itemsPerColorAllB)]            
    entropyB = []
    for p in itemsPerColorRunsB:
        entropyB.append(-np.sum(p * np.log2(p)))
        
    itemsPerColorRunsR = [list(v) for v in zip(*itemsPerColorAllR)]            
    entropyR = []
    for p in itemsPerColorRunsR:
        entropyR.append(-np.sum(p * np.log2(p)))    
        
    
    itemsPerColorAllAlternate = []
    for i in range(len(itemsPerColorAllB)):
        itemsPerColorAllAlternate.append(itemsPerColorAllB[i])
        itemsPerColorAllAlternate.append(itemsPerColorAllR[i])        

    pvaluesColorDistr = []     
    for i in range(int(len(itemsPerColorAllAlternate)/2)):
        pvaluesColorDistr.append(stats.mannwhitneyu(itemsPerColorAllAlternate[i * 2], itemsPerColorAllAlternate[i * 2 + 1]))
    print("Pvalues mannwhit color proportions")
    for i,pval in enumerate(pvaluesColorDistr):
        print(str(i) + " color proportion, Best vs. Random Mann-Whitney pval: " + str(pval[1]))
    #TODO test if proportion around a value     (12.5) =>entropy see below
    positionsByPairs = []
    for i in range(len(itemsPerColorAllB)):
        positionsByPairs.append(1.25 + i * 2)
        positionsByPairs.append(1.75 + i * 2)
    violin_parts = axisColorItems.violinplot(itemsPerColorAllAlternate, positions=positionsByPairs,
                                             points=2 * len(itemsPerColorAllB[0]), 
                                             widths=0.4, showmeans=False, 
                                             showextrema=True, showmedians=True) # , color = colorArr)    
    
    for i,pc in enumerate(violin_parts['bodies']):
        pc.set_facecolor(list(colorArr[int(i / 2)]) + [0.5])#'red')
        pc.set_edgecolor(colorArr[int(i / 2)]) #'black')           
        
    #labels = [item.get_text() for item in axisColorItems.get_xticklabels()] 
    
    colorLabels = np.arange(-1,1.0,0.25)    
    
    #for xpos in pos:                
    #   labels[xpos] = colorLabels[xpos % len(itemsPerColorAllB)]
    positions = np.arange(1.5, 16, 2)
    labels = []
    for i,v in enumerate(colorLabels):
        if i < (len(colorLabels) - 1):
            labels.append("[" + str(v) + "," + str(colorLabels[i+1]) + ")" )
        else:
            labels.append("[" + str(v) + ",1.0)")
    #axisColorItems.set_xticks(positions)
    #axisColorItems.set_xticklabels(labels)
    for i,colorDistr in enumerate(itemsPerColorAllAlternate):
        axisColorItems.plot([positionsByPairs[i] for j in colorDistr],colorDistr, 'x', color=colorArr[int(i / 2)], mew=3, ms=3)
                
    axisColorItems.set_ylim((0.0,1.0))
    #plt.suptitle("Proportion of collected items per color",fontsize=labelFontSize)
    axisColorItems.set_xlabel(#"Value of i" + \
                               "Item color (Best and Random)",fontsize=labelFontSize)
    #axisColorItems.set_ylabel("Proportion of collected items",fontsize=labelFontSize)  
    #axisColorItems.grid(axis='y',color="#FFFFFF", linestyle='-', linewidth=1)
    axisColorItems.xaxis.set_ticks_position('none') 
    axisColorItems.set_xticklabels([])
    #plt.tick_params(axis='x',          # changes apply to the x-axis
    #                which='both',      # both major and minor ticks are affected
    #                bottom='off',      # ticks along the bottom edge are off
    #                top='off',         # ticks along the top edge are off
    #                labelbottom='off')
    #for tick in axisColorItems.xaxis.get_major_ticks():
    #    tick.label.set_fontsize(7)
    #plt.show()    
    #if isToPng:        
    #    print("Items per color violins saved to ", outfileBest)
    #    time.sleep(2)
    #    figColorItems.savefig(outfileBest + "ViolinItemsPerColorBoth.png"  , dpi=dpi*3)        
    #    plt.close(figColorItems)       

    mannWhitEntropy = stats.mannwhitneyu(entropyB,entropyR)
    print("P-val mann whitney U Entropy proportion colors best vs. random" + str(mannWhitEntropy[1]))
 