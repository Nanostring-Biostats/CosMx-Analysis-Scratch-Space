"""
Created on  29 Feb 2024 for RnD use only
@author: mkorukonda
"""

import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import sys
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from typing import Tuple
paper_size_w = 11.7 # A4 paper width

                                                       
# # #  https://stackoverflow.com/questions/71517614/set-absolute-size-of-subplots
def generate_figsize_based_subplots(num_subplots, xFovs, yFovs):

    yscale = yFovs/xFovs # Ydimension scaled as function Xdimension

    nRows = num_subplots[0]
    nCols = num_subplots[1]


    fig_width = paper_size_w # Fixed Page width
    subplot_abs_spacing_width = 0.1*paper_size_w # The width of the spacing between subplots
    subplot_abs_excess_width = 0.2*paper_size_w # The width of the excess space on the left and right of the subplots
    subplot_abs_width = (fig_width-subplot_abs_excess_width-(nCols-1)*subplot_abs_spacing_width)/nCols

    subplot_abs_height = yscale * subplot_abs_width # Scaled by height of the sample
    subplot_abs_spacing_height = subplot_abs_spacing_width # The height of the spacing between subplots
    subplot_abs_excess_height = subplot_abs_excess_width # The height of the excess space on the top and bottom of the subplots    
    fig_height = (nRows * subplot_abs_height) + (nRows-1) * subplot_abs_spacing_height +subplot_abs_excess_height

    msize = max(5, round(500*subplot_abs_width/xFovs))

    #print(subplot_abs_width, xFovs, subplot_abs_width/xFovs)
    return fig_width, fig_height, msize
    

def ScatterSubPlot1D(fig, axs, sp, x,y,z,msize, title,isColorbar = True):
    p = axs[sp].scatter(x, y, s=msize, c=z, marker='s', facecolors='none')
    axs[sp].set_title(title)
    axs[sp].set_aspect('equal')
    axs[sp].invert_xaxis()
    axs[sp].invert_yaxis()
    divider = make_axes_locatable(axs[sp])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if(isColorbar):  
        plt.colorbar(p, cax=cax)


def ScatterSubPlot2D(fig, axs,spi,spj,x,y,z,msize, title,  clim = None, colormap = None, isColorbar = True):

    if(clim==None):
        p = axs[spi,spj].scatter(x, y, s=msize, c=z, marker='s', facecolors='none')
    else:
        p = axs[spi,spj].scatter(x, y, s=msize, c=z,vmin = clim[0], vmax = clim[1], marker='s', facecolors='none')

    if not (colormap is None):
        p.set_cmap(colormap)

    axs[spi,spj].set_title(title)
    axs[spi,spj].set_aspect('equal')
    axs[spi,spj].invert_xaxis()
    axs[spi,spj].invert_yaxis()
    axs[spi,spj].set_xticks([])
    axs[spi,spj].set_yticks([])
    
    if(isColorbar):
        divider = make_axes_locatable(axs[spi,spj])
        cax = divider.append_axes("right", size="5%", pad=0.05)  
        plt.colorbar(p, cax=cax)
    return p
    
    # Decide the rows and cols of subplots
def SubPlotConfig(nSpots):
    #1K
    if (nSpots ==16): 
        sr = 4
        sc = 4
    #6K
    elif(nSpots ==27):
        sr = 4
        sc = 7
    #WTx
    elif (nSpots ==39):
        sr = 5
        sc = 8
    else:
        sr = 4     
        sc = math.ceil(nSpots/sr)
    return sr, sc

def SubPlotSize(x,y):
    dx = max(x) - min(x)
    dy = max(y) - min(y)
    fov_size = 0.5 # 0.5mm x0.5mm
    xFovs = math.floor(dx/fov_size)+1
    yFovs = math.floor(dy/fov_size)+1
   
    return xFovs, yFovs

def SavePDFFile(fig, axs, plotFile,  globalcb=False, im=None):   
    fig.tight_layout()       
    if(globalcb): # Add a global colorbar 
        # fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.99)     
        fig.subplots_adjust(right=0.9)
        cbar_ax = fig.add_axes([0.95, 0.02, 0.02, 0.95])
        fig.colorbar(im, cax=cbar_ax)
    plt.gcf().savefig(plotFile, dpi=150, bbox_inches='tight', pad_inches=0.2)


def Plot_FovLoc_AcqOrder_FocalZ(FOVLocFile,PlotDir):

    #  Read FOV info 
    df = pd.read_csv(FOVLocFile)
    x = df['X_mm']
    y = df['Y_mm']
    z = df['Z_um']
    Fov = df['FOV']
    if 'Order' in df:
        Order = df['Order']
    else:
        Order = df['FOV']
    [xFovs, yFovs] = SubPlotSize(x,y)


    # Estimate Size of the plots         
    num_plots = [1,3]
    [fig_width, fig_height, msize] = generate_figsize_based_subplots(num_plots, xFovs,yFovs)


    # Plot FOV info
    plotFile = os.path.join(PlotDir,"FOV_Locations_AcqOrder_FocalZ.pdf")
    fig, axs = plt.subplots(1, 3, figsize=(fig_width, fig_height))
    fig.suptitle('Experiment Acquisition Stats')
   
    ScatterSubPlot1D(fig, axs, 0, x, y, Fov, msize, 'FOVs Order')
    ScatterSubPlot1D(fig, axs, 1, x, y, Order, msize, 'Acquisition Order')
    ScatterSubPlot1D(fig, axs, 2, x, y, z, msize, 'Focal Z (um)')
 
    
    # #plt.show()
    # fig1 = plt.gcf()
    # fig1.savefig(plotFile, dpi=150, bbox_inches='tight', pad_inches=0.2)

    SavePDFFile(fig, axs, plotFile,  False)    



def Plot_RegQ_perFOV_Metrics(x,y,RegQFile, PlotDir):
    df = pd.read_csv(RegQFile, skiprows=6); 
    score = df['Score']
    dx_sd = df['X-StdDev']
    dy_sd = df['Y-StdDev']
    dz_sd = df['Z-StdDev']
    ncc_mu = df['NCC-Mean']
    ncc_sd = df['NCC-StdDev']

    # Estimate Size of the plots         
    [xFovs, yFovs] = SubPlotSize(x,y)

   # Estimate Size of the plots         
    num_plots = [2,3]
    [fig_width, fig_height, msize] = generate_figsize_based_subplots(num_plots, xFovs,yFovs)


    # Plot FOV info
    plotFile = os.path.join(PlotDir,"FOV_RegistrationStats.pdf")

    # Default Fig size [6.4, 4.8]
    fig, axs = plt.subplots(2, 3, figsize=(fig_width, fig_height))
    fig.suptitle('Experiment Registration Stats')

    ScatterSubPlot2D(fig, axs, 0,0, x, y, dx_sd, msize, 'X-shift Std Dev')
    ScatterSubPlot2D(fig, axs, 0,1, x, y, dy_sd, msize, 'Y-shift Std Dev')
    ScatterSubPlot2D(fig, axs, 0,2, x, y, dz_sd, msize, 'Z-shift Std Dev')
    ScatterSubPlot2D(fig, axs, 1,0, x, y, score, msize, 'Final Score',[0,1],plt.cm.RdYlGn)
    ScatterSubPlot2D(fig, axs, 1,1, x, y, ncc_mu, msize, 'NCC Average',[0,100])
    ScatterSubPlot2D(fig, axs, 1,2, x, y, ncc_sd, msize, 'NCC Std Dev')

    SavePDFFile(fig, axs,  plotFile,  False)    

def Plot_XYReg_PerFOV(x,y,RegQFile, PlotDir):

        # Estimate Size of the plots         
    [xFovs, yFovs] = SubPlotSize(x,y)

    # Open the file for reading
    file = open(RegQFile, 'r')

    # Read the first line of the file
    first_line = file.readline().split(",")
    nSpots = int(first_line[1])

    # Read the first line of the file
    second_line = file.readline().split(",")
    nCycles = int(second_line[1])
    
    # Close the file
    file.close()
    [sr,sc] = SubPlotConfig(nSpots)

   # Estimate Size of the plots         
    num_plots = [sr,sc]
    [fig_width, fig_height, msize] = generate_figsize_based_subplots(num_plots, xFovs,yFovs)


    df = pd.read_csv(RegQFile, skiprows=6); 

    n_index = df.columns.get_loc('PerReporter')+1
    c_index = df.columns.get_loc('PerCycle')+1
    
    # Plot Score per Reporter
    plotFile = os.path.join(PlotDir,"XYRegistration_PerReporter.pdf")
    fig, axs = plt.subplots(sr, sc, figsize=(fig_width, fig_height))
    fig.suptitle('XY Registration Per Reporter')
    
    counter=0
    for cycle in range(nSpots):
        nr = math.floor(cycle/sc)
        nc = cycle % sc
        spot = 2*cycle+1
        N = 'N' + str(spot).zfill(2)
        dxy = df.iloc[:,cycle+n_index]
        im = ScatterSubPlot2D(fig,axs,nr,nc,x,y,dxy,msize, N,[0,1],plt.cm.RdYlGn, False)
        counter=counter+1

    # Blank the remaining axes    
    for n in range(counter,sr*sc):
        nr = math.floor(n/sc)
        nc = n % sc
        axs[nr,nc].set_visible(False)  

    # fig.tight_layout()    
    # fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.99)
    # #plt.show()
    # fig1 = plt.gcf()
    # fig1.savefig(plotFile, dpi=150, bbox_inches='tight', pad_inches=0.2)
    SavePDFFile(fig, axs, plotFile,  True, im)    



    # Plot Score per Cycle
    sr = 2
    sc = math.ceil(nCycles/sr)

    # Estimate Size of the plots         
    num_plots = [sr,sc]
    [fig_width, fig_height, msize] = generate_figsize_based_subplots(num_plots, xFovs,yFovs)

    plotFile = os.path.join(PlotDir,"XYRegistration_PerCycle.pdf")
    fig, axs = plt.subplots(sr, sc, figsize=(fig_width, fig_height))
    fig.suptitle('XY Registration Per Cycle')
    
    counter=0
    for cycle in range(nCycles):
        nr = math.floor(cycle/sc)
        nc = cycle % sc
        N = 'C' + str(cycle+1).zfill(1)
        dxy = df.iloc[:,cycle+c_index]
        im = ScatterSubPlot2D(fig,axs,nr,nc,x,y,dxy,msize, N,[0,1],plt.cm.RdYlGn, False)
        counter=counter+1

    # Blank the remaining axes    
    for n in range(counter,sr*sc):
        nr = math.floor(n/sc)
        nc = n % sc
        axs[nr,nc].set_visible(False)  

    SavePDFFile(fig, axs, plotFile, True, im)



def Plot_DZ_perSpot(x,y,dzFile, PlotDir):
    plotFile = os.path.join(PlotDir,"ZShift_PerReporter.pdf")
    df = pd.read_csv(dzFile)
    df = df.drop('DZ', axis=1)
    df = df.drop('X_mm', axis=1)
    df = df.drop('Y_mm', axis=1)
    df = df.abs()
    nSpots = df.shape[1]
    [sr,sc] = SubPlotConfig(nSpots)    

    # Estimate Size of scanned sample    
    [xFovs, yFovs] = SubPlotSize(x,y)


    # Estimate Size of the plots         
    num_plots = [sr,sc]
    [fig_width, fig_height, msize] = generate_figsize_based_subplots(num_plots, xFovs,yFovs)

    dzmax = df.to_numpy().max()

    fig, axs = plt.subplots(sr, sc, figsize=(fig_width, fig_height))
    fig.suptitle('ZShift Distribution Per Reporter; Max Shift:' + str(dzmax))

    counter=0
    for n in range(nSpots):
        nr = math.floor(n/sc)
        nc = n % sc
        spot = 2*n+1
        N = 'N' + str(spot).zfill(2)
        dz = df.iloc[:,n]
        im = ScatterSubPlot2D(fig,axs,nr,nc,x,y,dz,msize, N,[0,dzmax], None, False)
        counter=counter+1

    # Blank the remaining axes    
    for n in range(counter,sr*sc):
        nr = math.floor(n/sc)
        nc = n % sc
        axs[nr,nc].set_visible(False)  
 
    SavePDFFile(fig, axs, plotFile, True, im)


def Plot_Counts_perSpot(x,y,spotFile, PlotDir):
    cols = list(pd.read_csv(spotFile, nrows=1))
    noCols = len(cols)
    df = pd.read_csv(spotFile,usecols=[i for i in range(3,noCols)]) # Skip FOV,X_mm,Y_mm Columns
    
  #  df = df.drop('FOV vs Cnt (M)', axis=1)
  #  df = df.drop('X_mm', axis=1)
  #  df = df.drop('Y_mm', axis=1)

    nSpots = int(df.shape[1]/4) # Assume BGYR  
    noFov = df.shape[0]
    [sr,sc] = SubPlotConfig(nSpots)    

    # Estimate Size of scanned sample    
    [xFovs, yFovs] = SubPlotSize(x,y)


    # Estimate Size of the plots         
    num_plots = [sr,sc]
    [fig_width, fig_height, msize] = generate_figsize_based_subplots(num_plots, xFovs,yFovs)

    # # # Calculate total spots for each reporter
    CntPerSpot = np.zeros((noFov, nSpots))
    for n in range(nSpots):
        for ch in range(4):
            CntPerSpot[:,n]+= df.iloc[:,n+ch]

    spotMax = CntPerSpot.max()
    spotMin = CntPerSpot.min()
    plotFile = os.path.join(PlotDir,"SpotCounts_All_Per_Reporter.pdf")
    fig, axs = plt.subplots(sr, sc, figsize=(fig_width, fig_height))
    fig.suptitle('Counts Per Reporter (Millions); Max Count:' + "{:.2f}".format(spotMax) + "M")

    counter=0
    for n in range(nSpots):
        nr = math.floor(n/sc)
        nc = n % sc
        spot = 2*n+1
        N = 'N' + str(spot).zfill(2)
        im = ScatterSubPlot2D(fig,axs,nr,nc,x,y,CntPerSpot[:,n],msize, N,[spotMin,spotMax],plt.cm.jet, False)
        counter=counter+1

    # Blank the remaining axes    
    for n in range(counter,sr*sc):
        nr = math.floor(n/sc)
        nc = n % sc
        axs[nr,nc].set_visible(False)  
 
    SavePDFFile(fig, axs, plotFile, True, im)
    
    # # # Calculate spots per Channel for each reporter
    rrcode=["BB","GG","YY","RR"]
    cbar_choice = [plt.cm.Blues, plt.cm.Greens, plt.cm.YlOrBr, plt.cm.Reds]
    for ch in range(4):
        plotFile = os.path.join(PlotDir,"SpotCounts_"+rrcode[ch] +"_Per_Reporter.pdf")
        CntPerSpot = np.zeros((noFov, nSpots))
        for n in range(nSpots):
            CntPerSpot[:,n] = df.iloc[:,n*4+ch]

        spotMax = CntPerSpot.max()
        spotMin = CntPerSpot.min()
        fig, axs = plt.subplots(sr, sc, figsize=(fig_width, fig_height))
        fig.suptitle(rrcode[ch] +' Counts Per Reporter (Millions); Max Count:' +"{:.2f}".format(spotMax) + "M")

        counter=0
        for n in range(nSpots):
            nr = math.floor(n/sc)
            nc = n % sc
            spot = 2*n+1
            N = 'N' + str(spot).zfill(2)
            #im = axs[nr,nc].scatter(x, y, c = CntPerSpot[:,n], vmin = spotMin, vmax = spotMax, marker = 's', facecolors ='none', cmap=cbar_choice[ch])
            im = ScatterSubPlot2D(fig,axs,nr,nc,x,y,CntPerSpot[:,n],msize, N,[spotMin,spotMax],cbar_choice[ch],False)
            counter=counter+1

        # Blank the remaining axes    
        for n in range(counter,sr*sc):
            nr = math.floor(n/sc)
            nc = n % sc
            axs[nr,nc].set_visible(False)  
        
        SavePDFFile(fig, axs, plotFile, True, im)

def Plot_RawInt_perSpot(x,y,intFile, PlotDir):
    cols = list(pd.read_csv(intFile, nrows=1))
    noCols = len(cols)
    df = pd.read_csv(intFile,usecols=[i for i in range(3,noCols)]) # Skip FOV,X_mm,Y_mm Columns


    nSpots = int(df.shape[1]/4) # Assume BGYR  
    noFov = df.shape[0]
    [sr,sc] = SubPlotConfig(nSpots)    

    [xFovs, yFovs] = SubPlotSize(x,y)
    # Estimate Size of the plots         
    num_plots = [sr,sc]
    [fig_width, fig_height, msize] = generate_figsize_based_subplots(num_plots, xFovs,yFovs)

    # # # Plot average spot intensity per Channel for each reporter
    rrcode=["BB","GG","YY","RR"]
    cbar_choice = [plt.cm.Blues, plt.cm.Greens, plt.cm.YlOrBr, plt.cm.Reds]
    for ch in range(4):
        plotFile = os.path.join(PlotDir,"SpotChIntensity_"+rrcode[ch] +"_Per_Reporter.pdf")
        IntPerSpot = np.zeros((noFov, nSpots))
        for n in range(nSpots):
            IntPerSpot[:,n] = df.iloc[:,n*4+ch]

        intMax = IntPerSpot.max()
        intMin = IntPerSpot.min()
        fig, axs = plt.subplots(sr, sc, figsize=(fig_width, fig_height))
        fig.suptitle(rrcode[ch] +' Average Spot Channel Intensity; Max Int:' +"{:.2f}".format(intMax))
        fig.tight_layout()

        counter=0
        for n in range(nSpots):
            nr = math.floor(n/sc)
            nc = n % sc
            spot = 2*n+1
            N = 'N' + str(spot).zfill(2)
            im = ScatterSubPlot2D(fig,axs,nr,nc,x,y,IntPerSpot[:,n],msize, N,[intMin,intMax],cbar_choice[ch], False)
            counter=counter+1

        # Blank the remaining axes    
        for n in range(counter,sr*sc):
            nr = math.floor(n/sc)
            nc = n % sc
            axs[nr,nc].set_visible(False)  

        SavePDFFile(fig, axs, plotFile, True, im)


def Plot_All_RunQC_Metrics(RunSumDir):

    PlotDir = os.path.join(RunSumDir,'QC_Plots')
    if not os.path.exists(PlotDir):
        os.makedirs(PlotDir)

    # Create the list of all your folder content with os.listdir()
    folder_content = os.listdir(RunSumDir) 

    # # # Search for the FOV Locations file in the folder content
    FOVLocFile =""
    for file in folder_content:
        if file.endswith('_FOV_Locations.csv'):
            FOVLocFile = os.path.join(RunSumDir,file)
            break
    
    if(FOVLocFile==""):
        print("No FOV Locations file found, returning...")
        return

    # # #  Read FOV info 
    df = pd.read_csv(FOVLocFile)
    x = df['X_mm']
    y = df['Y_mm']
    Fov = df['FOV']
    print('Load FOV locations')

    Plot_FovLoc_AcqOrder_FocalZ(FOVLocFile,PlotDir)
    print('Completed FOV Loc Plots')



    # # # # Search for the regQ per FOV in the folder content
    regQFile =""
    for file in folder_content:
        if file.endswith('_Quality_perFOV.csv'):
            regQFile = os.path.join(RunSumDir,file)
            break
                         
    if(regQFile==""):
        print("No Quality_perFOV file found, skpping RegStats Plots...")           
    else:
        Plot_RegQ_perFOV_Metrics(x,y,regQFile, PlotDir)
        print('Completed RegStats Plots')

        # # # Plot the RegScore function of reporter and cycle
        Plot_XYReg_PerFOV(x,y,regQFile, PlotDir)
        print('Completed RegScore/Reporter/Cycle Plots')



    # # # Plot the variance in z-shift as a function of reporter
    dzFile =""
    for file in folder_content:
        if file.endswith('_Z_perSpotperFOV.csv'):
            dzFile = os.path.join(RunSumDir,file)
            break

    if(dzFile==""):
        print("No _Z_perSpotperFOV file found, skpping DzStats per Reporter Plots...")          
    else:
        Plot_DZ_perSpot(x,y,dzFile, PlotDir)
        print('Completed DzStats per Reporter Plots')

    # # # # Plot the variance in Spot Counts as a function of reporter
    spotFile =""
    for file in folder_content:
        if file.endswith('_SpotCount_perFOV.csv'):
            spotFile = os.path.join(RunSumDir,file)
            break
    if(spotFile==""):
        print("No _SpotCount_perFOV file found, skpping Spot Counts per Reporter Plots...")          
    else:
        Plot_Counts_perSpot(x,y,spotFile, PlotDir)
        print('Completed Spot Counts per Reporter Plots')
    
    # # # # Plot the variance in Average Channel Intensity as a function of reporter
    intFile =""
    for file in folder_content:
        if file.endswith('_ChIntensity_Spot_perFOV.csv'):
            intFile = os.path.join(RunSumDir,file)
            break
    if(intFile==""):
        print("No _ChIntensity_Spot_perFOV file found, skipping Average Spot Channel Intensity per Reporter Plots...")          
    else:
        Plot_RawInt_perSpot(x,y,intFile, PlotDir)
        print('Completed Average Channel Intensity per Reporter Plots')

def main():
    RunSumDir = sys.argv[1]
    Plot_All_RunQC_Metrics(RunSumDir)

if __name__ == "__main__":
    main()                       

# RunSumDir=r'E:\Data\RunQC\Run_cee38acc-ce36-4450-a190-fc3e545b7ed2_20240320_030317_S1_2208HB017\RunSummary'
# Plot_All_RunQC_Metrics(RunSumDir)