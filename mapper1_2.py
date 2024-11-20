# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 12:03:18 2021

@author: Piotr
"""

import sys
import os
import subprocess
import numpy as np
import math as m
from shutil import copyfile
import re
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import more_itertools as mit


#SORTAV functions


def sortav(hklspath, sortav_path, sortav_inp_path):
    '''
    Run a sortav routine on all hkls in the selected directory.
    path hklspath :: path to all hkl files
    path sortav_path :: path to sortav.exe file
    path sortav_inp_path :: path to template sortav.inp file
    '''
    create_dirs(hklspath)
    try:
        os.mkdir(os.path.join(hklspath, 'output'))
    except:
        pass
    for i in os.listdir(hklspath):
        if os.path.isdir(os.path.join(hklspath,i)) and i != 'output':
            wpath = os.path.join(hklspath,i)
            copyfile(sortav_path, os.path.join(wpath, 'sortav.exe'))
            f = open(sortav_inp_path)
            lines = f.readlines()
            lines[5] = i + '.hkl\n'
            lines[8] = i + '.sortav\n'
            f.close()
            f = open(os.path.join(wpath, 'sortav.inp'), '+w')
            for j in lines:
                f.write(j)
            f.close()
            curdir = os.getcwd()
            os.chdir(wpath)
            try:
                subprocess.run([os.path.join(wpath, 'sortav.exe'), os.path.join(wpath, 'sortav.inp')], timeout = 10)
            except:
                print(str(i) + '   HAD PROBLEMS!')
            os.chdir(curdir)
            f = open(os.path.join(wpath, i + '.sortav'))
            lines = f.readlines()
            f.close()
            lines = lines[7:]
            splitlines = []
            for j in lines:
                splitlines.append(list(j.split()))
            f = open(os.path.join(wpath, i + '_sortaved.hkl'), '+w')
            for j in splitlines:
                f.write(str(int(j[0])) + ' ')
                f.write(str(int(j[1])) + ' ')
                f.write(str(int(j[2])) + ' ')
                f.write(str(float(j[3])) + ' ')
                f.write(str(float(j[4])))
                f.write('\n')
            f.close()
            copyfile(os.path.join(wpath, i + '_sortaved.hkl'), os.path.join(os.path.join(hklspath, 'output'), i + '_sortaved.hkl'))
            copyfile(os.path.join(wpath, 'sortav.lst'), os.path.join(os.path.join(hklspath, 'output'), i + '.lst'))

#MAKEMAP functions:

def readhkl(workpath):
    '''
    finds and imports a hkl file in the workpath
	Usage::
    path workpath :: path to the folder containing .hkl file
    '''
    for path in os.listdir(workpath):
        FullPath = os.path.join(workpath,path)
        if os.path.isfile(FullPath) and FullPath[-4:] == '.hkl':
            hkl = open(FullPath, 'r')
            hkllist = np.loadtxt(FullPath, dtype = np.float_)
            hkl.close()
    return hkllist

def readm80(workpath):
    '''
    finds and imports an .m80 file in the workpath
	Usage::
    path workpath :: path to the folder containing .m80 file
    '''
    janalist = []
    for path in os.listdir(workpath):
        FullPath = os.path.join(workpath, path)
        if os.path.isfile(FullPath) and FullPath[-4:] == '.m80':
            jana = open(FullPath, 'r')
            while line:= jana.readline():
                spl = list(mit.split_into(line,[4, 4, 4, 4, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]))
                nspl = list(float(''.join(i)) for i in spl)
                janalist.append(nspl)
            jana.close()
    return np.array(janalist)

def readm50(workpath):
    '''
    Reads and returns cell parameters from the m50 file
	Usage::
    path workpath :: path to the folder containing .m50 file
    '''
    cellparams = []
    line = []
    for path in os.listdir(workpath):
        FullPath = os.path.join(workpath,path)
        if os.path.isfile(FullPath) and FullPath[-4:] == '.m50':
            m50 = open(FullPath, 'r')
            lines = m50.readlines()
            line = lines[2]
            cellparams = list(line[5:].split(' '))
            m50.close()
    return np.array(cellparams, dtype=np.float_)

def sinthl(cell, h, k, l):
    '''
    Calculates and returns sin(theta)/lambda.
	Usage::
    int cell :: list of cell parameters
    int h,k,l :: hkl indices
    '''
    a = cell[0]
    b = cell[1]
    c = cell[2]
    alpha = m.radians(cell[3])
    beta = m.radians(cell[4])
    gamma = m.radians(cell[5])
    delta=1-m.cos(alpha)**2-m.cos(beta)**2-m.cos(gamma)**2+2*m.cos(alpha)*m.cos(beta)*m.cos(gamma)
    vol=a*b*c*m.sqrt(delta)
    TR = np.zeros(6)
    TR[0]=h**2*b**2*c**2*m.sin(alpha)**2
    TR[1]=k**2*a**2*c**2*m.sin(beta)**2
    TR[2]=l**2*a**2*b**2*m.sin(gamma)**2
    TR[3]=2*h*k*a*b*c**2*(m.cos(alpha)*m.cos(beta)-m.cos(gamma))
    TR[4]=2*h*l*a*c*b**2*(m.cos(alpha)*m.cos(gamma)-m.cos(beta))
    TR[5]=2*k*l*b*c*a**2*(m.cos(beta)*m.cos(gamma)-m.cos(alpha))
    normhsq=np.sum(TR)/(vol**2)
    normh=m.sqrt(normhsq)
    return normh/2

#Constructs an array of data for a mopro input
def mopro(templates_path, workpath):
    '''
    returns an array of data for mopro output
	Usage::
    path templates_path :: path containing .m50 and .m80 files
    path workpath :: path containing .hkl file
    '''
    moprolist = []
    hkl = readhkl(workpath)
    m80 = readm80(templates_path)
    cell = readm50(templates_path)
    for i in range (hkl.shape[0]):
        for j in range (m80.shape[0]):
            if hkl[i][0] == int(m80[j][0]) and hkl[i][1] == int(m80[j][1]) and hkl[i][2] == int(m80[j][2]):
                refl = []
                h = hkl[i][0]
                k = hkl[i][1]
                l = hkl[i][2]
                sin = sinthl(cell, h, k, l)
                FobsOFF = m80[j][4]
                FestON = m.sqrt(hkl[i][3])*FobsOFF
                phase = m.degrees(m.atan2(m80[j][8],m80[j][7]))
                refl.extend([h,k,l,sin,FestON,FobsOFF,phase,1.0,1.0,1.0,FestON,0])
                moprolist.append(refl)
    return np.array(moprolist, dtype = np.float_)

def finalize_fourier (templates_path, workpath):
    '''
    creates .four file in the workpath, based on files in templates_path.
	Usage::
    path templates_path :: path containing .m50 and .m80 files
    path workpath :: path containing .hkl file
    '''
    result = mopro(templates_path, workpath)
    pdif = open(os.path.join(workpath, str(os.path.split(workpath)[1]) + '.four'), 'w+')
    for i in range (result.shape[0]):
            pdif.write('{:>3}'.format('%.0f'% result[i,0]) + '{:>4}'.format('%.0f'% result[i,1]) + '{:>4}'.format('%.0f'% result[i,2]) + '{:>6}'.format('%.3f'% result[i,3]) + '{:>10}'.format('%.3f'% result[i,4]) + '{:>10}'.format('%.3f'% result[i,5]) + '{:>7}'.format('%.1f'% result[i,6]) + '{:>7}'.format('%.3f'% result[i,7]) + '{:>6}'.format('%.3f'% result[i,8]) + '{:>7}'.format('%.3f'% result[i,9]) + '{:>9}'.format('%.2f'% result[i,10]) + '{:>2}'.format('%.0f'% result[i,11]) + '\n')
    pdif.close()

def change_mp_cwd(workpath, mp_ini_path):
    '''
    This function changes the mopro.ini file, sets the working directory of mopro to a desired path specified in the input path
    Usage::
    change_mp_cwd(path 'C:\\processing\\AgCu_02', path 'C:\\mopro\\mopro.ini')
    '''
    f = open(mp_ini_path) 
    lines = f.readlines()
    f.close()
    o = open(mp_ini_path, '+w')
    o.write('FILE  DIR  ')
    o.write(workpath)
    o.write('\n')
    for i in lines[1:]:
        o.write(i)
    o.close()

def prep_mp_input(workpath, temp_path):
    '''
    creates a vmopro_mv.inp file with info taken from template located at temp_path and corrected for the directory name at workpath
    Usage::
    prep_mp_input(path workpath, path temp_path)
    '''
    parpath = []
    fourpath = []
    setname = os.path.split(workpath)[1]
    for i in os.listdir(workpath):
        if i[-4:] == '.par': 
            parpath.append(i)
        if i[-5:] == '.four': 
            fourpath.append(i)
    parpath = os.path.join(workpath,parpath[0])
    fourpath = os.path.join(workpath,fourpath[0])
    
    f = open(temp_path) 
    lines = f.readlines()
    f.close()
    
    lines[0] = str(parpath) + '\n'
    lines[2] = str(fourpath) + '\n'
    lines[-3] = str(setname) + '_PDplot.ps\n'
    lines[-4] = str(setname) + '_PDgrid\n'
    
    f = open(os.path.join(workpath,'vmopro_mv.inp'), 'w+')
    for i in lines:
        f.write(i)
    f.close()
    return (str(setname) + '_PDplot.ps', str(setname) + '_PDgrid')
    
def prep_mp_input2(workpath, temp_path, c1, c2, c3):
    '''
    creates a vmopro_mv2.inp file with info taken from template located at temp_path and corrected for the directory name at workpath. The desired contour levels can be ajdusted as follows:
    c1,c2,c3: Contour increment, minimal and maximal
    Usage::
    prep_mp_input(path workpath, path temp_path, double c1, double c2, double c3)
    '''
    parpath = []
    setname = os.path.split(workpath)[1]
    for i in os.listdir(workpath):
        if i[-4:] == '.par': 
            parpath.append(i)
    parpath = os.path.join(workpath,parpath[0])
    f = open(temp_path) 
    lines = f.readlines()
    f.close()
    lines[0] = str(parpath) + '\n'
    lines[4] = str(c1) + '\n'
    lines[5] = str(c2) + '\n'
    lines[6] = str(c3) + '\n'
    
    filename = str(setname) + 'PDplot_cont_' + str(c1) + '_' + str(c2) + '_'+ str(c3) 
    filename = filename.replace('-', 'n')
    filename = filename.replace('.', '')
    filename = filename + '.ps'
    
    lines[-3] = filename + '\n'
    lines[-4] = str(setname) + '_PDgrid\n'
    
    f = open(os.path.join(workpath,'vmopro_mv2.inp'), 'w+')
    for i in lines:
        f.write(i)
    f.close()
    return filename
    
def mp_run(mp_ini_path, mp_path, workpath, temp_path, temp_path2, c_inc, c_min, c_max):
    '''
    runs compiled VMoPro-win32.exe program from the mp_path location, creates maps based on templates
	Usage::
    path mp_ini_path :: path to the mopro.ini file - usually located at Users/user/AppData/Roaming/mopro
    path mp_path :: path to the compiled VMoPro-win32.exe file
    path workpath :: path to the folder containing .four file
    path temp_path, temp_path2 :: paths to template vmopro inputs
    double c_int, c_min, c_max :: photodifference contour paramaters for the map
    '''
    change_mp_cwd(workpath,  mp_ini_path)
    plotname, gridname = prep_mp_input(workpath, temp_path)
    filename = prep_mp_input2(workpath, temp_path2, c_inc, c_min, c_max)
    if not any(i[-3:] == '.ps' for i in os.listdir(workpath)):
        print (os.path.join(workpath, 'vmopro_mv.inp'))
        subprocess.run([mp_path, os.path.join(workpath, 'vmopro_mv.inp')], stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
    subprocess.run([mp_path, os.path.join(workpath, 'vmopro_mv2.inp')], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return (filename, gridname, plotname)
    
def create_dirs(hklspath):
    '''
    create a separate directory for each .hkl file in hklspath and move it there
    path hklspath :: path to all hkl files
    '''
    for i in os.listdir(hklspath):
        if i[-4:] == '.hkl':
            try:
                raw = os.path.splitext(i)[0]
                os.mkdir(os.path.join(hklspath, raw))
                copyfile(os.path.join(hklspath, i), os.path.join((os.path.join(hklspath, raw)), i))
            except:
                continue

def ProduceMaps(hklspath, mp_ini_path, template_path, temp_path, temp_path2, mp_path):
    '''
    main program creating maps, and copying them to Images folder, as well as copying grids to Grids folder
	Usage::
    path hklspath :: path to the folder containing ALL hkl files
    path mp_ini_path :: path to mopro.ini file, which is usually located at C:\\Users\\user\\AppData\\Roaming\\mopro\\mopro.ini
    path templates_path :: path containing .m50 and .m80 files and .par
    path temp_path :: path to the primary template MoPro input file, containing basic plot and grid creationMP routine
    path temp_path2 :: path to the secndary template MoPro input file, containing specific contour plot creation MP routine
    path mp_path :: path to the VMoPro-win32.exe file
    '''
    if len(sys.argv) == 4:
        c_inc = sys.argv[1]
        c_min = sys.argv[2]
        c_max = sys.argv[3]
    else:
        c_inc = 0.30
        c_min = -0.30
        c_max = 0.30
        print('Using default contours' + ' c_inc = ' + str(c_inc) + ' c_min = ' + str(c_min) + ' c_max = ' + str(c_max) + '\n')
    
    create_dirs(hklspath)
    try:
        os.mkdir(os.path.join(hklspath, 'Images'))
    except:
        pass
    try:
        os.mkdir(os.path.join(hklspath, 'Grids'))
    except:
        pass
    
    for i in os.listdir(hklspath):
        if i[-4:] == '.hkl':
            hkldir = os.path.splitext(i)[0]
            print (hkldir + '...')
            workpath = os.path.join(hklspath, hkldir)
            for i in os.listdir(template_path):
                if i[-4:] == '.par':
                    try:
                        copyfile(os.path.join(template_path,i),os.path.join(workpath,i))
                    except:
                        pass
            if (os.path.exists(os.path.join(hklspath, str(hkldir) + '.four' ))):
                copyfile(os.path.join(hklspath, str(hkldir) + '.four'), os.path.join(workpath, str(hkldir) + '.four' ))
            if (os.path.exists(os.path.join(workpath, str(hkldir) + '.four' )) != 1) :
                finalize_fourier(template_path, workpath)
                try:
                    finalize_fourier(template_path, workpath)
                except:
                    print(i + ' FAILED AT FOURIER')
                    continue
            filename, gridname, plotname = mp_run(mp_ini_path, mp_path, workpath, temp_path, temp_path2, c_inc, c_min, c_max)
            print('Created:   ' + filename)
            try:
                copyfile(os.path.join(workpath, filename), os.path.join(os.path.join(hklspath, 'Images'), os.path.split(filename)[1]))
                copyfile(os.path.join(workpath, gridname), os.path.join(os.path.join(hklspath, 'Grids'), os.path.split(gridname)[1]))
                copyfile(os.path.join(workpath, plotname), os.path.join(os.path.join(hklspath, 'Images'), os.path.split(plotname)[1]))
            except:
                continue

def ImportGrid(gridPath, pngPath, epsi = 0.0, showMax = 1):
    '''
    Imports the grid.txt file located at os.path GridPath and produces a 2D photodifference map at pngPath. Uses epsi as countour cutoff
	Usage::
    path gridPath :: path to grid.txt
    path pngPath :: path to desired output image
    double epsi :: contour cutoff 
	bool showMax :: point to the max point or not
    '''    
    #Import data into a data array of correct shape, and atom data into separate array
    f = open(gridPath)
    lines = f.readlines()
    f.close()
    lines = lines[2:-1]
    atomLines = []
    cordLines = []
    paramLines = []
    for i in lines:
        if i[1].isalpha():
            atomLines.append(i)
        elif(len(paramLines) == 0):
            paramLines.append(i)
        else:
            cordLines.append(i)
    x = int(paramLines[0].split('  ')[2])
    y = int(paramLines[0].split('  ')[1])
    data = []
    for i in cordLines:
        line = re.split('\s+', i)
        for j in line:
            if j != '':
                if abs(float(j)) < epsi:
                    data.append('0')
                else:
                    data.append(j)
    atmData = []
    for i in atomLines:
        line = re.split('\s+', i)
        atomFS = line[1][0]
        atomSS = line[1][1]
        if (atomFS == 'C' or atomFS == 'H') and atomSS.isdigit():
            continue
        atmData.append((line[1],line[2],line[3]))
    dataFloat = []
    for i in data:
        dataFloat.append(float(i))
        
    ##Check maximum value for later colormap normalization
    maxAbs= 0
    if max(dataFloat) > abs(max(dataFloat)):
        maxAbs = max(dataFloat)
    else:
        maxAbs = abs(min(dataFloat))
        
    ##Check and write SD
    f = open(str(pngPath[:-4]) + '_SD.txt', '+w')
    f.write(str(np.std(dataFloat))[:5])
    f.close()    
    
    ##Reshape and redefine data to be compatible with plt
    data = np.reshape(data, [x,y])
    data = data.astype('float64')
    
    ##Check maximum value and its position for later maximum identification
    cordsMax = (np.unravel_index(np.argmax(data), np.shape(data)))
    valMax = data[cordsMax]
	
    
    plt.box(False)
    cmap = plt.get_cmap('bwr')
    
    normalize =  matplotlib.colors.Normalize(vmin=-maxAbs, vmax=maxAbs)
    
    xData = []
    yData = []
    atmLabels = []
    for i in atmData:
        atmLabels.append(i[0])
        xData.append(int(float(i[1])*y))
        yData.append(int(float(i[2])*x))
    
    fig = plt.contour(data, 30, colors = 'black', linewidths = 0.1)
    fig4 = plt.imshow(data, interpolation = 'none', cmap=cmap, norm = normalize)
    fig2 = plt.scatter(xData, yData, s=10, c='black')
    
    x_axis = fig4.axes.get_xaxis()
    x_axis.set_visible(False)

    y_axis = fig4.axes.get_yaxis()
    y_axis.set_visible(False)
    
    for i, ad in enumerate(atmLabels):
        plt.annotate(ad, (xData[i]+2, yData[i]-2))
    
    ##Add maximum point
    if (showMax):
        xMax = cordsMax[1]
        yMax = cordsMax[0]
        fig3 = plt.scatter(xMax, yMax, s=0.15, c='green')
        plt.annotate(str(valMax), (xMax+2, yMax-2))
    
    plt.savefig(pngPath,  bbox_inches='tight', dpi=300)
    plt.close()

def PlotDirectory(epsi = 0.0, showMax = 1, pathGrids = os.getcwd()):
    '''
    Plots every grid in the selected directory
	Usage::
    path pathGrids :: path to all grids that are supposed to be sketched
    double epsi :: contour cutoff 
	bool showMax :: point to the max point or not
    '''
    for i in os.listdir(pathGrids):
            if i[-4:] == 'grid':           
                gridPath = os.path.join(pathGrids, i)
                epsiStr = str(epsi)
                plotName = '_'
                for i in epsiStr:
                    if i.isdigit():
                        plotName = plotName + i
                    else:
                        plotName = plotName + '_'
                plotPath = (gridPath + plotName + '.png')
                ImportGrid(gridPath, plotPath, epsi, showMax)

def rescale(hklPath, cutoff_low = 0.5, cutoff_high = 1.5):
    
    
    #Finding .hkl files and extracting all data into a list of numpy arrays
    AllSetsHkl = []
    AllFileNames = []
    for path in os.listdir(hklPath):
        FullPath = os.path.join(hklPath, path)
        if os.path.isfile(FullPath) and FullPath[-10:] == 'ratios.hkl':
            CurrentHkl = np.loadtxt(FullPath)
            AllSetsHkl.append(CurrentHkl)
            AllFileNames.append(FullPath)    
            
    #Extracting average data from numpy arrays one by one
    
    SetAvgEta = []
    for j in range(len(AllSetsHkl)):
        #single hkl file:
        hkl = AllSetsHkl[j]
        TotalSetEta = 0
        for i in range (hkl.shape[0]):
            if hkl[i][3] < cutoff_low or hkl[i][3] > cutoff_high:
                continue
            TotalSetEta+=np.abs(hkl[i][3] - 1)
        SetAvgEta.append(TotalSetEta/hkl.shape[0])
    AllAvgEta = np.average(SetAvgEta)
    
    #Write average eta info:
    avgfile = open (os.path.join(hklPath,'averages.txt'), 'w+')
    for i in range (len(SetAvgEta)):
        avgfile.write(str(AllFileNames[i]) + '\n' + str('%.4f'%SetAvgEta[i]) + '\n')
    avgfile.close()
    
    #Writing new rescaled .hkl files one by one
    Iteration = 0
    for path in os.listdir(hklPath):
        FullPath = os.path.join(hklPath, path)
        if os.path.isfile(FullPath) and FullPath[-10:] == 'ratios.hkl':
    #       Creating appropriate new file and opening stream
            FileName = os.path.split(FullPath)[1]
            NewFileName = FileName[:-4] + '_rescaled.hkl'
            NewPath = os.path.join(os.path.split(FullPath)[0], NewFileName)
            f = open(NewPath, 'w+')
    #        Rewriting the .hkl file with corrected ratio and sigma values
            CurrentHkl = np.loadtxt(FullPath)
            for i in range (CurrentHkl.shape[0]):
                if CurrentHkl[i][3] < cutoff_low or CurrentHkl[i][3] > cutoff_high:
                    continue
                for j in range (3):
                    f.write('{:>5}'.format(str(int(CurrentHkl[i][j]))))
                    f.write(' ')
    #           Calculating Scaled Ratios 
                RNScaled = CurrentHkl[i][3]
                RScaled = abs(1 + (RNScaled - 1)*(AllAvgEta/SetAvgEta[Iteration]))
                f.write('{:>15}'.format(str('%.7f'%RScaled)))
                f.write(' ')
    #           Calculating Scaled Sigmas
                SScaled = CurrentHkl[i][4]*(AllAvgEta/SetAvgEta[Iteration])
                f.write('{:>15}'.format(str('%.7f'%SScaled)))
                f.write('{:>5}'.format(str(1)))
                f.write('\n')
            Iteration += 1
            f.close()

def glue(hklsPath):
    
    out = open(os.path.join(hklsPath, 'glued.hkl'), '+w')
    for i in os.listdir(hklsPath):
        if i[-4:] == '.hkl':
            f = open(os.path.join(hklsPath, i))
            lines = f.readlines()
            for j in lines:
                out.write(j)
            f.close()
    out.close()
    
def imageSweep(path):
    try:
        os.mkdir(os.path.join(path, 'PD_maps'))
    except:
        pass
    for root, dirs, files in os.walk(path):
        for file in files:
            if file[-4:] == '.png':
                 copyfile(os.path.join(root, file), os.path.join(path, 'PD_maps', file))
                
#EXECUTE

#sortav config:
template_path = os.path.normpath(r"C:\Users\piotr\Documents\VS_Code\working_dirs\agppps\AgPPPS_Maps\templates")
hklspath = os.path.normpath(r"C:\Users\piotr\Documents\VS_Code\working_dirs\agppps\AgPPPS_Maps\250ns")

mp_ini_path = os.path.normpath('C:\\Users\\piotr\\AppData\\Roaming\\mopro\\mopro.ini')

sortav_path = os.path.join(template_path, 'sortav.exe')
sortav_inp_path = os.path.join(template_path, 'sortav.inp')

temp_path = os.path.join(template_path, 'temp_inp.inp')
temp_path2 = os.path.join(template_path, 'temp_inp2.inp')
mp_path = os.path.join(template_path, 'VMoPro-win32.exe')

#sortav, mopro:

sortMoprPath = os.path.join(hklspath, 'sortav_mopro')
os.mkdir(sortMoprPath)
for i in os.listdir(hklspath):
    if i[-4:] == '.hkl':
        path = os.path.join(hklspath, i)
        copyfile (path, os.path.join(sortMoprPath, i))

sortav(sortMoprPath, sortav_path, sortav_inp_path)

ProduceMaps(os.path.join(sortMoprPath, 'output'), mp_ini_path, template_path, temp_path, temp_path2, mp_path)
PlotDirectory(pathGrids = os.path.join(os.path.join(sortMoprPath, 'output', 'Grids')))

# rescale, sortav, mopro

sortMoprPath = os.path.join(hklspath, 'rescale_sortav_mopro')
os.mkdir(sortMoprPath)
for i in os.listdir(hklspath):
    if i[-4:] == '.hkl':
        path = os.path.join(hklspath, i)
        copyfile (path, os.path.join(sortMoprPath, i))

rescale(sortMoprPath, cutoff_low = 0.5, cutoff_high = 1.5)

for i in os.listdir(sortMoprPath):
    if i[-12:] != 'rescaled.hkl':
        os.remove(os.path.join(sortMoprPath, i))


sortav(sortMoprPath, sortav_path, sortav_inp_path)

ProduceMaps(os.path.join(sortMoprPath, 'output'), mp_ini_path, template_path, temp_path, temp_path2, mp_path)
PlotDirectory(pathGrids = os.path.join(os.path.join(sortMoprPath, 'output', 'Grids')))

# rescale, glue, sortav, mopro

sortMoprPath = os.path.join(hklspath, 'rescale_sortav_glue_mopro')
os.mkdir(sortMoprPath)
for i in os.listdir(hklspath):
    if i[-4:] == '.hkl':
        path = os.path.join(hklspath, i)
        copyfile (path, os.path.join(sortMoprPath, i))
        
rescale(sortMoprPath, cutoff_low = 0.5, cutoff_high = 1.5)

for i in os.listdir(sortMoprPath):
    if i[-12:] != 'rescaled.hkl':
        os.remove(os.path.join(sortMoprPath, i))
        
glue(sortMoprPath)
os.mkdir(os.path.join(sortMoprPath, 'glued'))
copyfile(os.path.join(sortMoprPath, 'glued.hkl'), os.path.join(sortMoprPath, 'glued', 'glued.hkl'))
sortMoprPath = os.path.join(sortMoprPath, 'glued')

sortav(sortMoprPath, sortav_path, sortav_inp_path)

ProduceMaps(os.path.join(sortMoprPath, 'output'), mp_ini_path, template_path, temp_path, temp_path2, mp_path)
PlotDirectory(pathGrids = os.path.join(os.path.join(sortMoprPath, 'output', 'Grids')))
        
# copy images
imageSweep(hklspath)     



        
        