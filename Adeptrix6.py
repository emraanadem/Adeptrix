import math
import csv
import json
from tkinter import Canvas, PhotoImage, StringVar
import csv
import time
import pandas
import zlib
import sys
import pandas as pd
import numpy as np
import os
import threading
from bokeh.plotting import *
from bokeh.models import *
from sklearn.metrics import *
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.utils import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import text
from numpy import poly1d
from numpy import polyfit
import numpy as np
from scipy import integrate
from scipy.stats import f_oneway, mannwhitneyu, kruskal
import imgkit
from PyInstaller.utils.hooks import collect_data_files, eval_statement

datas = collect_data_files('tkinterDnD2')

class Adeptrix:

    datamini = []
    PCAlist = {}
    peaks = []
    rawdata = []
    gendata = []
    maxintens = 0
    allpeaks = []
    filteredpeaks = []
    negpeaks = []
    finalallpeaks = []
    removepeaks = []
    negdata = []
    filepathh = ''
    possrempeaks = []
    red = []
    orange = []
    datas = {}
    finalclearpeaks = []
    negcontrolfile = ''
    datafile = ''
    datafiles = []
    masses = []
    intensities = []
    ratios = []
    filename = ''
    peakpointcount = 5
    alignmentfile = ''
    totaldatas = {}
    tag = []

    @staticmethod
    def importer():
        from rpy2 import robjects
        filelisttwo = []
        masses = []
        intens = []
        robjects.r('library(MALDIquant)')
        robjects.r('library(readBrukerFlexData)')
        robjects.r("peaks <- readBrukerFlexDir('" + Adeptrix.negcontrolfile + "', removeCalibrationScans = FALSE, verbose = TRUE, useSpectraNames = FALSE)")
        robjects.r("masses <- peaks[[1]]$spectrum$mass")
        robjects.r("intensities <- peaks[[1]]$spectrum$intensity")
        robjects.r('library(MALDIquantForeign)')
        robjects.r("peaks <- importBrukerFlex('"+ Adeptrix.negcontrolfile + "')")
        robjects.r('peaks <- alignSpectra(peaks)')
        robjects.r("bSnip <- removeBaseline(peaks, method='TopHat')")
        #robjects.r("bSnip <- alignSpectra(bSnip, halfWindowSize = 20, SNR = 2, tolerance = .002, warpingMethod = 'lowess'")
        robjects.r("setClass('spec', slots = c(x='numeric', y='numeric'))")
        robjects.r("myspec <- new('spec', x = (bSnip[[1]]@mass), y = (bSnip[[1]]@intensity))")
        robjects.r('str(myspec)')
        robjects.r('''S4_to_dataframe <- function(s4obj) {nms <- slotNames(s4obj)
        lst <- lapply(nms, function(nm) slot(s4obj, nm))
        as.data.frame(setNames(lst, nms))}''')
        robjects.r("myspec <- S4_to_dataframe(myspec)")
        robjects.r("write.table(myspec, file = '" + Adeptrix.negcontrolfile + ".txt" +"', sep = ' ', row.names = FALSE, col.names = FALSE)")
        Adeptrix.negcontrolfile = str(Adeptrix.negcontrolfile) + ".txt"
        Adeptrix.negcontrolfilter()
        for folders in Adeptrix.datafiles:
            intens = []
            masses = []
            robjects.r("peaks <- readBrukerFlexDir('" + str(folders) + "', removeCalibrationScans = FALSE, verbose = TRUE, useSpectraNames = FALSE)")
            robjects.r("masses <- peaks[[1]]$spectrum$mass")
            robjects.r("intensities <- peaks[[1]]$spectrum$intensity")
            robjects.r('library(MALDIquantForeign)')
            robjects.r("peaks <- importBrukerFlex('"+ str(folders) + "')")
            robjects.r('peaks <- alignSpectra(peaks)')
            robjects.r("bSnip <- removeBaseline(peaks, method='TopHat')")
            #robjects.r("bSnip <- alignSpectra(bSnip, halfWindowSize = 20, SNR = 2, tolerance = .002, warpingMethod = 'lowess'")
            robjects.r("setClass('spec', slots = c(x='numeric', y='numeric'))")
            robjects.r("myspec <- new('spec', x = (bSnip[[1]]@mass), y = (bSnip[[1]]@intensity))")
            robjects.r('str(myspec)')
            robjects.r('''S4_to_dataframe <- function(s4obj) {nms <- slotNames(s4obj)
            lst <- lapply(nms, function(nm) slot(s4obj, nm))
            as.data.frame(setNames(lst, nms))}''')
            robjects.r("myspec <- S4_to_dataframe(myspec)")
            robjects.r("write.table(myspec, file = '" + str(folders) + ".txt" +"', sep = ' ', row.names = FALSE, col.names = FALSE)")
            Adeptrix.datafile = str(folders) + '.txt'
            filelisttwo.append(Adeptrix.datafile)
        Adeptrix.datafiles = filelisttwo
        for item in Adeptrix.datafiles:
            Adeptrix.datafile = item
            Adeptrix.datasplitter()

    @staticmethod
    def negcontrolfilter():
        import os
        # assign directory
        rowcount = 0
        data = []
        with open(Adeptrix.negcontrolfile, 'r' ) as adep:
            datas = list(csv.reader(adep, delimiter = ' '))
            for row in datas:
                row[0] = float(row[0])
                row[1] = int(float(row[1]))
                data.append(row)
                Adeptrix.negdata.append(row)
                Adeptrix.rawdata = data
                rowcount += 1
        rowcount2 = int(rowcount/2000 + 1)
        intensities = []
        for row in Adeptrix.rawdata:
            intensities.append(float(row[1]))
        Adeptrix.maxintens = int(max(intensities))
        numthree = 0
        for num in range(0, rowcount2):
            data = []
            for numtwo in range(0, 2000):
                if(numthree < rowcount):
                    data.append([float(datas[numthree][0]), float(datas[numthree][1])])
                    numthree += 1
            t = threading.Thread(target=Adeptrix.controllarge, args = [data]).start()
        for t in threading.enumerate():
            if t != threading.enumerate()[0]:
                if t != threading.current_thread():
                    t.join()
        Adeptrix.allpeaks.sort()
        templist = []
        for peak in Adeptrix.allpeaks:
            if peak not in templist:
                templist.append(peak)
        templist.sort()
        Adeptrix.negpeaks.extend(templist)
        Adeptrix.datamini = []
        Adeptrix.peaks = []
        Adeptrix.rawdata = []
        Adeptrix.gendata = []
        Adeptrix.maxintens = 0
        Adeptrix.allpeaks = []
        Adeptrix.filteredpeaks = []     
        Adeptrix.finalallpeaks = []
        Adeptrix.removepeaks = [] 
        Adeptrix.negpeaks.sort()
        templist = []
        for peak in Adeptrix.negpeaks:
            if peak not in templist:
                templist.append(peak)
            templist.sort()
        Adeptrix.negpeaks = templist
    
    @staticmethod
    def compare():
        Adeptrix.finalallpeaks = []
        for peak in Adeptrix.allpeaks:
            for peaks in Adeptrix.negdata:
                if abs(peak[0] - peaks[0]) < .5:
                    if peaks[1]/peak[1] > .6:
                        Adeptrix.removepeaks.append(peak)
                    elif peaks[1]/peak[1] <= .6 and peaks[1]/peak[1] >= .3:
                        Adeptrix.possrempeaks.append(peak)
        Adeptrix.rawdata.sort()
        maxmass = Adeptrix.rawdata[-1][0]
        iterations = int(maxmass/2000) + 1
        buffer = 50
        indexcount = 0
        for vals in range (0, iterations):
            indexcount2 = vals*2000
            indexcount = indexcount2 + buffer
            if ((vals+1)*2000 + 20 < len(Adeptrix.rawdata)):
                threading.Thread(target = Adeptrix.peakarea, args = [Adeptrix.possrempeaks, 
                    Adeptrix.rawdata[indexcount-50:(vals+1)*2000 + 20], Adeptrix.negdata[indexcount:(vals+1)*2000 + 20]]).start()
                indexcount = vals*2000
            else:
                threading.Thread(target = Adeptrix.peakarea, args = [Adeptrix.possrempeaks, 
                    Adeptrix.rawdata[indexcount-50:-1], Adeptrix.negdata[indexcount:-1]]).start()
        for thread in threading.enumerate():
            if thread != threading.enumerate()[0]:
                if thread != threading.current_thread():
                    thread.join()
        indexcount = 0
        leftoverpeaks = []
        for peak in Adeptrix.allpeaks:
            if peak not in Adeptrix.possrempeaks:
                leftoverpeaks.append(peak)
        for vals in range (0, iterations):
            indexcount2 = vals*2000
            indexcount = indexcount2 + buffer
            if ((vals+1)*2000 + 20 < len(Adeptrix.rawdata)):
                threading.Thread(target = Adeptrix.leftoverpeakarea, args = [leftoverpeaks, 
                    Adeptrix.rawdata[indexcount-50:(vals+1)*2000 + 20], Adeptrix.negdata[indexcount:(vals+1)*2000 + 20]]).start()
                indexcount = vals*2000
            else:
                threading.Thread(target = Adeptrix.leftoverpeakarea, args = [leftoverpeaks, 
                    Adeptrix.rawdata[indexcount-50:-1], Adeptrix.negdata[indexcount:-1]]).start()
        for thread in threading.enumerate():
            if thread != threading.enumerate()[0]:
                if thread != threading.current_thread():
                    thread.join()
        indexcount = 0
        for vals in range (0, iterations):
            indexcount2 = vals*2000
            indexcount = indexcount2 + buffer
            if ((vals+1)*2000 + 20 < len(Adeptrix.rawdata)):
                threading.Thread(target = Adeptrix.peakrebound, args = [Adeptrix.possrempeaks, 
                    Adeptrix.rawdata[indexcount-50:(vals+1)*2000 + 20], Adeptrix.negdata[indexcount:(vals+1)*2000 + 20]]).start()
                indexcount = vals*2000
            else:
                threading.Thread(target = Adeptrix.peakrebound, args = [Adeptrix.possrempeaks, 
                    Adeptrix.rawdata[indexcount-50:-1], Adeptrix.negdata[indexcount:-1]]).start()
        for thread in threading.enumerate():
            if thread != threading.enumerate()[0]:
                if thread != threading.current_thread():
                    thread.join()
        for peak in Adeptrix.allpeaks:
            for peaks in Adeptrix.allpeaks:
                if peak[0] == peaks[0]:
                    if len(peak) < len(peaks):
                        Adeptrix.removepeaks.append(peak)
                    if len(peaks) < len(peak):
                        Adeptrix.removepeaks.append(peaks)
        for peak in Adeptrix.ratios:
            for peaks in Adeptrix.ratios:
                if peak[0] == peaks[0]:
                    if len(peak) < len(peaks):
                        Adeptrix.removepeaks.append(peak)
                    if len(peaks) < len(peak):
                        Adeptrix.removepeaks.append(peaks)
        for peak in Adeptrix.orange:
            for peaks in Adeptrix.orange:
                if peak[0] == peaks[0]:
                    if len(peak) < len(peaks):
                        Adeptrix.removepeaks.append(peak)
                    if len(peaks) < len(peak):
                        Adeptrix.removepeaks.append(peaks)
        for peak in Adeptrix.red:
            for peaks in Adeptrix.red:
                if peak[0] == peaks[0]:
                    if len(peak) < len(peaks):
                        Adeptrix.removepeaks.append(peak)
                    if len(peaks) < len(peak):
                        Adeptrix.removepeaks.append(peaks)
        Adeptrix.red = [list(t) for t in set(tuple(element) for element in Adeptrix.red)]
        Adeptrix.orange = [list(t) for t in set(tuple(element) for element in Adeptrix.orange)]
        Adeptrix.finalclearpeaks = [list(t) for t in set(tuple(element) for element in Adeptrix.finalclearpeaks)]
        Adeptrix.finalallpeaks = [list(t) for t in set(tuple(element) for element in Adeptrix.finalallpeaks)]
        Adeptrix.red.sort()
        Adeptrix.orange.sort()
        Adeptrix.finalclearpeaks.sort()
        Adeptrix.finalallpeaks.sort()
        Adeptrix.allpeaks.sort()

    @staticmethod
    def peakrebound(peaks, data, negdata):
        peaks.sort()
        data.sort()
        negdata.sort()
        for peak in peaks:
            if peak[0] > data[0][0] and peak[0] < data[-1][0]:
                breakerpossur = False
                breakernegsur = False
                surroundgreater = []
                surroundless = []
                surroundvalues = []
                surroundvalues.append(peak)
                ref = peak
                intens = []
                peak[1] = int(peak[1])
                if((data.index(peak[0:2]) > 15) and (data.index(peak[0:2]) < len(data)-16)):
                    for num in range(1, 11):
                        poscount = 0
                        negcount = 0
                        posnum = data[data.index(peak[0:2])+num][1]
                        negnum = data[data.index(peak[0:2])-num][1]
                        for counter in range(1, Adeptrix.peakpointcount+1):
                            if(data[data.index(peak[0:2])+num+counter][1] > posnum):
                                poscount += 1
                            if(data[data.index(peak[0:2])-num-counter][1] > negnum):
                                negcount += 1
                            posnum = data[data.index(peak[0:2])+num+counter][1]
                            negnum = data[data.index(peak[0:2])-num-counter][1]
                        if(data[data.index(peak[0:2])+num][1] < 0):
                            breakerpossur = True
                        if(data[data.index(peak[0:2])-num][1] < 0):
                            breakernegsur = True
                        if((breakerpossur == False) and (poscount >= 5)):
                            breakerpossur = True
                        if(breakerpossur == False):
                            surroundvalues.append(data[data.index(peak[0:2])+num])
                        if((breakernegsur == False) and (negcount >= 5)):
                            breakernegsur = True
                        if(breakernegsur == False):
                            surroundvalues.append(data[data.index(peak[0:2])-num])
                    for val in surroundvalues:
                        if val[0] > peak[0]:
                            surroundgreater.append(val)
                        if val[0] < peak[0]:
                            surroundless.append(val)
                    surroundvalues.sort()
                    for val in surroundvalues:
                        intens.append(val[1])
                    distance = ref[1]-min(intens)
                    #if((peak[0] < 2708 and peak[0] > 2707) or (peak[0] < 1710 and peak[0] > 1709)):
                        #print(min(intens))
                        #print(peak)
                        #print(surroundgreater)
                        #print(surroundless)
                        #print(distance)
                    slopegreater = []
                    slopeless = []
                    for num in surroundgreater:
                        if surroundgreater.index(num) < len(surroundgreater)-1:
                            for numb in range(1,2):
                                if(surroundgreater[surroundgreater.index(num)+numb][1] - num[1] >= .3*distance):
                                    Adeptrix.removepeaks.append(peak)
                    for num in surroundless:
                        if surroundless.index(num) < len(surroundless)-1:
                            for numb in range(1,2):
                                if(surroundless[surroundless.index(num)+numb][1] - num[1] >= .3*distance):
                                    Adeptrix.removepeaks.append(peak)

    @staticmethod
    def signalnoise(peak, data):
        data.sort()
        if peak[0] > data[0][0] and peak[0] < data[-1][0]:
            breakerpossur = False
            breakernegsur = False
            surroundgreater = []
            surroundless = []
            surroundvalues = []
            surroundvalues.append(peak)
            intens = []
            peak[1] = int(peak[1])
            if((data.index(peak[0:2]) > 15) and (data.index(peak[0:2]) < len(data)-16)):
                for num in range(1, 11):
                    poscount = 0
                    negcount = 0
                    posnum = data[data.index(peak[0:2])+num][1]
                    negnum = data[data.index(peak[0:2])-num][1]
                    for counter in range(1, Adeptrix.peakpointcount+1):
                        if(data[data.index(peak[0:2])+num+counter][1] > posnum):
                            poscount += 1
                        if(data[data.index(peak[0:2])-num-counter][1] > negnum):
                            negcount += 1
                        posnum = data[data.index(peak[0:2])+num+counter][1]
                        negnum = data[data.index(peak[0:2])-num-counter][1]
                    if(data[data.index(peak[0:2])+num][1] < 0):
                        breakerpossur = True
                    if(data[data.index(peak[0:2])-num][1] < 0):
                        breakernegsur = True
                    if((breakerpossur == False) and (poscount >= 5)):
                        breakerpossur = True
                    if(breakerpossur == False):
                        surroundvalues.append(data[data.index(peak[0:2])+num])
                        if data[data.index(peak[0:2])+num] > peak[0:2]:
                            surroundgreater.append(data[data.index(peak[0:2])+num])
                    if((breakernegsur == False) and (negcount >= 5)):
                        breakernegsur = True
                    if(breakernegsur == False):
                        surroundvalues.append(data[data.index(peak[0:2])-num])
                        if(data[data.index(peak[0:2])-num] < peak[0:2]):
                            surroundless.append(data[data.index(peak[0:2])-num])
                surroundvalues.sort()
                for val in surroundvalues:
                    if val[0] > peak[0]:
                        surroundgreater.append(val)
                    if val[0] < peak[0]:
                        surroundless.append(val)
                for val in surroundvalues:
                    intens.append(val[1])
        slopegreater = []
        surroundgreater.append(peak)
        surroundless.append(peak)
        for num in surroundgreater:
            if num != surroundgreater[-1]:
                if(surroundgreater[surroundgreater.index(num)+1][1] > 0 and num[1] > 0):
                    slopegreater.append(num[1]/surroundgreater[surroundgreater.index(num)+1][1])
                else:
                    greaterrev = slopegreater
                    greaterrev.sort()
                    for nums in greaterrev:
                        if nums > 0:
                            slopegreater.append(abs(peak[1] - nums))
                            break
        slopeless = []
        for num in surroundless:
            if num != surroundless[-1]:
                if(surroundless[surroundless.index(num)+1][1] > 0 and num[1] > 0):
                    slopeless.append(num[1]/surroundless[surroundless.index(num)+1][1])
                else:
                    lessrev = slopeless
                    lessrev.sort()
                    for nums in lessrev:
                        if nums > 0:
                            slopeless.append(abs(peak[1] - nums))
                            break
        slopegreater.reverse()
        slopeless.reverse()
        negsn = 0
        possn = 0
        for val in slopegreater:
            if val > 2:
                possn = val
                break
        if possn == 0:
            for val in slopegreater:
                if val > 0:
                    possn = abs(peak[1] - val)
                    break
        if possn == 0:
            possn = peak[1]
        for val in slopeless:
            if val > 2:
                negsn = val
                break
        if negsn == 0:
            for val in slopeless:
                if val > 0:
                    negsn = abs(peak[1] - val)
                    break
        if negsn == 0:
            negsn = peak[1]
        peaker = list(peak)
        if possn > negsn:
            peaker.append(round(negsn, 2))
        if negsn > possn:
            peaker.append(round(possn, 2))
        if possn == negsn:
            peaker.append(round(possn, 2))
        return peaker

    @staticmethod
    def peakwidth(peaks, data):
        for peak in peaks:
            nearvals = []
            peakwratio = []
            peak[1] = int(peak[1])
            for num in range(0,3):
                nearvals.append(data[data.index(peak[0:2])+num])
                nearvals.append(data[data.index(peak[0:2])-num])
            nearvals.sort()
            ratio = abs(math.log((max(nearvals)[0] - min(nearvals)[0])/2.5)*1.2)
            peakwratio.append([peak, ratio])
        return peakwratio

    @staticmethod
    def leftoverpeakarea(peaks, data, negdata):
        safe = []
        matchfound = False
        #ratios = Adeptrix.peakwidth(peaks, data)
        ratios = []
        for numbeee in peaks:
            ratios.append([numbeee, 1])
        peaks.sort()
        data.sort()
        negdata.sort()
        Adeptrix.negpeaks.sort()
        for peak in peaks:
            if peak[0] > data[0][0] and peak[0] < data[-1][0]:
                breakerpossur = False
                breakernegsur = False
                matchfound = False
                surroundvalues = []
                peak[1] = int(peak[1])
                surroundvalues.append(peak)
                if((data.index(peak[0:2]) > 15) and (data.index(peak[0:2]) < len(data)-16)):
                    for num in range(1, 11):
                        poscount = 0
                        negcount = 0
                        posnum = data[data.index(peak[0:2])+num][1]
                        negnum = data[data.index(peak[0:2])-num][1]
                        for counter in range(1, Adeptrix.peakpointcount+1):
                            if(data[data.index(peak[0:2])+num+counter][1] > posnum):
                                poscount += 1
                            if(data[data.index(peak[0:2])-num-counter][1] > negnum):
                                negcount += 1
                            posnum = data[data.index(peak[0:2])+num+counter][1]
                            negnum = data[data.index(peak[0:2])-num-counter][1]
                        if(data[data.index(peak[0:2])+num][1] < 0):
                            breakerpossur = True
                        if(data[data.index(peak[0:2])-num][1] < 0):
                            breakernegsur = True
                        if((breakerpossur == False) and (poscount >= 5)):
                            breakerpossur = True
                        if(breakerpossur == False):
                            surroundvalues.append(data[data.index(peak[0:2])+num])
                        if((breakernegsur == False) and (negcount >= 5)):
                            breakernegsur = True
                        if(breakernegsur == False):
                            surroundvalues.append(data[data.index(peak[0:2])-num])
                    surroundvalues.sort()
                    surroundx = []
                    surroundy = []
                    for numb in surroundvalues:
                        surroundx.append(numb[0])
                        surroundy.append(numb[1])
                    for info in Adeptrix.negdata:
                        breakerposneg = False
                        breakernegneg = False
                        if matchfound == True:
                            break
                        if info[0] > negdata[0][0] and info[0] < negdata[-1][0]:
                            if((negdata.index(info) > 15) and (negdata.index(info) < len(negdata)-16)):
                                if(abs(peak[0] - info[0]) < .5):
                                    negsurround = []
                                    negsurround.append(info)
                                    for numbe in range(1, 11):
                                        backposcount = 0
                                        backnegcount = 0
                                        backposnum = negdata[negdata.index(info)+numbe][1]
                                        backnegnum = negdata[negdata.index(info)-numbe][1]
                                        for counterr in range(1, Adeptrix.peakpointcount+1):
                                            if(negdata[negdata.index(info)+numbe+counterr][1] > backposnum):
                                                backposcount += 1
                                            if(negdata[negdata.index(info)-numbe-counterr][1] > backnegnum):
                                                backnegcount += 1
                                            backposnum = negdata[negdata.index(info)+numbe+counterr][1]
                                            backnegnum = negdata[negdata.index(info)-numbe-counterr][1]
                                        if(negdata[negdata.index(info)+numbe][1] < 0):
                                            breakerposneg = True
                                        if(negdata[negdata.index(info)-numbe][1] < 0):
                                            breakernegneg = True
                                        if((breakerposneg == False) and (backposcount >= 5)):
                                            breakerposneg = True
                                        if(breakerposneg == False):
                                            negsurround.append(negdata[negdata.index(info)+numbe])
                                        if((breakernegneg == False) and (backnegcount >= 5)):
                                            breakernegneg = True
                                        if(breakernegneg == False):
                                            negsurround.append(negdata[negdata.index(info)-numbe])
                                    negsurround.sort()
                                    negsurroundx = []
                                    negsurroundy = []
                                    matchfound = True
                                    for numb in negsurround:
                                        negsurroundx.append(numb[0])
                                        negsurroundy.append(numb[1])
                                    peakintegral = integrate.cumulative_trapezoid(surroundy, surroundx)
                                    negintegrals = integrate.cumulative_trapezoid(negsurroundy, negsurroundx)
                                    negintegral = []
                                    for peakss in ratios:
                                        if peakss[0][0] == peak[0]:
                                            for val in negintegrals:
                                                negintegral.append(val*peakss[1])
                                    peaker = list(Adeptrix.signalnoise(peak, data))
                                    peaker.append(round(peakintegral[-1]/negintegral[-1], 2))
                                    if(len(peakintegral) < 1):
                                        #print(peak, info)
                                        #print(surroundvalues)
                                        #print(negsurround)
                                        Adeptrix.removepeaks.append(peaker)
                                    elif(len(negintegral) == 0 and len(peakintegral > 0)):
                                        Adeptrix.finalclearpeaks.append(peaker)
                                    elif(len(negintegral) > 0 and len(peakintegral) > 0):
                                        if(info[1] != 0):
                                            if(peak[1]/info[1] > 2.222222):
                                                if((abs(((negintegral[-1]/peakintegral[-1]) < .50)))):
                                                    if peakintegral[-1]/negintegral[-1] < 2.5:
                                                        if peak not in Adeptrix.red:
                                                            peaker.append('Highly Questionable')
                                                            Adeptrix.red.append(peaker)
                                                            Adeptrix.finalallpeaks.append(peaker)
                                                    elif peakintegral[-1]/negintegral[-1] < 3:
                                                        if peak not in Adeptrix.orange:
                                                            peaker.append('Questionable')
                                                            Adeptrix.orange.append(peaker)
                                                            Adeptrix.finalallpeaks.append(peaker)
                                                    else:
                                                        if peak not in Adeptrix.finalclearpeaks:
                                                            peaker.append('Clear Peak')
                                                            Adeptrix.finalclearpeaks.append(peaker)
                                                            Adeptrix.finalallpeaks.append(peaker)


    @staticmethod
    def peakarea(peaks, data, negdata):
        matchfound = False
        safe = []
        #ratios = Adeptrix.peakwidth(peaks, data)
        ratios = []
        for numbeee in peaks:
            ratios.append([numbeee, 1])
        peaks.sort()
        data.sort()
        negdata.sort()
        Adeptrix.negpeaks.sort()
        for peak in peaks:
            if peak[0] > data[0][0] and peak[0] < data[-1][0]:
                breakerpossur = False
                breakernegsur = False
                matchfound = False
                surroundvalues = []
                peak[1] = int(peak[1])
                surroundvalues.append(peak)
                if((data.index(peak[0:2]) > 15) and (data.index(peak[0:2]) < len(data)-16)):
                    for num in range(1, 11):
                        poscount = 0
                        negcount = 0
                        posnum = data[data.index(peak[0:2])+num][1]
                        negnum = data[data.index(peak[0:2])-num][1]
                        for counter in range(1, Adeptrix.peakpointcount+1):
                            if(data[data.index(peak[0:2])+num+counter][1] > posnum):
                                poscount += 1
                            if(data[data.index(peak[0:2])-num-counter][1] > negnum):
                                negcount += 1
                            posnum = data[data.index(peak[0:2])+num+counter][1]
                            negnum = data[data.index(peak[0:2])-num-counter][1]
                        if(data[data.index(peak[0:2])+num][1] < 0):
                            breakerpossur = True
                        if(data[data.index(peak[0:2])-num][1] < 0):
                            breakernegsur = True
                        if((breakerpossur == False) and (poscount >= 5)):
                            breakerpossur = True
                        if(breakerpossur == False):
                            surroundvalues.append(data[data.index(peak[0:2])+num])
                        if((breakernegsur == False) and (negcount >= 5)):
                            breakernegsur = True
                        if(breakernegsur == False):
                            surroundvalues.append(data[data.index(peak[0:2])-num])
                    surroundvalues.sort()
                    surroundx = []
                    surroundy = []
                    for numb in surroundvalues:
                        surroundx.append(numb[0])
                        surroundy.append(numb[1])
                    for info in Adeptrix.negdata:
                        breakerposneg = False
                        breakernegneg = False
                        if matchfound == True:
                            break
                        if info[0] > negdata[0][0] and info[0] < negdata[-1][0]:
                            if((negdata.index(info) > 15) and (negdata.index(info) < len(negdata)-16)):
                                if(abs(peak[0] - info[0]) < .5):
                                    negsurround = []
                                    negsurround.append(info)
                                    for numbe in range(1, 11):
                                        backposcount = 0
                                        backnegcount = 0
                                        backposnum = negdata[negdata.index(info)+numbe][1]
                                        backnegnum = negdata[negdata.index(info)-numbe][1]
                                        for counterr in range(1, Adeptrix.peakpointcount+1):
                                            if(negdata[negdata.index(info)+numbe+counterr][1] > backposnum):
                                                backposcount += 1
                                            if(negdata[negdata.index(info)-numbe-counterr][1] > backnegnum):
                                                backnegcount += 1
                                            backposnum = negdata[negdata.index(info)+numbe+counterr][1]
                                            backnegnum = negdata[negdata.index(info)-numbe-counterr][1]
                                        if(negdata[negdata.index(info)+numbe][1] < 0):
                                            breakerposneg = True
                                        if(negdata[negdata.index(info)-numbe][1] < 0):
                                            breakernegneg = True
                                        if((breakerposneg == False) and (backposcount >= 5)):
                                            breakerposneg = True
                                        if(breakerposneg == False):
                                            negsurround.append(negdata[negdata.index(info)+numbe])
                                        if((breakernegneg == False) and (backnegcount >= 5)):
                                            breakernegneg = True
                                        if(breakernegneg == False):
                                            negsurround.append(negdata[negdata.index(info)-numbe])
                                    negsurround.sort()
                                    negsurroundx = []
                                    negsurroundy = []
                                    for numb in negsurround:
                                        negsurroundx.append(numb[0])
                                        negsurroundy.append(numb[1])
                                    peakintegral = integrate.cumulative_trapezoid(surroundy, surroundx)
                                    negintegrals = integrate.cumulative_trapezoid(negsurroundy, negsurroundx)
                                    negintegral = []
                                    matchfound = True
                                    for peakss in ratios:
                                        if peakss[0][0] == peak[0]:
                                            for val in negintegrals:
                                                negintegral.append(val*peakss[1])
                                    peaker = list(Adeptrix.signalnoise(peak, data))
                                    peaker.append(round(peakintegral[-1]/negintegral[-1], 2))
                                    if(len(peakintegral) < 1):
                                        #print(peak, info)
                                        #print(surroundvalues)
                                        #print(negsurround)
                                        Adeptrix.removepeaks.append(peaker)
                                    elif(len(negintegral) == 0 and len(peakintegral > 0)):
                                        Adeptrix.finalclearpeaks.append(peaker)
                                    elif(len(negintegral) > 0 and len(peakintegral) > 0):
                                        if(info[1] != 0):
                                            if(peak[1]/info[1] > 2.222222):
                                                if((abs(((negintegral[-1]/peakintegral[-1]) < .50)) 
                                                and abs(peakintegral[-1]/negintegral[-1] < 3))):
                                                    if peakintegral[-1]/negintegral[-1] < 2.5:
                                                        if peak not in Adeptrix.red:
                                                            peaker.append('Highly Questionable')
                                                            Adeptrix.red.append(peaker)
                                                            Adeptrix.finalallpeaks.append(peaker)
                                                    elif peakintegral[-1]/negintegral[-1] < 3:
                                                        if peak not in Adeptrix.orange:
                                                            peaker.append('Questionable')
                                                            Adeptrix.orange.append(peaker)
                                                            Adeptrix.finalallpeaks.append(peaker)
                                                    else:
                                                        if peak not in Adeptrix.finalclearpeaks:
                                                            peaker.append('Clear Peak')
                                                            Adeptrix.finalclearpeaks.append(peaker)
                                                            Adeptrix.finalallpeaks.append(peaker)


    @staticmethod
    def filter(peak):
        peaks = Adeptrix.allpeaks
        data = Adeptrix.rawdata
        negdata = Adeptrix.negdata
        numb = 1
        ratios = Adeptrix.datafitone(peak)
        for peakk in peaks:
            if peakk == peak:
                for info in negdata:
                    if abs(info[0] - peak[0]) < .5:
                        for item in ratios:
                            if item[0] == peak:
                                numb = .8
        return numb

    @staticmethod
    def datasplitter():
        Adeptrix.filename = Adeptrix.datafile.split("/")[-1]
        Adeptrix.filename = str(Adeptrix.filename.split(".")[0])
        rowcount = 0
        directories = str(os.getcwd()) + '/Radx_data_8_24/Mutants/Sample Stuff/Peak Data'
        with open(Adeptrix.datafile, 'r' ) as adep:
            datas = list(csv.reader(adep, delimiter = ' '))
            for row in datas:
                rowcount += 1
        rowcount2 = int(rowcount/2000 + 1)
        numthree = 0
        Adeptrix.loader()
        intensities = []
        for row in Adeptrix.rawdata:
            intensities.append(float(row[1]))
        Adeptrix.maxintens = int(max(intensities))
        for num in range(0, rowcount2):
            data = []
            for numtwo in range(0, 2000):
                if(numthree < rowcount):
                    data.append([round(float(datas[numthree][0]), 2), float(datas[numthree][1])])
                    numthree += 1
            t = threading.Thread(target=Adeptrix.controllarge, args = [data]).start()
        for t in threading.enumerate():
            if t != threading.enumerate()[0]:
                if t != threading.current_thread():
                    t.join()
        Adeptrix.allpeaks.sort()
        templist = []
        for peak in Adeptrix.allpeaks:
            if peak not in templist:
                templist.append(peak)
        templist.sort()
        Adeptrix.allpeaks = templist
        Adeptrix.compare()
        dict = {}
        dict['Peaks for ' + Adeptrix.datafile] = {}
        newlistclear = []
        newlistorange = []
        newlistred = []
        for item in Adeptrix.finalclearpeaks:
            newlistclear.append(list(dict.fromkeys(item)))
        for item in Adeptrix.orange:
            newlistorange.append(list(dict.fromkeys(item)))
        for item in Adeptrix.red:
            newlistred.append(list(dict.fromkeys(item)))
        Adeptrix.finalclearpeaks = newlistclear
        Adeptrix.orange = newlistorange
        Adeptrix.red = newlistred
        dict['Peaks for ' + Adeptrix.datafile]['Clear Peaks'] = Adeptrix.finalclearpeaks
        dict['Peaks for ' + Adeptrix.datafile]['Questionable Peaks'] = Adeptrix.orange
        dict['Peaks for ' + Adeptrix.datafile]['Highly Questionable Peaks'] = Adeptrix.red
        newlisttot = []
        newlisttot.extend(newlistclear)
        newlisttot.extend(newlistorange)
        newlisttot.extend(newlistred)
        for peak in newlisttot:
            for peaks in newlisttot:
                if peak[0] == peaks[0]:
                    if len(peak) < len(peaks):
                        newlisttot.remove(peak)
                    if len(peaks) < len(peak):
                        newlisttot.remove(peaks)
        newlisttot.sort(key=lambda x: x[0], reverse=False)
        dict5 = {}
        dict5['All Peaks'] = newlisttot
        dict6 = dict5['All Peaks']
        cols = ['Mass/Charge Ratio', 'Intensity', 'Signal/Noise Ratio', 'Intensity Ratio', 'Peak Confidence']
        if os.path.exists("./Peak Images/" + str(Adeptrix.filename)):
            csv_file = "./Peak Images/" + str(Adeptrix.filename) + "/PeakDataTable.csv"            
        else:
            os.makedirs("./Peak Images/" + str(Adeptrix.filename))
            csv_file = "./Peak Images/" + str(Adeptrix.filename) + "/PeakDataTable.csv"
        with open(csv_file, 'w') as adep:
            try:
                write = csv.writer(adep)
                write.writerow(cols)
                write.writerows(dict6)
            except IOError:
                print("I/O error")
        if(len(Adeptrix.finalallpeaks) > 0):
            print("Peaks have been uploaded to the appropriate subfolder in the Peak Images folder.")
            # Develop dictionary to get more specific to mutations and amino acid sequences as program develops further
        else:
            print('No Peaks have been found.')
        Adeptrix.plot()
        for row in Adeptrix.rawdata:
            Adeptrix.masses.append(row[0])
            Adeptrix.intensities.append(row[1])
            Adeptrix.tag.append(Adeptrix.datafile.split("/")[-1].split('.')[0])
        Adeptrix.datas['Masses'] = Adeptrix.masses
        Adeptrix.datas['Intensities'] = Adeptrix.intensities
        Adeptrix.datas['Target'] = Adeptrix.tag
        Adeptrix.totaldatas = Adeptrix.datas
        Adeptrix.datamini = []
        Adeptrix.rawdata = []
        Adeptrix.red = []
        Adeptrix.orange = []
        Adeptrix.finalclearpeaks = []
        Adeptrix.possrempeaks = []
        Adeptrix.peaks = []
        Adeptrix.gendata = []
        Adeptrix.maxintens = 0
        Adeptrix.allpeaks = []
        Adeptrix.filteredpeaks = []     
        Adeptrix.finalallpeaks = []
        Adeptrix.removepeaks = [] 

    @staticmethod
    def statistics():
        dataset = []
        datasets = []
        subsets = []
        filename = './Peak Images/'
        if '_MEIPASS2' in os.environ:
            filename = os.path.join(os.environ['_MEIPASS2'], filename)
        for item in os.listdir(filename):
            subsets = []
            if os.path.isdir(filename + item):
                for items in os.listdir(filename + item):
                    if items.endswith('.csv'):
                        with open(filename + item + '/' + items, 'r') as file:
                            reader = csv.reader(file)
                            for row in reader:
                                if row[0] != "Mass/Charge Ratio":
                                    subsets.append([float(row[0]), float(row[1])])
                            for row in subsets:
                                datasets.append([row[0], row[1], item.split("_")[len(item.split("_"))-1] + " " + item.split("_")[len(item.split("_"))-6]])
                                dataset.append([row[0], row[1], item.split("_")[len(item.split("_"))-1] + " " + item.split("_")[len(item.split("_"))-6]])
        anovas = []
        kruskals = []
        mannwhit = []
        peaks = {}
        test = pd.DataFrame(datasets)
        test.columns =['Masses', 'Intensities', 'Target']
        lisss = datasets
        targerrr = test.loc[:,['Target']].values
        targerrrtwo = list(np.unique(targerrr))
        listofstuff = ['Masses']
        for val in targerrrtwo:
            listofstuff.append(str(val))
        testflipped = pd.DataFrame(columns=listofstuff)
        for index, row in test.iterrows():
            for columnName, columnData in testflipped.iteritems():
                if row['Target'] == columnName:
                    mass = row["Masses"]
                    intensity = row["Intensities"]
                    colnum = listofstuff.index(columnName)
                    rowid = len(testflipped.index)
                    testflipped.loc[len(testflipped.index)] = [row["Masses"], 0, 0, 0, 0]
                    testflipped.iloc[rowid][columnName] = intensity
        for index, row in testflipped.iterrows():
            for indextwo, rowtwo in testflipped.iterrows():
                if (row["Masses"] < rowtwo["Masses"] + .5) and (row["Masses"] > rowtwo["Masses"] - .5):
                    for col in row.iteritems():
                        if col[0] != "Masses":
                            if col[1] == 0 and testflipped.iloc[indextwo][col[0]] != 0:
                                testflipped.iloc[index][col[0]] = testflipped.iloc[indextwo][col[0]]
                                testflipped.iloc[index]["Masses"] = testflipped.iloc[indextwo]["Masses"]
        testflipped = testflipped.drop_duplicates()
        for row in subsets:
            num = round(row[0])
            if str(num) not in peaks:
                peaks[str(num)] = []
                for rows in subsets:
                    if((rows[0] > row[0] - 1) and (rows[0] < row[0] + 1)):
                        peaks[str(num)].append(rows[1])
                if(len(peaks[str(num)]) > 1):
                    #anovas.append(f_oneway(peaks[str(num)]))
                    #kruskals.append(kruskal(peaks[str(num)]))
                    #mannwhit.append(mannwhitneyu(peaks[str(num)], method="auto"))
                    pass
        peaklist = list(peaks.keys())
        """
        anovadata = {}
        anovadata['Peaks'] = peaklist
        anovadata['Anova'] = anovas
        kruskalsdata = {}
        kruskalsdata['Peaks'] = peaklist
        kruskalsdata['Kruskal Wallis'] = kruskals
        mannwhitdata = {}
        mannwhitdata['Peaks'] = peaklist
        mannwhitdata['Mann Whitney'] = mannwhit
        anovasdf = pd.DataFrame.from_dict(anovadata)
        kruskalsdf = pd.DataFrame.from_dict(kruskalsdata)
        mannwhitdf = pd.DataFrame.from_dict(mannwhitdata)
        """
        pcas = PCA(n_components=.95)
        features = listofstuff
        totaldatass = testflipped
        x = totaldatass.loc[:, features].values
        targerrr = listofstuff
        totaldatass.to_csv('Final Data Before PCA.csv', sep = '\t')
        x = StandardScaler(with_mean=False).fit_transform(x)
        principalComponents = pcas.fit_transform(x)
        pcas.fit(x)
        principalDf = pd.DataFrame(data = principalComponents, columns = targerrr)
        fig = plt.figure(figsize = (8,8))
        figtwo = plt.figure(figsize = (8,8))
        figthree = plt.figure(figsize = (8,8))
        elbowplot = figthree.add_subplot(1,1,1)
        scoreplot = figtwo.add_subplot(1,1,1)
        scoreplot.set_xlabel('Principal Component 1', fontsize = 15)
        scoreplot.set_ylabel('Principal Component 2', fontsize = 15)
        scoreplot.set_title('2 Component PCA', fontsize = 20)
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('Principal Component 1', fontsize = 15)
        ax.set_ylabel('Principal Component 2', fontsize = 15)
        ax.set_title('2 Component PCA', fontsize = 20)
        dictstuff = {}
        listofpoints = []
        targetcounter = []
        for target in targerrr:
            listofpointers = []
            if target != "Masses":
                for index, row in principalDf.iterrows():
                    for col in row.iteritems():
                        if col[0] == target:
                            listofpointers.append([row["Masses"], col[1]])
                x = []
                y = []
                xval = 0
                yval = 0
                for row in listofpointers:
                    x.append(row[0])
                    y.append(row[1])
                for item in x:
                    xval += item
                for item in y:
                    yval += item
                xval = xval/len(x)
                yval = yval/len(y)
                listofpoints.append([xval, yval])
                targetcounter.append([[xval, yval], target])
        import matplotlib.colors as mcolors
        colorlist = list(mcolors.BASE_COLORS)
        colornum = len(targerrr)
        colors = []
        for num in range(0, colornum):
            colors.append(colorlist[num])
        """for target, colorval in zip(targerrr, colors):
            x = []
            y = []
            if target != "Masses":
                for coord in dictstuff[target]:
                    x.append(coord[0])
                    y.append(coord[1])
                scoreplot.scatter(x,y,color = colorval)
        """
        kmeans = KMeans(n_clusters=2)
        listofpointss = pd.DataFrame(listofpoints)
        q = kmeans.fit_predict(listofpointss)
        #filter rows of original data
        filtered_label0 = listofpointss[q == 1]
        filtered_label1 = listofpointss[q == 0]
        filtered_label0.columns = ["x", "y"]
        filtered_label1.columns = ["x", "y"]
        targlist1 = []
        targlist2 = []
        for index,row in filtered_label0.iterrows():
            mass = row[0]
            intensity = row[1]
            for rowtwo in targetcounter:
                    if rowtwo[1] != "Masses":
                        if rowtwo[0][0] == mass and rowtwo[0][1] == intensity:
                            targlist1.append(rowtwo[1])
        for index,row in filtered_label1.iterrows():
            mass = row[0]
            intensity = row[1]
            for rowtwo in targetcounter:
                    if rowtwo[1] != "Masses":
                        if rowtwo[0][0] == mass and rowtwo[0][1] == intensity:
                            targlist2.append(rowtwo[1])
        filtered_label1["Target"] = targlist2
        filtered_label0["Target"] = targlist1
        #plotting the results
        plottt = plt.figure(figsize = (8,8))
        scorer = plottt.add_subplot(1,1,1)
        print(filtered_label0)
        print(filtered_label1)
        for index,row in filtered_label0.iterrows():
            scorer.scatter(row["x"], row["y"], color = "red")
            scorer.annotate(row["Target"], xy = (row["x"], row["y"]))
        for index,row in filtered_label1.iterrows():
            scorer.scatter(row["x"], row["y"], color = "blue")
            scorer.annotate(row["Target"], xy = (row["x"], row["y"]))
        scorer.set_xlabel('Principal Component 1', fontsize = 15)
        scorer.set_ylabel('Principal Component 2', fontsize = 15)
        scorer.set_title('2 Component PCA', fontsize = 20)
        plottt.savefig("test.png")
        """
        targets = list((totaldatass.loc[:,['Target']].values).tolist())
        targs = []
        for targe in targets:
            if targe[0] not in targs:
                targs.append(targe[0])
        colornum = len(targs)
        colors = []
        finalDf = pd.concat([principalDf, totaldatass[['Target']]], axis = 1)
        import matplotlib.colors as mcolors
        colorlist = list(mcolors.BASE_COLORS)
        points = []
        targsss = []
        finalDf = principalDf
        for targer in targerrr:
            x = []
            y = []
            targ = ''
            for i in range(len(finalDf)):
                if finalDf.loc[i, 'Target'] == targer:
                    x.append(finalDf.loc[i, 'principal component 1'])
                    y.append(finalDf.loc[i, 'principal component 2'])
                    targ = targer
            pc1 = 0
            pc2 = 0
            for num in x:
                pc1 += num
            for num in y:
                pc2 += num
            pc1 = pc1/len(x)
            pc2 = pc2/len(y)
            points.append([pc1, pc2])
            targsss.append(targ)
        inertias = []
        newtargs = list(np.unique(targerrr))
        print(newtargs)
        q = len(newtargs) + 1
        for i in range(1,q):
            kmeans = KMeans(n_clusters=i)
            kmeans.fit(points)
            inertias.append(kmeans.inertia_)
        elbowplot.plot(range(1,q), inertias, marker='o')
        elbowplot.set_title('Elbow method')
        elbowplot.set_xlabel('Number of clusters')
        elbowplot.set_ylabel('Inertia')
        figthree.savefig('Elbow Grouping.png')
        pointss = pd.DataFrame(points)
        print(pointss)
        pointss.columns =['Masses', 'Intensities']
        kmeans = KMeans(n_clusters=2)
        q = kmeans.fit_predict(pointss)
        #filter rows of original data
        filtered_label0 = pointss[q == 1]
        filtered_label1 = pointss[q == 0]
        #plotting the results
        plottt = plt.figure(figsize = (8,8))
        scorer = plottt.add_subplot(1,1,1)
        scorer.scatter(filtered_label0['Masses'], filtered_label0['Intensities'], color = "red")
        scorer.scatter(filtered_label1['Masses'], filtered_label1['Intensities'], color = "blue")
        plottt.savefig('test.png')
        x = []
        y = []
        kmeans.fit(points)
        for point in points:
            x.append(point[0])
            y.append(point[1])
        colorss = kmeans.labels_
        for point in points:
            if colorss[points.index(point)] == 0:
                scoreplot.scatter(point[0], point[1], c = 'b', s = 50)
                scoreplot.annotate(targsss[points.index(point)], xy = (point[0], point[1]), xytext = (point[0], point[1]))
            if colorss[points.index(point)] == 1:
                scoreplot.scatter(point[0], point[1], c = 'g', s = 50)
                scoreplot.annotate(targsss[points.index(point)], xy = (point[0], point[1]), xytext = (point[0], point[1]))
            if colorss[points.index(point)] == 2:
                scoreplot.scatter(point[0], point[1], c = 'r', s = 50)
                scoreplot.annotate(targsss[points.index(point)], xy = (point[0], point[1]), xytext = (point[0], point[1]))
        for target, color in zip(targs,colors):
            indicesToKeep = finalDf['Target'] == target
            ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
                    , finalDf.loc[indicesToKeep, 'principal component 2']
                    , c = color
                    , s = 50)
        ax.legend(['Delta', 'Omicron'])
        ax.grid()
        """
        figtwo.savefig("PCA Score Plot.png")
        fig.savefig('PCA Loading Plot.png')

    @staticmethod
    def loader():
        data = []
        # with open('./Radx_data_8_24/Mutant1/AF_220823_220014_0354_VT_TR_Mix  1_0_C11_1.txt', 'r', ) as adep:
        with open(Adeptrix.datafile, 'r' ) as adep:
            datas = csv.reader(adep, delimiter = ' ')
            for row in datas:
                row[0] = round(float(row[0]), 2)
                row[1] = int(float(row[1]))
                data.append(row)
        Adeptrix.rawdata = data
        return data

    @staticmethod
    def plot():
        rawdatas = Adeptrix.rawdata
        negdatas = Adeptrix.negdata
        negpeak = []
        for peak in Adeptrix.finalallpeaks:
            peak[1] = int(peak[1])
            if(rawdatas.index(peak[0:2]) > 24):
                rawdata = rawdatas[rawdatas.index(peak[0:2])-25:rawdatas.index(peak[0:2])+25]
                for info in negdatas:
                    if abs(peak[0] - info[0]) < .5:
                        negpeak = info
                negdata = negdatas[negdatas.index(negpeak)-25:negdatas.index(negpeak)+25]
            else:
                rawdata = rawdatas[0:rawdatas.index(peak[0:2])+25]
                for info in negdatas:
                    if abs(peak[0] - info[0]) < .5:
                        negpeak = info
                negdata = negdatas[0:negdatas.index(negpeak)+25]
            plotter = figure(width = 1000, height = 500)
            negmass = []
            negintens = []
            for row in negdata:
                negmass.append(row[0])
                negintens.append(row[1])
            masses = []
            intensities = []
            for row in rawdata:
                masses.append(row[0])
                intensities.append(row[1])
            plotter.line(masses, intensities, line_width = 2, color = 'orange')
            plotter.line(negmass, negintens, line_width = 2, color = 'gray')
            plotter.yaxis.axis_label = 'Intensity'
            plotter.xaxis.axis_label = 'm/z'
            plotter.title = 'Intensity to Mass-Charge Ratio' + Adeptrix.filename
            from bokeh.models import Label as labe
            from bokeh.models import Span as spanner
            if peak in Adeptrix.finalclearpeaks:
                plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'black', line_dash = 'dashed', line_width = 1))
                plotter.add_layout(labe(x = peak[0], y = 350, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
            if peak in Adeptrix.orange:
                plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'orange', line_dash = 'dashed', line_width = 1))
                plotter.add_layout(labe(x = peak[0], y = 375, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
            if peak in Adeptrix.red:
                plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'red', line_dash = 'dashed', line_width = 1))
                plotter.add_layout(labe(x = peak[0], y = 400, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
            if os.path.exists("./Peak Images/" + str(Adeptrix.filename)):
                output_file("./Peak Images/" + str(Adeptrix.filename) + "/" + str(peak) + ".html")
            else:
                os.makedirs("./Peak Images/" + str(Adeptrix.filename))
                output_file("./Peak Images/" + str(Adeptrix.filename) + "/" + str(peak) + ".html")
            save(plotter)
        rawdata = Adeptrix.rawdata
        negdata = Adeptrix.negdata
        plotter = figure(width = 1000, height = 500)
        negmass = []
        negintens = []
        for row in negdata:
            negmass.append(row[0])
            negintens.append(row[1])
        masses = []
        intensities = []
        for row in rawdata:
            masses.append(row[0])
            intensities.append(row[1])
        plotter.line(masses, intensities, line_width = 2, color = 'orange')
        plotter.line(negmass, negintens, line_width = 2, color = 'gray')
        plotter.yaxis.axis_label = 'Intensity'
        plotter.xaxis.axis_label = 'm/z'
        plotter.title = 'Intensity to Mass-Charge Ratio' + Adeptrix.filename
        from bokeh.models import Label as labe
        from bokeh.models import Span as spanner
        from bokeh.models import Arrow, VeeHead
        for peak in Adeptrix.finalclearpeaks:
            plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'black', line_dash = 'dashed', line_width = 1))
            plotter.add_layout(labe(x = peak[0], y = 350, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
        for peak in Adeptrix.orange:
            plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'orange', line_dash = 'dashed', line_width = 1))
            plotter.add_layout(labe(x = peak[0], y = 375, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
        for peak in Adeptrix.red:
            plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'red', line_dash = 'dashed', line_width = 1))
            plotter.add_layout(labe(x = peak[0], y = 400, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
        if os.path.exists("./Peak Images/" + str(Adeptrix.filename)):
            output_file("./Peak Images/" + str(Adeptrix.filename) + "/" + "Allpeaks.html")
        else:
            os.makedirs("./Peak Images/" + str(Adeptrix.filename))
            output_file("./Peak Images/" + str(Adeptrix.filename) + "/" + "Allpeaks.html")
        save(plotter)

    @staticmethod
    def condenser(data):
        for row in data:
            pass
        return data

    @staticmethod
    def datafitter():
        peaks = Adeptrix.allpeaks
        for stuff in peaks:
            Adeptrix.datafitone(stuff)

    @staticmethod
    def datafitone(point):
        ratios = []
        data = Adeptrix.rawdata
        negdata = Adeptrix.negdata
        indexx = data.index(point)
        more_index = []
        less_index = []
        negmore = []
        negless = []
        if indexx < len(data) - 11:
            for neg in negdata:
                if abs(neg[0] - point[0]) < .5:
                    for num in range(0, 10):
                        negmore.append(negdata[negdata.index(neg)+num])
                        negless.append(negdata[negdata.index(neg)-num])
            negless.extend(negmore)
            negless.sort()
            negmass = []
            negintens = []
            for vals in negless:
                negmass.append(vals[0])
                negintens.append(vals[1])
            negcurves = []
            for nums in range(1, 2):
                negcurve = poly1d(polyfit(negmass, negintens, nums))
                negrscore = r2_score(negintens, negcurve(negmass))
                negcurves.append([negcurve, negrscore])
            scorecount = 0
            for item in negcurves:
                if(item[1] > .8):
                    scorecount += 1
                    ratios.append([point, .8])
        return ratios        

    @staticmethod
    def datacondenser(data):
        obvval = []
        splicepoints = []
        newdata = []
        for row in data:                
            if(row[1] > 500):
                obvval.append([row, data.index(row)])
        for row in obvval:
            for rows in obvval:
                if(abs(row[1] - rows[1]) > 99):
                    splicepoints.append(int(abs(row[1] - rows[1])/2))
        for val in splicepoints:
            newdata.append(data.pop[0:val])
        return newdata

    @staticmethod
    def controllarge(data):
        integral = int(Adeptrix.maxintens/(2000))
        subdata = []
        pot_peak = []
        for integ in range(Adeptrix.maxintens, 500, -integral):
            for stuff in data:
                if(float(stuff[1]) >= integ-integral and float(stuff[1]) <= integ):
                    subdata.append(stuff)
        testmain = []
        testmain = Adeptrix.localmax(subdata, data)        
        test = testmain[1]
        if(test): 
            pot_peak = testmain[0].copy()
        dict = {}
        dict["Peaks"] = pot_peak
        for peak in pot_peak:
            Adeptrix.allpeaks.append(peak)

    @staticmethod
    def controlsmall(data):
        # Logic for smaller peaks that cannot be detected with 'controllarge()' logic. Must plot, smoothen, and find peaks that way
        pass

    @staticmethod
    def maxrevamp(data):
        maindata = Adeptrix.loader()
        newdata = []
        newpeaks = []
        pos_count = 0
        neg_count = 0
        for row in range(0, len(data)):
            pos_count = 0
            neg_count = 0
            for rows in range(0, len(maindata)):
                if(rows > 9 and rows < len(maindata)-10):
                    if(data[row][0] == maindata[rows][0]):
                        for counting in range(1, 10):
                            if(float(maindata[rows][1]) > float(maindata[rows+counting][1])):
                                pos_count += 1
                            if(float(maindata[rows][1]) < float(maindata[rows-counting][1])):
                                neg_count += 1
                        if(pos_count + neg_count >= 18):
                            newpeaks.append(data[row])
        Adeptrix.peaks = newpeaks

    @staticmethod
    def localmax(data, maindata):
        # A bunch of logic to find a local max in the dataset
        newpeaks = []
        finalpeaks = []
        for row in range(0, len(data)):
            pos_count = 0
            neg_count = 0
            for rows in range(0, len(maindata)):
                if(rows > 10 and rows < len(maindata)-11):
                    if(data[row][0] == maindata[rows][0]):
                        for counting in range(1, Adeptrix.peakpointcount+1):
                            if(float(maindata[rows][1]) > float(maindata[rows+counting][1])):
                                pos_count += 1
                            if(float(maindata[rows][1]) > float(maindata[rows-counting][1])):
                                neg_count += 1
                        if(pos_count + neg_count >= 10):
                            newpeaks.append(data[row])
        newpeaks.sort()
        peaktransfer = newpeaks
        for item in peaktransfer:
            finalpeaks.append(item)
        peaks = False
        if len(finalpeaks) > 0:
            peaks = True
            for row in finalpeaks:
                Adeptrix.peaks.append(row)
        for row in newpeaks:
            Adeptrix.gendata.append(row)
        return [Adeptrix.peaks, peaks]

    @staticmethod
    def start():
        Adeptrix.importer()
        Adeptrix.statistics()

class Gui:
    from tkinterdnd2 import Tk
    window = Tk()
    var = StringVar()
    vartwo = StringVar()


    def browseFileBackground(stringvar):
        i = 0
        def test():
            while i == 0:
                if stringvar.get() != '':
                    Adeptrix.negcontrolfile = stringvar.get()
                    break
        test()
        
    def browseFileSample(stringvar):
        i = 0
        def test():
            while i == 0:
                if stringvar.get() != '':
                    data = stringvar.get()
                    data = data.replace("}{",",").replace("{", "").replace("}", "").replace(" /", ",/")
                    Adeptrix.datafiles = data.split(",")
                    break
        test()

    def openNewWindow(filename, address):
            from tkinter import Toplevel
            from tkinter import Label
            from tkinter import Button
            from tkinter import Canvas
            from tkinter import PhotoImage
            from tkinter import Image   
            from tkinter import Grid   
            from tkinter import Listbox
            from tkinterweb import HtmlFrame
            import tkinter as tk
            # Toplevel object which will
            # be treated as a new window
            newWindow = Toplevel(Gui.window)
            # sets the title of the
            # Toplevel widget
            newWindow.title(filename)
        
            # sets the geometry of toplevel
            newWindow.geometry("750x400")
            # A Label widget to show in toplevel
            Label(newWindow,
                text ="Analysis Results of: " + filename).pack()
            Label(newWindow, 
                text = 'Loading Information ...').pack()
            lists = Listbox(newWindow)
            lists.pack()
            files = os.listdir(address)
            while len(files) == 0:
                if len(files) > 0:
                    break
            for items in files:
                lists.insert(tk.END, items)
            def html():
                import webview
                name = lists.get(tk.ANCHOR)
                names = os.getcwd() + '/' + 'Peak Images/' + filename + '/' + name
                webview.create_window('Data', names, fullscreen = False, min_size=(2000, 1000))
                webview.start()
            Button(newWindow, text = 'Open File for Analysis', command = html).pack()

    def openSummary():
        from tkinter import Toplevel
        from tkinter import Label
        from tkinter import Listbox
        from tkinter import Button
        import tkinter as tk
        # Toplevel object which will
        # be treated as a new window
        newWindow = Toplevel(Gui.window)
        # sets the title of the
        # Toplevel widget
        newWindow.title('Data Summary')
        Label(newWindow, text ="This window contains data summarizing all of your analyses.").pack()
        # sets the geometry of toplevel
        newWindow.geometry("1000x533")
        datas = str(Adeptrix.statistics())
        # A Label widget to show in toplevel
        lists = Listbox(newWindow)
        lists.pack()
        files = os.listdir('./')
        while len(files) == 0:
            if len(files) > 0:
                break
        for items in files:
            if items.endswith(".png"):
                lists.insert(tk.END, items)
        def html():
            import webview
            name = lists.get(tk.ANCHOR)
            names = os.getcwd() + '/' + name
            webview.create_window('Data', names, fullscreen = False, min_size=(2000, 1000))
            webview.start()
        Button(newWindow, text = 'Open Analysis', command = html).pack()
        Label(newWindow, text = datas).pack()

    def begin():        
        from tkinter import Label
        from tkinter import Button
        Gui.window.title('Adeptrix Peak Analyzer')
        Gui.window.geometry("1500x800")

        label_file_explorer = Label(Gui.window,
                                    text = "Search Files for Analysis",
                                    width = 100, height = 4,
                                    fg = "blue")
        label_file_explorer.pack()

        from tkinter import Entry
        from tkinterdnd2 import DND_FILES

        def drop(event):
            Gui.var.set(event.data.replace("{", "").replace("}", ""))
        textvarone = StringVar()
        textvartwo = StringVar()
        textvarone.set("Drop Background Files Here")
        textvartwo.set("Drop Sample Files Here")
        e_box = Entry(Gui.window, textvar=textvarone, width=80)
        e_box.pack()
        e_box.drop_target_register(DND_FILES)
        e_box.dnd_bind('<<Drop>>', drop)

        def droptwo(event):
            value = Gui.vartwo.get()
            value += (event.data)
            Gui.vartwo.set(value)
        e_box = Entry(Gui.window, textvariable=textvartwo, width=80)
        e_box.pack()
        e_box.drop_target_register(DND_FILES)
        e_box.dnd_bind('<<Drop>>', droptwo)

        def starter():
            Gui.browseFileBackground(Gui.var)
            Gui.browseFileSample(Gui.vartwo)
            if Adeptrix.negcontrolfile != '' and len(Adeptrix.datafiles) > 0:
                label_analyzing = Label(Gui.window,
                                    text = "The results of your analysis will pop up in a new Window.",
                                    width = 100, height = 4,
                                    fg = "blue")
                label_analyzing.pack()
                Adeptrix.datafiles.sort()
                t = threading.Thread(target=Adeptrix.start)
                t.start()
                t.join()
                for items in Adeptrix.datafiles:
                    filename = (items.split("/")[-1].split('.')[0])
                    if os.path.isdir('./Peak Images/' + items.split("/")[-1].split('.')[0]):
                        path = ('./Peak Images/' + items.split("/")[-1].split('.')[0])
                        threading.Thread(target=Gui.openNewWindow, args = [filename, path]).start()
                threading.Thread(target=Gui.openSummary).start()
        
        button_starter = Button(Gui.window,
                                text = 'Begin Peak Analysis',
                                command = starter)
        button_starter.pack()
        Button(Gui.window, text = "Open Summary Window", command = Gui.openSummary).pack()
        Gui.window.mainloop()

Gui.begin()
# make cut off variable for lowest intensity of peak to detect based on user input
# make cut off variable for lowest mass ratio to detect based on user input
# change ratio to some factor of peakshape so that it can single out more specific/borderline peaks such as 1709 and 2707