import math
import csv
import json
from tkinter import Canvas, PhotoImage
import csv
import time
import pandas
import zlib
import sys
import os
import threading
from bokeh.plotting import *
from bokeh.models import *
from sklearn.metrics import *
from numpy import poly1d
from numpy import polyfit
from scipy import integrate

class Adeptrix:

    datamini = []
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
    finalclearpeaks = []
    negcontrolfile = ''
    datafile = ''
    datafiles = []
    ratios = []
    filename = ''
    peakpointcount = 5
    alignmentfile = ''

    @staticmethod
    def alignment():
        maxintens = 0
        maxintensxval = 0
        adjustment = 0
        found = False
        data = []
        Adeptrix.datafiles = list(dict.fromkeys(Adeptrix.datafiles))
        for file in Adeptrix.datafiles:
            intens = []
            with open(file, 'r' ) as adep:
                datas = list(csv.reader(adep, delimiter = ' '))
                for row in datas:
                    row[0] = float(row[0])
                    row[1] = float(row[1])
                    intens.append(row[1])
                    data.append(row)
        maxintens = max(intens)
        for row in data:
            if maxintens in row:
                maxintensxval = row[0]
        for file in Adeptrix.datafiles:
            with open(file, 'r' ) as adep:
                datas = list(csv.reader(adep, delimiter = ' '))
                for row in datas:
                    row[0] = float(row[0])
                    row[1] = float(row[1])
                    if maxintens == row[1]:
                        Adeptrix.alignmentfile = file
                        found = True
                        break
            if found == True:
                break
        datastoanalyze = []
        for file in Adeptrix.datafiles:
            data2 = []
            with open(file, 'r' ) as adep:
                datas = list(csv.reader(adep, delimiter = ' '))
                for row in datas:
                    row[0] = float(row[0])
                    row[1] = float(row[1])
                    data2.append(row)
            adep.close()
            diff = []
            for row in data2:
                diff.append([abs(row[0] - maxintensxval), row[0], maxintensxval])
            adjustment = min(diff)[0]
            if(min(diff)[2] > min(diff)[1]):
                adjustment = min(diff)[0]
            elif(min(diff)[2] < min(diff)[1]):
                adjustment = -min(diff)[0]
            else:
                adjustment = 0
            for row in data2:
                row[0] = row[0] + adjustment
            with open (file + "edited.txt", 'w+') as filewriter:
                for row in data2:
                    filewriter.write(str(row[0]) + ' ' + str(row[1]) + '\n')
            datastoanalyze.append(file + 'edited.txt')
        datastoanalyze = list(dict.fromkeys(datastoanalyze))
        for file in datastoanalyze:
            Adeptrix.datafile = file
            Adeptrix.negcontrolfilter()
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
            if(int(peak[1]) < 1000):
                Adeptrix.removepeaks.append(peak)
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
                thread.join()
        for peak in Adeptrix.allpeaks:
            peaker = peak
            for peakss in Adeptrix.ratios:
                if peak not in Adeptrix.orange and peak not in Adeptrix.red:
                    if peakss[0] == peaker[0]:
                        peaker.append(peakss[2])
            if peak not in Adeptrix.removepeaks:
                Adeptrix.finalallpeaks.append(peaker)
                if peak not in Adeptrix.orange and peak not in Adeptrix.red:
                    peaker.append('Clear Peak')
                    Adeptrix.finalclearpeaks.append(peaker)
        Adeptrix.finalallpeaks.sort()
        Adeptrix.red = [list(t) for t in set(tuple(element) for element in Adeptrix.red)]
        Adeptrix.orange = [list(t) for t in set(tuple(element) for element in Adeptrix.orange)]
        Adeptrix.finalclearpeaks = [list(t) for t in set(tuple(element) for element in Adeptrix.finalclearpeaks)]
        Adeptrix.finalallpeaks = [list(t) for t in set(tuple(element) for element in Adeptrix.finalallpeaks)]
        Adeptrix.red.sort()
        Adeptrix.orange.sort()
        Adeptrix.finalclearpeaks.sort()

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
                    for num in surroundless:
                        if surroundless.index(num) < len(surroundless)-1:
                            for numb in range(1,2):
                                if(surroundless[surroundless.index(num)+numb][1] - num[1] >= .3*distance):
                                    Adeptrix.removepeaks.append(peak)
                    for num in surroundgreater:
                        if surroundgreater.index(num) < len(surroundgreater)-1:
                            for numb in range(1,2):
                                if(surroundgreater[surroundgreater.index(num)+numb][1] - num[1] >= .3*distance):
                                    Adeptrix.removepeaks.append(peak)


    @staticmethod
    def peakwidth(peaks, data):
        for peak in peaks:
            nearvals = []
            peakwratio = []
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
                                    #if((peak[0] < 1710 and peak[0] > 1709)):
                                    #    print(peak, info)
                                    #    print(surroundvalues)
                                    #    print(peakintegral[-1])
                                    #    print(negsurround)
                                    #    print(negintegral[-1])
                                    #    print(((negintegral[-1]/peakintegral[-1])))
                                    if(len(peakintegral) < 1):
                                        #print(peak, info)
                                        #print(surroundvalues)
                                        #print(negsurround)
                                        Adeptrix.removepeaks.append(peak)
                                    if(float(peak[1]) < 1000):
                                        Adeptrix.removepeaks.append(peak)
                                    if((len(negintegral) > 0) and (len(peakintegral) > 0)):
                                        #print(peak, info)
                                        #print(surroundvalues)
                                        #print(peakintegral[-1])
                                        #print(negsurround)
                                        #print(negintegral[-1])
                                        #print(((negintegral[-1]/peakintegral[-1]) >= .3))
                                        peaker = peak
                                        if len(peaker) < 3:
                                            peaker.append(peakintegral[-1]/negintegral[-1])
                                    elif((len(negintegral) == 0) and (len(peakintegral) > 0)):
                                        peaker = peak
                                        if len(peaker) < 3:
                                            peaker.append('Null')
                                            peaker.append('Highly Questionable')
                                            Adeptrix.red.append(peak)
                                    Adeptrix.ratios.append(peaker)
                                    if((len(negintegral) > 0) and (len(peakintegral) > 0)):
                                        if((peak[1] > 5000) and peak[1]/info[1] > 2.222 and (abs(((negintegral[-1]/peakintegral[-1]) < .50)) 
                                            and abs(negintegral[-1]/peakintegral[-1] > .3))):
                                            if peak not in safe:
                                                safe.append(peak)
                                            if peakintegral[-1]/negintegral[-1] < 2.5:
                                                if peak not in Adeptrix.red:
                                                    peaker.append('Highly Questionable')
                                                    Adeptrix.red.append(peaker)
                                            elif peakintegral[-1]/negintegral[-1] < 3:
                                                if peak not in Adeptrix.orange:
                                                    peaker.append('Questionable')
                                                    Adeptrix.orange.append(peaker)
                                        elif(abs(((negintegral[-1]/peakintegral[-1]) > .3))):
                                            if peak not in safe:
                                                Adeptrix.removepeaks.append(peak)
                                    elif((len(negintegral) > 0) and len(peakintegral) == 0):
                                        Adeptrix.removepeaks.append(peak)
                                    else:
                                        peaker.append("Null")
                                        peaker.append('Highly Questionable')
                                        Adeptrix.red.append(peak)

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
                                    #if((peak[0] < 1710 and peak[0] > 1709)):
                                    #    print(peak, info)
                                    #    print(surroundvalues)
                                    #    print(peakintegral[-1])
                                    #    print(negsurround)
                                    #    print(negintegral[-1])
                                    #    print(((negintegral[-1]/peakintegral[-1])))
                                    if(len(peakintegral) < 1):
                                        #print(peak, info)
                                        #print(surroundvalues)
                                        #print(negsurround)
                                        Adeptrix.removepeaks.append(peak)
                                    if(float(peak[1]) < 1000):
                                        Adeptrix.removepeaks.append(peak)
                                    if((len(negintegral) > 0) and (len(peakintegral) > 0)):
                                        #print(peak, info)
                                        #print(surroundvalues)
                                        #print(peakintegral[-1])
                                        #print(negsurround)
                                        #print(negintegral[-1])
                                        #print(((negintegral[-1]/peakintegral[-1]) >= .3))
                                        peaker = peak
                                        if len(peaker) < 3:
                                            peaker.append(peakintegral[-1]/negintegral[-1])
                                    if((len(negintegral) == 0) and (len(peakintegral) > 0)):
                                        peaker = peak
                                        if len(peaker) < 3:
                                            peaker.append('Null')
                                            peaker.append('Highly Questionable')
                                    Adeptrix.ratios.append(peaker)
                                    if((len(negintegral) > 0) and (len(peakintegral) > 0)):
                                        if((peak[1] > 5000) and peak[1]/info[1] > 2.222 and (abs(((negintegral[-1]/peakintegral[-1]) < .50)) 
                                            and abs(negintegral[-1]/peakintegral[-1] > .3))):
                                            if peak not in safe:
                                                safe.append(peak)
                                            if peakintegral[-1]/negintegral[-1] < 2.5:
                                                if peak not in Adeptrix.red:
                                                    peaker.append('Highly Questionable')
                                                    Adeptrix.red.append(peaker)
                                            elif peakintegral[-1]/negintegral[-1] < 3:
                                                if peak not in Adeptrix.orange:
                                                    peaker.append('Questionable')
                                                    Adeptrix.orange.append(peaker)
                                        elif(abs(((negintegral[-1]/peakintegral[-1]) > .3))):
                                            if peak not in safe:
                                                Adeptrix.removepeaks.append(peak)
                                    elif((len(negintegral) > 0) and len(peakintegral) == 0):
                                        Adeptrix.removepeaks.append(peak)
                                    else:
                                        peaker.append("Null")
                                        peaker.append('Highly Questionable')
                                        Adeptrix.red.append(peak)




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
                    data.append([float(datas[numthree][0]), float(datas[numthree][1])])
                    numthree += 1
            t = threading.Thread(target=Adeptrix.controllarge, args = [data]).start()
        for t in threading.enumerate():
            if t != threading.enumerate()[0]:
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
        dict['Peaks for ' + Adeptrix.datafile] = []
        if(os.path.exists(directories)):
            with open(directories+'/Potential Peaks.json', 'r') as inst:
                dict = json.loads(inst.read())
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
        newlisttot.sort()
        dict5 = {}
        dict5['All Peaks'] = newlisttot
        dict6 = dict5['All Peaks']
        cols = ['Mass/Charge Ratio', 'Intensity', 'Intensity Ratio', 'Peak Confidence']
        if os.path.exists("./Peak Images/" + str(Adeptrix.filename)):
            csv_file = "./Peak Images/" + str(Adeptrix.filename) + "/PeakDataTable.csv"            
        else:
            os.mkdir("./Peak Images/" + str(Adeptrix.filename))
            csv_file = "./Peak Images/" + str(Adeptrix.filename) + "/PeakDataTable.csv"
        with open(csv_file, 'w') as adep:
            try:
                write = csv.writer(adep)
                write.writerow(cols)
                write.writerows(dict6)
            except IOError:
                print("I/O error")
        if(len(Adeptrix.finalallpeaks) > 0):
            with open(directories+'/Potential Peaks.json', 'w+') as file:
                json.dump(dict, file, indent = 2)
                print("Potential Peaks have been uploaded to the File named 'Potential Peaks.json' inside of the 'Peak Data' folder.")
                # Develop dictionary to get more specific to mutations and amino acid sequences as program develops further
        else:
            with open(directories+'/Potential Peaks.json', 'w+') as file:
                json.dump(dict, file, indent = 2)
                print('No Peaks have been found.')
        Adeptrix.plot()
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
    def loader():
        data = []
        # with open('./Radx_data_8_24/Mutant1/AF_220823_220014_0354_VT_TR_Mix  1_0_C11_1.txt', 'r', ) as adep:
        with open(Adeptrix.datafile, 'r' ) as adep:
            datas = csv.reader(adep, delimiter = ' ')
            for row in datas:
                row[0] = float(row[0])
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
            plotter = figure(width = 2000, height = 1000)
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
            if peak in Adeptrix.finalallpeaks:
                plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'black', line_dash = 'dashed', line_width = 1))
                plotter.add_layout(labe(x = peak[0], y = 700, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
            if peak in Adeptrix.orange:
                plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'orange', line_dash = 'dashed', line_width = 1))
                plotter.add_layout(labe(x = peak[0], y = 750, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
            if peak in Adeptrix.red:
                plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'red', line_dash = 'dashed', line_width = 1))
                plotter.add_layout(labe(x = peak[0], y = 800, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
            from bokeh.io import export_svg
            plotter.output_backend = "svg"
            if os.path.exists("./Peak Images/" + str(Adeptrix.filename)):
                output_file("./Peak Images/" + str(Adeptrix.filename) + "/" + str(peak) + ".html")
            else:
                os.mkdir("./Peak Images/" + str(Adeptrix.filename))
                output_file("./Peak Images/" + str(Adeptrix.filename) + "/" + str(peak) + ".html")
            save(plotter)
        rawdata = Adeptrix.rawdata
        negdata = Adeptrix.negdata
        plotter = figure(width = 2000, height = 1000)
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
        for peak in Adeptrix.finalallpeaks:
            plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'black', line_dash = 'dashed', line_width = 1))
            plotter.add_layout(labe(x = peak[0], y = 700, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
        for peak in Adeptrix.orange:
            plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'orange', line_dash = 'dashed', line_width = 1))
            plotter.add_layout(labe(x = peak[0], y = 750, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
        for peak in Adeptrix.red:
            plotter.add_layout(spanner(location = peak[0], dimension = 'height', line_color = 'red', line_dash = 'dashed', line_width = 1))
            plotter.add_layout(labe(x = peak[0], y = 800, y_units = 'screen', text=str(peak[0]), text_font_size = '8pt'))
        if os.path.exists("./Peak Images/" + str(Adeptrix.filename)):
            output_file("./Peak Images/" + str(Adeptrix.filename) + "/" + "Allpeaks.html")
        else:
            os.mkdir("./Peak Images/" + str(Adeptrix.filename))
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
        Adeptrix.alignment()

class Gui:
    def browseFileBackground():
        from tkinter import filedialog  
        filename = str(filedialog.askopenfilename())
        Adeptrix.negcontrolfile = filename
        
    def browseFileSample():
        from tkinter import filedialog
        filename = str(filedialog.askopenfilename())
        Adeptrix.datafiles.append(filename)

    def begin():        
        from tkinter import Tk
        from tkinter import Label
        from tkinter import Button
        from tkinter import Canvas
        from tkinter import PhotoImage
        from tkinter import Image   
        from tkinter import Grid   
        from PIL import ImageTk                            
        window = Tk()
        window.title('Adeptrix Peak Analyzer')
        window.canvas = Canvas(width = 300, height = 300)
        window.canvas.place(x = 500, y = 250)
        window.geometry("1500x800")
        canvas = Canvas(width=600, height=800, bg='blue')
        image = ImageTk.PhotoImage(file="background.gif")
        canvas.create_image(10, 10, image=image)

        label_file_explorer = Label(window,
                                    text = "Search Files for Analysis",
                                    width = 100, height = 4,
                                    fg = "blue")

        button_explorer = Button(window,
                                text = "Select Background File",
                                command = Gui.browseFileBackground)

        button_explorer2 = Button(window,
                            text = "Select Sample File",
                            command = Gui.browseFileSample)

        label_file_explorer.grid()
        button_explorer.grid()
        button_explorer2.grid()

        def starter():
            if Adeptrix.negcontrolfile != '' and len(Adeptrix.datafiles) > 0:
                label_analyzing = Label(window,
                                    text = "Analyzing Peaks, Please Wait...",
                                    width = 100, height = 4,
                                    fg = "blue")
                label_analyzing.grid()
                Adeptrix.start()
        button_starter = Button(window,
                                text = 'Begin Peak Analysis',
                                command = starter)
        button_starter.grid()
        window.mainloop()


Gui.begin()
# make cut off variable for lowest intensity of peak to detect based on user input
# make cut off variable for lowest mass ratio to detect based on user input
# change ratio to some factor of peakshape so that it can single out more specific/borderline peaks such as 1709 and 2707