import math
import json
import matplotlib as plot
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
from numpy import *
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

    @staticmethod
    def negcontrolfilter():
        import os
        # assign directory
        directory = './Radx_data_8_24/Background'
        for filename in os.listdir(directory):
            f = os.path.join(directory, filename)
            if f.endswith(".txt"):
                rowcount = 0
                data = []
                with open(f, 'r' ) as adep:
                    datas = list(csv.reader(adep, delimiter = ' '))
                    for row in datas:
                        row[0] = float(row[0])
                        row[1] = int(row[1])
                        data.append(row)
                        Adeptrix.negdata.append(row)
                        Adeptrix.rawdata = data
                        rowcount += 1
                rowcount2 = int(rowcount/2000 + 1)
                intensities = []
                for row in Adeptrix.rawdata:
                    intensities.append(int(row[1]))
                Adeptrix.maxintens = max(intensities)
                numthree = 0
                for num in range(0, rowcount2):
                    data = []
                    for numtwo in range(0, 2000):
                        if(numthree < rowcount):
                            data.append([float(datas[numthree][0]), int(datas[numthree][1])])
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
            for peaks in Adeptrix.negdata:
                if abs(peak[0] - peaks[0]) < .5:
                    if peaks[1]/peak[1] >= .3333333333 and peaks[1]/peak[1] <= .618:
                        Adeptrix.possrempeaks.append(peak)
                    elif peaks[1]/peak[1] >= .333333333:
                        Adeptrix.removepeaks.append(peak)
        Adeptrix.peakarea(Adeptrix.possrempeaks, Adeptrix.rawdata, Adeptrix.negdata)
        for peak in Adeptrix.allpeaks:
            if peak not in Adeptrix.removepeaks:
                Adeptrix.finalallpeaks.append(peak)
        Adeptrix.finalallpeaks.sort()

    @staticmethod
    def peakarea(peaks, data, negdata):
        peaks.sort()
        data.sort()
        negdata.sort()
        Adeptrix.negpeaks.sort()
        for peak in peaks:
            breakerpossur = False
            breakernegsur = False
            surroundvalues = []
            surroundvalues.append(peak)
            if((data.index(peak) > 15) and (data.index(peak) < len(data)-16)):
                for num in range(1, 11):
                    poscount = 0
                    negcount = 0
                    posnum = data[data.index(peak)+num][1]
                    negnum = data[data.index(peak)-num][1]
                    for counter in range(1, 6):
                        if(data[data.index(peak)+num+counter][1] > posnum):
                            poscount += 1
                        if(data[data.index(peak)-num-counter][1] > negnum):
                            negcount += 1
                        posnum = data[data.index(peak)+num+counter][1]
                        negnum = data[data.index(peak)-num-counter][1]
                    if((breakerpossur == False) and (poscount >= 5)):
                        breakerpossur = True
                    if(breakerpossur == False):
                        surroundvalues.append(data[data.index(peak)+num])
                    if((breakernegsur == False) and (negcount >= 5)):
                        breakernegsur = True
                    if(breakernegsur == False):
                        surroundvalues.append(data[data.index(peak)-num])
                surroundvalues.sort()
                surroundx = []
                surroundy = []
                for numb in surroundvalues:
                    surroundx.append(numb[0])
                    surroundy.append(numb[1])
                for info in Adeptrix.negdata:
                    breakerposneg = False
                    breakernegneg = False
                    if((negdata.index(info) > 15) and (negdata.index(info) < len(negdata)-16)):
                        if(abs(peak[0] - info[0]) < .5):
                            negsurround = []
                            negsurround.append(info)
                            for numbe in range(1, 11):
                                backposcount = 0
                                backnegcount = 0
                                backposnum = negdata[negdata.index(info)+numbe][1]
                                backnegnum = negdata[negdata.index(info)-numbe][1]
                                for counterr in range(1, 6):
                                    if(negdata[negdata.index(info)+numbe+counterr][1] > backposnum):
                                        backposcount += 1
                                    if(negdata[negdata.index(info)-numbe-counterr][1] > backnegnum):
                                        backnegcount += 1
                                    backposnum = negdata[negdata.index(info)+numbe+counterr][1]
                                    backnegnum = negdata[negdata.index(info)-numbe-counterr][1]
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
                            negintegral = integrate.cumulative_trapezoid(negsurroundy, negsurroundx)
                            if(len(peakintegral) < 1):
                                #print(peak, info)
                                #print(surroundvalues)
                                #print(negsurround)
                                Adeptrix.removepeaks.append(peak)
                            if((len(negintegral) > 0) and (len(peakintegral) > 0)):
                                #print(peak, info)
                                #print(surroundvalues)
                                #print(peakintegral[-1])
                                #print(negsurround)
                                #print(negintegral[-1])
                                if(((negintegral[-1]/peakintegral[-1]) >= .3)):
                                    Adeptrix.removepeaks.append(peak)

    @staticmethod
    def datasplitter():
        directory = './Radx_data_8_24/Mutants'
        for filename in os.listdir(directory):
            f = os.path.join(directory, filename)
            directories = f
            if filename == 'Sample Stuff':
                for filenames in os.listdir(directories):
                    q = os.path.join(directories, filenames)
                    Adeptrix.filepathh = q
                    if q.endswith(".txt"):
                        rowcount = 0
                        data = []
                        with open(q, 'r' ) as adep:
                            datas = list(csv.reader(adep, delimiter = ' '))
                            for row in datas:
                                rowcount += 1
                        rowcount2 = int(rowcount/2000 + 1)
                        numthree = 0
                        Adeptrix.loader()
                        intensities = []
                        for row in Adeptrix.rawdata:
                            intensities.append(int(row[1]))
                        Adeptrix.maxintens = max(intensities)
                        for num in range(0, rowcount2):
                            data = []
                            for numtwo in range(0, 2000):
                                if(numthree < rowcount):
                                    data.append([float(datas[numthree][0]), int(datas[numthree][1])])
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
                        Adeptrix.datafitter()
                        Adeptrix.compare()
                        dict = {}
                        dict['Peaks for ' + filenames] = []
                        if(os.path.exists(directories+'/Peak Data/Potential Peaks.json')):
                            with open(directories+'/Peak Data/Potential Peaks.json', 'r') as inst:
                                dict = json.loads(inst.read()) 
                        dict['Peaks for ' + filenames] = Adeptrix.finalallpeaks
                        if(len(Adeptrix.finalallpeaks) > 0):
                            with open(directories+'/Peak Data/Potential Peaks.json', 'w+') as file:
                                file.write(json.dumps(dict))
                                print("Potential Peaks have been uploaded to the File named 'Potential Peaks.json' inside of the 'Peak Data' folder.")
                                # Develop dictionary to get more specific to mutations and amino acid sequences as program develops further
                        else:
                            with open(directories+'/Peak Data/Potential Peaks.json', 'w+') as file:
                                file.write(json.dumps(dict))
                                print('No Peaks have been found.')
                        Adeptrix.plot()
                        Adeptrix.datamini = []
                        Adeptrix.peaks = []
                        Adeptrix.rawdata = []
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
        with open(Adeptrix.filepathh, 'r' ) as adep:
            datas = csv.reader(adep, delimiter = ' ')
            for row in datas:
                row[0] = float(row[0])
                row[1] = int(row[1])
                data.append(row)
        Adeptrix.rawdata = data
        return data

    @staticmethod
    def plot():
        plotter = figure(width = 2000, height = 1000)
        rawdata = Adeptrix.rawdata
        negdata = Adeptrix.negdata
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
        plotter.title = 'Intensity to Mass-Charge Ratio' + Adeptrix.filepathh
        for val in Adeptrix.finalallpeaks:
            plotter.add_layout(Span(location = val[0], dimension = 'height', line_color = 'black', line_dash = 'dashed', line_width = 1))
            plotter.add_layout(Label(x = val[0], y = 800, y_units = 'screen', text='m/z: ' + str(val[0])))
        output_file(Adeptrix.filepathh + ".html")
        show(plotter)

    @staticmethod
    def condenser(data):
        for row in data:
            pass
        return data

    @staticmethod
    def datafitter():
        peaks = Adeptrix.peaks
        for stuff in peaks:
            Adeptrix.datafitone(stuff)

    @staticmethod
    def datafitone(point):
        data = Adeptrix.rawdata
        indexx = data.index(point)
        more_index = []
        less_index = []
        if indexx < len(data) - 50:
            for num in range(0, 50):
                more_index.append(data[indexx+num])
                less_index.append(data[indexx-num])
            less_index.extend(more_index)
            less_index.sort()
            masses = []
            intensities = []
            for val in less_index:
                masses.append(val[0])
                intensities.append(val[1])
            curve = poly1d(polyfit(masses, intensities, 3))
            rscore = r2_score(intensities, curve(masses))
            if(rscore > .8):
                try:
                    Adeptrix.peaks.remove(point)
                except:
                    print("peak already removed")            

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
                if(int(stuff[1]) >= integ-integral and int(stuff[1]) <= integ):
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
                            if(int(maindata[rows][1]) > int(maindata[rows+counting][1])):
                                pos_count += 1
                            if(int(maindata[rows][1]) < int(maindata[rows-counting][1])):
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
                        for counting in range(1, 11):
                            if(int(maindata[rows][1]) > int(maindata[rows+counting][1])):
                                pos_count += 1
                            if(int(maindata[rows][1]) < int(maindata[rows-counting][1])):
                                neg_count += 1
                        if(pos_count + neg_count >= 10):
                            newpeaks.append(data[row])
        newpeaks.sort()
        peaktransfer = newpeaks
        pointer = 0
        oldindex = 0
        listcounter = 0
        listlength = 0
        for index in range(0, len(peaktransfer)):
            for index in range(0, len(newpeaks)-1):
                if(pointer < len(newpeaks)-1):
                    if((int(newpeaks[pointer+1][1]) - int(newpeaks[pointer][1])) > 300):
                        newlist = newpeaks[oldindex:pointer+1]
                        peakintensities = []
                        for val in newlist:
                            peakintensities.append(val[1])
                        peakintensities.sort()
                        last = max(peakintensities)
                        for val in newlist:
                            if(int(val[1]) == last):
                                finalpeaks.append(val)
                        oldindex = pointer
                        listlength += len(newlist)
                        newlist = []
                    listcounter += 1
                    pointer += 1
        if(listlength < len(newpeaks)):
            newlist = newpeaks[oldindex:len(newpeaks)]
            peakintensities = []
            for val in newlist:
                peakintensities.append(val[1])
            peakintensities.sort()
            last = max(peakintensities)
            for val in newlist:
                if(int(val[1]) == last):
                    finalpeaks.append(val)
        peaks = False
        if len(finalpeaks) > 0:
            peaks = True
            for row in finalpeaks:
                Adeptrix.peaks.append(row)
        for row in newpeaks:
            Adeptrix.gendata.append(row)
        return [Adeptrix.peaks, peaks]

Adeptrix.negcontrolfilter()
Adeptrix.datasplitter()



# see shape of peaks, maybe change peakarea allowance from 5 points down to 3 points