"""
File for reading Caesium lab data
"""

from math import sin
from math import pi
from math import exp
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import arange

DoubleDip = "C:\\Users\\maxim\OneDrive\\Caesium Lab Data 2\\Double Peak Data.csv"
secondDip = "C:\\Users\\maxim\OneDrive\\Caesium Lab Data 2\\Second Dip.csv"
smallPeak = "C:\\Users\\maxim\OneDrive\\Caesium Lab Data 2\\Small Peak.csv"
firstDip = "C:\\Users\\maxim\OneDrive\\Caesium Lab Data 2\\First Dip.csv"

def readData(path):
    timeData = []
    sinData = []
    absData = []
    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if (line_count >= 2):
                timeData.append(float(row[0]))
                sinData.append(float(row[2]))
                absData.append(float(row[1]))
            line_count += 1
    return timeData, sinData, absData

def prunedData(file):
    timeData, sinData , absData = readData(file)
    ##return timeData, sinData, absData
    fileLength = len(timeData)

    if (file == DoubleDip):
        return timeData[int(fileLength*0.10):int(fileLength*0.90)], sinData[int(fileLength*0.10):int(fileLength*0.90)],absData[int(fileLength*0.10):int(fileLength*0.90)]
    if (file == secondDip):
        return timeData[int(fileLength*0.10):int(fileLength*0.90)], sinData[int(fileLength*0.10):int(fileLength*0.90)],absData[int(fileLength*0.10):int(fileLength*0.90)]

    if (file == smallPeak):
        return timeData[int(fileLength*0.10):int(fileLength*0.90)], sinData[int(fileLength*0.10):int(fileLength*0.90)],absData[int(fileLength*0.10):int(fileLength*0.90)]

    if (file == firstDip):
        return timeData[int(fileLength*0.10):int(fileLength*0.90)], sinData[int(fileLength*0.10):int(fileLength*0.90)],absData[int(fileLength*0.10):int(fileLength*0.90)]


def convertTimeToFrequency(file, converter):

    freqData = []
    time, _, _  = prunedData(file)
    offSet = time[0]
    for t in time:
        freqData.append((t-offSet)*converter)
    return freqData

def normalizeData(dataList):
    maxValue = max(dataList)
    newList = [entry/maxValue for entry in dataList]
    return newList

def plotRawData(file):
    time, sin, abs = readData(file)
    time, sin, abs = time, normalizeData(sin), normalizeData(abs)
    plt.plot(time, sin)
    plt.plot(time, abs)
    plt.xlabel("Time (s)")
    plt.ylabel("Light Intensity (arb. units)")
    plt.legend(("Interferometer Intensity", "Transmitted Intensity"))
    plt.savefig("smallPeakRaw.png")
    plt.show()

plotRawData(smallPeak)

def plotData(file, period):
    time, sin, abs = prunedData(file)
    abs = normalizeData(abs)

    converter = 393416133.5/period

    newTime = convertTimeToFrequency(file, converter )
    plt.plot(newTime, abs)
    plt.xlabel("Frequency Deviation (Hz)")
    plt.ylabel("Transmittance Through Caesium (arb. units)")
    plt.savefig(file[0:len(file)-4]+".png")
    plt.clf()
    ##plt.show()
DoubleDipPeriod = 0.002208004
smallPeakPeriod = 0.002167203
secondDipPeriod = 0.001002348
firstDipPeriod = 0.001049945
plotData(DoubleDip, DoubleDipPeriod)
plotData(smallPeak, smallPeakPeriod)
plotData(secondDip, secondDipPeriod)
plotData(firstDip, firstDipPeriod)


def sinObjective(x, amp, phase, omega, linegrad, lineint):
    ##return x * linegrad + lineint

    ##return amp * phase*x + linegrad * x**2 + lineint
    return [(amp * sin(xi*omega + phase) + xi * linegrad + lineint) for xi in x]

def funcOmegaUncert(file):
    timeData, sinData, absData = readData(file)
    ##popt, thingo = curve_fit(sinObjective, timeData, sinData, p0=[10,pi/2-2.84564048e+03*(-0.04698),2.84564048e+03, 0.01, 5])
    popt, thingo = curve_fit(sinObjective, timeData, sinData, p0=[10,0,2*pi*1000, 0.01, 5])
    ##popt, thingo = curve_fit(sinObjective, timeData, sinData,maxfev=10000)
    amp, phase, omega, linegrad, lineint = popt
    print(popt)
    print(thingo)
    x_line = timeData
    print(x_line[0])
    print("omega", omega)
    print("omega varience", thingo[2][2]**2)
    print("percent omega uncert", thingo[2][2]**2/omega * 100)

    print("here is len timeData", len(timeData))
    # calculate the output for the range
    y_line = sinObjective(x_line, amp, phase, omega, linegrad, lineint)
    plt.plot(x_line, y_line)
    plt.plot(timeData, sinData)
    plt.show()
    ##plt.scatter(timeData, sinData)
    ##plt.show()

def gaussianObjective(x, amp, center, std, verticalShift):
    return [verticalShift*0.0001 + abs(amp)*exp(-((xi-center)**2)/(2*std**2)) for xi in x]

def gaussianFit(timeStart, timeEnd, file):
    timeData, sinData, absData = readData(file)
    approximateCenter = (timeStart + timeEnd)/2
    print("approximateCenter", approximateCenter)
    startIndex = None
    endIndex = None
    counter = 0
    while True:
        if (timeData[counter] <= timeStart and timeData[counter+1] >= timeStart):
            startIndex=counter
        if (timeData[counter] >= timeEnd and timeData[counter - 1] <= timeEnd):
            endIndex = counter
        counter += 1
        if (counter == len(timeData)-2):
            print("captain we have a problem")
            startIndex = 2
            endIndex = len(timeData) - 2
        if (startIndex != None and endIndex != None):
            break
    popt, thingo = curve_fit(gaussianObjective, timeData[startIndex:endIndex], absData[startIndex:endIndex], p0=[3.5*10**(-3), approximateCenter, 8*10**(-5), 10**2],maxfev=100000)
    amp, center, std, verticalShift = popt
    print(popt)
    print("gaussian center of", popt[1])
    print("uncertainty in center is", 2**0.5 * abs(std) )
    print("percent uncertainty in gaussian center is", 2**0.5 * std/center * 100)
    print("gaussian standardDeviation of", thingo[1][1]**2)
    x_line = timeData[startIndex:endIndex]
    y_line=gaussianObjective(x_line, amp, center, std, verticalShift)
    plt.plot(x_line, y_line)
    plt.plot(timeData[startIndex:endIndex], absData[startIndex:endIndex])
    plt.show()


def getAllGaussians(intervals, file):
    counter = 1
    for interval in intervals:
        print("peak number", counter)
        counter += 1
        gaussianFit(interval[0], interval[1], file)


## peak locations used in calibrating
f33=[1.870*10**(-2),1.914*10**(-2)]
f34=[-3.288*10**(-2), -3.253*10**(-2)]

## reading small data peak
f32 = [0.01735,0.017550]
f32WITHf33 = [0.018050,0.01840]
f33 = [0.01869,0.019010]

## second dip data
f33s = [0.000024+1.433*10**(-1),0.000037+1.434*10**(-1)]
f34s = [0.000058 +1.435*10**(-1),0.000008+1.4371*10**(-1)]
f34Plusf34s = [0.000066+1.438*10**(-1),0.000068 + 1.439*10**(-1)]
f44=[0.000045+1.441*10**(-1),0.144211]

## first Dip

f43f = [0.00008+1.186*10**(-1), 0.00004+1.187*10**(-1)]
f43Plusf44f = [0.118910, 0.0000355+1.19*10**(-1)]
f44f = [0.000035+1.191*10**(-1), 0.000090+1.192*10**(-1)]
f44Plusf45f = [0.1193920,0.000030+1.195*10**(-1)]
f45f = [0.00003+1.196*10**(-1),0.11979]






#funcOmegaUncert(firstDip)
##getAllGaussians([f45f], firstDip)

##f33 got
