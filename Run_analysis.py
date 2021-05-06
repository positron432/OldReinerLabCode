# STEPS
# 1. open abf file
# 2. find open current
# 3. scan for shift to particle
# 4. delete everything prior to particle capture


import pyabf
import pyabf.filter
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from colorama import Fore


def main():
    print(Fore.WHITE)
    printHeader()
    trace, traceTime, startTime = loadTrace("4mba.abf")



    plt.figure(figsize=(8, 5))
    plt.plot(traceTime, trace)
    plt.show()

    openCurrent = findOpenCurrent(trace, startTime)
    print(f"Open current is {Fore.GREEN}{openCurrent} pA{Fore.WHITE}.")
    print()


    particleStart = findParticle(trace, startTime, openCurrent)
    print(f"Particle start occured at {Fore.GREEN}{particleStart / 50000} seconds{Fore.WHITE}.")
    print()

    # Delete all data prior to capture of particle and make this a new array
    trace2 = trace[particleStart + 50000: trace.size]
    traceTime2 = traceTime[particleStart + 50000: traceTime.size]

    plt.plot(traceTime2, trace2)
    plt.show()

    particleExit = findParticleExit(trace2, openCurrent)
    print(f"Particle exit occured at {Fore.GREEN}{particleExit / 50000} seconds{Fore.WHITE}.")
    print()

    trace3 = trace2[0:particleExit]
    traceTime3 = traceTime2[0:particleExit]

    plt.plot(traceTime3, trace3)
    plt.show()

    deltas = buildDeltas(trace3)

    plt.plot(traceTime3, deltas)
    plt.show()

    deltasExpanded = expandDeltas(deltas)

    traceTime4 = np.zeros(deltasExpanded.size)
    for x in range (0, traceTime4.size):
        traceTime4.itemset(x, x)

    plt.plot(traceTime4, deltasExpanded)
    plt.show()

    traceExpanded = integrateDeltas(deltasExpanded, trace3.item(0))

    plt.plot(traceTime4, traceExpanded)
    plt.show()
    print()
    print("Code finished without catastrophic errors.")


# Loads the abf file and returns a ndArray (numpy) of the amperages at each point
def loadTrace(fileName):
    print(f"Loading file {Fore.GREEN}{fileName}{Fore.WHITE}...")
    print()
    abf = pyabf.ABF(fileName)
    print(Fore.GREEN)
    print(abf)
    print(Fore.WHITE)
    abf.setSweep(0, absoluteTime=True)
    print(f"Loaded a trace of {Fore.CYAN}{abf.sweepY.size}{Fore.WHITE} points.")
    print()
    print(
        f"Voltage was turned on at {Fore.GREEN}{abf.tagTimesSec[0]}{Fore.WHITE} seconds, which is point {Fore.GREEN}{int(abf.tagTimesSec[0] * 50000)}{Fore.WHITE}.")
    print()
    print("Now running a 100 Hz Gaussian filter over the trace please wait. This is likely the longest step...")
    print()

    pyabf.filter.gaussian(abf, 10)
    return abf.sweepY, abf.sweepX, abf.tagTimesSec[0] * 50000


def findOpenCurrent(trace, startTime):
    # scan average of next 5 seconds after voltage applied
    # define this average as the open current and return it
    print("Now estimating open current...")
    print()
    return np.average(trace[int(startTime + 5000):int(startTime + 5000 * 5)])


def findParticle(trace, startTime, openCurrent):
    # look for when current is below 70%
    # check half a second later to confirm this isn't a spike
    # return particle start index
    print("Looking for when the particle entered the pore. This may take some time....")
    print()
    particleStart = 0
    for x in range(int(startTime + 5000), trace.size):
        if trace.item(x) / openCurrent < 0.7:
            if trace.item(x + 1 * 25000) / openCurrent < 0.85:
                print(
                    f"Found possible particle at point {Fore.GREEN}{x}{Fore.WHITE} with value {Fore.GREEN}{trace.item(x)} pA{Fore.WHITE}.")
                print()
                particleStart = x
                break
    return particleStart


def findParticleExit(trace2, openCurrent):
    print("Looking for when the particle exited the pore. This may take some time....")
    print()
    particleExit = 0
    for x in range(0, trace2.size):
        if trace2.item(x) / openCurrent > 0.9:
            if trace2.item(x + 1 * 25000) / openCurrent > 0.9:
                print(
                    f"Found possible particle exit at {Fore.GREEN}{x}{Fore.WHITE}with value{Fore.GREEN}{trace2.item(x)}{Fore.WHITE} pA.")
                print()
                particleExit = x - 2500
                break
    return particleExit


def buildDeltas(trace3):
    deltas = np.zeros(trace3.size)
    print("Now making deltas table. (taking the derivative)")
    print()
    for x in range(0, trace3.size - 1):
        deltas.itemset(x, trace3.item(x + 1) - trace3.item(x))
    return deltas


def expandDeltas(deltas):
    print("Now expanding deltas... (cubing the derivative)")
    print()
    deltas = shrinkTrace(deltas)
    deltasExpanded = np.zeros(deltas.size)
    for x in range(0, deltas.size):
        deltasExpanded.itemset(x, deltas.item(x) * deltas.item(x) * deltas.item(x))
    # deltasExpanded = cleanExpanded(deltasExpanded)
    return deltasExpanded


def integrateDeltas(deltasExpanded, startingPoint):
    print("Now integrating expanded deltas...")
    print()
    traceExpanded = np.zeros(deltasExpanded.size)
    traceExpanded.itemset(0, startingPoint)
    for x in range(1, deltasExpanded.size):
        traceExpanded.itemset(x, traceExpanded.item(x - 1) + deltasExpanded.item(x))
    return traceExpanded


def cleanExpanded(expanded):
    print("Now attempting to clean trace. This is still broken.")
    print()
    average = np.var(expanded)
    for x in range(0, expanded.size):
        if expanded.item(x) > 3 * average:
            expanded.itemset(x, 3 * average)
    return expanded


def printHeader():
    print("------------------------------")
    print("    Cluster Pore Analysis     ")
    print("         Version 0.5          ")
    print("  Made for Reiner Lab by BDC  ")
    print("------------------------------")
    print()
    print()
    print("Now executing code...")
    return True


def shrinkTrace(trace1):
    counter = 0
    newTraceCounter = 0
    tempArray = np.zeros(100)
    newTrace = np.zeros(int(trace1.size / 100))
    for x in range(0, trace1.size):
        if counter < 100:
            tempArray.itemset(counter, trace1.item(x))
        else:
            newTrace.itemset(newTraceCounter, np.average(tempArray))
            counter = 0
            newTraceCounter += 1
    return newTrace

if __name__ == '__main__':
    main()
