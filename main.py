import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import numpy
import functools
from functools import reduce
import math

centerNum = -1
averagesArray = []
centerArray = []


def averages(accumulator, current):
    currentIndex = len(accumulator) - 1
    if currentIndex == -1:
        accumulator.append((current[0], current[2]))
    else:
        if accumulator[currentIndex][0] == current[0]:
            tupleList = list(accumulator[currentIndex])
            tupleList[1] += current[2]
            accumulator[currentIndex] = tuple(tupleList)
        else:
            tupleList = list(accumulator[currentIndex])
            tupleList[1] /= 5068
            accumulator[currentIndex] = tuple(tupleList)
            accumulator.append((current[0], current[2]))
    return accumulator


def center(curElement):
    global centerNum
    global averagesArray
    centerNum += 1
    return curElement[0], curElement[1], curElement[2] - averagesArray[int(centerNum / 5068)][1]


def variance(accumulator, current):
    currentIndex = len(accumulator) - 1
    if currentIndex == -1:
        accumulator.append((current[0], pow(current[2], 2)))
    else:
        if accumulator[currentIndex][0] == current[0]:
            tupleList = list(accumulator[currentIndex])
            tupleList[1] += pow(current[2], 2)
            accumulator[currentIndex] = tuple(tupleList)
        else:
            tupleList = list(accumulator[currentIndex])
            tupleList[1] /= 5068
            accumulator[currentIndex] = tuple(tupleList)
            accumulator.append((current[0], current[2]))
    return accumulator


if __name__ == '__main__':
    # startingData

    df = pd.read_table('ekspresije.tsv', index_col=0)

    data = [(cell, gene, value) for cell in df.columns
            for gene, value in df[cell].items()]

    # 1.1
    averagesArray = reduce(averages, data, [])
    # 1.2
    centerArray = list(map(center, data))
    centerNum = -1
    # 1.3
    varianceArray = reduce(variance, centerArray, [])
    # 1.4
    standardDeviationArray = list(map(lambda x: (x[0], math.sqrt(x[1])), varianceArray))
    standardDeviationDict = dict(map(lambda x: (x[0], math.sqrt(x[1])), varianceArray))
    # 1.5
    standardisedValuesArray = list(map(lambda x: (x[0], x[1], x[2] / standardDeviationDict[x[0]]), centerArray))

    # 2.1

