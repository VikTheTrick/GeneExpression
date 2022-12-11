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
varianceArray = []
standardisedValuesArray = []
geneAveragesDict = dict()
geneCenterArray = []
mostVariableGenesDict = dict()


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
            tupleList[1] /= 534
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
            tupleList[1] /= 534
            accumulator[currentIndex] = tuple(tupleList)
            accumulator.append((current[0], current[2]))
    return accumulator


def geneAverage(accumulator, current):
    if not (current[1] in accumulator):
        accumulator[current[1]] = [current[2], 1]
    else:
        accumulator[current[1]] = [accumulator[current[1]][0] + current[2], accumulator[current[1]][1] + 1]
        if accumulator[current[1]][1] == 5068:
            accumulator[current[1]][0] = accumulator[current[1]][0] / 5068

    return accumulator


def topG(accumulator, current):
    if accumulator[1] < 500:
        accumulator[0][current[0]] = current[1]
        accumulator[1] += 1
    return accumulator


def geneFilter(accumulator, current):
    if current[1] in mostVariableGenesDict:
        accumulator.append(current)
    return accumulator


def rankGenes(accumulator, current):
    if len(accumulator[1]) == 0:
        accumulator[1].append((current[0], current[1], current[2]))
    else:
        if accumulator[1][len(accumulator[1]) - 1][1] == current[1]:
            accumulator[1].append((current[0], current[1], current[2]))
        else:
            accumulator[1].sort(key=lambda x: x[2], reverse=True)
            for i in range(0, len(accumulator[1])):
                accumulator[0].append((accumulator[1][i][0], accumulator[1][i][1], accumulator[1][i][2], i + 1))
            accumulator[1] = []
    return accumulator


def groupCells(accumulator, current):
    if not (current[0] in accumulator):
        accumulator[current[0]] = [current[2]]
    else:
        accumulator[current[0]].append(current[2])
    return accumulator


if __name__ == '__main__':
    # startingData

    df = pd.read_table('ekspresije.tsv', index_col=0)
    print(df.shape)
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
    geneAveragesDict = reduce(geneAverage, standardisedValuesArray, dict())
    geneCenterArray = list(map(lambda x: (x[1], x[0], x[2] - geneAveragesDict[x[1]][0]), standardisedValuesArray))
    geneCenterArray.sort(key=lambda x: x[0])
    geneVarianceArray = reduce(variance, geneCenterArray, [])

    # 2.2
    mostVariableGenesArray = sorted(geneVarianceArray, key=lambda x: abs(x[1]), reverse=True)
    mostVariableGenesDict = reduce(topG, geneVarianceArray, [dict(), 0])[0]

    # 2.3
    filteredGenes = sorted(reduce(geneFilter, standardisedValuesArray, []), key=lambda x: x[1])

    # 2.5
    rankedGenes = reduce(rankGenes, filteredGenes, [[], []])[0]

    # 3.1
    groupCellsArray = reduce(groupCells, rankedGenes, dict())
    print(groupCellsArray)
