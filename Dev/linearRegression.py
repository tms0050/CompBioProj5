'''
Created on Apr 17, 2019

@author: Tyler
'''
import os
import decisionTree as dt
import random
import math
import time

def bayesExtract():
    dirPath = os.path.abspath(os.path.dirname(__file__))
    bayesPath = os.path.join(dirPath, 'HEC_guess.txt')
    bayes = open(bayesPath)
    lines = bayes.readlines()
    returnDict = dict()
    fastaPath = os.path.join(dirPath, '..', '5970_6970_SP_19_PROJECT_5', 'fasta')
    fileList = os.listdir(fastaPath)
    colPosition = 0
    rowPosition = 0
    for fileName in fileList:
        strValue = ""
        fastaFile = open(os.path.join(fastaPath, fileName))
        fastaLines = fastaFile.readlines()
        lengthToSearch = len(fastaLines[len(fastaLines) - 1])
        for _ in range(0, lengthToSearch):
            if colPosition == len(lines[rowPosition]):
                colPosition = 0
                rowPosition += 1
            strValue += lines[rowPosition][colPosition]
            colPosition += 1
        returnDict[fileName[:4]] = strValue
    return returnDict

def pssmExtract(filename):   
    dirPath = os.path.abspath(os.path.dirname(__file__))
    pssmPath = os.path.join(dirPath, '..', '5970_6970_SP_19_PROJECT_5', 'pssm', filename + '.pssm')
    fp = open(pssmPath)
    lines = fp.readlines()
    pssmSum = [0.0]*20  
    counter = 0     
    for y in lines[3:]:            
        if y == '\n':
            break            
        trimmedline = y[9:]
        trimmedline = trimmedline[trimmedline.index('   '):]
        trimmedline = trimmedline.split(' ')
        trimmedline = trimmedline[:len(trimmedline)-2]
        counter = -1
        for i in range(0, len(trimmedline)):
            if len(trimmedline[i]) == 0:
                continue
            if trimmedline[i][0] == ' ':
                trimmedline[i] = trimmedline[i][1:]
            if trimmedline[i][len(trimmedline[i]) - 1] == ' ':
                trimmedline[i] = trimmedline[i][:(len(trimmedline[i]) - 1)]
            counter += 1
            pssmSum[counter] += (float(trimmedline[i])/100.0)
    for i in range(0, 20):
        pssmSum[i] = round(pssmSum[i]/150,5)
    return pssmSum

def combineFeatures():
    featureList = [0.0] * 25
    featureDict = dict()
    fastapath = os.path.abspath(os.path.join(os.path.abspath(__file__), '..', '..', '5970_6970_SP_19_PROJECT_5', 'fasta'))
    fileList = os.listdir(fastapath)
    dtArray = dt.implementDecisionTree()
    for fileName in fileList:
        if not(os.path.isfile(os.path.join(fastapath, fileName))):
            continue
        fileName = fileName[:4]
        featureVals = pssmExtract(fileName)[:]
        for i in range(0, 20):
            featureList[i] = featureVals[i]
        featureList[20] = round(dtArray[fileName][0],5)
        featureList[21] = round(dtArray[fileName][1],5)
        featureDict[fileName] = featureList[:]
    return featureDict
    
def gradientDescent(featureDict, currentWeights, currentStep, testTMVals, timeSinceLastStep):
    dirPath = os.path.abspath(os.path.dirname(__file__))
    weightsToTest = currentWeights[:]
    pssmPath = os.path.join(dirPath, '..', '5970_6970_SP_19_PROJECT_5', 'pssm')
    fileList = os.listdir(pssmPath)
    for _ in range(0, 200):
        val1 = int(random.random() * 100)
        val2 = int(random.random() * 150)
        while(val1 == val2):
            val2 = int(random.random() * 110)
        newWeights = [0.0] * len(currentWeights)
        for i in range(0, 25):
            newWeights[i] = featureDict[fileList[val1][:4]][i]
        for i in range(0, len(featureDict[fileList[val2][:4]])):            
            newWeights[i + 25] = featureDict[fileList[val2][:4]][i]
        newWeights[50] = getTMValue(fileList[val1], fileList[val2])
        for j in range(0, len(currentWeights) - 2):
            try:
                slope = (weightsToTest[len(currentWeights)-1] - newWeights[len(newWeights)-1])/(weightsToTest[j] - newWeights[j])
            except:
                slope = 1
            weightsToTest[j] = weightsToTest[j] + currentStep*(slope - weightsToTest[j])
    newVal = getAvgSqError(featureDict, weightsToTest, testTMVals)
    oldVal = getAvgSqError(featureDict, currentWeights, testTMVals)
    newTime = timeSinceLastStep + 1
    if (oldVal < newVal) & (timeSinceLastStep > 200):
        return currentWeights, (currentStep/10), 0
    elif newVal < oldVal:
        return weightsToTest, currentStep, 0
    else:
        return currentWeights, currentStep, newTime

def getTMValue(protein1, protein2):
    tmalign = os.path.abspath(os.path.join(os.path.abspath(__file__), '..', '..', '5970_6970_SP_19_PROJECT_5', 'tmalign', protein1[:4] + '_' + protein2[:4] + '_tmalign'))
    tmalignfile = open(tmalign)
    lines = tmalignfile.readlines()
    correctRow = lines[17]
    tmvalue = float(correctRow.split(' ')[1])
    return tmvalue

def getAvgSqError(featureDict, weightList, testTMVals):
    counter = 0
    totSqError = 0.0
    for key in testTMVals:
        correctAnswer = testTMVals[key]
        first25 = featureDict[key[:4]]
        second25 = featureDict[key[4:]]
        predTMVal = 0.0
        for i in range(0, 22):
            predTMVal += first25[i] * weightList[i]
        for j in range(0, 22):
            predTMVal += second25[j] * weightList[j + 25]
        totSqError += math.pow(correctAnswer - predTMVal, 2)
        counter += 1
    return totSqError/counter
    
def getTMAlignVals(test=True):
    retDict = dict()
    dirPath = os.path.abspath(os.path.dirname(__file__))
    pssmPath = os.path.join(dirPath, '..', '5970_6970_SP_19_PROJECT_5', 'pssm')
    fileList = os.listdir(pssmPath)
    if test:
        for i in range(100, 150):
            for j in range(0, 150):
                if i==j:
                    continue
                alignVal = getTMValue(fileList[i][:4], fileList[j][:4])
                keyName = fileList[i][:4] + '' + fileList[j][:4]
                retDict[keyName] = alignVal
    else:
        for i in range(0, 100):
            for j in range(0, 150):
                if i==j:
                    continue
                alignVal = getTMValue(fileList[i][:4], fileList[j][:4])
                keyName = fileList[i][:4] + '' + fileList[j][:4]
                retDict[keyName] = alignVal
    return retDict

def implementLinearRegression():
    myfeatures = combineFeatures()
    currentWeights = [0.0] * 51
    currentStep = 1.0
    trainTMVals = getTMAlignVals(False)
    testTMVals = getTMAlignVals()
    timeSinceLastStep = 0
    while currentStep > 0.00000001:
        currentWeights, currentStep, timeSinceLastStep = gradientDescent(myfeatures, currentWeights, currentStep, trainTMVals, timeSinceLastStep)
        if(timeSinceLastStep == 0):
            print 'Current Weights: ', currentWeights
            print 'Current Step Size: ', currentStep
            print 'Current Avg Square Error on Train Data: ', getAvgSqError(myfeatures, currentWeights, trainTMVals)
            print 'Current Avg Square Error on Test Data: ', getAvgSqError(myfeatures, currentWeights, testTMVals), '\n'
    
if __name__ == '__main__':
    start = time.time()
    implementLinearRegression()
    end = time.time()
    print('Took %d seconds to run' % (end - start))