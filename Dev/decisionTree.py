'''
Created on Feb 25, 2019

@author: Tyler
'''
import math
import os
import DecisionTreeObj as dtobj
import sys

def calcInfoGain(infoArray, totPos, totNeg):
    negativeVals = 0
    positiveVals = 0
    infoGainTotal = 0.0
    for i in range(0, 2):
        positiveVals = infoArray[i][1]
        negativeVals = infoArray[i][0]
        infoGainTotal += ((positiveVals + negativeVals)/(totPos + totNeg)) * calcEntropy(positiveVals, negativeVals)
    return calcEntropy(totPos, totNeg) - infoGainTotal

def calcEntropy(positives, negatives):
    if ((positives == 0) | (negatives == 0)):
        return 0
    positiveRatio = float(positives/(positives+negatives))
    negativeRatio = 1 - positiveRatio
    entropyVal = -1 * positiveRatio * math.log(positiveRatio,2) - negativeRatio * math.log(negativeRatio,2)
    return entropyVal

# How to iterate through directory source:
# https://stackoverflow.com/questions/10377998/how-can-i-iterate-over-files-in-a-given-directory
# How to get a singular file:
# https://stackoverflow.com/questions/13223737/how-to-read-a-file-in-other-directory-in-python
def calculateDecisionTree(currentRoot):
    bestAttSoFar = ""
    bestIGSoFar = -2.0
    currentValues = currentRoot.getValues()
    rowLen = len(currentValues[0]) - 1
    posNegMatrix = [[0,0],[0,0]]
    for j in range(0, rowLen):
        for i in range(0, len(currentValues)):
            if(currentValues[i][j] == 0):
                if(currentValues[i][rowLen] == 0):
                    posNegMatrix[0][0] += 1
                else:
                    posNegMatrix[0][1] += 1
            else:
                if(currentValues[i][rowLen] == 0):
                    posNegMatrix[1][0] += 1
                else:
                    posNegMatrix[1][1] += 1
        thisValue = calcInfoGain(posNegMatrix, posNegMatrix[1][0] + posNegMatrix[1][1], posNegMatrix[0][0] + posNegMatrix[0][1])
        if thisValue > bestIGSoFar:
            bestIGSoFar = thisValue
            bestAttSoFar = currentRoot.getDictVal(j)
        posNegMatrix = [[0,0],[0,0]]
    if bestIGSoFar <= 0:
        counter = 0
        for i in range(0, len(currentValues)):
            if currentValues[i][len(currentValues[0]) - 1] == 1:
                counter+=1
            else:
                counter-=1
        if counter >= 0:
            currentRoot.paramName = "True"
        else:
            currentRoot.paramName = "False"
    elif (len(currentValues) == 2):
        if posNegMatrix[0][0] > posNegMatrix[0][1]:
            if posNegMatrix[1][0] > posNegMatrix[1][1]:
                currentRoot.paramName = "False"
            else:
                currentRoot.leftVal.paramName = "False"
                currentRoot.rightVal.paramName = "True"
        else:
            if posNegMatrix[1][1] > posNegMatrix[1][0]:
                currentRoot.paramName = "True"
            else:
                currentRoot.leftVal.paramName = "True"
                currentRoot.rightVal.paramName = "False"
    else:
        currentRoot.insertChild(bestAttSoFar, currentValues)
        currentRoot.insertChild(bestAttSoFar, currentValues)
        calculateDecisionTree(currentRoot.leftVal)
        calculateDecisionTree(currentRoot.rightVal)

def getFeatureList(aminoName):
    aminoProperties = {
        # [Hydrophobic, Polar, Small, Proline, Tiny, Aliphatic,
        # Aromatic, Positive, Negative, Charged, Placeholder for B/E Val]
        'A': [1,0,1,0,1,0,0,0,0,0,0],
        'C': [1,0,1,0,0,0,0,0,0,0,0],
        'D': [0,1,1,0,0,0,0,0,1,1,0],
        'E': [0,1,0,0,0,0,0,0,1,1,0],
        'F': [1,0,0,0,0,0,1,0,0,0,0],
        'G': [1,0,1,0,1,0,0,0,0,0,0],
        'H': [0,1,0,0,0,0,1,1,0,1,0],
        'I': [1,0,0,0,0,1,0,0,0,0,0],
        'K': [0,1,0,0,0,0,0,1,0,1,0],
        'L': [1,0,0,0,0,1,0,0,0,0,0],
        'M': [1,0,0,0,0,0,0,0,0,0,0],
        'N': [0,1,1,0,0,0,0,0,0,0,0],
        'P': [1,0,1,1,0,0,0,0,0,0,0],
        'Q': [0,1,0,0,0,0,0,0,0,0,0],
        'R': [0,1,0,0,0,0,0,1,0,1,0],
        'S': [0,1,1,0,1,0,0,0,0,0,0],
        'T': [1,1,1,0,0,0,0,0,0,0,0],
        'V': [1,0,1,0,0,1,0,0,0,0,0],
        'W': [1,0,0,0,0,0,1,0,0,0,0],
        'Y': [1,1,0,0,0,0,1,0,0,0,0]}
    return aminoProperties[aminoName]
            
def getBaseArrayDict(isReverse=False):
    if not(isReverse):
        arrayDict = {
            0:'Hydrophobic', 
            1:'Polar', 
            2:'Small', 
            3:'Proline',
            4:'Tiny', 
            5:'Aliphatic', 
            6:'Aromatic', 
            7:'Positive', 
            8:'Negative', 
            9:'Charged'}
    else:
        arrayDict = {
            'Hydrophobic':0, 
            'Polar':1, 
            'Small':2, 
            'Proline':3,
            'Tiny':4, 
            'Aliphatic':5, 
            'Aromatic':6, 
            'Positive':7, 
            'Negative':8, 
            'Charged':9}
    return arrayDict     

def getBEValues(aminoString, buriedExposedString):
    BEValueArray = [[0] for _ in range(0, len(aminoString))]
    for i in range(0, len(aminoString)):
        if not(aminoString[i] == '\n'):
            BEValueArray[i] = getFeatureList(aminoString[i])
            if(buriedExposedString[i] == 'B'):
                BEValueArray[i][10] = 0
            else:
                BEValueArray[i][10] = 1
    return BEValueArray

def getTestValues(dtRoot, testAminoSequence):
    BEExpectedVals = ""
    arrayDict = getBaseArrayDict(True)
    for i in range(0, len(testAminoSequence)):
        searching = True
        currentNode = dtRoot
        while(searching):
            currentProperties = getFeatureList(testAminoSequence[i])
            index = arrayDict[currentNode.getName()]
            if currentProperties[index] == 1:
                currentNode = currentNode.rightVal
            else:
                currentNode = currentNode.leftVal
            if currentNode.getName() == 'False':
                BEExpectedVals += 'B'
                searching = False
            elif currentNode.getName() == 'True':
                BEExpectedVals += 'E'
                searching = False
    return BEExpectedVals
        
        
# How to get a singular file:
# https://stackoverflow.com/questions/13223737/how-to-read-a-file-in-other-directory-in-python

if __name__ == "__main__":
    trainAminoString = ""
    trainBuriedExposedString = ""
    testAminoString = ""
    trueBuriedExposedString = ""
    professorAminoString = []
    professorSAFiles = []
    fileCounter = 0
    for file in os.listdir(sys.argv[1]):
        filename = os.path.basename(file)
        fastaFile = open(sys.argv[1] + filename, 'r')
        fastaLine = ""
        saLine = ""
        for i, line in enumerate(fastaFile):
            fastaLine = line
        saFilename = filename[:4] + '.sa'
        if saFilename != '.DS_.sa':
            if(fileCounter < 150):
                saFile = open(sys.argv[2]+saFilename, 'r')
                if saFile == '.DS_':
                    continue
                for i, line in enumerate(saFile):
                    saLine = line
                    if saLine[0] == '>':
                        saLine = saLine[5:len(saLine)]
            else:
                professorSAFiles.append(saFilename)
            fastaLine = fastaLine[0:len(fastaLine) - 1]
            if fastaLine[0] == '>':
                fastaLine = fastaLine[5:len(fastaLine)]
            if(fileCounter < 100):
                trainAminoString += fastaLine
                trainBuriedExposedString += saLine
                fileCounter+=1
            elif(fileCounter < 150):
                testAminoString += fastaLine
                trueBuriedExposedString += saLine
                fileCounter+=1
            else:
                professorAminoString.append(fastaLine)
    trainingMatrix = getBEValues(trainAminoString, trainBuriedExposedString)
    arrayDict = getBaseArrayDict()
    rootObject = dtobj.DecisionTreeObj("root", trainingMatrix, arrayDict)
    calculateDecisionTree(rootObject)
    expectedBuriedExposedString = getTestValues(rootObject, testAminoString)
    truePos = 0.0
    trueNeg = 0.0
    falsePos = 0.0
    falseNeg = 0.0
    for i in range(0, len(testAminoString)):
        if expectedBuriedExposedString[i] == 'B':
            if trueBuriedExposedString[i] == 'B':
                trueNeg+=1.0
            else:
                falseNeg+=1.0
        else:
            if trueBuriedExposedString[i] == 'B':
                falsePos+=1.0
            else:
                truePos+=1.0
    accuracy = (truePos + trueNeg)/len(expectedBuriedExposedString)
    precision = truePos/(truePos + falsePos)
    recall = truePos/(truePos+falseNeg)
    f1score = 2*(precision*recall)/(precision+recall)
    mccscore = (truePos * trueNeg - falsePos * falseNeg)/math.sqrt((truePos + falsePos) * (truePos + falseNeg) * (trueNeg + falsePos) * (trueNeg + falseNeg))
    print("True Pos: %d \n True Neg: %d \n False Pos: %d \n False Neg: %d \n Accuracy: %f \n Precision: %f \n Recall: %f \n F1: %f \n MCC: %f" % (truePos, trueNeg, falsePos, falseNeg, accuracy, precision, recall, f1score, mccscore))
    for i in range(0, len(professorAminoString)):
        testResult = getTestValues(rootObject, professorAminoString[i])
        print("Provided String: %s\nDT Output      : %s\n" %(professorAminoString[i], testResult))
        if testResult.count('B') < 3*testResult.count('E'):
            print("The provided amino acid string is classified as a primarily exposed Protein")
        else:
            print("The provided amino acid string is classified as a primarily buried Protein")
        saFile = open(sys.argv[2]+saFilename, 'w+')
        saFile.write('>%s%s' %(saFilename, getTestValues(rootObject, professorAminoString[i])))
        saFile.close()