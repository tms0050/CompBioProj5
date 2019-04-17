'''
Created on Feb 28, 2019

@author: Tyler
'''

class DecisionTreeObj(object):
    '''
    classdocs
    '''


    def __init__(self, objectName, valueMatrix=None, arrayDict=None):
        self.leftVal = None
        self.rightVal = None
        self.paramName = objectName
        self.values = valueMatrix
        self.valueDict = arrayDict
        '''
        Constructor
        '''
# http://code.activestate.com/recipes/192401-quickly-remove-or-order-columns-in-a-list-of-lists/
# Source used for deleting a single column out of an array
  
    def insertChild(self, objectName, valueMatrix=None):
        arrayDict = self.valueDict.copy()
        self.paramName = objectName
        newValueMatrix = []
        finalValueMatrix = valueMatrix.copy()
        if not(valueMatrix is None):  
            if self.leftVal is None:
                valToCheck = 0
            else:
                valToCheck = 1
            colNum = 0
            removed = False
            keyVal = 0
            for key in self.valueDict:
                keyVal = key - 1
                if removed:
                    arrayDict[keyVal] = arrayDict[key]
                if arrayDict[key] == objectName:
                    removed = True
                    colNum = key
                    if colNum != 0:
                        newValueMatrix = [[col[0:key] + col[key+1:len(valueMatrix[0])]] for col in valueMatrix]
                    else:
                        newValueMatrix = [[0] for _ in range(0,len(valueMatrix))]
                        for i in range(0, len(valueMatrix)):
                            newLine = valueMatrix[i][1:len(valueMatrix[0])]
                            newValueMatrix[i] = newLine
            del arrayDict[keyVal + 1]
            locator = 0
            for i in range(0, len(valueMatrix)):
                if valueMatrix[i][colNum] == valToCheck:
                    if colNum != 0:
                        finalValueMatrix[locator] = newValueMatrix[i][0]
                    else:
                        finalValueMatrix[locator] = newValueMatrix[i]
                    locator += 1
            finalValueMatrix = finalValueMatrix[0:locator]
        # 0 values are always in leftVal, 1 values are always in rightVal            
        if self.leftVal is None:
            self.leftVal = DecisionTreeObj("", finalValueMatrix, arrayDict)
        else:
            self.rightVal = DecisionTreeObj("", finalValueMatrix, arrayDict)
            
    def getName(self):
        return self.paramName
    
    def getValues(self):
        return self.values
    
    def getDictVal(self, key):
        return self.valueDict[key]