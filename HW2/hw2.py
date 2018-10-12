# import libraries
import sys
import os

# global variables

global optimalScore, numSolutions, solutions
numSolutions = 0
solutions = list()

# helper methods

def parseOptions():
    """System Argument #1: path to the input file"""
    path = str(sys.argv[1])
    
    # read from input file
    inputFile = open(os.getcwd() + "/" + path, "r")
    if inputFile.mode == 'r':
        arguments = inputFile.read().splitlines()
        
    # determine input arguments
    """
    Arguments:
        alignment:     0 if global alignment, 1 if local alignment
        matchScore:    score for match
        mismatchScore: score for mismatch
        indelScore:    score for insertion and deletion
        firstSeq:      first sequence to be aligned
        secondSeq:     second sequence to be aligned
    """
    global alignment, matchScore, mismatchScore, indelScore, firstSeq, secondSeq
    
    # first line of input file
    if arguments[0] == 'g':
        alignment = 0
    elif arguments[0] == 'l':
        alignment = 1
    else:
        print("ERROR: Alignment must be either 'g' or 'l'")
        
    # second line of input file
    scores = arguments[1].split()
    matchScore = int(scores[0])
    mismatchScore = int(scores[1])
    indelScore = int(scores[2])
    
    # third line of input file
    firstSeq = arguments[2]
    
    # fourth line of input file
    secondSeq = arguments[3]
    
    # close the input file
    inputFile.close()

"""
Global alignment via the Needleman-Wunsch Algorithm.
"""
def globalAlign():
    # construct the grid
    scoreTable = []
    tracebackTable = []
    
    for j in range(len(secondSeq) + 2):
        sRow = []
        tRow = []
        for i in range(len(firstSeq) + 2):
            sRow.append(0)
            tRow.append('')
        
        scoreTable.append(sRow)
        tracebackTable.append(tRow)
    
    # fill in the grid with base case
    scoreTable[1][1] = 0                                               # base case
    tracebackTable[1][1] = list('S')                                   # stop
    
    for i in range(len(firstSeq)):
        scoreTable[0][2 + i] = firstSeq[i]
        scoreTable[1][2 + i] = scoreTable[1][1 + i] + indelScore       # indel pairing 
        tracebackTable[0][2 + i] = firstSeq[i]
        tracebackTable[1][2 + i] = list('L')                           # left
        
    for j in range(len(secondSeq)):
        scoreTable[2 + j][0] = secondSeq[j]
        scoreTable[2 + j][1] = scoreTable[1 + j][1] + indelScore       # indel pairing
        tracebackTable[2 + j][0] = secondSeq[j]
        tracebackTable[2 + j][1] = list('U')                           # up
    
    # fill in the grid row by row
    for j in range(len(secondSeq)):
        for i in range(len(firstSeq)):
            scoreTable[2 + j][2 + i], tracebackTable[2 + j][2 + i] = gDetMaxScore(scoreTable, 2 + j, 2 + i)
    
    # trace back to origin
    traceback(tracebackTable, int(len(secondSeq) + 1), int(len(firstSeq) + 1), "", "")
    
    # update output results
    global optimalScore
    optimalScore = scoreTable[len(secondSeq) + 1][len(firstSeq) + 1]
    
"""
Helper method for global alignment that determines the direction that best optimizes the score.
"""
def gDetMaxScore(scoreTable, row, col):
    # determine previous score
    diagonalCandidate = scoreTable[row - 1][col - 1]        # diagonal (top-left) neighbor
    topCandidate = scoreTable[row - 1][col]                 # top neighbor
    leftCandidate = scoreTable[row][col - 1]                # left neighbor
    
    # consider action to make
    ## diagonal path represents a match/mismatch
    if scoreTable[0][col] == scoreTable[row][0]:
        diagonalCandidate += matchScore
    else:
        diagonalCandidate += mismatchScore
    
    ## path from top or left represents an indel pairing
    topCandidate += indelScore
    leftCandidate += indelScore
    
    # return the maximizing score
    maxScore = diagonalCandidate
    maxPath = list('D')
    
    if topCandidate > maxScore:
        maxScore = topCandidate
        maxPath.clear()
        maxPath.append('U')
    elif topCandidate == maxScore:
        maxPath.append('U')
    
    if leftCandidate > maxScore:
        maxScore = leftCandidate
        maxPath.clear()
        maxPath.append('L')
    elif leftCandidate == maxScore:
        maxPath.append('L')
    
    return maxScore, maxPath

"""
Local alignment via the Smith-Waterman Algorithm.
"""
def localAlign():
    # construct the grid
    scoreTable = []
    tracebackTable = []
    
    for j in range(len(secondSeq) + 2):
        sRow = []
        tRow = []
        for i in range(len(firstSeq) + 2):
            sRow.append(0)
            tRow.append('')
        
        scoreTable.append(sRow)
        tracebackTable.append(tRow)
    
    # fill in the grid with base case
    scoreTable[1][1] = 0                                  # base case
    tracebackTable[1][1] = list('S')                      # stop
    
    for i in range(len(firstSeq)):
        scoreTable[0][2 + i] = firstSeq[i]
        scoreTable[1][2 + i] = 0
        tracebackTable[0][2 + i] = firstSeq[i]
        tracebackTable[1][2 + i] = list('S')              # stop
        
    for j in range(len(secondSeq)):
        scoreTable[2 + j][0] = secondSeq[j]
        scoreTable[2 + j][1] = 0
        tracebackTable[2 + j][0] = secondSeq[j]
        tracebackTable[2 + j][1] = list('S')              # stop
    
    # fill in the grid row by row
    # determine the global maximum in scores
    maxScore = 0
    maxScoreRow = list()
    maxScoreCol = list()
    for j in range(len(secondSeq)):
        for i in range(len(firstSeq)):
            scoreTable[2 + j][2 + i], tracebackTable[2 + j][2 + i] = lDetMaxScore(scoreTable, 2 + j, 2 + i)
            
            if scoreTable[2 + j][2 + i] > maxScore:
                maxScore = scoreTable[2 + j][2 + i]
                maxScoreRow.clear()
                maxScoreRow.append(2 + j)
                maxScoreCol.clear()
                maxScoreCol.append(2 + i)
            elif scoreTable[2 + j][2 + i] == maxScore:
                maxScoreRow.append(2 + j)
                maxScoreCol.append(2 + i)  

    # trace back to origin
    for i in range(len(maxScoreRow)):
        traceback(tracebackTable, maxScoreRow[i], maxScoreCol[i], "", "")
    
    # update output results
    global optimalScore
    optimalScore = maxScore

"""
Helper method for local alignment that determines the direction that best optimizes the score.
(Assumes linear gap penalty where linear gap penalty = indel pairing score)
"""
def lDetMaxScore(scoreTable, row, col):
    # determine previous score
    diagonalCandidate = scoreTable[row - 1][col - 1]        # diagonal (top-left) neighbor
    topCandidate = scoreTable[row - 1][col]                 # top neighbor
    leftCandidate = scoreTable[row][col - 1]                # left neighbor
    
    # consider action to make
    ## diagonal path represents a match/mismatch
    if scoreTable[0][col] == scoreTable[row][0]:
        diagonalCandidate += matchScore
    else:
        diagonalCandidate += mismatchScore
    
    ## path from top or left represents an indel pairing
    topCandidate += indelScore
    leftCandidate += indelScore
    
    # return the maximizing score
    maxScore = diagonalCandidate
    maxPath = list('D')
    
    if topCandidate > maxScore:
        maxScore = topCandidate
        maxPath.clear()
        maxPath.append('U')
    elif topCandidate == maxScore:
        maxPath.append('U')
    
    if leftCandidate > maxScore:
        maxScore = leftCandidate
        maxPath.clear()
        maxPath.append('L')
    elif leftCandidate == maxScore:
        maxPath.append('L')
        
    if 0 > maxScore:
        maxScore = 0
        maxPath.clear()
        maxPath.append('S')
    elif 0 == maxScore:
        maxPath.append('S')
    
    return maxScore, maxPath

"""
Helper method that traces arrows back to the origin.
"""
def traceback(tracebackTable, row, col, firstPath, secondPath):   
    if 'S' in tracebackTable[row][col]:
        global numSolutions
        numSolutions += 1
        
        global solutions
        solutions.append(firstPath + "\n" + secondPath + "\n")
    
    if 'D' in tracebackTable[row][col]:
        dFirstPath = tracebackTable[0][col] + firstPath
        dSecondPath = tracebackTable[row][0] + secondPath
        traceback(tracebackTable, row - 1, col - 1, dFirstPath, dSecondPath)
        
    if 'U' in tracebackTable[row][col]:
        uFirstPath = '-' + firstPath
        uSecondPath = tracebackTable[row][0] + secondPath
        traceback(tracebackTable, row - 1, col, uFirstPath, uSecondPath)
    
    if 'L' in tracebackTable[row][col]: 
        lFirstPath = tracebackTable[0][col] + firstPath
        lSecondPath = '-' + secondPath
        traceback(tracebackTable, row, col - 1, lFirstPath, lSecondPath)

"""
Prints the following output:
    a. The optimal score
    b. The number of solutions that achieve the optimal score
    c. The actual alignments with the optimal score
"""
def output():
    print("Optimal Score: " + str(optimalScore))
    print("Number of Solutions that Achieve the Optimal Score: " + str(numSolutions))
    print("Actual Alignments with the Optimal Score:\n")
    for solution in solutions:
        print(solution)

"""Main method."""
def main():
    parseOptions()
    
    if alignment == 0:
        globalAlign()
    else:
        localAlign()
    
    output()

if __name__ == "__main__": main()