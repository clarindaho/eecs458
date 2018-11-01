# import libraries
import sys
import os

# global variables

global states, observations
states = ('F', 'B')                     # "Fair" and "Biased"
observations = ('H', 'T')               # "Head" and "Tail"

global startProb, transProb, emitProb
startProb = (0.5, 0.5)
transProb = (0.9, 0.1, 0.1, 0.9)        # "Fair to Fair", "Fair to Biased", "Biased to Fair", "Biased to Biased"
emitProb = (0.5, 0.5, 0.75, 0.25)       # "Fair and Head", "Fair and Tail", "Biased and Head", "Biased and Tail"

global viterbiDecode, posterioriProb
viterbiDecode = list()
posterioriProb = 1

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
        sequence: the sequence of coin tosses
    """
    global sequence, position
    
    # first line of input file
    sequence = arguments[0]
    
    # second line of input file
    position = int(arguments[1])
    
    # close the input file
    inputFile.close()

"""
Determining optimal path using Viterbi Algorithm.
"""
def viterbi():
    # construct the grid
    probTable = []
    tracebackTable = []
    
    for j in range(len(states) + 1):
        pRow = []
        tRow = []
        for i in range(len(sequence) + 1):
            pRow.append(0)
            tRow.append('')
        
        probTable.append(pRow)
        tracebackTable.append(tRow)
        
    # fill in the grid with base case
    probTable[0][0] = 'States'
    tracebackTable[0][0] = 'States'
    
    for i in range(len(sequence)):
        probTable[0][1 + i] = sequence[i]
        tracebackTable[0][1 + i] = sequence[i]
        
    for j in range(len(states)):
        probTable[1 + j][0] = states[j]
        tracebackTable[1 + j][0] = states[j]
    
    # fill in the grid column by column
    for i in range(len(sequence)):
        for j in range(len(states)):
            probTable[1 + j][1 + i], tracebackTable[1 + j][1 + i] = detViterbiProb(probTable, 1 + j, 1 + i)

    # trace back to origin
    maxProb = 0
    maxState = list()
    for j in range(len(states)):
        if probTable[1 + j][len(sequence)] > maxProb:
            maxProb = probTable[1 + j][len(sequence)]
            maxState.clear()
            maxState.append(j)
        elif probTable[1 + j][len(sequence)] is maxProb:
            maxState.append(j)
    
    for state in maxState:
        traceback(tracebackTable, 1 + state, len(sequence), "")
    
    # update output results
    global posterioriProb 
    posterioriProb = probTable[2][position]

"""
Helper method for Viterbi algorithm that determines the probability of the current state given the observation
and determines the previous state that best optimizes the probability of the current state.
"""
def detViterbiProb(probTable, row, col):    
    # determine emission
    fairEmission = -1
    biasedEmission = -1
    if probTable[0][col] is 'H':            # head
        fairEmission = emitProb[0]
        biasedEmission = emitProb[2]
    elif probTable[0][col] is 'T':          # tail
        fairEmission = emitProb[1]
        biasedEmission = emitProb[3]
    
    if col is not 1:
        # determine previous probabilities
        fairCandidate = probTable[1][col - 1]                   # probability of fair at previous position
        biasedCandidate = probTable[2][col - 1]                 # probability of biased at previous position
        
        # determine all state transition probabilities
        # determine the transition with the maximum probability
        maxProb = 0
        nextState = list()
        
        if probTable[row][0] is 'F':
            fairToFair = fairCandidate * transProb[0] * fairEmission
            if fairToFair > maxProb:
                maxProb = fairToFair
                nextState.clear()
                nextState.append('F')
            elif fairToFair is maxProb:
                nextState.append('F') 
                
            biasedToFair = biasedCandidate * transProb[2] * fairEmission
            if biasedToFair > maxProb:
                maxProb = biasedToFair
                nextState.clear()
                nextState.append('F')
            elif biasedToFair is maxProb:
                nextState.append('F') 
                
        elif probTable[row][0] is 'B':
            fairToBiased = fairCandidate * transProb[1] * biasedEmission
            if fairToBiased > maxProb:
                maxProb = fairToBiased
                nextState.clear()
                nextState.append('B')
            elif fairToBiased is maxProb:
                nextState.append('B') 
            
            biasedToBiased = biasedCandidate * transProb[3] * biasedEmission
            if biasedToBiased > maxProb:
                maxProb = biasedToBiased
                nextState.clear()
                nextState.append('B')
            elif biasedToBiased is maxProb:
                nextState.append('B') 
            
        return maxProb, nextState
    
    elif col is 1:
        if probTable[row][0] is 'F':
            return startProb[row - 1] * fairEmission, list('S')
        elif probTable[row][0] is 'B':
            return startProb[row - 1] * biasedEmission, list('S')   

"""
Helper method that traces arrows back to the origin.
"""
def traceback(tracebackTable, row, col, path):   
    if 'S' in tracebackTable[row][col]:
        global viterbiDecode
        viterbiDecode.append(path)
    
    if 'F' in tracebackTable[row][col]:
        fPath = 'F' + path
        traceback(tracebackTable, 1, col - 1, fPath)
        
    if 'B' in tracebackTable[row][col]:
        bPath = 'B' + path
        traceback(tracebackTable, 2, col - 1, bPath)

"""
Prints the following output:
    a. The Viterbi path
    b. The posteriori probability generated by the biased coin at specified position
"""
def output():       
    print("Viterbi Path:")
    for path in viterbiDecode:
        print(path)
        
    print()
    
    print("Posteriori Probability at Position " + str(position) + " Generated by Biased Coin: " + str(posterioriProb))

"""Main method."""
def main():
    parseOptions()
    viterbi()
    output()

if __name__ == "__main__": main()