/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/
/* $Revision: 1.1.10.4 $ */

#include "mex.h"
#include <stdlib.h>

//Global vars

short *cfgs;
static const int  DEFAULT_CFG_SIZE = 50;
int numCfgs;
int cfgArrSize;
int numRxns;
int numPatterns;
int totalNumConfigsCut;
short *unfeasiblePatterns;
short *dxns;
short *lastNonZero;



void recursiveCall(short *dirTemp, int level){
    //Base case 1: hit the bottom
    int i, j, k, foundPattern, numNonZero, numConfigsCut, numDirs;
    short patternDir;
    short thisDir;
    /* Report on the number of configs found and cut */
    if (numCfgs % 100 == 0) {
        printf("Number of configs: %d\n", numCfgs);
    }        
    
    level += 1;    
    //Base case 1: Pattern is in list of unfeasible patterns    
    for (i = 0; i < numPatterns; i++){
        //Check that we are capable of matching the pattern at this level
        numNonZero = lastNonZero[i];
        if (level < numNonZero+1) {
            continue;
        }
        foundPattern = 1;    
        for (j = 0; j < numNonZero; j++) {
            patternDir = unfeasiblePatterns[numPatterns * j + i];
            thisDir = dirTemp[j];
            if (patternDir != 2 && thisDir != patternDir){
                // Failed pattern matching. Move to the next pattern
                foundPattern = 0;
                break;
            }
        }
        if (foundPattern == 0){
            continue;
        }
        else {
            //We have a match for an unfeasible pattern. Don't continue this branch
            
            //Report number of configurations cut
            numConfigsCut = 1;
            for (j = level; j < numRxns; j++){
                numDirs = dxns[numRxns * 0 + j] + dxns[numRxns * 1 + j] + dxns[numRxns * 2 + j];
                numConfigsCut *= numDirs;
            }
            printf("Cut %d configurations", numConfigsCut);
            totalNumConfigsCut += numConfigsCut;
            return;
        }
    }  
    //Base case 2: This is a feasible configuration
    if (level > numRxns) {        
        //Check if the array of feasible patterns needs to be extended
        if (numCfgs == cfgArrSize){
            //We need to reallocate the array
            printf("reallocating\n");
            cfgArrSize *= 2;
            cfgs = realloc(cfgs, sizeof(short) * numRxns * cfgArrSize);
            if (cfgs == NULL) {
                // Assigning more memory failed
                fprintf(stdout, "Unable to allocate memory");
                free(cfgs);
                free(dirTemp);
                printf("Exiting");
                exit(0);
            }
        }
        for (i = 0; i < numRxns; i++){
            //Maybe we could increment the pointer instead. This could improve caching (probably not significantly)
            cfgs[numCfgs * numRxns + i] = dirTemp[i];
        }
        numCfgs = numCfgs + 1;
        return;
    }  
    //Make the next recursive call   
    for (i = 0; i < 3; i++){
        if (dxns[numRxns * i + (level-1)] == 1){
            //This direction is allowed. Assign it to the array
            dirTemp[level-1] = i-1; //Convert to -1, 0, 1 notation            
            recursiveCall(dirTemp, level);
        }
    }       
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    short *outputData;
    short *dirTemp;
    int i, j, numCfgsSaved;
    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("JuicyError:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("JuicyError:nlhs","One output required.");
    }
    
    /* make sure the first input argument is an array of shorts */
    if( !mxIsInt16(prhs[0]) || !(mxGetM(prhs[0]) > 1) || !(mxGetN(prhs[0]) > 1)) {
        mexErrMsgIdAndTxt("JuicyError:notScalar","Unfeasible pattern list must be a matrix.");
    }
    
    /* make sure the second input argument is an array of shorts */
    if( !mxIsInt16(prhs[1]) || !(mxGetM(prhs[1]) > 1) || !(mxGetN(prhs[1]) > 1)) {
        mexErrMsgIdAndTxt("JuicyError:notScalar","Directionalities must be a matrix.");
    }
    /* make sure the third input argument is an array of shorts */
    if( !mxIsInt16(prhs[2]) || !(mxGetM(prhs[2]) > 1) || !(mxGetN(prhs[2]) == 1)) {
        mexErrMsgIdAndTxt("JuicyError:notScalar","Last non-zero pattern element list must be a vector.");
    }
    /* Make sure that lastNonZero and unfeasiblePatterns have the same N */
    if (mxGetM(prhs[0]) != mxGetM(prhs[2])){
        mexErrMsgIdAndTxt("JuicyError:invalidArguments","lastNonZero and unfeasiblePatterns must be the same height.");
    }
    
    /* Set the counters */
    numCfgs = 0;
    cfgArrSize = 50;    
    totalNumConfigsCut = 0;
    /* Read in our inputs  */
    numRxns = mxGetN(prhs[0]);
    numPatterns = mxGetM(prhs[0]);
    unfeasiblePatterns = mxGetData(prhs[0]);
    dxns = (short *) mxGetData(prhs[1]);
    lastNonZero = (short *) mxGetData(prhs[2]);
    
    /* Run */
    cfgs = calloc(sizeof(short), numRxns * DEFAULT_CFG_SIZE); //Must be cleared by caller
    dirTemp = calloc(sizeof(short), numRxns);
    recursiveCall(dirTemp, 0);
    //Free dirtemp
    free(dirTemp);
    
    //Transfer the output
    
    //Maybe we could create and assign for matlab directly with mxCalloc
//     plhs[0] = mxCreateNumericMatrix(0, 0, mxINT16_CLASS, mxREAL);
    /* Point mxArray to dynamicData */
//     mxSetData(plhs[0], dynamicData);
//     mxSetM(plhs[0], ROWS);
//     mxSetN(plhs[0], COLUMNS);
    // (From example)
    
    /* Allocate the output array */
    plhs[0] = mxCreateNumericMatrix(numRxns,numCfgs, mxINT16_CLASS, mxREAL);
    outputData = (short *) mxGetData(plhs[0]);
    printf("Num configs cut with branches: %d", totalNumConfigsCut);
    
    
    numCfgsSaved = 0;
    for (j = 0; j < numCfgs; j++){        
        for (i = 0; i < numRxns; i++){
            numCfgsSaved += 1;
//             printf("%d,",cfgs[i * numCfgs + j]);
            outputData[i * numCfgs + j] = cfgs[i * numCfgs + j];
        }
//         printf("\n");
    }    
    free(cfgs);
}
