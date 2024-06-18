#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
//Global vars

unsigned char *cfgs;
int fileNum;
// static const int  DEFAULT_CFG_SIZE = 75000000;
static const int  DEFAULT_CFG_SIZE = 10000000;
// static const int  DEFAULT_CFG_SIZE = 500000;
int numCfgs;
unsigned long long totalNumCfgs;
int cfgArrSize;
int numRxns;
int numPatterns;
int totalNumConfigsCut;
short *unfeasiblePatterns;
short *lastNonZero;
int dispFlag;
int clusterNum;
int numCharPerRow;



void saveFile() {
    //Save the current configurations to a file
    char filename[17];
    FILE *fileID;
    clock_t startTime;
    int i;
    //Make the filename
    fileNum += 1;
    sprintf(filename, "TFM_%03d_%04d.bin", clusterNum, fileNum);
    fileID = fopen(filename, "wb");
    startTime = clock();
//     for (i = 0; i < numCfgs; i++){
//         printf("%u:  %d\n", i, cfgs[i]);        
//     }    
    fwrite(cfgs, sizeof(char), numCharPerRow * numCfgs, fileID); 
    printf("Total configs: %llu\n", totalNumCfgs);
    printf("time taken: %d\n", clock() - startTime);
    fclose(fileID);
}


void recursiveCall(short *dirTemp, int level){
    //Base case 1: hit the bottom    
    //The stack seems to be increasing quite steadily - maybe if we reduce the number of variables here we can reduce this effect
    int i, j, foundPattern, bitBufferPos, rowPos;
    float newArrSize;
    unsigned char bitBuffer;
    bitBufferPos = 0;
    /* Report on the number of configs found and cut */
    if (totalNumCfgs % 1000000 == 0 && dispFlag == 0) {
        mexPrintf("Number of configs: %llu\n", totalNumCfgs);
        mexEvalString("drawnow");
        dispFlag = 1;
    } else if (numCfgs % 100000 != 0) {
        dispFlag = 0;
    }    
    level += 1;    
    //Base case 1: Pattern is in list of unfeasible patterns    
    for (i = 0; i < numPatterns; i++){
        //Check that we are capable of matching the pattern at this level
        if (level < lastNonZero[i]+1) {
            continue;
        }
        foundPattern = 1;    
        for (j = 0; j < lastNonZero[i]; j++) {
            //Replacing these for less memory consumption
            if (unfeasiblePatterns[numPatterns * j + i] != 0 && dirTemp[j] != unfeasiblePatterns[numPatterns * j + i]){
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
            return;
        }
    }  
    //Base case 2: This is a feasible configuration
    if (level > numRxns) {        
        //Check if the array of feasible patterns needs to be extended
        if (numCfgs == cfgArrSize){
            // Save the configurations to file and start saving from the beginning
           saveFile();
           //Reset the number of configurations to 0 
           numCfgs = 0;        
        }
        //Save the new configuration to the configurations array
        bitBuffer = 0;
        rowPos = 0;
        bitBufferPos = 0;
        for (i = 0; i < numRxns; i++){
            //Fill the buffer with our number            
            if (dirTemp[i] == 1){
                //Add a 1 to the bit buffer
                bitBuffer = bitBuffer + 1;                
            }
           //If the buffer is full, write it to the current position in the list of configurations
           if (bitBufferPos == 7){                         
              cfgs[numCfgs * numCharPerRow + rowPos] = bitBuffer;
              bitBufferPos = -1; //this gets corrected by subsequent increment
              bitBuffer = 0;
              rowPos += 1;
           }
            //Shift the bits along
           bitBuffer <<= 1;                
           bitBufferPos += 1;
        }

        //Write the final buffer position
        bitBuffer >>= 1;
        cfgs[numCfgs * numCharPerRow + rowPos] = bitBuffer;       
        numCfgs += 1;
        totalNumCfgs += 1;
        return;
    }     
    //Make the next recursive calls
    dirTemp[level-1] = -1; //Reverse            
    recursiveCall(dirTemp, level);
    dirTemp[level-1] = 1; //Forward           
    recursiveCall(dirTemp, level); 
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {

    short *outputData;
    short *dirTemp;
    long *outputNumRxns;
    int i, j, numCfgsSaved;
    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("JuicyError:nrhs","Two inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("JuicyError:nlhs","One output required.");
    }
    
    /* make sure the first input argument is an array of shorts */
    if( !mxIsInt16(prhs[0]) || !(mxGetM(prhs[0]) > 1) || !(mxGetN(prhs[0]) > 1)) {
        mexErrMsgIdAndTxt("JuicyError:notScalar","Unfeasible pattern list must be a matrix.");
    }
    
    /* make sure the third input argument is an array of shorts */
    if( !mxIsInt16(prhs[1]) || !(mxGetM(prhs[1]) > 1) || !(mxGetN(prhs[1]) == 1)) {
        mexErrMsgIdAndTxt("JuicyError:notScalar","Last non-zero pattern element list must be a vector.");
    }
    /* Set the counters */
    numCfgs = 0;
    totalNumCfgs = 0;
    fileNum = 0;
    cfgArrSize = DEFAULT_CFG_SIZE;   
    totalNumConfigsCut = 0;
    dispFlag = 1;
    /* Read in our inputs  */
    numRxns = mxGetN(prhs[0]);
    numPatterns = mxGetM(prhs[0]);
    unfeasiblePatterns = mxGetData(prhs[0]);
    lastNonZero = (short *) mxGetData(prhs[1]);
    clusterNum = (int) mxGetScalar(prhs[2]);   
    //Calculate the number of chars (8 bits) required per line. 
    numCharPerRow = ceil(numRxns / (double) 8);    
    cfgs = calloc(sizeof(short), numCharPerRow * DEFAULT_CFG_SIZE); //Must be cleared by caller
    dirTemp = calloc(sizeof(short), numRxns);
//     saveFile();
    printf("Starting recursion\n");
    /* Run */
    recursiveCall(dirTemp, 0);
    printf("Total num configs: %llu\n", totalNumCfgs);
    printf("Done\n");
    //Save the final file
    saveFile();    
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
    printf("Num configs cut with branches: %d\n", totalNumConfigsCut);
    plhs[1] = mxCreateNumericMatrix(1,1, mxUINT64_CLASS, mxREAL);
    outputNumRxns = (unsigned long long *) mxGetData(plhs[1]);
    *outputNumRxns = totalNumCfgs;
    
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
