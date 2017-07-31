#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "string.h"
#include "time.h"

/**/
typedef struct iter iter;
struct iter {
    int index;
    iter *next;
};

// Constructor
iter *iter_Init(iter *iterHead) {
	iterHead = NULL;
	return iterHead;
}

// Destructor
iter *iter_Finalize(iter *iterHead) {
    iter* iterTemp;
    while(iterHead != NULL) {
	iterTemp = iterHead;
        iterHead = iterHead->next;
        iterTemp->next = NULL;
        free(iterTemp);
    }
}

// Add a iter to the list of iters
iter *iter_Add(iter *iterHead, int index) {
	if(iterHead == NULL) {
		iterHead = malloc(sizeof(iter));
		iterHead->index = index;
		iterHead->next = NULL;
	}
	else {
		iter *iterTemp = iterHead;
		while(iterTemp->next != NULL) 
			iterTemp = iterTemp->next;
		iterTemp->next = malloc(sizeof(iter));
        iterTemp->next->index = index;
        iterTemp->next->next = NULL;
	} 
	return iterHead;
}
////

/**/
typedef struct {
    int numSend;
    int numRecv;
    int sumSrc;
    int sumDest;
    int sumRankSend;
    int sumRankRecv;
    int xorSend;
    int xorRecv;
} pattern;

//Constructor
void pattern_Init(pattern *Pattern, int size) {
    int i;
    for(i = 0; i < size; i++) {
    	Pattern[i].numSend = 0;
    	Pattern[i].numRecv = 0;
   		Pattern[i].sumSrc = 0;
    	Pattern[i].sumDest = 0;
    	Pattern[i].sumRankSend = 0;
    	Pattern[i].sumRankRecv = 0;
    	Pattern[i].xorSend = 0;
    	Pattern[i].xorRecv = 0;
    }
}

//Update patterns
void pattern_Update(pattern *Pattern, int index, int procIndex, int numSend, int numRecv,
                            int sumSrc, int sumDest, int xorSend, int xorRecv) {
    Pattern[index].numSend += numSend;
    Pattern[index].numRecv += numRecv;
    Pattern[index].sumSrc += sumSrc;
    Pattern[index].sumDest += sumDest;
    Pattern[index].sumRankSend += numSend * procIndex;
    Pattern[index].sumRankRecv += numRecv * procIndex;
    Pattern[index].xorSend ^= xorSend;
    Pattern[index].xorRecv ^= xorRecv;
}
////

/**/
typedef struct loop loop;
struct loop {
	int numIters;
	pattern* patternArray;
	loop* next;
};

//Constructor
loop *loop_Init(loop *loopHead) {
	loopHead = NULL;
	return loopHead;
}

//Destructor
loop *loop_Finalize(loop *loopHead) {
	loop *loopTemp;
	while(loopHead != NULL) {
		loopTemp = loopHead;
		loopHead = loopHead->next;
		loopTemp->next = NULL;
		free(loopTemp->patternArray);
		free(loopTemp);
	}
	return loopHead;
}

//Add a loop to the list of loops
loop *loop_Add(loop *loopHead, int numIters) {
	if(loopHead == NULL) {
		loopHead = malloc(sizeof(loop));
		loopHead->numIters = numIters;
		loopHead->patternArray = malloc(numIters * sizeof(pattern));
		pattern_Init(loopHead->patternArray, numIters);
		loopHead->next = NULL;
	}
	else {
		loop *loopTemp = loopHead;
		while(loopTemp->next != NULL) 
			loopTemp = loopTemp->next;
		loopTemp->next = malloc(sizeof(loop));
		loopTemp->next->numIters = numIters;
		loopTemp->next->patternArray = malloc(numIters * sizeof(pattern));
		pattern_Init(loopHead->next->patternArray, numIters);
		loopTemp->next->next = NULL;
	}
	return loopHead;
}

////

int main(int argc, char* argv[]) {
    int numProcs, numIters, loopIndex, i, j, begin, end, lineIndex;
    int numSendTemp, numRecvTemp, sumSrcTemp, sumDestTemp, xorSendTemp, xorRecvTemp;
    char fileName[8];
    char *buffer, *strTemp;
    iter *iterHead, *iterTemp;
    long start_time, end_time, elapsed;
	loop *loopHead, *loopTemp;
	
	loopHead = loop_Init(loopHead);

    start_time = clock();
    if(argc != 2) {
        printf("Please enter number of process as first argument\n");
        return 0;
    } else {
        numProcs = atoi(argv[1]);
    }

    /*FILE * fp = fopen("iter", "r");
    char * line = NULL;
    size_t len = 0;
    getline(&line, &len, fp);
    numIters = atoi(line);
    fclose(fp);*/
    
    for(i = 0; i < numProcs; i++) {
		//!**caution: don't use "i" at anywhere in this for if you don't understand it's function**!
        sprintf(fileName, "%d", i);
        FILE *pFile = fopen(fileName, "r");
        if(pFile == NULL)
            perror("No such file exists");
        else {
            buffer = (char*) malloc(10000000);
            lineIndex = 0;
			loopIndex = 0;
            while(fgets(buffer, 10000000, pFile) != NULL) {
				if(buffer[0] == '[') {
					loopIndex++;
					if(i == 0) {
						strTemp = (char*) malloc((strlen(buffer) - 2) * sizeof(char));
						for(j = 1; j < strlen(buffer) - 2; j++)
							strTemp[j - 1] = buffer[j]; 
						strTemp[strlen(buffer) - 3] = '\0';
						numIters = atoi(strTemp);
						free(strTemp);
						loopHead = loop_Add(loopHead, numIters);
					}
					else {}
				}
				else if(buffer[0] == '!') {}
				else {
					lineIndex++;
					switch(lineIndex % 8) {
					    case 0:
					        iterTemp = iterHead;
					        while(iterTemp != NULL) {
								loopTemp = loopHead;
								for(j = 0; j < loopIndex - 1; j++)
									loopTemp = loopTemp->next;
					            pattern_Update(loopTemp->patternArray, iterTemp->index, i, numSendTemp, numRecvTemp,
														 sumSrcTemp, sumDestTemp, xorSendTemp, xorRecvTemp);
					            iterTemp = iterTemp->next;
					        }
					        iterHead = iter_Finalize(iterHead);
					        break;
					    case 1:
					        iterHead = iter_Init(iterHead);
					        for(begin = 0; begin < strlen(buffer); begin++) {
					            int temp = buffer[begin] - '0';
					            if(temp >= 0 && temp <= 9) {
					                end = begin + 1;
					                while(buffer[end] != ' ')
					                    end++;
					                strTemp = (char*) malloc((end - begin + 1) * sizeof(char));
					                for(j = begin; j < end; j++) {
					                    strTemp[j - begin] = buffer[j];
					                    //printf("%c", buffer[j]);
					                }
					                //printf("\n");
					                strTemp[end - begin] = '\0';
					                iterHead = iter_Add(iterHead, atoi(strTemp));
					                free(strTemp);
					                begin = end;
					            }
					        }	
					        break;
					    case 2:
					        buffer[strlen(buffer) - 1] = '\0';
					        numSendTemp = atoi(buffer);
					        break;
					    case 3:
					        buffer[strlen(buffer) - 1] = '\0';
					        numRecvTemp = atoi(buffer);
					    break;
					    case 4:
					        buffer[strlen(buffer) - 1] = '\0';
					        sumSrcTemp = atoi(buffer);
					    break;
					    case 5:
					        buffer[strlen(buffer) - 1] = '\0';
					        sumDestTemp = atoi(buffer);
					    break;
					    case 6:
					        buffer[strlen(buffer) - 1] = '\0';
					        xorSendTemp = atoi(buffer);
					    break;
					    case 7:
					        buffer[strlen(buffer) - 1] = '\0';
					        xorRecvTemp = atoi(buffer);
					        break;
					    default:
					        break;
					}
				}
				free(buffer);
				buffer = (char*) malloc(10000000);
            }
            free(buffer);
        }
		fclose(pFile);
    }
    
    /* Testing */
	loopTemp = loopHead;

    int sNumSend = 0;
    int sNumRecv = 0;
    int sSumSrc = 0;
    int sSumDest = 0;
    int sSumRankSend = 0;
    int sSumRankRecv = 0;
    int sXorSend = 0;
    int sXorRecv = 0;
	int loop_index = 1;
	FILE *fo = fopen("output","w");
    if (fo == NULL) {
       	printf("Error opening file !\n");
		exit(1);
    }	

	while(loopTemp != NULL) {	
		printf("Loop index = %d ; numIters = %d\n", loop_index, loopTemp->numIters);
		fprintf(fo, "Loop index = %d ; numIters = %d\n", loop_index, loopTemp->numIters);
		for(j = 0; j < loopTemp->numIters; j++) {
			fprintf(fo, "Iters[%d] : ", j);
			printf("Iters[%d] : ", j);

			sNumSend += loopTemp->patternArray[j].numSend;
        	sNumRecv += loopTemp->patternArray[j].numRecv;
        	sSumSrc += loopTemp->patternArray[j].sumSrc;
        	sSumDest += loopTemp->patternArray[j].sumDest;
        	sSumRankSend += loopTemp->patternArray[j].sumRankSend;
        	sSumRankRecv += loopTemp->patternArray[j].sumRankRecv;
        	sXorSend ^= loopTemp->patternArray[j].xorSend;
        	sXorRecv ^= loopTemp->patternArray[j].xorRecv;

			if (sNumSend != sNumRecv) {
				fprintf(fo, "Message Leaks !!!\n");
	            printf("Message Leaks !!!\n");
    	    } else {
    	        if (sSumSrc != sSumRankSend || sSumDest != sSumRankRecv) {
					fprintf(fo, "Message Leaks !!!\n");
            	    printf("Message Leaks !!!\n");
            	} else {
                	if (sXorSend != sXorRecv) {
						fprintf(fo, "Message Leaks !!!\n");
		    			printf("Message Leaks !!!\n");
                	} else {
						fprintf(fo, "No Message leaks !!! \n");
                    	printf("No Message leaks !!! \n");
                	}
            	}
        	}

		    fprintf(fo, "numSend = %d\n", loopTemp->patternArray[j].numSend);
		    fprintf(fo, "numRecv = %d\n", loopTemp->patternArray[j].numRecv);
		    fprintf(fo, "sumSrc = %d\n", loopTemp->patternArray[j].sumSrc);
		    fprintf(fo, "sumDest = %d\n", loopTemp->patternArray[j].sumDest);
		    fprintf(fo, "sumRankSend = %d\n", loopTemp->patternArray[j].sumRankSend);
		    fprintf(fo, "sumRankRecv = %d\n", loopTemp->patternArray[j].sumRankRecv);
		    fprintf(fo, "xorSend = %d\n", loopTemp->patternArray[j].xorSend);
		    fprintf(fo, "xorRecv = %d\n", loopTemp->patternArray[j].xorRecv);
		}
		loopTemp = loopTemp->next;
		loop_index++;
		fprintf(fo, "\n");
		printf("\n");
	}
    /*int sNumSend = 0;
    int sNumRecv = 0;
    int sSumSrc = 0;
    int sSumDest = 0;
    int sSumRankSend = 0;
    int sSumRankRecv = 0;
    int sXorSend = 0;
    int sXorRecv = 0;
    FILE *fo = fopen("output","w");
    if (fo == NULL) {
        printf("Error opening file !\n");
	exit(1);
    }	
    for(i = 0; i < numIters; i++) {
        sNumSend += patternArray[i].numSend;
        sNumRecv += patternArray[i].numRecv;
        sSumSrc += patternArray[i].sumSrc;
        sSumDest += patternArray[i].sumDest;
        sSumRankSend += patternArray[i].sumRankSend;
        sSumRankRecv += patternArray[i].sumRankRecv;
        sXorSend ^= patternArray[i].xorSend;
        sXorRecv ^= patternArray[i].xorRecv;
	fprintf(fo, "Iter[%d] : \n", i);
        printf("Iter[%d]:\n", i);
        if (sNumSend != sNumRecv) {
	    fprintf(fo, "Message Leaks !!!\n");
            printf("Message Leaks !!!\n");
        } else {
            if (sSumSrc != sSumRankSend || sSumDest != sSumRankRecv) {
		fprintf(fo, "Message Leaks !!!\n");
                printf("Message Leaks !!!\n");
            } else {
                if (sXorSend != sXorRecv) {
		    fprintf(fo, "Message Leaks !!!\n");
                    printf("Message Leaks !!!\n");
                } else {
		    fprintf(fo, "No message leaks !!!\n");
                    printf("No Message leaks !!! \n");
                }
            }
        }

        fprintf(fo, "numSend = %d\n", patternArray[i].numSend);
        fprintf(fo, "numRecv = %d\n", patternArray[i].numRecv);
        fprintf(fo, "sumSrc = %d\n", patternArray[i].sumSrc);
        fprintf(fo, "sumDest = %d\n", patternArray[i].sumDest);
        fprintf(fo, "sumRankSend = %d\n", patternArray[i].sumRankSend);
        fprintf(fo, "sumRankRecv = %d\n", patternArray[i].sumRankRecv);
        fprintf(fo, "xorSend = %d\n", patternArray[i].xorSend);
        fprintf(fo, "xorRecv = %d\n", patternArray[i].xorRecv);
    } */

    end_time = clock();
    elapsed = (end_time - start_time) / CLOCKS_PER_SEC * 1000;
    fprintf(fo, "Elapsed Time : %ld", elapsed);
    printf("Elapsed Time : %ld\n", elapsed);
    fclose(fo);

    

    loopHead = loop_Finalize(loopHead);
    return 0;
}
