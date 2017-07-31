#ifndef __MEMORY_H__
#define __MEMORY_H__

#include <stdio.h>
#include <stdlib.h>
#include "string.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int parseLine(char* line);

int getMemory();

void printMemory();

#ifdef __cplusplus 
}
#endif /* __cplusplus */

#endif /* __MEMORY_H__ */ 
