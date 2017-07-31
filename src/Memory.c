#include "Memory.h"

/*struct sysinfo memInfo;

sysinfo (&memInfo);

long long totalVirtualMem = memInfo.totalram;
totalVirtualMem += memInfo.totalswap;
totalVirtualMem *= memInfo.mem_unit;*/

int parseLine(char* line){
	int i = strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
	i = atoi(line);
	return i;
}

int getMemory(){ //Note: this value is in KB!
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];

	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmRSS:", 6) == 0){
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
}

void printMemory() {
	printf("%i KB \n", getMemory());
}

/*int main(int argc, char* argv[]) {
	printf("%i KB \n", getValue());
	return 0;
}*/