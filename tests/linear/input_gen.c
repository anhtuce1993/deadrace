#include <stdio.h>
#include <stdlib.h>

double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

int main(int argc, char* argv[]) {

	int size;
	if(argc != 2) {
        	printf("Please enter size as first argument\n");
        	return 0;
 	} else {
        	size = atoi(argv[1]);
 	}

	int i, j;
	FILE *fi = fopen("mdatgaus.inp","w");
	if (fi == NULL) {
        	printf("Error opening file !\n");
		exit(1);
    	}
	fprintf(fi, "%d %d", size, size);
	for(i = 0; i < size; i++){
		fprintf(fi, "\n");
		for(j = 0; j < size; j++){
			fprintf(fi, "%f ", randfrom(1.0, 100000.0));
		}	
	}
	fclose(fi);	

	FILE *fo = fopen("vdatgaus.inp","w");
	if (fo == NULL) {
        	printf("Error opening file !\n");
		exit(1);
    	}
	fprintf(fo, "%d\n", size);
	for(i = 0; i < size; i++){
		fprintf(fo, "%f ", randfrom(1.0, 1000.0));		
	}
	fclose(fo);	
	return 0;

}
