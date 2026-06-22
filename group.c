#include "group.h"

#include <stdio.h>
#include <stdlib.h>

unsigned P;
unsigned C;
double ReTr[PMAX];
double ImTr[PMAX];
double ReChar[CMAX][CMAX];
double ImChar[CMAX][CMAX];
group_t mult[PMAX][PMAX];
group_t inv[PMAX];
group_t conclass[PMAX];
group_t id;

void load_group(const char *fn) {
	FILE *fin = fopen(fn, "r");
	/* File format for group specification is as follows:
	 *
	 *  <order>
	 *  ReChar00 ReChar01 ReChar02 ...
	 *  ...
	 *  ImChar00 ImChar01 ImChar02 ...
	 *  ...
	 *  0x0 0x1 0x2 ...
	 *  1x0 1x1 1x2 ...
	 *  ...
	 *  ConClass0 ConClass1 ConClass2 ...
	 * For a group of order P, there are thus 1 + P + P^2 entries. Whitespace is ignored.
	 */
	fscanf(fin, "%d", &P);
	if (P > PMAX) {
		fprintf(stderr, "Order of group too large: %d > %d\n", P, PMAX);
		abort();
	}
        fscanf(fin, "%d", &C);
        if (C > CMAX) {
                fprintf(stderr, "Number of Classes too large: %d > %d\n", C, CMAX);
                abort();
        }

	for (unsigned n = 0; n < C; n++){
		for (unsigned m = 0; m < C; m++){
			fscanf(fin, "%lf", &ReChar[n][m]);
	//		printf("%e ",ReChar[n][m]);
		}
	//	printf("\n");	
	}
	//printf("\n");	
	for (unsigned n = 0; n < C; n++){
		for (unsigned m = 0; m < C; m++){
                	fscanf(fin, "%lf", &ImChar[n][m]);
        //                printf("%e ",ImChar[n][m]);
		}
        //        printf("\n");
	}
	for (unsigned n = 0; n < P; n++){
		for (unsigned m = 0; m < P; m++){
			fscanf(fin, "%d", &mult[n][m]);
        //               printf("%d ",mult[n][m]);
		}
        //        printf("\n");
	}
	for (unsigned n = 0; n < P; n++){
		fscanf(fin, "%d", &conclass[n]);
        //        printf("%d ",conclass[n]);
	}
        //        printf("\n");
	for (unsigned n = 0; n < P; n++){
                //The row of ReChar put here depends on the group s1080=1, s108=4!!!
                ReTr[n] = ReChar[1][conclass[n]];
                ImTr[n] = ImChar[1][conclass[n]];
		//The row of ReChar put here depends on the group!!!
//        	ReTr[n] = ReChar[4][conclass[n]];	
//        	ImTr[n] = ImChar[4][conclass[n]];
//		printf("%e + i*%e \n", ReTr[n], ImTr[n]);
	}
	//printf("\n %e + i*%e ",ReTr[0],ImTr[0]);
	fclose(fin);

	// Find the identity.
	char id_found = 0;
	for (unsigned n = 0; n < P; n++) {
		char is_id = 1;
		for (unsigned m = 0; m < P; m++)
			if (mult[n][m] != m || mult[m][n] != m)
				is_id = 0;
		if (is_id) {
			id_found = 1;
			id = n;
			break;
		}
	}
	if (!id_found) {
		fprintf(stderr, "Group does not have identity element\n");
		abort();
	}

	// Find inverses.
	for (unsigned n = 0; n < P; n++) {
		char inv_found = 0;
		for (unsigned m = 0; m < P; m++)
			if (mult[n][m] == id && mult[m][n] == id) {
				inv[n] = m;
				inv_found = 1;
				break;
			}
		if (!inv_found) {
			fprintf(stderr, "Group does not have inverses\n");
			abort();
		}
	}
}
