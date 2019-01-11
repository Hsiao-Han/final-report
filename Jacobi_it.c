#include "function.h"

double* Jacobi_it(double** A, int ASize, double* b, double* u_f, double w, double TOL, int itmax) {
	int i,j;
	int nit=1;  //控制Jacobi執行次數 
	double *New_u_f;  //調整後的u_f 
	New_u_f = (double*)malloc(sizeof(double)*ASize);
	
	while(nit <= itmax) {
		for(i = 0; i < ASize; i++) {
			double temp = 0;
			for(j = 0; j < ASize; j++) {
				if(i != j)
					temp += A[i][j] * u_f[j];
			}
			New_u_f[i] = w * (b[i] - temp) / A[i][i] + (1.0-w) * u_f[i];  //不知原因 
		}
		
		printf("Jacobit loop %d\n", nit);  //輸出每次做Jacobi後的u_f 

		double *temp1 = matrix_vector_multiply(A, New_u_f, ASize, ASize);  //temp1 = A * New_u_f
		double *temp2 = vector_sub(b, temp1, ASize);  //temp2 = F - A*New_u_f  b為F在這個函式的local name 
		
		if( (norm(temp2, ASize) / norm(b, ASize) ) < TOL) {  //計算是否已到最小誤差 
			printf("success at %d -th step!\n", nit);
			return New_u_f;
		}
		else {  //將新的u_f copy 到u_f上繼續做下一輪的Jacobi 
			nit++;
			for(i = 0; i < ASize; i++)
				u_f[i] = New_u_f[i];
		}
		free(temp1);
		free(temp2);
	}
	
	printf("fail...\n\n");
	return New_u_f;
}

