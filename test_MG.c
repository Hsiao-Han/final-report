#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "function.h"

#define Pstart 0
#define Pend 1
#define part_num 16

int main() {
	int i, j;
	int input_size = part_num;  //切等分 
	int array_length = (input_size - 1) * (input_size - 1);  //依據等分 確定矩陣大小 

	double **A;  //宣告A, F 
	double *F;

	A = (double**)malloc(sizeof(double*) * array_length);  //分配A, F的記憶體空間 
	for (i = 0; i < array_length; i++)
		A[i] = (double*)malloc(sizeof(double) * array_length);
	F = (double*)malloc(sizeof(double) * array_length);
	
	//給定 
	for (i = 0; i < array_length; i++) {
		for (j = 0; j < array_length; j++) {  //初始化A 
			if (i == j)
				A[i][j] = 1;
			else if(j - i == 1)
				A[i][j] = -1;
			else
				A[i][j] = 0;
			
		}
	}

	for (i = 0; i < array_length; i++)  //初始化F 
		F[i] = 0.5;

	double *u_f;
	u_f = (double*)malloc(sizeof(double) * array_length);
	srand(time(NULL));
	for (i = 0; i < array_length; i++)
		u_f[i] = (double)rand() / (RAND_MAX + 1.0);  //宣告, 分配空間, 初始化u_f, u_f為Ax=F的解 

	int m = 3;  //V-cycle深度 
	double TOL = pow(10, -8);  //最小誤差 
	printf("TOL = %e\n\n",TOL);
	double *u;
	u = (double*)malloc(sizeof(double) * array_length);  //u複製u_f, u,u_f可刪一個?? 
	for (i = 0; i < array_length; i++)
		u[i] = u_f[i];
	double res = 1;  //u_f與F的誤差 
	int k = 1;  //迴圈控制條件,避免無窮迴圈 
	while (res > TOL && k < 500) {  //當誤差小於最小誤差或迴圈條件到 跳出迴圈 
		printf("\nMG loop %d\n\n", k);
		u = MG(A, array_length, F, u, m, TOL, 3);  //對u_f做multigrid 
		k++;
		double *temp1, *temp2;  //計算u_f與F的誤差  
		temp1 = matrix_vector_multiply(A, u, array_length, array_length);  //temp1 = A*u_f 
		temp2 = vector_sub(temp1, F, array_length);  //temp2 = F-A*u_f
		res = norm(temp2, array_length);
		res /= norm(F, array_length);    //res = ||F-A*u_f|| / ||F|| 
		printf("residual = %e\n\n", res);
		
		free(temp1);  //回收記憶體空間 
		free(temp2);
	}

	if (k < 500)
		printf("success in loop %d\n", k);  //輸出成功或失敗 
	else
		printf("fail within %d loop", k);
	
	return 0;
}
