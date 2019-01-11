#include "function.h"
#include "math.h"

//matlab支援的指令用C做出 因為會重複呼叫 所以額外寫成function 
//這裡的function都沒有額外檢查輸入的矩陣及向量能不能相乘 有可能出錯 
//這裡的每個function的output都有分配記憶體空間 
double* vector_sub(double* v1, double* v2, int size) {  //向量減法 
	int i;
	double *output;
	output = (double*)malloc(sizeof(double) * size);
	
	for(i = 0; i < size; i++)
		output[i] = v1[i] - v2[i];
	
	return output;
}

double* matrix_vector_multiply(double** A, double* v, int Arow, int Acol) {  //矩陣*向量 
	double* output;
	output = (double*)malloc(sizeof(double) * Arow);
	
	int i,j;
	for(i = 0; i < Arow; i++)
		output[i] = 0;
	for(i = 0; i < Arow; i++) {
		for(j = 0; j < Acol; j++)
			output[i] += A[i][j] * v[j];
	}
	
	return output;
}

double norm(double* v, int size) {  //計算向量v的長度 即||v|| 
	int i;
	double output = 0;
	for(i = 0; i < size; i++)
		output += v[i]*v[i];
	
	return sqrt(output);
}

//matlab的指令
//做法為先宣告一個剛好大小的矩陣 將它切成維度與A相同的區塊矩陣 每個區塊放入B*A[i][j]的矩陣 
double** kron(double** A, int Arow, int Acol, double** B, int Brow, int Bcol) {
	
	int i,j,i_Block,j_Block;
	int idx_na,idx_ma,idx_nb,idx_mb;
	int new_row = Arow * Brow;
	int new_col = Acol * Bcol;
	double **output;
	output = (double**)malloc(sizeof(double*) * new_row);
	for(i = 0; i < new_row; i++)
		output[i] = (double*)malloc(sizeof(double) * new_col);
	
	for( i = 0; i < new_row; i++ )
	{
		for( j = 0; j < new_col; j++ )
		{

			idx_na = (i+1)/Brow + 1;
			idx_ma = (j+1)/Bcol + 1;
			idx_nb = (i+1)%Brow;
			idx_mb = (j+1)%Bcol;

			if(idx_nb == 0)
			{
				idx_na = idx_na - 1;
				idx_nb = Brow;
			}
			if(idx_mb == 0)
			{
				idx_ma = idx_ma - 1;
				idx_mb = Bcol;
			}

			output[i][j] = A[idx_na-1][idx_ma-1] * B[idx_nb-1][idx_mb-1];
		}
	}
	return output;
}

double** two_matrix_multiply(double** A, double** B, int rowA, int colA, int colB) {  //矩陣乘法 
	int i,j,k;
	double **output;
	output = (double**)malloc(sizeof(double*) * rowA);
	for(i = 0; i < rowA; i++)
		output[i] = (double*)malloc(sizeof(double) * colB);
	
	for(i = 0; i < rowA; i++) {
		for(j = 0; j < colB; j++) {
			output[i][j] = 0;
			for(k = 0; k < colA; k++) 
				output[i][j] += A[i][k] * B[k][j];
		}
	}
	
	return output;
} 

double** inverse(double** A, int ASize) {  //求矩陣的反矩陣 用最基本的高斯消去法 [A|I] -> [I|A'] 
	int i,j,k;
	double** output;
	output = (double**)malloc(sizeof(double*) * ASize);
	for(i = 0; i < ASize; i++) 
		output[i] = (double*)malloc(sizeof(double) * ASize);
		
	for(i = 0; i < ASize; i++) 
		for(j = 0; j < ASize; j++)
			output[i][j] = 0;
	
	for(i = 0; i < ASize; i++)
		output[i][i] = 1;  //output = I
	
	for(i = 0; i < ASize-1; i++) {  //將A的下三角清為0 
		for(j = i+1; j < ASize; j++) {
			if(A[i][i] != 0) {
				double div = A[j][i] / A[i][i] * -1;
				for(k = 0; k < ASize; k++) {
					A[j][k] += A[i][k] * div;
					output[j][k] += output[i][k] * div;
				}
			}
		}
	}
	
	for(i = 0; i < ASize; i++) {  //A的對角線除為1 
		double temp = A[i][i];
		for(j = 0; j < ASize; j++) {
			A[i][j] /= temp;
			output[i][j] /= temp;
		}
	}
	
	for(i = ASize-1; i > 0; i--) {  //將A的上三角清為0 
		for(j = i-1; j >= 0; j--) {
			if(A[i][i] != 0) {
				double div = A[j][i] / A[i][i] * -1;
				for(k = 0; k < ASize; k++) {
					A[j][k] += A[i][k] * div;
					output[j][k] += output[i][k] * div;
				}
			}
		}
	}
	
	return output;
}

double* vector_add(double* v1, double alpha, double* v2,int size) {  //向量加法 
	int i;
	double *output;
	output = (double*)malloc(sizeof(double) * size);
	
	for(i = 0; i < size; i++)
		output[i] = v1[i] + alpha*v2[i];
	
	return output;
}
