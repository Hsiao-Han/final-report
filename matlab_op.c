#include "function.h"
#include "math.h"

//matlab�䴩�����O��C���X �]���|���ƩI�s �ҥH�B�~�g��function 
//�o�̪�function���S���B�~�ˬd��J���x�}�ΦV�q�ण��ۭ� ���i��X�� 
//�o�̪��C��function��output�������t�O����Ŷ� 
double* vector_sub(double* v1, double* v2, int size) {  //�V�q��k 
	int i;
	double *output;
	output = (double*)malloc(sizeof(double) * size);
	
	for(i = 0; i < size; i++)
		output[i] = v1[i] - v2[i];
	
	return output;
}

double* matrix_vector_multiply(double** A, double* v, int Arow, int Acol) {  //�x�}*�V�q 
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

double norm(double* v, int size) {  //�p��V�qv������ �Y||v|| 
	int i;
	double output = 0;
	for(i = 0; i < size; i++)
		output += v[i]*v[i];
	
	return sqrt(output);
}

//matlab�����O
//���k�����ŧi�@�ӭ�n�j�p���x�} �N���������׻PA�ۦP���϶��x�} �C�Ӱ϶���JB*A[i][j]���x�} 
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

double** two_matrix_multiply(double** A, double** B, int rowA, int colA, int colB) {  //�x�}���k 
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

double** inverse(double** A, int ASize) {  //�D�x�}���ϯx�} �γ̰򥻪��������h�k [A|I] -> [I|A'] 
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
	
	for(i = 0; i < ASize-1; i++) {  //�NA���U�T���M��0 
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
	
	for(i = 0; i < ASize; i++) {  //A���﨤�u����1 
		double temp = A[i][i];
		for(j = 0; j < ASize; j++) {
			A[i][j] /= temp;
			output[i][j] /= temp;
		}
	}
	
	for(i = ASize-1; i > 0; i--) {  //�NA���W�T���M��0 
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

double* vector_add(double* v1, double alpha, double* v2,int size) {  //�V�q�[�k 
	int i;
	double *output;
	output = (double*)malloc(sizeof(double) * size);
	
	for(i = 0; i < size; i++)
		output[i] = v1[i] + alpha*v2[i];
	
	return output;
}
