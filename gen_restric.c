#include "function.h"

struct two_matrix* gen_restric(int size_fine, int size_coarse, int op, int flag_cut) {
	struct two_matrix* output;  //宣告回傳值 
	output = (struct two_matrix*)malloc(sizeof(struct two_matrix));  //分配output空間 
	 
	double **R_f_c_1d, **R_f_c_2d;
	int i,j;
	R_f_c_1d = (double**)malloc(sizeof(double*) * size_coarse);
	for(i = 0; i < size_coarse; i++)
		R_f_c_1d[i] = (double*)malloc(sizeof(double) * size_fine);

	int square_size_coarse = size_coarse * size_coarse;
	int square_size_fine = size_fine * size_fine;
	R_f_c_2d = (double**)malloc(sizeof(double*) * square_size_coarse);
	for(i = 0; i < square_size_coarse; i++)
		R_f_c_2d[i] = (double*)malloc(sizeof(double) * square_size_fine);
	
	//generate 1d restriction. ex. [x1,x2,x3,x4,x5,x6,x7] -> [x2,x4,x6]
	for(i = 0; i < size_coarse; i++) {
		for(j = 0; j < size_fine; j++)
			R_f_c_1d[i][j] = 0;
	}
	for(i = 0; 2*i+1 < size_fine && (i<size_coarse); i++)
		R_f_c_1d[i][2*i+1] = 1;
	//generate 2d restriction
	R_f_c_2d = kron(R_f_c_1d, size_coarse, size_fine, R_f_c_1d, size_coarse, size_fine);

	//generate filter matrix
	double **Ef, **Ef_temp, **E_f;
	if( flag_cut == 0 )
	{
		Ef = (double**)malloc(sizeof(double*) * size_fine);
		for(i = 0; i < size_fine; i++)
			Ef[i] = (double*)malloc(sizeof(double) * size_fine);	
	}
	else if( flag_cut == 1 )
	{
		Ef = (double**)malloc(sizeof(double*) * size_fine);
		for(i = 0; i < size_fine; i++)
			Ef[i] = (double*)malloc(sizeof(double) * (size_fine-2));	
	}	

	Ef_temp = (double**)malloc(sizeof(double*) * size_fine);
	for(i = 0; i < size_fine; i++)
		Ef_temp[i] = (double*)malloc(sizeof(double) * size_fine);
	
	for(i = 0; i < size_fine; i++) {
		for(j = 0; j < size_fine; j++) {
			if(i == j)
				Ef_temp[i][j] = 0.5;
			else if(i-j ==1 || j-i ==1)
				Ef_temp[i][j] = 0.25;
			else
				Ef_temp[i][j] = 0;
		}
	}
    if( flag_cut == 0 )
	{
		for( i = 0; i < size_fine; i++ ){
			for( j = 0; j < size_fine; j++ )
				Ef[i][j] = Ef_temp[i][j];
		}
		E_f = kron(Ef_temp, size_fine, size_fine, Ef_temp, size_fine, size_fine);
		output->Restrict = two_matrix_multiply(R_f_c_2d, E_f, square_size_coarse, square_size_fine, square_size_fine);
		output->Interpolate = (double**)malloc(sizeof(double*) * square_size_fine);
		
		for(i = 0; i < square_size_fine; i++)
			output->Interpolate[i] = (double*)malloc(sizeof(double) * square_size_coarse);

		for(i = 0; i < square_size_fine; i++) {
			for(j = 0; j < square_size_coarse; j++)
				output->Interpolate[i][j] = 4 * output->Restrict[j][i];
		}

	}
	else if( flag_cut == 1 )
	{
		for( i = 0; i < size_fine; i++ ){
			for( j = 0; j < size_fine-2; j++ ){
				Ef[i][j] = Ef_temp[i][j+1];
			}
		}
		E_f = kron(Ef, size_fine, size_fine-2, Ef, size_fine, size_fine-2);

		output->Restrict = two_matrix_multiply(R_f_c_2d, E_f, square_size_coarse, square_size_fine, (size_fine-2)*(size_fine-2));
		output->Interpolate = (double**)malloc(sizeof(double*) * (size_fine-2)*(size_fine-2));
		
		int c_square_size_fine = (size_fine-2)*(size_fine-2);
		for(i = 0; i < c_square_size_fine; i++)
			output->Interpolate[i] = (double*)malloc(sizeof(double) * square_size_coarse);
		for(i = 0; i < (size_fine-2)*(size_fine-2); i++){
			for(j = 0; j < square_size_coarse; j++)
				output->Interpolate[i][j] = 4 * output->Restrict[j][i];
		}

	}		

	for(i = 0; i < size_coarse; i++)
		free(R_f_c_1d[i]);
	free(R_f_c_1d);
	
	for(i = 0; i < square_size_coarse; i++)
		free(R_f_c_2d[i]);
	free(R_f_c_2d);
	
	for(i = 0; i < size_fine; i++)
		free(Ef_temp[i]);
	free(Ef_temp);
	
	for(i = 0; i < square_size_fine; i++)
		free(E_f[i]);
	free(E_f);
	
	return output;
}
