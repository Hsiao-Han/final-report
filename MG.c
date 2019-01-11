#include "function.h"
#include <math.h>
#include <pthread.h>
#include <assert.h>

pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
double* e_f;

void *threadFunc(void *thread_data) {
	unsigned int thread_id   = ((struct THREAD_DATA*)thread_data) -> thread_id;
	unsigned int flag_cut    = ((struct THREAD_DATA*)thread_data) -> flag_cut;
	unsigned int size_fine   = ((struct THREAD_DATA*)thread_data) -> size_fine;
	unsigned int size_coarse = ((struct THREAD_DATA*)thread_data) -> size_coarse;
	unsigned int square_size_fine   = ((struct THREAD_DATA*)thread_data) -> square_size_fine;
	unsigned int square_size_coarse = ((struct THREAD_DATA*)thread_data) -> square_size_coarse;

	double* r_f = ((struct THREAD_DATA*)thread_data) -> r_f;
	double** A  = ((struct THREAD_DATA*)thread_data) -> A;
	unsigned int ASize;

	struct two_matrix* Rst_Itp = gen_restric(size_fine  , size_coarse, (int) 1, flag_cut);  //ºâ¥XRestrict matrix ¤Î Interpolate matrix 
	//struct two_matrix* Rst_Itp_red  = gen_restric(size_fine+2, size_coarse+1, op, 1);  //ºâ¥XRestrict matrix ¤Î Interpolate matrix 

	double *r_c = matrix_vector_multiply(Rst_Itp->Restrict, r_f, square_size_coarse, square_size_fine);

	//­pºâ²Êºô®æªº¯x°} ¤é«á­n§ï 
	double **temp2 = two_matrix_multiply(Rst_Itp->Restrict, A, square_size_coarse, square_size_fine, square_size_fine);
	double **A_c = two_matrix_multiply(temp2, Rst_Itp->Interpolate, square_size_coarse, square_size_fine, square_size_coarse);

	double *e_c = matrix_vector_multiply(inverse(A_c, square_size_coarse), r_c, square_size_coarse, square_size_coarse);

	double *e_f_local = matrix_vector_multiply(Rst_Itp->Interpolate, e_c, square_size_fine, square_size_coarse);  //±Ncoarse-grid Âà¦¨ fine-grid­×¥¿¤è¦V 
    
	pthread_mutex_lock( &mutex1 ); // 銝?
	for(int i = 0; i < ASize; i++)
	{
		e_f[i] = e_f[i] + 0.5*e_f_local[i];
	}
	pthread_mutex_unlock( &mutex1 ); // 閫??
	
	free(r_c);
	free(e_c);
	free(e_f_local);
	
	int i;
	for (i = 0; i < square_size_coarse; i++)
		free(temp2[i]);
	free(temp2);
	
	for (i = 0; i < square_size_coarse; i++)
		free(A_c[i]);
	free(A_c);
	
	int c_square_size_fine = (size_fine-2)*(size_fine-2);
	for(i = 0; i < c_square_size_fine; i++)
		free(Rst_Itp->Interpolate[i]);
	free(Rst_Itp->Interpolate);
	
	for(i = 0; i < square_size_coarse; i++)
		free(Rst_Itp->Restrict[i]);
	free(Rst_Itp->Restrict);
	
	free(Rst_Itp);

	pthread_exit(NULL);
}

double* MG(double** A, int ASize, double* F, double* u_f, int m, double TOL, int itmax) {
	unsigned int thread_num  = 2;
	int err;
	pthread_t *thread_array = malloc(thread_num * sizeof(pthread_t));
	pthread_mutex_init(&mutex1, NULL);
	struct THREAD_DATA* thread_data = malloc(thread_num * sizeof(struct THREAD_DATA));

	int i,j;
	double w = 0.5;  //op w ¼È®Éµ¹©w 
	double *output;  //³Ì«áªº¿é¥X 

	unsigned int size_fine = (int)sqrt((double)ASize);  //²Óºô®æ¤j¤p 
	unsigned int size_coarse = (size_fine-1) / 2;       //²Êºô®æ¤j¤p 
	unsigned int square_size_fine = size_fine * size_fine;
	
	e_f = (double*) calloc(ASize, sizeof(double));

	u_f = Jacobi_it(A, ASize, F, u_f, w, TOL, itmax);  //¹ïu_f°µpre-smooth 
	
	double *temp1 = matrix_vector_multiply(A, u_f, ASize, ASize);  //ºâ¥Xfine-grid residual  r_f = F-A*u_f 
	double *r_f = vector_sub(F, temp1, ASize);  
	
	for(unsigned int i = 0; i < thread_num; i++)
  	{
  		thread_data[i].thread_id = i;
	    thread_data[i].flag_cut = i;
	    thread_data[i].square_size_fine = square_size_fine;
	    thread_data[i].r_f = r_f;
	    thread_data[i].A = A;
	    thread_data[i].ASize = ASize;
  		if(i == 0){
	    	thread_data[i].size_fine = size_fine;
	    	thread_data[i].size_coarse = size_coarse;
	    	thread_data[i].square_size_coarse = size_coarse * size_coarse;
	    }
	    else if( i == 1 ){
	    	thread_data[i].size_fine = size_fine + 2;
	    	thread_data[i].size_coarse = size_coarse + 1;
	    	thread_data[i].square_size_coarse = (size_coarse+1) * (size_coarse+1);
	    }
	    err = pthread_create(&thread_array[i], NULL, threadFunc, (void*) &thread_data[i]);
	    assert(err == 0);
	}
	for (int i = 0; i < thread_num; i++)
	{
		err = pthread_join(thread_array[i], NULL);
		assert(err == 0);
	}
	pthread_mutex_destroy(&mutex1);
  	free(thread_array);

  	for( i = 0; i < ASize; i++)
  		u_f[i] = u_f[i] + e_f[i];

	output = Jacobi_it(A, ASize, F, u_f, w, TOL, itmax);  //¹ïu_f°µpost-smooth 
	
	free(e_f);
	free(temp1);  //¦¬¦^°O¾ÐÅéªÅ¶¡ ¦]¬°matlab_op§Ú³£¦³¤À°t°O¾ÐÅé ©Ò¥H¤j¦h¼Æªº°}¦C³£­n¦^¦¬ 
	free(u_f);    //©Ò¦³.cÀÉ³£¥u¦³²Ê²¤¦^¦¬ ¥i¥H§ó¶i¤@¨Bªº´î¤Ö¨Ï¥ÎªÅ¶¡ 
		
	return output;
}
