n: gen_restric.c Jacobi_it.c matlab_op.c MG.c test_MG.c
	gcc -pthread -std=gnu99 -O2 gen_restric.c Jacobi_it.c matlab_op.c MG.c test_MG.c -o test_MG -lm

clean: 
	-rm test_MG
