#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mackiepheap.h"

//#define DEBUG

void *matrix_gentrification(void *ms);
void printMatrix(matrix l);

double matrix_multiplication(matrix *kapMatrix, pthread_t *threads){
  matrix_package matrix_state[NUM_THREADS];
  int jac;      //just another counter(jac)
  int rv;       //check for thread creation

  //Thread dispatch functionality - dispatch a thread
  //for every thread that is capable of doing work.
  for(jac = 0;jac < NUM_THREADS; jac++){

    matrix_state[jac].currentMatrix = kapMatrix;
    matrix_state[jac].thread_max = NUM_THREADS;
    matrix_state[jac].thread_id = jac;
#ifdef DEBUG
    printf("DISPATCHING THREAD %d\n", jac);
#endif
    
    rv = pthread_create(&threads[jac],
			NULL,
			matrix_gentrification,
			(void*)&matrix_state[jac]);
    if(rv<0) printf ("thread creation error");
  }

  for(jac = 0;jac<NUM_THREADS;jac++){
    pthread_join(threads[jac], NULL);
  }

  //printMatrix(*kapMatrix);

  //return 23.4;
}

void *matrix_gentrification(void *ms){
  matrix_package *m_p = (matrix_package *)ms;
  double *A, *b, *t, *t1;
  double sum = 0; //sum through matrix addition
  int rc; //matrix index
  int entry_row; //first row vector.
  int leave_row; //end of for loop
  int dimension = m_p->currentMatrix->dimension;
  A = m_p->currentMatrix->A;
  b = m_p->currentMatrix->b;
  t = m_p->currentMatrix->t;
  t1 = m_p->currentMatrix->t1;
  
    //Here, we shall break up each thread to do its own
    //share of work.  
    //(dimension * threadID)/nofthreads
  entry_row = (dimension*m_p->thread_id) / NUM_THREADS;
  leave_row = entry_row + (dimension/NUM_THREADS);
#ifdef DEBUG
  printf("entry row: %d leave_row: %d dimension %d nt %d\n", entry_row,
	 leave_row, dimension, NUM_THREADS);
#endif
  //t1[0] = 1;

  for(;entry_row<leave_row;entry_row++){
    for(rc = 0;rc<dimension;rc++){
      sum += A[entry_row*dimension+rc]*b[rc];
    } //end for(rc)
    t1[entry_row] = sum;
#ifdef DEBUG
    printf("SUM = %f\n", sum);
#endif
    sum = 0;
  } //end for(entry_row)

}

void printMatrix(matrix x){
  int dimension = x.dimension;
  int j = 0;
    for(j = 0; j<dimension; j++){
      printf("%d:%f\n", j, x.t1[j]);
    }
}
