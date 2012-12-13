#include <pthread.h>

#define NUM_THREADS 4

typedef struct{
  double *A;   //An [n][n] matrix
  double *b;   //B 
  double *t;   //t
  double *t1;  //t1
  int dimension; //dimension of the rowspace of t.
}matrix;

/*Packages the previous matrix struct, and includes other */
/*metainformation. */
typedef struct {
  matrix *currentMatrix;
  int thread_max;  //maximum number of threads
  int thread_id;  //current thread number 
  double *return_storage_p;
}matrix_package;

double matrix_multiplication(matrix *system, pthread_t * );
