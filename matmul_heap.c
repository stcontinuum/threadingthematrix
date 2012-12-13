//#define DEBUG
// matmul_heap.c
// CS4540 Fall 2010
// kapenga
//
// This program is one of a set of three programs that show a
// matrix itteration. The difference between the three is the way
// space for the arrays is allocated:
//   matmul_stack uses local (automatic) variables for the arrays,
//        which use space on the stack for the arrays. The array
//        sizes are dynamic.
//   matmul_static uses global (external) variables for the arrays,
//        which uses space on the stack for the arrays. The maximum
//        array sizes are compiled in.
//   matmul_heap uses dynamic (malloced) space for the arrays,
//        which uses space on the heap for the arrays. The array
//        sizes are dynamic.
//
// The following itteretion can be used to solve linear systems
//   t_{i+1} = A t_i + b
// If the itteration converges to t, then t == t_{i+1} == t_i
// So t = A t + b
//   or  (I-a) t = b
//   where, I is the n*n idenity matrix
// There are several important applied problems where convergence 
// will take place. One such case is when for
// each row of A ( rows 0 <= i < n)
//             sum(j=0 ... n-1) abs(a[i][j])  < 1.0    
// Then the itteration will converge, assuming no roundoff or overflow.
// Example
// % ./matmul_heap 4 10 5
//
//  a=
//  0.189331   0.147829  -0.009582   0.012830
// -0.020409   0.222627   0.073037   0.042701
//  0.069882   0.228326  -0.001161   0.024936
//  0.116375  -0.100117   0.229832   0.022235
//
//  b=
//  2.411774   9.837874   6.251698   6.576916
//
//  itt  error
//    0   2.878398e+00
//    1   8.266521e-01
//    2   2.688652e-01
//    3   8.817662e-02
//    4   2.832084e-02
//    5   9.015857e-03
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>
#include "mackiepheap.h"
// These two function are not ansi C so they do not appear from the
// libstd.h  header if the gcc option -std=c99 option is used.
// I should see if there is a safe way to include them from stdlib.h
// and not place them explicitly here, which is bad style.
void srand48(long int seedval);
double drand48(void);


int main(int argc, char *argv[]) {
  int	n=400;	// problenm size
double	*a;	// pointer to transfromation array, space to be malloced
double	*b;	// pointer to transfromation vector, space to be malloced
int	seed=10;// seed for srand48() / drand48()
double  *t;     // pointer to solution vector, space to be malloced 
double  *t1;    // pointer to next itteraton of solution vector, 
                // space to be malloced
//double	*ttemp;	// used to swap t1 and t at each itteration
 int	itt_max=4000;// number of itterations to preform
int	itt;	// current itteration 
int	i, j;   // indices into arrays
//double	sum;	// computes the inner products for A * t
double 	error;  // max | t1[i] - t[i] |
//double 	errori; // | t1[i] - t[i] |
char	ch;	// for error checking on command line args.
 matrix kapMatrix;
 pthread_t threads[NUM_THREADS];
 long elapsed;
 int jac = 0;

if( argc == 4 ) {
  if( (sscanf(argv[1],"%d %[^ /t]", &n, &ch) != 1) ||
      (sscanf(argv[2],"%d %[^ /t]", &seed, &ch) != 1) ||
      (sscanf(argv[3],"%d %[^ /t]", &itt_max, &ch) != 1) ) {
    fprintf(stderr," ERROR : useage: %s [ <n> <seed> <itt_max>]\n", argv[0]); 
    return(1);
  }
} else if(argc != 1 ) {
  fprintf(stderr," ERROR : useage: %s [ <n> <seed> <itt_max>]\n", argv[0]); 
  return(1);
} 
if( n<1 ) {
  fprintf(stderr," ERROR : n must be positive\n");
  return(1);
}
if( (a=(double *)malloc(sizeof(double)*n*n)) == NULL) {
  fprintf(stderr," ERROR : malloc for a failed\n");
  return(1);
}
if( (b=(double *)malloc(sizeof(double)*n)) == NULL) {
  fprintf(stderr," ERROR : malloc for b failed\n");
  return(1);
}
if( (kapMatrix.t=(double *)malloc(sizeof(double)*n)) == NULL) {
  fprintf(stderr," ERROR : malloc for t failed\n");
  return(1);
}
if( (kapMatrix.t1=(double *)malloc(sizeof(double)*n)) == NULL) {
  fprintf(stderr," ERROR : malloc for t1 failed\n");
  return(1);
}
 kapMatrix.A = (double*)malloc(sizeof(double*)*n*n);
 kapMatrix.b = (double*)malloc(sizeof(double*)*n);
 kapMatrix.A = a;
 kapMatrix.b = b;
 kapMatrix.t = t;
 kapMatrix.dimension = n;
// Generate matrix a with | eigenvalues | < 1
srand48((long int)seed);
printf("\n  a=\n");
for(i=0; i< n; i++) {
  for(j=0; j< n; j++) {
    *(a+n*i+j) = 1.999 * (drand48() - 0.5) / n;
    //printf("%10.6f ", *(a+n*i+j) );
  }
  //printf("\n");
}
printf("\n  b=\n");
// Generate vector b 
for(i=0; i< n; i++) {
  b[i] = 10.0 * drand48();
  //printf("%10.6f ", b[i]);
}
printf("\n");
// Initialize t
for(i=0; i< n; i++) {
  t[i] = b[i];
}
 printf("\n  itt  error\n");

 struct timeval tv1, tv2;
 double current;
 gettimeofday(&tv1, NULL);

 for(itt=0; itt<=itt_max; itt++) {
   error=0.0;
   matrix_multiplication(&kapMatrix, (pthread_t *)&threads);
   //printf("%5d %14.6e\n", itt, error); 
 }
 for(jac = 0;jac<NUM_THREADS;jac++){
   //free(threads[jac]);
 }
 free(a);
 free(b);
 free(kapMatrix.t1);
 //free(kapMatrix.t); not used
 
 gettimeofday(&tv2, NULL);
 //printf("%f=l\n", (long)((tv2.tv_sec-tv1.tv_sec)*1000000
		       //+tv2.tv_usec-tv1.tv_usec));
 //double k = 5332.325;
 //elapsed = (tv2.tv_sec-tv1.tv_sec)*1000000 + tv2.tv_usec-tv1.tv_usec;
 //printf("Time elapsed: %f Threads: %d", elapsed, NUM_THREADS);
return(0);
}
