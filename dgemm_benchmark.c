#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "cblas.h"

#define FLOAT double

#ifndef MAX
#define MAX(a,b)   (a<b? b:a)
#endif

#define COMPSIZE 1
#define TOUPPER(a) {if ((a) > 0x60) (a) -= 0x20;}

int main(int argc, char *argv[]){

  FLOAT *a, *b, *c;
  FLOAT alpha[] = {1.0, 0.0};
  FLOAT beta [] = {0.0, 0.0};
  char transa = 'N';
  char transb = 'N';
  blasint m, n, k, i, j, lda, ldb, ldc;
  int loops = 1;
  int has_param_m = 0;
  int has_param_n = 0;
  int has_param_k = 0;
  char *p;

  int from =   1;
  int to   = 200;
  int step =   1;

  struct timeval start, stop;
  double time1, timeg;

  argc--;argv++;

  if (argc > 0) { from = atol(*argv);            argc--; argv++; }
  if (argc > 0) { to   = MAX(atol(*argv), from); argc--; argv++; }
  if (argc > 0) { step = atol(*argv);            argc--; argv++; }

  //if ((p = getenv("OPENBLAS_TRANS"))) {
  //  transa=*p;
  //  transb=*p;
  //}
  //if ((p = getenv("OPENBLAS_TRANSA"))) {
  //  transa=*p;
  //}
  //if ((p = getenv("OPENBLAS_TRANSB"))) {
  //  transb=*p;
  //}
  TOUPPER(transa);
  TOUPPER(transb);

  fprintf(stderr, "From : %3d  To : %3d Step=%d : Transa=%c : Transb=%c\n", from, to, step, transa, transb);

  //p = getenv("OPENBLAS_LOOPS");
  //if ( p != NULL ) {
  //  loops = atoi(p);
  //}

  //if ((p = getenv("OPENBLAS_PARAM_M"))) {
  //  m = atoi(p);
  //  has_param_m=1;
  //} else {
    m = to;
  //}
  //if ((p = getenv("OPENBLAS_PARAM_N"))) {
  //  n = atoi(p);
  //  has_param_n=1;
  //} else {
    n = to;
  //}
  //if ((p = getenv("OPENBLAS_PARAM_K"))) {
  //  k = atoi(p);
  //  has_param_k=1;
  //} else {
    k = to;
  //}

  if (( a = (FLOAT *)malloc(sizeof(FLOAT) * m * k * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  if (( b = (FLOAT *)malloc(sizeof(FLOAT) * k * n * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  if (( c = (FLOAT *)malloc(sizeof(FLOAT) * m * n * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

#ifdef linux
  srandom(getpid());
#endif

  for (i = 0; i < m * k * COMPSIZE; i++) {
    a[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
  }
  for (i = 0; i < k * n * COMPSIZE; i++) {
    b[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
  }
  for (i = 0; i < m * n * COMPSIZE; i++) {
    c[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
  }
 
  fprintf(stderr, "          SIZE                   Flops             Time\n");

  for (i = from; i <= to; i += step) {
    
    timeg=0;

    if (!has_param_m) { m = i; }
    if (!has_param_n) { n = i; }
    if (!has_param_k) { k = i; }

    if (transa == 'N') { lda = m; }
    else { lda = k; }
    if (transb == 'N') { ldb = k; }
    else { ldb = n; }
    ldc = m;

    fprintf(stderr, " M=%4d, N=%4d, K=%4d : ", (int)m, (int)n, (int)k);
    gettimeofday( &start, (struct timezone *)0);

    for (j=0; j<loops; j++) {
      cblas_dgemm (CblasColMajor, CblasTrans, CblasTrans, m, n, k, 1.0, a, lda, b, ldb, 0.0, c, ldc);
    }

    gettimeofday( &stop, (struct timezone *)0);
    time1 = (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_usec - start.tv_usec)) * 1.e-6;

    timeg = time1/loops;
    fprintf(stderr,
	    " %10.2f MFlops %10.6f sec\n",
	    COMPSIZE * COMPSIZE * 2. * (double)k * (double)m * (double)n / timeg * 1.e-6, time1);
    
  }

  return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
