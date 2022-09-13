#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include "cblas.h"

int main() {
    float ret1 = 0.0;
    double ret2 = 0.0;
    blasint N=4,incx=1;
    double xd[]={1.0,3.0,5.0,7.0, 1.0, 1.0, 1.0, 1.0};
    float xs[]={1.0,3.0,5.0,7.0, 1.0, 1.0, 1.0, 1.0};
    //double xd[]={1.0,3.0,5.0,7.0,9.0,11.0,13.0};
    //float xs[]={1.0,3.0,5.0,7.0,9.0,11.0,13.0};

    ret1 = cblas_scnrm2(N, xs, incx);
    printf("cblas_scnrm2 ret1 = %.14f\n", ret1);
    ret2 = cblas_dznrm2(N, xd, incx);
    printf("cblas_dznrm2 ret2 = %.14f\n", ret2);
    return 0;
}
