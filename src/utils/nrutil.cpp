#include "nrutil.h"

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
    fprintf(stderr, "Numerical Recipes run-time error...\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now exiting to system...\n");
    exit(1);
}


double* darray(long n)
/* allocate a float vector with subscript range v[nl..nh] */
{
    double* da;

    da = new double[n];//(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));

    if (!da)
    {
        nrerror("allocation failure in vector()");
    }

    return da;
}

void free_darray(float* da)
/* free a double vector allocated with vector() */
{
    delete[] da;
}