#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <iostream>

#include "morton.h"
extern "C"
{
#include "myfunc.h"
}

using namespace std;

int main()
{
    const int n = 16;
    const int N = n*n;
    uint_fast32_t coordinates_int[N][2];
    // Column oriented
    for(int i = 0; i<n; i++)
        for(int j = 0; j<n; j++)
        {
            coordinates_int[j*n+i][0] = i;
            coordinates_int[j*n+i][1] = j;
        }
    uint_fast64_t encodedVal[N];
    for(int i = 0; i<N; i++)
        encodedVal[i] = morton2D_64_encode(coordinates_int[i][0],coordinates_int[i][1]);
    sort(&encodedVal[0], &encodedVal[N]);
    double coordinates_double[N][2];
    for(int i = 0; i<N; i++)
    {
        morton2D_64_decode(encodedVal[i],coordinates_int[i][0],coordinates_int[i][1]);
        coordinates_double[i][0] = (double) coordinates_int[i][0];
        coordinates_double[i][1] = (double) coordinates_int[i][1];
    }
    output_exp_hie_cov(coordinates_double,N,2,32,1,1e-3);
    return 0;
}
