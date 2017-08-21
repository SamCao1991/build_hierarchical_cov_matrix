#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "morton.h"

using namespace std;

int main()
{
    const int n = 4;
    const int N = n*n;
    uint32_t coordinates[N][2];
    // Column oriented
    for(int i = 0; i<n; i++)
        for(int j = 0; j<n; j++)
        {
            coordinates[j*n+i][0] = i;
            coordinates[j*n+i][1] = j;
        }
    uint64_t encodedVal[N];
    for(int i = 0; i<N; i++)
        encodedVal[i] = morton2D_64_encode(coordinates[i][0],coordinates[i][1]);
    vector<size_t> ind(N);
    iota(ind.begin(),ind.end(),0);
    sort(ind.begin(),ind.end(),[&encodedVal](size_t i1, size_t i2){return encodedVal[i1] < encodedVal[i2];});
    for(int i = 0; i<N; i++)
        cout << coordinates[ind[i]][0] << " " << coordinates[ind[i]][1] << endl;
    return 0;
}
