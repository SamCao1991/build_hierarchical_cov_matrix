#include <algorithm>
#include <iostream>
#include <fstream>
#include <chrono>
#include "morton.h"
extern "C"
{
#include "myfunc.h"
}

typedef std::chrono::time_point<std::chrono::steady_clock> TimeStamp;

using namespace std;


int main()
{
    const int n = 256;
    const int N = n*n;
    const int bsz = 256;
    const double beta = 1.0;
    const double epsl = 1e-1;
    TimeStamp myStart, myEnd;
    uint_fast32_t coordinates_int[N][2];
    cout << "Begin building geometry" << endl;
    myStart = chrono::steady_clock::now();
    // Column oriented
    for(int i = 0; i<n; i++)
        for(int j = 0; j<n; j++)
        {
            coordinates_int[j*n+i][0] = i;
            coordinates_int[j*n+i][1] = j;
        }
    myEnd = chrono::steady_clock::now();
    cout << "Building geomtry took " << std::chrono::duration<double> (myEnd - myStart).count() << "s" << endl;
    uint_fast64_t encodedVal[N];
    for(int i = 0; i<N; i++)
        encodedVal[i] = morton2D_64_encode(coordinates_int[i][0],coordinates_int[i][1]);
    cout << "Begin sorting" << endl;
    myStart = chrono::steady_clock::now();
    sort(&encodedVal[0], &encodedVal[N]);
    myEnd = chrono::steady_clock::now();
    cout << "Sorting took " << std::chrono::duration<double> (myEnd - myStart).count() << "s" << endl;
    double coordinates_double[N][2];
    for(int i = 0; i<N; i++)
    {
        morton2D_64_decode(encodedVal[i],coordinates_int[i][0],coordinates_int[i][1]);
        coordinates_double[i][0] = (double) coordinates_int[i][0];
        coordinates_double[i][1] = (double) coordinates_int[i][1];
    }
    cout << "Begin building hmatrix" << endl;
    myStart = chrono::steady_clock::now();
    output_exp_hie_cov(coordinates_double,N,2,bsz,beta,epsl);
    myEnd = chrono::steady_clock::now();
    cout << "Building hmatrix took " << std::chrono::duration<double> (myEnd - myStart).count() << "s" << endl;

    ofstream myFile;
    myFile.open("geom.txt");
    for(int i = 0; i<N; i++)
        myFile << coordinates_double[i][0] << " " << coordinates_double[i][1] << endl;
    myFile.close();
    return 0;
}
