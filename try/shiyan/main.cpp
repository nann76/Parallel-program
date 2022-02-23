#include<iostream>
#include<windows.h>
#include <stdlib.h>
#include<ctime>
using namespace std;



const int n = 1000;

static  int** matrix = new int*[n];
static int* vector = new int[n];
static int inner_production1[n];
static int inner_production2[n];





int main() {

    for (int i = 0; i < n; i++) {
        matrix[i] = new int[n];
        vector[i] = rand();
    }

    for (int i = 0; i < n; i++) {

        for (int j = 0; j < n; j++) {
            matrix[i][j] = rand();

        }

    }



    //trivial algorithm
    long long head1, tail1, freq1;

    QueryPerformanceFrequency((LARGE_INTEGER*)&freq1);
    QueryPerformanceCounter((LARGE_INTEGER*)&head1);

    for (int i = 0; i < n; i++) {
        inner_production1[i] = 0;
    }

    for (int k = 0; k < 100; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inner_production1[i] += matrix[j][i] * vector[j];
            }
        }
    }
    QueryPerformanceCounter((LARGE_INTEGER*)&tail1);	// end time
    cout << "trivial algorithm: " << ((tail1 - head1) * 1000.0 / freq1 )/100<< "ms" << endl;




    // cache algorithm

    long long head2, tail2, freq2;

    QueryPerformanceFrequency((LARGE_INTEGER*)&freq2);	// similar to CLOCKS_PER_SEC
    QueryPerformanceCounter((LARGE_INTEGER*)&head2);// start time

    for (int i = 0; i < n; i++) {
        inner_production1[i] = 0;
    }

    for (int k = 0; k < 100; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inner_production2[j] = matrix[i][j] * vector[i];
            }
        }
    }
    QueryPerformanceCounter((LARGE_INTEGER*)&tail2);	// end time
    cout << "cache algorithm: " << ((tail2 - head2) * 1000.0 / freq2) / 100 << "ms" << endl;



}






