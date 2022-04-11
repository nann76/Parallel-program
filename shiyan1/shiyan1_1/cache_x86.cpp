#include <iostream>
#include<windows.h>
#include<random>

using namespace std;


/*const int n = 1000;

static  float** matrix = new float* [n];
static float* vec = new float[n];
static float* sum = new float[n];*/


void trivial(int n,float* sum, float** matrix, float* vec) {

    for (int i = 0; i < n; i++) {
        sum[i] = 0.0;
        for (int j = 0; j < n; j++) {

            sum[i] += matrix[j][i] * vec[j];

        }
    }

};

void optimize(int n,float* sum, float** matrix, float* vec) {

    for (int i = 0; i < n; i++) {

        sum[i] = 0.0;
    }

    for (int i = 0; i < n; i++) {

        for (int j = 0; j < n; j++) {

            sum[j] += matrix[i][j] * vec[i];
        }
    }

};




int main()
{
    int n = 2000;
    float** matrix = new float* [n];
    float* vec = new float[n];
    float* sum = new float[n];

    for (int i = 0; i < n; i++) {
        matrix[i] = new float[n];
        vec[i] = rand() % 100;
    }

    for (int i = 0; i < n; i++) {

        for (int j = 0; j < n; j++) {
            matrix[i][j] = rand() % 100;

        }

    }





    //trivial algorithm

    long long head1, tail1, freq1;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq1);
    QueryPerformanceCounter((LARGE_INTEGER*)&head1);

  for (int i = 0; i < 100; i++) 
    trivial(n, sum, matrix, vec);


    
    QueryPerformanceCounter((LARGE_INTEGER*)&tail1); // end time
    cout << "trivial algorithm: " << ((tail1 - head1) * 1000.0 / freq1) /100 << "ms" << endl;


    // cache optimize algorithm

    long long head2, tail2, freq2;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq2);
    QueryPerformanceCounter((LARGE_INTEGER*)&head2);

    for (int i = 0; i < 100; i++)
    optimize(n,sum, matrix, vec);



    

    QueryPerformanceCounter((LARGE_INTEGER*)&tail2);
    cout << "cache algorithm: " << ((tail2 - head2) * 1000.0 / freq2)/100  << "ms" << endl;

    int m;
    cin >> m;
}
