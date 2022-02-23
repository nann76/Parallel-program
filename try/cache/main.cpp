#include <iostream>

using namespace std;

int main()
{


ofstream materix("Matrix.txt");
if(!materix.is_open()){
    cout<<"Can't open";

}
random_device rd;



  //trivial algorithm
    long long head1, tail1, freq1;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq1);
    QueryPerformanceCounter((LARGE_INTEGER*)&head1);
    for (int k = 0; k < 100; k++) {
        for (int i = 0; i < n; i++) {
            inner_production1[i] = 0.0;
        }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    inner_production1[i] += matrix[j][i] * vector[j];
                }
            }
        }
        QueryPerformanceCounter((LARGE_INTEGER*)&tail1);	// end time
        cout << "trivial algorithm: " << ((tail1 - head1) * 1000.0 / freq1 )/100<< "ms" << endl;

// cache optimize algorithm
    long long head2, tail2, freq2;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq2);
    QueryPerformanceCounter((LARGE_INTEGER*)&head2);
    for (int k = 0; k < 100; k++) {
        for (int i = 0; i < n; i++) {
            inner_production1[i] = 0.0;
        }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    inner_production2[j] += matrix[i][j] * vector[i];
                }
            }
        }
        QueryPerformanceCounter((LARGE_INTEGER*)&tail2);
        cout << "cache algorithm: " << ((tail2 - head2) * 1000.0 / freq2) / 100 << "ms" << endl;


}
