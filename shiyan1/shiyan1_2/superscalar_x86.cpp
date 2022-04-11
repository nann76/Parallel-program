#include <iostream>
#include<windows.h>
#include<random>
#include<cmath>
using namespace std;







    void tri(long long int n, long long int* a) {

    int sum = 0;
    for (long long int i = 0; i < n; i++) {
        sum += a[i];
    }
    
}


void multilink(long long int n, long long int* a) {
    int sum = 0;
    int sum1 = 0;
    int sum2 = 0;
    for (long long int i = 0; i < n; i += 2) {

        sum1 += a[i];
        sum2 += a[i + 1];
    }
    sum = sum1 + sum2;
    
}



void recursion(long long int n, long long int* a) {
    if (n == 1) {
      
        return;
    }
    else {

        for (long long int i = 0; i < n / 2; i++)
            a[i] += a[n - i - 1];
        n = n / 2;
        
        recursion(n, a);
    }

}


void function2(long long int n, long long int* b) {
    
    for (int m = n; m > 1; m /= 2) {
        for (int i = 0; i < m / 2; i++)
            b[i] = b[i * 2] + b[i * 2 + 1];

    }
    
}




int main()
{



    long long int n = pow(2, 22);

    long long  int* a = new long long int[n];
    long long int* b = new long long int[n];

    for (long long int i = 0; i < n; i++) {

        a[i] = rand() % 100;
        b[i] = a[i];

    }



    //Trivial algorithm

    long long head1, tail1, freq1;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq1);
    QueryPerformanceCounter((LARGE_INTEGER*)&head1);

     //for(int i=0;i<10000;i++)
    tri(n, a);


    QueryPerformanceCounter((LARGE_INTEGER*)&tail1); // end time
    cout << "trivial algorithm: " << ((tail1 - head1) * 1000.0 / freq1)<< "ms" << endl;




    int sum1;
    int sum2;

    //多链路式

    long long head, tail, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);

    //for (int i = 0; i < 10000; i++)
    multilink(n, a);


    QueryPerformanceCounter((LARGE_INTEGER*)&tail); // end time
    cout << "多链路式: " << ((tail - head) * 1000.0 / freq) << "ms" << endl;




    //递归函数1
    long long head3, tail3, freq3;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq3);
    QueryPerformanceCounter((LARGE_INTEGER*)&head3);

    
    recursion(n, a);




    QueryPerformanceCounter((LARGE_INTEGER*)&tail3); // end time
    cout << "递归函数1: " << ((tail3 - head3) * 1000.0 / freq3) << "ms" << endl;



    //函数2
    long long head4, tail4, freq4;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq4);
    QueryPerformanceCounter((LARGE_INTEGER*)&head4);

    
    function2(n, b);


    QueryPerformanceCounter((LARGE_INTEGER*)&tail4); // end time
    cout << "函数2: " << ((tail4 - head4) * 1000.0 / freq4)<< "ms" << endl;



}
