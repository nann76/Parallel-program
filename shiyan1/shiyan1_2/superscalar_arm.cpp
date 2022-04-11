#include <iostream>
#include <sys/time.h>
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



    long long int n = pow(2, 5);

    long long  int* a = new long long int[n];
    long long int* b = new long long int[n];

    for (long long int i = 0; i < n; i++) {

        a[i] = rand() % 100;
        b[i] = a[i];

    }



    //Trivial algorithm

    timeval start, end;
    gettimeofday(&start, NULL);

     //for(int i=0;i<10000;i++)
    tri(n, a);


    gettimeofday(&end, NULL);
    cout << "trivial algorithm: " << 1000000* (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)<< "us" << endl;




    int sum1;
    int sum2;

    //多链路式

    timeval start1, end1;
    gettimeofday(&start1, NULL);

    //for (int i = 0; i < 10000; i++)
    multilink(n, a);


    gettimeofday(&end1, NULL);
    cout << "多链路式: " << 1000000 * (end1.tv_sec - start1.tv_sec) + (end1.tv_usec - start1.tv_usec)  << "us" << endl;




    //递归函数1
    timeval start2, end2;
    gettimeofday(&start2, NULL);

    
    recursion(n, a);




    gettimeofday(&end2, NULL);
    cout << "递归函数1: " << 1000000 * (end2.tv_sec - start2.tv_sec) + (end2.tv_usec - start2.tv_usec)  << "us" << endl;



    //函数2
    timeval start3, end3;
    gettimeofday(&start3, NULL);

    
    function2(n, b);


    gettimeofday(&end3, NULL);
    cout << "函数2: " << 1000000 * (end3.tv_sec - start3.tv_sec) + (end3.tv_usec - start3.tv_usec)  << "us" << endl;



}
