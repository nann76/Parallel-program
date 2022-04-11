#include <iostream>
#include <sys/time.h>
#include<random>

using namespace std;


const int n = 1000;

static  float** matrix = new float*[n];
static float* vec = new float[n];
static float*sum=new float[n];



int main()
{


 for (int i = 0; i < n; i++) {
        matrix[i] = new float[n];
        vec[i] = rand()%100;
    }

    for (int i = 0; i < n; i++) {

        for (int j = 0; j < n; j++) {
            matrix[i][j] = rand()%100;

        }

    }




 /*   for(int i=0;i<n;i++){
        vec[i]=rand()%100;
        for(int j=0;j<n;j++){
            matrix[i][j]=rand()%100;

        }
    }*/


        //trivial algorithm

timeval start, end;  
gettimeofday(&start, NULL);
    for(int k=0;k<100;k++){

      for(int i=0;i<n;i++){
        sum[i]=0.0;
         for(int j=0;j<n;j++){

        sum[i]+=matrix[j][i]*vec[j];

      }
      }

    }

gettimeofday(&end, NULL); 
    cout << "trivial algorithm: " << 1000*(end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000<< "ms" << endl;


    // cache optimize algorithm

timeval start2, end2;  
gettimeofday(&start2, NULL);

    for(int k=0;k<100;k++){

        for(int i=0;i<n;i++){

            sum[i]=0.0;
        }

        for(int i=0;i<n;i++){

            for(int j=0;j<n;j++){

                sum[j]+=matrix[i][j]*vec[i];
            }
        }



    }

gettimeofday(&end2, NULL); 
    cout << "cache algorithm: " <<1000*(end2.tv_sec - start2.tv_sec) + (end2.tv_usec - start2.tv_usec)/1000<< "ms" << endl;


}
