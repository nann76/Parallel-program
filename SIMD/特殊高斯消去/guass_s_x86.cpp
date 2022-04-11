#include<iostream>
#include<fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <iomanip>
#include <windows.h>
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX、AVX2
using namespace std;

const int col = 3799;//列
const int num_E = 1953;//被消元子

const string stre = "D:\\BaiduNetdiskDownload\\data\\测试样例6 矩阵列数3799，非零消元子2759，被消元行1953\\被消元行.txt";
const string strr = "D:\\BaiduNetdiskDownload\\data\\测试样例6 矩阵列数3799，非零消元子2759，被消元行1953\\消元子.txt";


void SplitString(const string& s, vector<string>& v, const string& c)
{
    string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while (string::npos != pos2)
    {
        v.push_back(s.substr(pos1, pos2 - pos1));

        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if (pos1 != s.length())
        v.push_back(s.substr(pos1));
}


bool Is_NULL(int **mm,int num ) {

    for (int k = col-1; k >=0; k--) {
           
        if (mm[num][k] == 1)
            return 0;
    
    }
    return 1;
}

int get_lp(int*m,int col) {

    for (int p = col - 1; p >= 0; p--) {
        if (m[p] == 1) {
            
            return p;
        }
    
    }
    return 0;
}



bool Is_NULL1(float** mm, int num) {

    for (int k = col - 1; k >= 0; k--) {

        if (mm[num][k] == 1)
            return 0;

    }
    return 1;
}

int get_lp1(float* m, int col) {

    for (int p = col - 1; p >= 0; p--) {
        if (m[p] == 1) {

            return p;
        }

    }
    return 0;
}


void guass_s_serial(int** E, int** R) {

    bool E_ab[num_E] ;

    for (int j = 0; j < num_E;j++) {
        E_ab[j] = 1;
    }


    for (int i = 0; i < num_E; i++) {

        while (E_ab[i]==1) {

            int lp_E = get_lp(E[i], col);
            if (Is_NULL(R, lp_E)==0) {

                for (int j = 0; j < col; j++) {

                    E[i][j] = E[i][j] ^ R[lp_E][j];
                }

                if (Is_NULL(E, i) == 1) {
                    E_ab[i] = 0;
                }

            }

            else {

                for (int j = 0; j < col; j++) {

                    R[lp_E][j] = E[i][j];
                    E_ab[i] = 0;
                }
            }
        }


    }

}


void guass_s_sse(float** E1, float** R1) {

    bool E_ab[num_E];

    for (int j = 0; j < num_E; j++) {
        E_ab[j] = 1;
    }


    for (int i = 0; i < num_E; i++) {

        while (E_ab[i] == 1) {

            int lp_E = get_lp1(E1[i], col);
            if (Is_NULL1(R1, lp_E) == 0) {

                __m128 t1, t2, t3;
                int j = 0;
                for (j ; j+4 <= col; j+=4) {

                    t1 = _mm_loadu_ps(E1[i]+j);
                    t2 = _mm_loadu_ps(R1[lp_E] + j);
                    t3 = _mm_xor_ps(t1, t2);
                    _mm_storeu_ps(E1[i] + j, t3);

                   
                }

                for (j; j < col; j++) {
                    E1[i][j] = (int)E1[i][j] ^ (int)R1[lp_E][j];
                }
                if (Is_NULL1(E1, i) == 1) {
                    E_ab[i] = 0;
                }

            }

            else {

                for (int j = 0; j < col; j++) {

                    R1[lp_E][j] = E1[i][j];
                    E_ab[i] = 0;
                }
            }
        }


    }

}


void guass_s_avx(float** E1, float** R1) {

    bool E_ab[num_E];

    for (int j = 0; j < num_E; j++) {
        E_ab[j] = 1;
    }


    for (int i = 0; i < num_E; i++) {

        while (E_ab[i] == 1) {

            int lp_E = get_lp1(E1[i], col);
            if (Is_NULL1(R1, lp_E) == 0) {

                __m256 t1, t2, t3;
                int j = 0;
                for (j; j + 8 <= col; j += 8) {

                    t1 = _mm256_loadu_ps(E1[i] + j);
                    t2 = _mm256_loadu_ps(R1[lp_E] + j);
                    t3 = _mm256_xor_ps(t1, t2);
                    _mm256_storeu_ps(E1[i] + j, t3);


                }

                for (j; j < col; j++) {
                    E1[i][j] = (int)E1[i][j] ^ (int)R1[lp_E][j];
                }
                if (Is_NULL1(E1, i) == 1) {
                    E_ab[i] = 0;
                }

            }

            else {

                for (int j = 0; j < col; j++) {

                    R1[lp_E][j] = E1[i][j];
                    E_ab[i] = 0;
                }
            }
        }


    }

}






int main() {


    int** E = new int* [col];
    int** R = new int* [col];
    float** E1 = new float* [col];
    float** R1 = new float* [col];
    float** E2 = new float* [col];
    float** R2 = new float* [col];

    for (int i = 0; i < col; i++) {
        E[i] = new int[col];
        E1[i] = new float[col];
        E2[i] = new float[col];
    }
    for (int i = 0; i < col; i++) {
        R[i] = new int[col];
        R1[i] = new float[col];
        R2[i] = new float[col];
    }

    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
            E[i][j] = 0;
            E1[i][j] = 0;
            E2[i][j] = 0;
        }
    }
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
            R[i][j] = 0;
            R1[i][j] = 0;
            R2[i][j] = 0;
        }
    }




    ifstream infileE;
    ifstream infileR;

    infileE.open(stre);
    infileR.open(strr);





    string lineE;
    if (infileE) // 有该文件
    {
        int count = 0;
        while (getline(infileE, lineE)) // line中不包括每行的换行符
        {
            // cout << line << endl;
            vector<string> v;
            SplitString(lineE, v, " "); //可按多个字符来分隔;
            for (vector<string>::size_type i = 0; i != v.size(); ++i)
                E[count][atoi(v[i].c_str())] = 1;
            count++;

        }
    }
    else // 没有该文件
    {
        cout << "no such fileE" << endl;
    }

    string lineR;
    if (infileR) // 有该文件
    {

        while (getline(infileR, lineR)) // line中不包括每行的换行符
        {
            // cout << line << endl;
            vector<string> v;
            SplitString(lineR, v, " "); //可按多个字符来分隔;
            for (vector<string>::size_type i = 0; i != v.size(); ++i)
                R[atoi(v[0].c_str())][atoi(v[i].c_str())] = 1;


        }
    }
    else // 没有该文件
    {
        cout << "no such fileR" << endl;
    }

    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
            R1[i][j] = R2[i][j] = R[i][j];
            E1[i][j] = E2[i][j] = E[i][j];
        }
    }




    //  for (int i = 0;i++)
      //g    cout << E1[0][i] << ",";

    long long head, tail, freq;
    double t1, t2, t3, t4, t5, t6;



    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    guass_s_serial(E, R);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);

    for (int i = 0; i < col; i++) {
    delete[] E[i];
    delete[] R[i];
}
    delete[] E;
    delete[] R;


    t1 = (double)(tail - head) * 1000 / freq;
    cout << "1串行算法:" << t1 << "ms" << endl;
    cout << endl;
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    guass_s_sse(E1, R1);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    t2 = (double)(tail - head) * 1000 / freq;

    for (int i = 0; i < col; i++) {
        delete[] E1[i];
        delete[] R1[i];
    }
    delete[] E1;
    delete[] R1;
    cout << "2SSE:" << t2 << "ms" << endl;
    cout << "提升率:" << ((t1 - t2) / t1) * 100 << "%" << endl << endl;
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    guass_s_avx(E2, R2);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    t3 = (double)(tail - head) * 1000 / freq;
    cout << "3AVX:" << t3 << "ms" << endl;
    cout << "提升率:" << ((t1 - t3) / t1) * 100 << "%" << endl << endl;

    cout.setf(ios::fixed);
    cout << "加速比1/2  " <<fixed << setprecision(2)<< t1 / t2 << endl;
    cout << "加速比1/3  " <<fixed << setprecision(2) <<t1/ t3 << endl;


    int m;
 
  
    cin >> m;
}