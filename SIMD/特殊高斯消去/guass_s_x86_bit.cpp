#include<iostream>
#include<fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <windows.h>
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX��AVX2

#include <bitset>


using namespace std;


const string stre1 = "��������1 ��������130��������Ԫ��22������Ԫ��8";
const string strr1 = "��������1 ��������130��������Ԫ��22������Ԫ��8";

const string stre2 = "��������2 ��������254��������Ԫ��106������Ԫ��53";
const string strr2 = "��������2 ��������254��������Ԫ��106������Ԫ��53";

const string stre3 = "��������3 ��������562��������Ԫ��170������Ԫ��53";
const string strr3 = "��������3 ��������562��������Ԫ��170������Ԫ��53";

const string stre4 = "��������4 ��������1011��������Ԫ��539������Ԫ��263";
const string strr4 = "��������4 ��������1011��������Ԫ��539������Ԫ��263";

const string stre5 = "��������5 ��������2362��������Ԫ��1226������Ԫ��453";
const string strr5 = "��������5 ��������2362��������Ԫ��1226������Ԫ��453";

const string stre6 = "��������6 ��������3799��������Ԫ��2759������Ԫ��1953";
const string strr6 = "��������6 ��������3799��������Ԫ��2759������Ԫ��1953";

const string stre7 = "��������7 ��������8399��������Ԫ��6375������Ԫ��4535";
const string strr7 = "��������7 ��������8399��������Ԫ��6375������Ԫ��4535";


const int col = 3799;//��
const int num_E = 1953;//����Ԫ��

const string str_e = "D:\\BaiduNetdiskDownload\\data\\"+stre6+"\\����Ԫ��.txt";
const string str_r = "D:\\BaiduNetdiskDownload\\data\\"+strr6+"\\����Ԫ��.txt";
//const int num_R = 22;//��Ԫ��


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


bool Is_NULL(bitset<col> *mm, int num) {

    for (int k = col - 1; k >= 0; k--) {

        if (mm[num][k] == 1)
            return 0;

    }
    return 1;
}

int get_lp(bitset<col> m, int col) {

    for (int p = col - 1; p >= 0; p--) {
        if (m[p] == 1) {

            return p;
        }

    }
    return 0;
}





void guass_s_serial(bitset<col> *E, bitset<col>* R) {

    bool E_ab[num_E];

    for (int j = 0; j < num_E; j++) {
        E_ab[j] = 1;
    }


    for (int i = 0; i < num_E; i++) {

        while (E_ab[i] == 1) {

            int lp_E = get_lp(E[i], col);
            if (Is_NULL(R, lp_E) == 0) {

              

                    E[i] = E[i] ^ R[lp_E];
                

                if (Is_NULL(E, i) == 1) {
                    E_ab[i] = 0;
                }

            }

            else {

                for (int j = 0; j < col; j++) {

                    R[lp_E] = E[i];
                    E_ab[i] = 0;
                }
            }
        }


    }

}


void guass_s_simd(bitset<col>* E, bitset<col>* R) {

    bool E_ab[num_E];

    for (int j = 0; j < num_E; j++) {
        E_ab[j] = 1;
    }


    for (int i = 0; i < num_E; i++) {

        while (E_ab[i] == 1) {

            int lp_E = get_lp(E[i], col);
            if (Is_NULL(R, lp_E) == 0) {



                E[i] = E[i] ^ R[lp_E];


                if (Is_NULL(E, i) == 1) {
                    E_ab[i] = 0;
                }

            }

            else {

                for (int j = 0; j < col; j++) {

                    R[lp_E] = E[i];
                    E_ab[i] = 0;
                }
            }
        }


    }

}






int main() {


 /*   int** E = new int* [col];
    int** R = new int* [col];
    int** E1 = new int* [col];
    int** R1 = new int* [col];
    int** E2 = new int* [col];
    int** R2 = new int* [col];

    for (int i = 0; i < col; i++) {
        E[i] = new int[col];
        E1[i] = new int[col];
        E2[i] = new int[col];
    }
    for (int i = 0; i < col; i++) {
        R[i] = new int[col];
        R1[i] = new int[col];
        R2[i] = new int[col];
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
    */

    bitset<col> *E;
    bitset<col> *R;

    E = new bitset<col>[num_E];
    R = new bitset<col>[col];
    ifstream infileE;
    ifstream infileR;
    infileE.open(str_e);
    infileR.open(str_r);
    //infileE.open("D:\\BaiduNetdiskDownload\\data\\��������9 ��������37960��������Ԫ��29304������Ԫ��14921\\����Ԫ��.txt");
    //infileR.open("D:\\BaiduNetdiskDownload\\data\\��������9 ��������37960��������Ԫ��29304������Ԫ��14921\\��Ԫ��.txt");


    


    string lineE;
    if (infileE) // �и��ļ�
    {
        int count = 0;
        while (getline(infileE, lineE)) // line�в�����ÿ�еĻ��з�
        {
            // cout << line << endl;
            vector<string> v;
            SplitString(lineE, v, " "); //�ɰ�����ַ����ָ�;
            for (vector<string>::size_type i = 0; i != v.size(); ++i) {
                //E[count][atoi(v[i].c_str())] = 1;
                
                E[count][atoi(v[i].c_str())] = 1;
            }
            count++;

        }
    }
    else // û�и��ļ�
    {
        cout << "no such fileE" << endl;
    }

    string lineR;
    if (infileR) // �и��ļ�
    {

        while (getline(infileR, lineR)) // line�в�����ÿ�еĻ��з�
        {
            // cout << line << endl;
            vector<string> v;
            SplitString(lineR, v, " "); //�ɰ�����ַ����ָ�;
            for (vector<string>::size_type i = 0; i != v.size(); ++i)
               // R[atoi(v[0].c_str())][atoi(v[i].c_str())] = 1;
            R[atoi(v[0].c_str())][atoi(v[i].c_str())] = 1;

        }
    }
    else // û�и��ļ�
    {
        cout << "no such fileR" << endl;
    }

  /*  for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
            R1[i][j] = R2[i][j] = R[i][j];
            E1[i][j] = E2[i][j] = E[i][j];
        }
    }*/

  //  guass_s_serial(E, R);

   /* for (int j = 0; j < num_E; j++) {
        for (int i = 0; i < col; i++)
    
            cout << E[j][i] << ",";
        cout << endl;
    }*/
   long long head, tail, freq;
   double t1, t2, t3, t4, t5, t6;



   QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    guass_s_serial(E, R);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);


    t1 = (double)(tail - head) * 1000 / freq;
    cout << "1�����㷨:" << t1 << "ms" << endl;
    cout << endl;
  /* QueryPerformanceCounter((LARGE_INTEGER*)&head);
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
    cout << "������:" << ((t1 - t2) / t1) * 100 << "%" << endl << endl;
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    guass_s_avx(E2, R2);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    t3 = (double)(tail - head) * 1000 / freq;
    cout << "3AVX:" << t3 << "ms" << endl;
    cout << "������:" << ((t1 - t3) / t1) * 100 << "%" << endl << endl;


    cout << "���ٱ�1/2  " << t1 / t2 << endl;
    cout << "���ٱ�1/3  " << t1 / t3 << endl;

    */

    // system("pause");

    //cout << sizeof(bitset<64>);
    int m;
    cin >> m;
}