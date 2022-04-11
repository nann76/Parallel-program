#include<iostream>
#include<fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <arm_neon.h>
#include <sys/time.h>
using namespace std;




const string stre1 = "1_130_22_8";
const string strr1 = "1_130_22_8";

const string stre2 = "2_254_106_53";
const string strr2 = "2_254_106_53";

const string stre3 = "3_562_170_53";
const string strr3 = "3_562_170_53";

const string stre4 = "4_1011_539_263";
const string strr4 = "4_1011_539_263";

const string stre5 = "5_2362_1226_453";
const string strr5 = "5_2362_1226_453";

const string stre6 = "6_3799_2759_1953";
const string strr6 = "6_3799_2759_1953";

const string stre7 = "7_8399_6375_4535";
const string strr7 = "7_8399_6375_4535";


const int col = 8399;//列
const int num_E = 4535;//被消元子

const string str_e = "//home//data//Groebner//"+stre7+"//2.txt";
const string str_r = "//home//data//Groebner//"+strr7+"//1.txt";


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


bool Is_NULL(unsigned int **mm,int num ) {

    for (int k = col-1; k >=0; k--) {
           
        if (mm[num][k] == 1)
            return 0;
    
    }
    return 1;
}

int get_lp(unsigned int*m,int col) {

    for (int p = col - 1; p >= 0; p--) {
        if (m[p] == 1) {
            
            return p;
        }
    
    }
    return 0;
}




void guass_s_serial(unsigned int** E, unsigned int** R) {

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


void guass_s_neon(unsigned int** E1,unsigned int** R1) {

    bool E_ab[num_E];

    for (int j = 0; j < num_E; j++) {
        E_ab[j] = 1;
    }


    for (int i = 0; i < num_E; i++) {

        while (E_ab[i] == 1) {

            int lp_E = get_lp(E1[i], col);
            if (Is_NULL(R1, lp_E) == 0) {

                uint32x4_t t1, t2, t3;
                int j = 0;
                for (j ; j+4 <= col; j+=4) {

                    t1 = vld1q_u32(E1[i]+j);
                    t2 = vld1q_u32(R1[lp_E] + j);
                    t3 = veorq_u32(t1, t2);
                    vst1q_u32(E1[i] + j, t3);

                   
                }

                for (j; j < col; j++) {
                    E1[i][j] = E1[i][j] ^ R1[lp_E][j];
                }
                if (Is_NULL(E1, i) == 1) {
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


    unsigned int** E = new unsigned int* [col];
    unsigned int** R = new unsigned int* [col];
    unsigned int** E1 = new unsigned int* [col];
    unsigned int** R1 = new unsigned int* [col];
   
    for (int i = 0; i < col; i++) {
        E[i] = new unsigned int[col];
        E1[i] = new unsigned int[col];
    }
    for (int i = 0; i < col; i++) {
        R[i] = new unsigned int[col];
        R1[i] = new unsigned int[col];
    }
    
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
        E[i][j] = 0;
        E1[i][j] = 0;
    }
    }
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
            R[i][j] = 0;
            R1[i][j] = 0;
        }
    }

   


	ifstream infileE;
    ifstream infileR;
	infileE.open(str_e);
    infileR.open(str_r);
    
    
    
    string lineE;
    if (infileE) // 有该文件
    {
        int count = 0;
        while (getline(infileE, lineE)) // line中不包括每行的换行符
        {
           // cout << line << endl;
            vector<string> v;
            SplitString(lineE, v, " "); //可按多个字符来分隔;
            for (vector<string>::size_type i = 0; i != v.size(); ++i){
           
             E[count][atoi(v[i].c_str())]=1;
	
}
            count++;
            
        }
    }
    else // 没有该文件
    {
        cout << "no such file" << endl;
    }

    string lineR;
    if (infileR) // 有该文件
    {
        
        while (getline(infileR, lineR)) // line中不包括每行的换行符
        {
            // cout << line << endl;
            vector<string> v;
            SplitString(lineR, v, " "); //可按多个字符来分隔;
            for (vector<string>::size_type i = 0; i != v.size(); ++i){
               
                R[atoi(v[0].c_str())][atoi(v[i].c_str())] = 1;
           }

        }
    }
    else // 没有该文件
    {
        cout << "no such file" << endl;
    }

    for (int i = 0; i < col; i++) {
        for (int j = 0; j < col; j++) {
	E[i][0]=0;
	R[i][0]=0;
            R1[i][j] = R[i][j];
            E1[i][j] = E[i][j];
        }
    }

    
 

   
    double t1, t2;



     timeval start, end;
     gettimeofday(&start, NULL);
    	
	guass_s_serial(E,R);


     gettimeofday(&end, NULL);

    t1 =1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000 ;
    cout << "1串行算法:" << t1 << "ms" << endl;
    cout << endl;


    
     


     gettimeofday(&start, NULL);

    guass_s_neon(E1, R1);
    gettimeofday(&end, NULL);
    t2 = 1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1000 ;
    cout << "2SSE:" << t2 << "ms" << endl;
    cout << "提升率:" << ((t1 - t2) / t1) * 100 << "%" << endl << endl;


    cout << "加速比1/2  " << t1 / t2 << endl;



  
   
}