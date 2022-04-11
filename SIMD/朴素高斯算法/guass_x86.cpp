#include<iostream>
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX、AVX2
#include<windows.h>
using namespace std;




const int n = 2000;

float m[n][n];

 void m_reset()
 {
  for (int i = 0; i < n; i++){
		 for (int j = 0; j < i; j++)
	     m[i][j] = 0;
		 m[i][i] = 1.0;
		 for (int j = i + 1; j < n; j++)
		 m[i][j] = rand();
		 }
  for (int k = 0; k < n; k++)
		 for (int i = k + 1; i < n; i++)
		 for (int j = 0; j < n; j++)
		 m[i][j] += m[k][j];
	 }



 void normal(float **m1) {
	 
	 for (int k = 0; k < n; k++) {
	 
		 for (int j = k+1; j < n; j++) {
			 m1[k][j] = m1[k][j] / m1[k][k];
		 
		 }
		 m1[k][k] = 1;
	 
		 for (int i = k + 1; i < n; i++) {
		 
			 for (int j = k + 1; j < n; j++) {
			 
				 m1[i][j] = m1[i][j] - m1[i][k] * m1[k][j];
			 
			 }
			 m1[i][k] = 0;
		 }


	 }
 
 }



 void sse_lu_part1(float **m1) {
 
	 __m128 t1, t2, t3;
	 for (int k = 0; k < n; k++) {
		 t1 = _mm_set1_ps(m1[k][k]);
		 int j = k + 1;
		 for (j; j + 4 <= n; j += 4) {
			 t2 = _mm_loadu_ps(m1[k] + j);
			 t3 = _mm_div_ps(t2, t1);

			 _mm_storeu_ps(m1[k] + j, t3);
		 }
		 for (; j < n; j++) {
			 m1[k][j] = m1[k][j] / m1[k][k];
		 }
		 m1[k][k] = 1;

		 for (int i = k + 1; i < n; i++) {

			 for (int j = k + 1; j < n; j++) {

				 m1[i][j] = m1[i][j] - m1[i][k] * m1[k][j];

			 }
			 m1[i][k] = 0;
		 }


	 }
 
 }



 void sse_lu_part2(float** m1) {

	 for (int k = 0; k < n; k++) {

		 for (int j = k + 1; j < n; j++) {
			 m1[k][j] = m1[k][j] / m1[k][k];

		 }
		 m1[k][k] = 1;
		 __m128 vaik, vakj, vaij, vx;
		 for (int i = k + 1; i < n; i++) {
			 vaik = _mm_set1_ps(m1[i][k]);
			 int j = k + 1;
			 for (j; j + 4 <= n; j += 4) {
				 vakj = _mm_loadu_ps(m1[k] + j);
				 vaij = _mm_loadu_ps(m1[i] + j);
				 vx = _mm_mul_ps(vakj, vaik);
				 vaij = _mm_sub_ps(vaij, vx);
				 _mm_storeu_ps(m1[i] + j, vaij);

			 }
			 for (j; j < n; j++) {
				 m1[i][j] = m1[i][j] - m1[i][k] * m1[k][j];
			 }

			 m1[i][k] = 0;


		 }

	 }
 }

 void sse_lu_both(float **b) {
	 __m128 t1, t2, t3;
	 for (int k = 0; k < n; k++) {
		 t1 = _mm_set1_ps(b[k][k]);
		 // float tmp[4] = { b[k][k], b[k][k], b[k][k], b[k][k] };
		 // t1 = _mm_loadu_ps(tmp);
		 int j = k + 1;
		 for (j; j + 4 <= n; j += 4) {
			 t2 = _mm_loadu_ps(b[k] + j);
			 t3 = _mm_div_ps(t2, t1);

			 _mm_storeu_ps(b[k] + j, t3);
		 }
		 for (; j < n; j++) {
			 b[k][j] = b[k][j] / b[k][k];
		 }
		 b[k][k] = 1;

		 __m128 vaik, vakj, vaij, vx;
		 for (int i = k + 1; i < n; i++) {
			 vaik = _mm_set1_ps(b[i][k]);
			 int j = k + 1;
			 for (j; j + 4 <= n; j += 4) {
				 vakj = _mm_loadu_ps(b[k] + j);
				 vaij = _mm_loadu_ps(b[i] + j);
				 vx = _mm_mul_ps(vakj, vaik);
				 vaij = _mm_sub_ps(vaij, vx);
				 _mm_storeu_ps(b[i] + j, vaij);

			 }
			 for (j; j < n; j++) {
				 b[i][j] = b[i][j] - b[i][k] * b[k][j];
			 }

			 b[i][k] = 0;
		 }

	 


	 }
 
 

 
 }


 void sse_lu_both_align(float** b) {
	

	 __m128 t1, t2, t3, t4;
	 int mod;
	 float temp;
	 for (int k = 0; k < n; k++) {
		 t1 = _mm_set1_ps(b[k][k]);
		 mod = 4 - (k & 3);
		 temp =b[k][k];
		 for (int j = k; j < k + mod; j++)
			 b[k][j] /= temp;
		 for (int j = k + mod; j < n; j += 4) {
			 t2 = _mm_load_ps(b[k] + j);
			 t2 = _mm_div_ps(t2, t1);
			 _mm_store_ps(b[k] + j, t2);
		 }
		 b[k][k] = 1;

		 for (int i = k + 1; i < n; i++) {
			 t1 = _mm_set1_ps(b[i][k]);
			 temp = b[i][k];
			 for (int j = k + 1; j < k + mod; j++) {
				 b[i][j] -= temp * b[k][j];
			 }
			 for (int j = k + mod; j < n; j += 4) {
				 t2 = _mm_load_ps(b[i] + j);
				 t3 = _mm_load_ps(b[k] + j);
				 t4 = _mm_sub_ps(t2, _mm_mul_ps(t1, t3));
				 _mm_store_ps(b[i] + j, t4);
			 }
			 b[i][k] = 0;
		 }
	 }

	 }



 void avx_both(float** matrix) {
	 __m256 t1, t2, t3, t4;
	 for (int k = 0; k < n; k++) {
		 //PART ONE
		 float tmp[8] = { matrix[k][k], matrix[k][k], matrix[k][k], matrix[k][k],
		 matrix[k][k], matrix[k][k], matrix[k][k], matrix[k][k] };
		 t1 = _mm256_loadu_ps(tmp);
		 for (int j = n - 8; j >= k; j -= 8) {
			 t2 = _mm256_loadu_ps(matrix[k] + j);
			 t3 = _mm256_div_ps(t2, t1);
			 _mm256_storeu_ps(matrix[k] + j, t3);
		 }

		 if (k & 7 != (n & 7)) {
			 for (int j = k; (j & 7) != (n & 7); j++) {
				 matrix[k][j] /= tmp[0];
			 }
		 }

		 for (int j = (n & 7) - 1; j >= 0; j--) {
			 matrix[k][j] /= tmp[0];
		 }
		 //PART TWO
		 for (int i = k + 1; i < n; i++) {
			 float tmp[8] = { matrix[i][k], matrix[i][k], matrix[i][k], matrix[i][k],
			 matrix[i][k], matrix[i][k], matrix[i][k], matrix[i][k] };
			 t1 = _mm256_loadu_ps(tmp);
			 for (int j = n - 8; j > k; j -= 8) {
				 t2 = _mm256_loadu_ps(matrix[i] + j);
				 t3 = _mm256_loadu_ps(matrix[k] + j);
				 t4 = _mm256_sub_ps(t2, _mm256_mul_ps(t1, t3));
				 _mm256_storeu_ps(matrix[i] + j, t4);
			 }

			 for (int j = k + 1; (j & 7) != (n & 7); j++) {
				 matrix[i][j] -= matrix[i][k] * matrix[k][j];
			 }
			 matrix[i][k] = 0;
		 }
	 }
 }

 int main() {
 

	 float** m1 = new float* [n];
	 float** m2 = new float* [n];
	 float** m3 = new float* [n];
	 float** m4 = new float* [n];
	 float** m5 = new float* [n];
	 float** m6 = new float* [n];
	 for (int i = 0; i < n; i++) {
		 m1[i] = new float[n];
		 m2[i] = new float[n];
		 m3[i] = new float[n];
		 m4[i] = new float[n];
		 m5[i] = new float[n];
		 m6[i] = new float[n];
	 }
	 m_reset();

	 
	 for (int i = 0; i < n; i++) {
		 for (int j = 0; j < n; j++) {
			 m1[i][j] = m2[i][j] = m3[i][j] = m4[i][j] = m5[i][j] = m[i][j];
			 
		 }
	 }
 

	
	 long long head, tail, freq;
	 double t1, t2, t3, t4, t5, t6;

	

	 QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	 QueryPerformanceCounter((LARGE_INTEGER*)&head);
	 normal(m1);
	 QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	 t1 = (double)(tail - head) * 1000 / freq;
	 cout << "1串行算法:" << t1 << "ms" << endl;
	 cout << endl;
	 QueryPerformanceCounter((LARGE_INTEGER*)&head);
	 sse_lu_part1(m2);
	 QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	 t2 = (double)(tail - head) * 1000 / freq;
	 cout << "2采用SSE并行只优化第一部分的串行算法:" << t2 << "ms" << endl;
	 cout << "采用SSE并行只优化第一部分对比无优化性能提升率:" << (t1 - t2) / t1 * 100 << "%" << endl << endl;

	 QueryPerformanceCounter((LARGE_INTEGER*)&head);
	 sse_lu_part2(m3);
	 QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	 t3 = (double)(tail - head) * 1000 / freq;
	 cout << "3采用SSE并行只优化第二部分的串行算法:" << t3 << "ms" << endl;
	 cout << "采用SSE并行只优化第二部分对比无优化性能提升率:" << (t1 - t3) / t1 * 100 << "%" << endl << endl;

	 QueryPerformanceCounter((LARGE_INTEGER*)&head);
	 sse_lu_both(m4);
	 QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	 t4 = (double)(tail - head) * 1000 / freq;
	 cout << "4采用SSE并行同时优化两部分的算法:" << t4 << "ms" << endl;
	 cout << "采用SSE并行同时优化两部分对比无优化性能提升率:" << (t1 - t4) / t1 * 100 << "%" << endl << endl;

	 QueryPerformanceCounter((LARGE_INTEGER*)&head);
	 avx_both(m5);
	 QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	 t5 = (double)(tail - head) * 1000 / freq;
	 cout << "5采用AVX256并行同时优化两部分的算法:" << t5 << "ms" << endl;
	 cout << "采用AVX256并行同时优化两部分对比无优化性能提升率:" << (t1 - t5) / t1 * 100 << "%" << endl;
	 cout << "采用AVX256并行同时优化两部分对比采用SSE性能提升率:" << (t4 - t5) / t4 * 100 << "%" << endl << endl;

	 QueryPerformanceCounter((LARGE_INTEGER*)&head);
	sse_lu_both_align(m6);
	 QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	 t6 = (double)(tail - head) * 1000 / freq;
	 cout << "6采用SSE并行采用对齐策略同时优化两部分的算法:" << t6 << "ms" << endl;
	 cout << "采用SSE并行采用对齐策略同时优化对比无优化性能提升率:" << (t1 - t6) / t1 * 100 << "%" << endl;
	 cout << "采用SSE并行采用对齐策略同时优化对比采用SSE不对齐策略性能提升率:" << (t4 - t6) / t4 * 100 << "%" << endl << endl;

	



	/*for (int i = 0; i < n; i++) {
		 for (int j = 0; j < n; j++) {
			 cout << m[i][j] << " ";
		 }
		 cout << endl;
	 }

	 cout << endl;
	 for (int i = 0; i < n; i++) {
		 for (int j = 0; j < n; j++) {
			 cout << m1[i][j] << " ";
		 }
		 cout << endl;
	 }
	 for (int i = 0; i < n; i++) {
		 for (int j = 0; j < n; j++) {
			 cout << m5[i][j] << " ";
		 }
		 cout << endl;
	 }*/


	 printf("12加速比: %.2f\n", t1 / t2);
	 printf("13加速比: %.2f\n", t1 / t3);
	// printf("24加速比: %.2f\n", t2 / t4);
	// printf("34加速比: %.2f\n", t3 / t4);
	 printf("14加速比: %.2f\n", t1 / t4);
	 printf("15加速比: %.2f\n", t1 / t5);
	// printf("45加速比: %.2f\n", t4 / t5);
	 printf("16加速比: %.2f\n", t1 / t6);
	// printf("46加速比: %.2f\n\n\n\n", t4 / t6);

	// g++ - g - march = corei7 guass.cpp - O0 - o g0
	 int m;
	 cin >> m;
 }