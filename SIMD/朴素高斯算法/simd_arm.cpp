#include<iostream>
#include <arm_neon.h>
#include <sys/time.h>
#include<random>
using namespace std;




const int n = 100;

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



 void neno_lu_part1(float **m1) {

	 float32x4_t t1, t2, t3;
	 for (int k = 0; k < n; k++) {
		 float temp[4] = { m1[k][k],m1[k][k],m1[k][k] ,m1[k][k] };
		 t1 = vld1q_f32(temp);

		 int j = k + 1;
		 for (j; j + 4 <= n; j += 4) {
			 t2 =  vld1q_f32(m1[k] + j);
			 t3 =  vdivq_f32(t2, t1);
			 
			 vst1q_f32(m1[k] + j, t3);
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



 void neno_lu_part2(float** m1) {

	 for (int k = 0; k < n; k++) {

		 for (int j = k + 1; j < n; j++) {
			 m1[k][j] = m1[k][j] / m1[k][k];

		 }
		 m1[k][k] = 1;
		 float32x4_t vaik, vakj, vaij, vx;
		 for (int i = k + 1; i < n; i++) {
			 float temp[4] = { m1[i][k],m1[i][k],m1[i][k] ,m1[i][k] };
			 vaik = vld1q_f32(temp);

			 //vaik =  vld1q_f32(m1[i][k]);
			 int j = k + 1;
			 for (j; j + 4 <= n; j += 4) {
				 vakj =  vld1q_f32(m1[k] + j);
				 vaij =  vld1q_f32(m1[i] + j);
				 vx =  vmulq_f32(vakj, vaik);
				 vaij =  vsubq_f32(vaij, vx);
				 vst1q_f32(m1[i] + j, vaij);

			 }
			 for (j; j < n; j++) {
				 m1[i][j] = m1[i][j] - m1[i][k] * m1[k][j];
			 }

			 m1[i][k] = 0;


		 }

	 }
 }

 void neno_lu_both(float **b) {
	 float32x4_t t1, t2, t3;
	 for (int k = 0; k < n; k++) {
		 float temp[4] = { b[k][k],b[k][k],b[k][k] ,b[k][k] };
		 t1 = vld1q_f32(temp);

		// t1 =  vld1q_f32(b[k][k]);
		 // float tmp[4] = { b[k][k], b[k][k], b[k][k], b[k][k] };
		 // t1 = _mm_loadu_ps(tmp);
		 int j = k + 1;
		 for (j; j + 4 <= n; j += 4) {
			 t2 =  vld1q_f32(b[k] + j);
			 t3 =  vdivq_f32(t2, t1);

			 vst1q_f32(b[k] + j, t3);
		 }
		 for (; j < n; j++) {
			 b[k][j] = b[k][j] / b[k][k];
		 }
		 b[k][k] = 1;

		 float32x4_t vaik, vakj, vaij, vx;
		 for (int i = k + 1; i < n; i++) {
			 float temp[4] = { b[i][k],b[i][k],b[i][k] ,b[i][k] };
			 vaik = vld1q_f32(temp);

			 //vaik =  vld1q_f32(b[i][k]);
			 int j = k + 1;
			 for (j; j + 4 <= n; j += 4) {
				 vakj =  vld1q_f32(b[k] + j);
				 vaij =  vld1q_f32(b[i] + j);
				 vx =  vmulq_f32(vakj, vaik);
				 vaij =  vsubq_f32(vaij, vx);
				 vst1q_f32(b[i] + j, vaij);

			 }
			 for (j; j < n; j++) {
				 b[i][j] = b[i][j] - b[i][k] * b[k][j];
			 }

			 b[i][k] = 0;
		 }




	 }




 }


 void neno_lu_both_align(float** b) {


	 float32x4_t t1, t2, t3, t4;
	 int mod;
	 float temp;
	 for (int k = 0; k < n; k++) {
		 float te[4] = { b[k][k],b[k][k],b[k][k] ,b[k][k] };
		 t1 = vld1q_f32(te);

		// t1 =  vld1q_f32(b[k][k]);
		 mod = 4 - (k & 3);
		 temp =b[k][k];
		 for (int j = k; j < k + mod; j++)
			 b[k][j] /= temp;
		 for (int j = k + mod; j < n; j += 4) {
			 t2 =  vld1q_f32(b[k] + j);
			 t2 =  vdivq_f32(t2, t1);
			 vst1q_f32(b[k] + j, t2);
		 }
		 b[k][k] = 1;

		 for (int i = k + 1; i < n; i++) {
			 float te[4] = { b[i][k],b[i][k],b[i][k] ,b[i][k] };
			 t1 = vld1q_f32(te);

			// t1 =  vld1q_f32(b[i][k]);
			 temp = b[i][k];
			 for (int j = k + 1; j < k + mod; j++) {
				 b[i][j] -= temp * b[k][j];
			 }
			 for (int j = k + mod; j < n; j += 4) {
				 t2 =  vld1q_f32(b[i] + j);
				 t3 =  vld1q_f32(b[k] + j);
				 t4 = vsubq_f32(t2,  vmulq_f32(t1, t3));
				 vst1q_f32(b[i] + j, t4);
			 }
			 b[i][k] = 0;
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




 






	  
	 double t1, t2, t3, t4, t5;

	 timeval start, end;
	 gettimeofday(&start, NULL);
	 normal(m1);
	 gettimeofday(&end, NULL);
  
	
	 t1 = 1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000;
	 cout << "1, normal:" << t1 << "ms" << endl;
	 cout << endl;

	 gettimeofday(&start, NULL);
	 neno_lu_part1(m2);
	 gettimeofday(&end, NULL);
     
	

	 t2 = 1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000;
	 cout << "2, neno_lu_part1:" << t2 << "ms" << endl;
     cout << "采用NENO并行只优化第一部分对比无优化性能提升率:" << (t1 - t2) / t1 * 100 << "%" << endl << endl;

	 gettimeofday(&start, NULL);
	 neno_lu_part2(m3);
	 gettimeofday(&end, NULL);

	 t3 = 1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000;
	 cout << "3, neno_lu_part2:" << t3 << "ms" << endl;
	 cout << "采用NENO并行只优化第二部分对比无优化性能提升率:" << (t1 - t3) / t1 * 100 << "%" << endl << endl;

	 gettimeofday(&start, NULL);
	 neno_lu_both(m4);
	 gettimeofday(&end, NULL);

	 t4 = 1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000 ;
	 cout << "4, neno_lu_both:" << t4 << "ms" << endl;
	 cout << "采用NENO并行同时优化两部分对比无优化性能提升率:" << (t1 - t4) / t1 * 100 << "%" << endl << endl;


	 gettimeofday(&start, NULL);
	neno_lu_both_align(m5);
	gettimeofday(&end, NULL);

	 t5 = 1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000;
	 cout << "5, neno_lu_both_align:" << t5 << "ms" << endl;
	 cout << "采用NENO并行采用对齐策略同时优化对比无优化性能提升率:" << (t1 - t5) / t1 * 100 << "%" << endl;
	 cout << "采用NENO并行采用对齐策略同时优化对比采用NENO不对齐策略性能提升率:" << (t4 - t5) / t4 * 100 << "%" << endl << endl;





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
	 printf("24加速比: %.2f\n", t2 / t4);
	 printf("34加速比: %.2f\n", t3 / t4);
	 printf("14加速比: %.2f\n", t1 / t4);
	 printf("15加速比: %.2f\n", t1 / t5);
	 printf("45加速比: %.2f\n\n\n\n", t4 / t5);
	 

	// g++ - g - march = native guass.cpp - O0 - o g0
	// g++ - g - march = native simd_arm.cpp - O0 - o g0
	 int m;
	 cin >> m;
 }
