#pragma once
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;


struct matrix {
   const double size = 7;                // ���������� ����������
   int N = 0, M = 0, maxiter = 0, k = 0; // ������ �������, ���-�� ������� ����������,
                                         // ������������ ���-�� ��������, ������� ��������
   vector<vector<double>> A;             // �������� �������
   vector<double> x_prev, x_next, b;     // ������ ���������� � ������ ����������� ��������,
                                         // ������ ������ �����
   double eps = 0.0;                     // ������������ ���������� �����������
   double w = 0.0;                       // �������� ����������
   double residual = 0.0;                // ������������� �������
   double norm = 0.0;                    // ����� ������� ������ �����
   void input_matrix();                  // ���� �������
   void input_vector();                  // ���� �������
   void iter_step(vector<double> &x0, vector<double> &x1, double omega); // ������������ ���
   void Seidel(double omega);                // ����� ������-�������
   void Jacobi(double omega);                // ����� �����
   void output();                        // �����
};
