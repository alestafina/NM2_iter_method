#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include <locale>

using namespace std;

struct matrix {
   const double size = 7;                // Количество диагоналей
   int N = 0, M = 0, maxiter = 0, k = 0; // Размер матрицы, кол-во нулевых диагоналей,
                                         // максимальное кол-во итераций, счетчик итераций
   vector<vector<double>> A;             // Исходная матрица
   vector<double> x_prev, x_next, b;     // Вектор предыдущей и вектор последующей итераций,
                                         // вектор правой части
   double eps = 0.0;                     // Максимальная абсолютная погрешность
   double w = 0.0;                       // Параметр релаксации
   double residual = 0.0;                // Относительная невзяка
   double norm = 0.0;                    // Норма вектора правой части
   void input_matrix();                  // Ввод матрицы
   void input_vector();                  // Ввод вектора
   void iter_step(vector<double> &x0, vector<double> &x1, double w); // Итерационный шаг
   void Seidel(double omega);            // Метод Гаусса-Зейделя
   void Jacobi(double omega);            // Метод Якоби
   void output();                        // Вывод
};

