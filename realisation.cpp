#include "lab2.h"

void matrix::input_matrix() {
   double tmp = 0;

   ifstream fin("input.txt");

   fin >> N >> M >> eps >> maxiter;

   vector<int> l = { N - M - 3, N - M - 2, N - 1, N,
                    N - 1, N - M - 2, N - M - 3 };

   for (int i = 0; i < size; i++) {
      vector<double> str;
      for (int j = 0; j < l[i]; j++) {
         fin >> tmp;
         str.push_back(tmp);
      }
      A.push_back(str);
   }
   A[4].insert(A[4].begin(), 0);
   for (int i = 0; i < M + 2; i++)
      A[5].insert(A[5].begin(), 0);
   for (int i = 0; i < M + 3; i++)
      A[6].insert(A[6].begin(), 0);

   for (int i = 0; i < N; i++) {
      fin >> tmp;
      b.push_back(tmp);
      norm += b[i] * b[i];
   }
   norm = sqrt(norm);
   fin.close();
}

void matrix::input_vector() {
   ifstream fin("input_v.txt");
   double tmp = 0;
   for (int i = 0; i < N; i++) {
      fin >> tmp;
      x_prev.push_back(tmp);
      x_next.push_back(x_prev[i]);
   }
   fin.close();
}

void matrix::iter_step(vector<double> &x0, vector<double> &x1, double w) {
   double sum = 0.0;
   int tmp = 0;
   vector<int> idx = { M + 3, M + 2, 1, 0, -1, -M - 2, -M - 3 };
     
   for (int i = 0; i < N; i++) {
      sum = 0.0;
      residual = 0.0;
      int j = 0;
      j = (i > 0 &&  N - idx[tmp] == i) ? ++tmp : tmp; 
      while (j < size && i + idx[j] >= 0) {
         sum += A[j][i] * x0[i + idx[j]];
         ++j;
      }
      x1[i] = x0[i] + w * (b[i] - sum) / A[3][i];
      residual += (b[i] - sum) * (b[i] - sum);
   }
   residual = sqrt(residual) / norm;
}

void matrix::Jacobi(double w) {
   input_vector();
   residual = 1.0;
   for (k = 0; k < maxiter && residual > eps; k++)
      iter_step(x_prev, x_next, w);
}

void matrix::Seidel(double w) {
   input_vector();
   residual = 10.0;
   for (k = 0; k < maxiter && residual > eps; k++)
      iter_step(x_next, x_next, w);

}

void matrix::output() {
   ofstream fout("output.txt");
   fout.precision(3);
   fout << w << ";" << endl;
   fout.precision(15);
   fout << k << ";" << endl;
   for (int i = 0; i < N; i++) {
      fout << x_next[i] << endl;
   }
  /* fout << x_next[0] << "; " << 1.0 - x_next[0] << ";" << endl;
   for (int i = 1; i < N - 1; i++)
      fout << x_next[i] << ";" << (i + 1) * 1.0 - x_next[i] << ";" << endl;
   fout << x_next[N - 1] << ";" << N * 1.0 - x_next[N - 1] << ";" << endl << residual << ";" << endl;*/

}
