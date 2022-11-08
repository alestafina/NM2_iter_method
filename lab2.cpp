#include "lab2.h"

int main() {
   matrix *A = new matrix;
   int count = 100000;
   double optimalW = 1.0;
   A->input_matrix();
   for (int i = 1; i < 200; i++) {
      A->Seidel(i/100.0);
      if (count > A->k) {
         count = A->k;
         optimalW = i / 100.0;
      }
      A->output();
      cout << "w = " << i / 100.0 << " Ready! residual = " << A->residual << " by " << A->k << " iterations" << endl;
   }
   cout << "Optimal w = " << optimalW << " by " << count << " iterations";
   return 0;
}