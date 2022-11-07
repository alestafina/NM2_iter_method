#include "lab2.h"

int main() {
   matrix A;
   int count = 100000;
   double optimalW = 1.0;
   A.input_matrix();
   A.Jacobi(optimalW);
   //for (int i = 1; i <= 100; i++) {
   //   if (count > A.k) {
   //      count = A.k;
   //      optimalW = i / 100.0;
   //   }
   //   //cout << "w = " << i / 100.0 << " Ready! residual = " << A.residual << " by " << A.k << " iterations" << endl;
   //}
   A.output();
   cout << "Optimal w = " << optimalW << " by " << count << " iterations";
   return 0;
}