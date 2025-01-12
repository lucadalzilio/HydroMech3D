#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
  {
    wall_clock timer;
    int dim = 100;

    cx_mat C = randu<cx_mat>(dim,dim);
    cx_mat D = C.t()*C;

    vec eigval2;
    cx_mat eigvec2;

    timer.tic(); // Initialize clock

    eig_sym(eigval2, eigvec2, D);

    double n = timer.toc();

    cout << "Elapsed time: " << n << " seconds" << endl;
    cout << eigval2 << endl;

  return 0;
  }
