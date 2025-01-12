#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
  {
  
  mat A(4, 5, fill::randu);
  mat B(4, 5, fill::randu);
  cout<<A.n_cols<<endl;
  cout<<A.n_rows<<endl;
  cout<<A.n_elem<<endl;
  cout<<dot(A,B)<<endl;
  cout << A*B.t() << endl;
  
  vec C(3,arma::fill::zeros);
  C[0]=1;
  C[1]=2;
  C[2]=3;
  cout<<norm(C, 1)<<endl;
  return 0;
  }

