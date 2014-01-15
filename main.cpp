#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void dispVector (const VectorXd, const unsigned int, const unsigned int);

int main () {
  const unsigned int N = 8;
  const unsigned int vec_size = 2 * N * (N - 1) + (N - 1) * (N - 1);

  VectorXd rhs (vec_size);
  Ref<VectorXd> uVect = rhs.segment (0, N * (N - 1));
  Ref<VectorXd> vVect = rhs.segment (N * (N - 1), N * (N - 1));
  Ref<VectorXd> pVect = rhs.segment (2 * N * (N - 1), (N - 1) * (N - 1));

  uVect[3] = 5.0;

  dispVector (uVect, N - 1, N);
  cout << endl;
  dispVector (vVect, N, N - 1);
  cout << endl;
  dispVector (pVect, N - 1, N - 1);
  cout << endl;

  return 0;
}

void dispVector (const VectorXd v, const unsigned int N, const unsigned int M) {
  for (unsigned int i = 0; i < N; ++i) {
    cout << v.segment (i * M, M).transpose () << endl;
  }
}
