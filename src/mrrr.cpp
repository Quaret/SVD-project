// #include "Eigen/src/Core/Matrix.h"
// #include "Eigen/src/Core/Matrix.h"
#include "mrrr_n.h"
// #include "dequdeas.h"
// #include "Eigen/src/Core/Matrix.h"
// #include "Eigen/src/Core/util/Constants.h"
// #include <cstddef>
// #include <iterator>
#include "random"
using namespace std;
using namespace Eigen;
int main()
{
    int n = 4;
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<double> distribution(0, 1);
    MatrixXd A (n,n);
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         A(i,j) = distribution(generator);
    //     }
    // }
    A(0,0) = 1; A(0,1) = 1; A(0,2) = 1; A(0,3) = 1;
    A(1,0) = 0; A(1,1) = 0; A(1,2) = 0; A(1,3) = 1;
    A(2,0) = 1; A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
    A(3,0) = 1; A(3,1) = 1; A(3,2) = 0; A(3,3) = 1;
    auto bid = Eigen::internal::UpperBidiagonalization(A); // Бидиагонализируем входную матрицу
    // auto bL = bid.householderU();
    // auto bR = bid.householderV();
    // MatrixXd bd = bid.bidiagonal();
    // cout << bid.householder() << endl;
    cout << "A: " << endl << A << endl;
    // cout << "A restored " << endl << bL*bd*bR.transpose() << endl;
    MRRR_SVD<int> Svd (A);
    // DQDS_SVD<int> Svdd(A);
    auto U = Svd.matrixU();
    auto V = Svd.matrixV();
    auto S = Svd.singularValues();
    cout << "U: " << endl << U << endl;
    cout << "V: " << endl << V << endl;
    cout << "S: " << endl << S << endl;
    auto B = U*S*V.transpose();
    cout << "B: " << endl << B << endl;
    // auto BB = U*S.diagonal(0)*V.transpose();
    // cout << "S.diagonal()" << endl << S.diagonal(0) << endl;
    // cout << "BB: " << endl << BB << endl;
}