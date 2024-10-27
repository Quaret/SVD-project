// Telegram: @Quareten
// Кочержук Александр, КМБО-03-22
// #include "Eigen/src/Core/Matrix.h"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <cstdint>
#include <lapacke.h>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <ostream>
#include <tuple>
#include <vector>
double sq = sqrt(2);
template<typename T>
class MRRR_SVD
{
private:
    Eigen::MatrixXd U;
    Eigen::MatrixXd S;
    Eigen::MatrixXd V;
    void Set_U(const Eigen::MatrixXd& A)
    {
        U = A;
    }

    void Set_V(const Eigen::MatrixXd& A)
    {
        V = A;
    }

    void Set_S(const Eigen::MatrixXd& A)
    {
        S = A;
    }

protected:
    MRRR_SVD<T> compute_bsvd (const Eigen::MatrixXd& matrix)
    {
        auto bid = Eigen::internal::UpperBidiagonalization(matrix); // Бидиагонализируем входную матрицу
        auto B = bid.bidiagonal(); // Обозначаем её как B
        auto L = bid.householderU();
        auto R = bid.householderV();
        // Eigen::MatrixXd LL = L;
        // Eigen::MatrixXd RR = R;
        // std::cout << "L: " << std::endl << LL << std::endl;
        // std::cout << "R: " << std::endl << RR << std::endl;
        int n = B.rows();        
        Eigen::MatrixXd TGK(2*n,2*n); // Создаём матрицу Голуба-Кахана TGK(B)
        TGK.setZero();
        // this->m_matrixU = bid.householderU(); //Инициализурем всяким мусором для того, чтобы потом туда записывать уже нормальные вектора с помощью обычных операций
        // this->m_matrixV = bid.householderV().transpose(); //без выноса мозга с Eigen
        Eigen::MatrixXd BB = B; //Переводим матрицу B к виду, где можно использовать взятие по индексу
        std::cout << "B: " << std::endl << B.diagonal(0) << std::endl << B.diagonal(1) << std::endl;
        std::cout << "BB: " << std::endl << BB << std::endl;
        // std::cout << "A original" << LL*BB*RR << std::endl;
        int nn = 2*n;
        for (int i = 0; i < n - 1; ++i) //Считаем TGK(B)
        {
            TGK(2*i,2*i+1) = BB(i,i); TGK(2*i+1,2*i+2) = BB(i,i+1);
            TGK(2*i+1,2*i) = BB(i,i); TGK(2*i+2,2*i+1) = BB(i,i+1);
        }
        TGK(nn-2,nn-1) = BB(n-1,n-1);
        TGK(nn-1,nn-2) = BB(n-1,n-1);
        std::cout << "TGK: " << std::endl << TGK << std::endl;
        int32_t nzc = std::max(1,nn); //Количество собственных векторов, которые будут храниться в z

        double *d = new double[2*n]; //Диагональ TGK(B)
        for (int i = 0; i < 2*n; i++)
        {
            d[i] = TGK.diagonal(0)[i];
        }
        
        double *e = new double [2*n];
        for (int i = 0; i < 2*n-1; i++) //Поддиагональ TGK(B)
        {
            e[i] = TGK.diagonal(-1)[i];
        }
        double *w = new double[nn]; //Собственные числа в возрастающем (ВАЖНО) порядке
        double *z = new double[nn*nn]; //Собственные вектора
        int32_t m; //Число найденных собственных значений - здесь не важно
        int32_t isuppz [2*nn]; //Индексы ненулевых векторов z 
        int32_t tryrac = 1; //Логическая переменная для проверки на точность, см. документацию
        int32_t info = LAPACKE_dstemr 
        (
            LAPACK_COL_MAJOR, // int matrix_layout, in
            'V', // char jobz, in
            'A', // char range, in
            nn, // int32_t n, in
            d,// double *d, in out
            e, // double *e, in out
            0, // double vl, in
            0, // double vu, in
            0, // int32_t il, in
            0, // int32_t iu, in
            &m, // int32_t *m, out
            w, // double in
            z, // double *z, out
            nn, // int32_t ldz, in
            nzc, // int32_t nzc, in
            isuppz,// int32_t *isuppz, out
            &tryrac// int32_t *tryrac, in outs
        );
        if (info != 0)
            throw std::runtime_error("LAPACK error: " + std::to_string(info));
        
        Eigen::VectorXd eigenvalues (nn);
        for (int i = 0; i < nn; ++i) 
        {
            eigenvalues(i) = w[nn-i-1];
        }
        std::cout << "eigenvalues: " << std::endl << eigenvalues << std::endl;
        Eigen::MatrixXd eigenvectors (nn,nn);
        for (int i = 0; i < nn; i++)
        {
            for (int j = 0; j < nn; j++)
            {
                eigenvectors(i,j) = i < n ? z[(nn-i-1)*nn+j]*sq : z[(nn-i-1)*nn+j] * (-1)*sq;
            }
        }
        for (int i = n-1; i < n+1; i++)
        {
            for (int j = 0; j < nn; j++) 
            {
                eigenvectors(i,j) *= -1;
            }
        }
        std::cout << "eigenvectors: " << std::endl << eigenvectors << std::endl;
        // eigenvectors = eigenvectors.transpose();
        // std::cout << "eigenvectors: " << std::endl << eigenvectors << std::endl;
        Eigen::MatrixXd matU (n,n);
        Eigen::MatrixXd matV (n,n);
        Eigen::MatrixXd singularValuess (n,n);
        int ind = 0;
        for (int i = 0; i < n; ++i) 
        {
            ind = i;
            double sigma = std::abs(eigenvalues(ind));
            if (sigma < 1e-15) 
            { 
                sigma = 0; ind = -i - 1;
            } 
            Eigen::VectorXd q (nn);
            for (int j = 0; j < nn; j++)
            {
                q[j] = eigenvectors(ind,j);
            }
            std::cout << "q: " << std::endl << q << std::endl;
            Eigen::VectorXd u (n);
            Eigen::VectorXd v (n);
            for (int j = 0; j < n; ++j)
            {
                v[j] = q[2*j];
                u[j] = q[2*j+1];
            }
            singularValuess(i,i) = sigma;
            matU.col(i) = u;
            matV.col(i) = v;
        }
        std::cout << "matU: " << std::endl << matU << std::endl;
        std::cout << "matV: " << std::endl << matV << std::endl;
        Set_S(singularValuess);
        Set_U(L*matU);
        Set_V(R*matV);
        delete[] d; delete[] e; delete[] w; delete[] z; 
        return *this;
    }
public:
    MRRR_SVD () {};
    MRRR_SVD (const Eigen::MatrixXd& A)
    {
        compute_bsvd(A);
    }
    Eigen::MatrixXd matrixV()
    {
        return V;
    }

    Eigen::MatrixXd matrixU()
    {
        return U;
    }

    Eigen::MatrixXd singularValues()
    {
        return S;
    }
};