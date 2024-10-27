// Telegram: @Quareten
// Кочержук Александр, КМБО-03-22
/*
    Написанное здесь является применением MR^3 для вычисления BSVD. По большей части это является имплементацией black-box подхода
    вычисления с помощью матрицы Голуба-Кахана. Читать -  
    On MR3-type Algorithms for the
    Tridiagonal Symmetric Eigenproblem
    and the Bidiagonal SVD, страница 120.
    В данной программе есть 2 проблемы: невозможность считать разложение для матрицы размера больше, чем 3x3 - связана с какими-то внутренними
    реализациями LAPACK и Eigen, которые убивают m_matrixU и m_matrixV после dstemr из-за того, что данные матрицы инициализируются заранее,
    (ПРЕДПОЛОЖЕНИЕ) однако выделяются на куче возле памяти для dstemr, который после себя вычищает всё и данные ( и не только данные) матрицы в том числе; 
    2-ая проблема же обнаруживается в mrrr_test_table.txt - последний столбец AVG relative err. sigma выглядит как параметр функции svd_test_fuct -1, что странно,
    ибо в m_singularvalues нормально передаются сингулярные числа.
    Неоптимальность кода даже не обсуждается, некоторые вещи сделаны посредственно (перенос значений из векторов).
*/
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
template<typename _MatrixType> class MRRR_SVD : public Eigen::SVDBase<MRRR_SVD<_MatrixType> >
{
private:
    typedef Eigen::SVDBase<MRRR_SVD<_MatrixType>> Base;
public:
    MRRR_SVD (const _MatrixType &matrix, unsigned int computationOptions) {
        compute_tgk(matrix);
    };
    MRRR_SVD<_MatrixType>& compute_tgk (const _MatrixType& matrix)
    {
        // std::cout << "Original matrix:" << std::endl << matrix << std::endl;
        auto bid = Eigen::internal::UpperBidiagonalization(matrix); // Бидиагонализируем входную матрицу
        auto B = bid.bidiagonal(); // Обозначаем её как B
        int n = B.rows();
        Eigen::MatrixXd TGK(2*n,2*n); // Создаём матрицу Голуба-Кахана TGK(B)
        TGK.setZero();
        auto L = bid.householderU();
        auto R = bid.householderV();
        Eigen::MatrixXd Base(n,n);
        Base.setZero();
        this->m_matrixU = Base;
        this->m_matrixV = Base;
        // this->m_matrixU = bid.householderU(); //Инициализурем всяким мусором для того, чтобы потом туда записывать уже нормальные вектора с помощью обычных операций
        // std::cout << "L:" << std::endl << L << std::endl;
        // this->m_matrixV = bid.householderV().transpose(); //без выноса мозга с Eigen
        // std::cout << "R:" << std::endl << R << std::endl;
        Eigen::MatrixXd BB = B; //Переводим матрицу B к виду, где можно использовать взятие по индексу

        int nn = 2*n;
        for (int i = 0; i < n - 1; ++i) //Считаем TGK(B)
        {
            TGK(2*i,2*i+1) = BB(i,i); TGK(2*i+1,2*i+2) = BB(i,i+1);
            TGK(2*i+1,2*i) = BB(i,i); TGK(2*i+2,2*i+1) = BB(i,i+1);
        }
        TGK(nn-2,nn-1) = BB(n-1,n-1);
        TGK(nn-1,nn-2) = BB(n-1,n-1);

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
        int ind = 0;
        Eigen::MatrixXd singularValuess (n,n);
        for (int i = 0; i < n; ++i) 
        {
            double sigma = std::abs(eigenvalues(i));
            Eigen::VectorXd q (nn);
            for (int j = 0; j < nn; j++)
            {
                q[j] = eigenvectors(i,j);
            }
            Eigen::VectorXd u (n);
            Eigen::VectorXd v (n);
            for (int j = 0; j < n; ++j)
            {
                v[j] = q[2*j];
                u[j] = q[2*j+1];
            }
            singularValuess(i,i) = sigma;
            this->m_matrixU.col(i) = u;
            this->m_matrixV.col(i) = v;
        }
        this->m_singularValues = singularValuess.diagonal(0);
        this->m_matrixU = L * this->m_matrixU;
        // this->m_matrixV = this->m_matrixV.transpose() * R.transpose();
        this->m_matrixV = R * this->m_matrixV;
        // std::cout << "S:" << std::endl << this->m_singularValues << std::endl;
        // std::cout << "U:" << std::endl << this->m_matrixU << std::endl;
        // std::cout << "V transposed:" << std::endl << this->m_matrixV << std::endl;
        // Eigen::MatrixXd Ans = this->m_matrixU*singularValuess*this->m_matrixV;
        // std::cout << "My matrix: " << std::endl << this->m_matrixU*singularValuess*this->m_matrixV << std::endl;
        this->m_isInitialized = true;
        this->m_computeFullU = true; 
        this->m_computeFullV = true;
        delete[] d; delete[] e; delete[] w; delete[] z; 
        return *this;
    }
    
};
template<typename _MatrixType> 
struct Eigen::internal::traits<MRRR_SVD<_MatrixType>>
        : Eigen::internal::traits<_MatrixType>
{
    typedef _MatrixType MatrixType;
};