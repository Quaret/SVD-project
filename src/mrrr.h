// Telegram: @Quareten
// Кочержук Александр, КМБО-03-22
/*
    Написанное здесь является применением MR^3 для вычисления BSVD. По большей части это является имплементацией black-box подхода
    вычисления с помощью матрицы Голуба-Кахана. Читать -  
    On MR3-type Algorithms for the
    Tridiagonal Symmetric Eigenproblem
    and the Bidiagonal SVD, страница 120.
    https://d-nb.info/1002687853/34/. 
    Дополнительная информация:
    Смысл данной имплементации состоит в том, чтобы показать неэффективность конкретно данного подхода, подробнее читать
    в указанном источнике. Также важно отметить, что данная имплементация не является дословным переводам формул из источника в C++, 
    ибо для верной работы в нашем фреймворке пришлось внести некоторые изменения, по типу неверного (с точки зрения статьи) разбиения 
    собственного вектора в 2 сингулярных. В добавок к этому, функция LAPACK dstemr выдаёт собственные числа и вектора в порядке возрастания,
    когда в нашем фреймворке требуется в порядке убывания. Это приводит к лишним вычислениям и ещё больше ухудшает результат.
    Возможные улучшения: 
    MR^3 знаменит тем, что он позволяет вычислять не полный спектр, а лишь его часть. Поэтому можно вычислять ровно половину собственных чисел,
    и этого будет достаточно.
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
        auto bid = Eigen::internal::UpperBidiagonalization(matrix); // Бидиагонализируем входную матрицу
        auto B = bid.bidiagonal(); // Обозначаем её как B
        int n = B.rows(); // Размерность матрицы
        Eigen::MatrixXd TGK(2*n,2*n); // Создаём матрицу Голуба-Кахана TGK(B)
        TGK.setZero();
        auto L = bid.householderU(); // Не забываем запомнить другие множители бидиагонализации
        auto R = bid.householderV();
        Eigen::MatrixXd Base(n,n); // Приходится создавать такую матрицу, ибо невозможно узнать, какой тип "скормить" m_matrixU и m_matrixV
        Base.setZero(); // А так становится ясно, что у них за тип и размерность
        this->m_matrixU = Base;
        this->m_matrixV = Base;
        Eigen::MatrixXd BB = B; //Переводим матрицу B к виду, где можно использовать взятие по индексу
        for (int i = 0; i < n - 1; ++i) //Считаем TGK(B)
        {
            TGK(2*i,2*i+1) = BB(i,i); TGK(2*i+1,2*i+2) = BB(i,i+1);
            TGK(2*i+1,2*i) = BB(i,i); TGK(2*i+2,2*i+1) = BB(i,i+1);
        }
        TGK(2*n-2,2*n-1) = BB(n-1,n-1);
        TGK(2*n-1,2*n-2) = BB(n-1,n-1);

        int32_t nzc = std::max(1,2*n); //Количество собственных векторов, которые будут храниться в z

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
        double *w = new double[2*n]; //Собственные числа в возрастающем (ВАЖНО) порядке
        double *z = new double[2*n*2*n]; //Собственные вектора
        int32_t m; //Число найденных собственных значений - здесь не важно
        int32_t isuppz [2*2*n]; //Индексы ненулевых векторов z 
        int32_t tryrac = 1; //Логическая переменная для проверки на точность, см. документацию
        int32_t info = LAPACKE_dstemr 
        (
            LAPACK_COL_MAJOR, // int matrix_layout, in
            'V', // char jobz, in
            'A', // char range, in
            2*n, // int32_t n, in
            d,// double *d, in out
            e, // double *e, in out
            0, // double vl, in
            0, // double vu, in
            0, // int32_t il, in
            0, // int32_t iu, in
            &m, // int32_t *m, out
            w, // double in
            z, // double *z, out
            2*n, // int32_t ldz, in
            nzc, // int32_t nzc, in
            isuppz,// int32_t *isuppz, out
            &tryrac// int32_t *tryrac, in outs
        );
        if (info != 0)
            throw std::runtime_error("LAPACK error: " + std::to_string(info));
        
        Eigen::VectorXd eigenvalues (2*n);
        for (int i = 0; i < 2*n; ++i) 
        {
            eigenvalues(i) = w[2*n-i-1];
        }
        
        Eigen::MatrixXd eigenvectors (2*n,2*n);
        for (int i = 0; i < 2*n; i++)
        {
            for (int j = 0; j < 2*n; j++)
            {
                eigenvectors(i,j) = i < n ? z[(2*n-i-1)*2*n+j]*sq : z[(2*n-i-1)*2*n+j] * (-1)*sq;
            }
        }
        for (int i = n-1; i < n+1; i++)
        {
            for (int j = 0; j < 2*n; j++) 
            {
                eigenvectors(i,j) *= -1;
            }
        }
        int ind = 0;
        Eigen::MatrixXd singularValuess (n,n);
        for (int i = 0; i < n; ++i) 
        {
            double sigma = std::abs(eigenvalues(i));
            Eigen::VectorXd q (2*n);
            for (int j = 0; j < 2*n; j++)
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
        this->m_matrixV = R * this->m_matrixV;
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