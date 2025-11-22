#ifndef HW1_H
#define HW1_H
#include<iostream>
#include <vector>
#include<random>
#include<iomanip>
#include<algorithm>

namespace algebra {
    using Matrix = std::vector<std::vector<double>>;

    // 基本矩阵创建函数
    Matrix zeros(size_t n, size_t m);
    Matrix ones(size_t n, size_t m);
    Matrix random(size_t n, size_t m, double min, double max);

    // 矩阵显示
    void show(const Matrix& matrix);

    // 矩阵运算
    bool multiplyable(const Matrix& matrix1, const Matrix& matrix2);
    Matrix multiply( Matrix& matrix, double c);
    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2);
    Matrix sum(const Matrix& matrix, double c);
    Matrix sum(const Matrix& matrix1, const Matrix& matrix2);
    Matrix transpose(const Matrix& matrix);

    // 矩阵操作
    Matrix minor(const Matrix& matrix, size_t n, size_t m);
    double determinant(const Matrix& matrix);
    Matrix inverse(const Matrix& matrix);
    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis = 0);

    // 初等行变换
    Matrix ero_swap( Matrix& matrix, size_t r1, size_t r2);
    Matrix ero_multiply( Matrix& matrix, size_t r, double c);
    Matrix ero_sum( Matrix& matrix, size_t r1, double c, size_t r2);

    // 上三角矩阵
    
    Matrix upper_triangular( Matrix& matrix);

} // namespace algebra

#endif