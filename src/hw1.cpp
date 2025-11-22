#include "hw1.h"

namespace algebra {

    Matrix zeros(size_t n, size_t m) {
        Matrix matrix(n, std::vector<double>(m, 0));
        return matrix;
    }
    Matrix ones(size_t n, size_t m) {
        Matrix matrix(n, std::vector<double>(m, 1));
        return matrix;
    }
    Matrix random(size_t n, size_t m, double min, double max) {
        Matrix matrix(n, std::vector<double>(m));
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(min, max);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                matrix[i][j] = dis(gen);
            }
        }
        return matrix;
    }

    void show(const Matrix& matrix) {
        size_t n = matrix.size();
        size_t m = matrix[0].size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                std::cout << std::setw(6) << std::fixed << std::setprecision(3) << matrix[i][j];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    Matrix multiply(Matrix& matrix, double c) {
        size_t n = matrix.size();
        size_t m = matrix[0].size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                matrix[i][j] *= c;
            }
        }
        return matrix;
    }

    bool multiplyable(const Matrix& matrix1, const Matrix& matrix2) {
        if (matrix1[0].size() == matrix2.size())return true;
        else return false;
    }
    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2) {
        if (!multiplyable(matrix1, matrix2)) {
            std::cerr << "Error the two matrixs cannot be multiplied!" << std::endl;
            return matrix1;
        }
        size_t n = matrix1.size();
        size_t m = matrix2[0].size();
        size_t s = matrix2.size();
        Matrix result = zeros(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                for (size_t k = 0; k < s; ++k) {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }

        }
        return result;
    }
    Matrix sum(Matrix& matrix, double c) {
        size_t n = matrix.size();
        size_t m = matrix[0].size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                matrix[i][j] += c;
            }
        }
        return matrix;
    }
    Matrix sum(const Matrix& matrix1, const Matrix& matrix2) {
        size_t n = matrix1.size();
        size_t m = matrix1[0].size();
        Matrix result = matrix1;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                result[i][j] += matrix2[i][j];
            }
        }
        return result;
    }

    Matrix transpose(const Matrix& matrix) {
        size_t n = matrix.size();
        size_t m = matrix[0].size();
        Matrix result(m, std::vector<double>(n));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }
    Matrix minor(const Matrix& matrix, size_t n, size_t m) {
        size_t row = matrix.size();
        size_t col = matrix[0].size();
        size_t r = 0;
        Matrix result(row - 1, std::vector<double>(col - 1));
        for (size_t i = 0; i < row; ++i) {
            if (i == n)continue;
            size_t c = 0;
            for (size_t j = 0; j < col; ++j) {
                if (j == m)continue;
                else {
                    result[r][c] = matrix[i][j];
                    c++;

                }

            }
            r++;
        }
        return result;
    }

    double determinant(const Matrix& matrix) {
        size_t n = matrix.size();
        if (n == 1)return matrix[0][0];
        else if (n == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }
        else {
            double sum = 0;
            for (size_t i = 0; i < n; ++i) {
                if (i % 2 == 0)
                    sum += matrix[0][i] * determinant(minor(matrix, 0, i));
                else
                    sum -= matrix[0][i] * determinant(minor(matrix, 0, i));
            }
            return sum;
        }
    }
    Matrix inverse(const Matrix& matrix) {
        double det = determinant(matrix);
        if (std::abs(det) < 1e-10) {
            std::cerr << "Error cannot be inversed!";
            return matrix;
        }
        size_t n = matrix.size();
        Matrix result(n, std::vector<double>(n, 0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if ((i + j) % 2 == 0)
                    result[i][j] = determinant(minor(matrix, j, i)) / det;
                else
                    result[i][j] = -determinant(minor(matrix, j, i)) / det;
            }
        }
        return result;
    }
    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis ) {
        if (axis == 0) {//put mat2 under mat1
            size_t n1 = matrix1.size();
            size_t m1 = matrix1[0].size();
            size_t n2 = matrix2.size();
            size_t m2 = matrix2[0].size();
            if (m1 != m2) {
                std::cerr << "Error:they cannot be concatenate due to different rows";
                Matrix tem;
                return tem;
            }
            else {
                Matrix result(n1 + n2, std::vector<double>(m1));
                for (size_t i = 0; i < n1; ++i) {
                    for (size_t j = 0; j < m1; ++j) {
                        result[i][j] = matrix1[i][j];
                    }
                }
                for (size_t i = 0; i < n2; ++i) {
                    for (size_t j = 0; j < m2; ++j) {
                        result[n1 + i][j] = matrix2[i][j];
                    }
                }
                return result;
            }
        }
        else if (axis == 1) {
            size_t n1 = matrix1.size();
            size_t m1 = matrix1[0].size();
            size_t n2 = matrix2.size();
            size_t m2 = matrix2[0].size();
            if (n1 != n2) {
                std::cerr << "Error:they cannot be concatenate due to different rows";
                Matrix tem;
                return tem;
            }
            else {
                Matrix result(n1, std::vector<double>(m1 + m2));
                for (size_t i = 0; i < n1; ++i) {
                    for (size_t j = 0; j < m1; ++j) {
                        result[i][j] = matrix1[i][j];
                    }
                }
                for (size_t i = 0; i < n2; ++i) {
                    for (size_t j = 0; j < m2; ++j) {
                        result[i][m1 + j] = matrix2[i][j];
                    }
                }
                return result;
            }
        }
    }

    Matrix ero_swap( Matrix& matrix, size_t r1, size_t r2) {
        swap(matrix[r1], matrix[r2]);
        return matrix;
    }
    Matrix ero_multiply( Matrix& matrix, size_t r, double c) {
        size_t m = matrix[0].size();
        for (int i = 0; i < m; ++i) {
            matrix[r][i] *= c;
        }
        return matrix;
    }
    Matrix ero_sum(Matrix& matrix, size_t r1, double c, size_t r2) {
        size_t m = matrix[0].size();
        for (int i = 0; i < m; ++i) {
            matrix[r1][i] += c * matrix[r2][i];
        }
        return matrix;
    }
    Matrix upper_triangular(const Matrix& matrix)
    {
        Matrix result = matrix;  // 创建副本
        size_t n = result.size();

        for (size_t k = 0; k < n - 1; ++k) {
            // 处理主对角线为0的情况（基础版本）
            if (result[k][k] == 0) {
                // 寻找下面行中同列非零元素进行行交换
                bool found = false;
                for (size_t i = k + 1; i < n; ++i) {
                    if (result[i][k] != 0) {
                        std::swap(result[k], result[i]);  // 行交换
                        found = true;
                        break;
                    }
                }
                if (!found) continue;  // 如果整列都为0，跳过
            }

            // 高斯消元
            for (size_t i = k + 1; i < n; ++i) {
                double multiplier = result[i][k] / result[k][k];
                for (size_t j = k; j < n; ++j) {
                    result[i][j] -= multiplier * result[k][j];
                }
            }
        }

        return result;
    }


}