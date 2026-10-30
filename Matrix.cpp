//
// Created by valleron on 10/15/26.
//

#include "Matrix.h"

#include <iostream>
#include <iomanip>
#include "MatrixView.h"
#include "ConstMatrixView.h"
#include "Matrix_operators.cpp"

Matrix::Matrix(const unsigned int rows, const unsigned int columns) {
    this->rows = rows;
    this->columns = columns;
    data.resize(rows * columns);
}

unsigned Matrix::n_rows() const {
    return this->rows;
}
unsigned Matrix::n_columns() const {
    return this->columns;
}


void Matrix::fillMatrix(double(*f)(unsigned int, unsigned int)) {
    for (unsigned int i = 0; i < rows; i++) {
        for (unsigned int j = 0; j < columns; j++) {
            data.at(i * rows + j) = f(i, j);
        }
    }
}

Matrix Matrix::binet_multiply(Matrix *a, Matrix *b) {
    if (a->rows != b->rows || a->columns != b->columns || a->rows != a->columns) {
        throw std::invalid_argument("Matrix multiplication matrix is not the same size");
    }
    int n = a->rows;
    Matrix out = Matrix(n, n);

    auto c = MatrixView(out, 0, 0, n, n);
    auto a_view = ConstMatrixView(*a, 0, 0, n, n);
    auto b_view = ConstMatrixView(*b, 0, 0, n, n);
    binet_recursion(a_view, b_view, c);
    return out;
}



void Matrix::binet_recursion(ConstMatrixView a, ConstMatrixView b, MatrixView c) {
    unsigned int n = a.n_rows(), m = b.n_rows(), p = b.n_cols();
    if (a.n_cols() != m) {

        throw std::invalid_argument("Matrix multiplication matrix is not the same size");
    }
    if (m*n*p == 0) {return;}
    unsigned int max = std::max(std::max(n, m), p);
    if (max < 2) {
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int j = 0; j < n; ++j) {
                for (unsigned int k = 0; k < n; ++k) {
                    c(i, j) += a(i, k) * b(k, j);
                }
            }
        }
        c.add_flops(n * n * n);
        return;
    }
    else if (max == n) {
        ConstMatrixView a1 = a.sub(0, 0, m, n / 2);
        ConstMatrixView a2 = a.sub(n / 2, 0, m, (n + 1)/ 2);
        MatrixView c1 = c.sub(0, 0, p, n / 2);
        MatrixView c2 = c.sub(n / 2, 0, p, (n + 1) / 2);
        binet_recursion(a1, b, c1);
        binet_recursion(a2, b, c2);
    }
    else if (max == p) {
        ConstMatrixView b1 = b.sub(0, 0, p / 2, m);
        ConstMatrixView b2 = b.sub(0, p / 2, (p + 1) / 2, m);
        MatrixView c1 = c.sub(0, 0, p / 2, n);
        MatrixView c2 = c.sub(0, p / 2, (p + 1) / 2, n);
        binet_recursion(a, b1, c1);
        binet_recursion(a, b2, c2);
    }
    else if (max == m) {
        ConstMatrixView a1 = a.sub(0, 0, m / 2, n);
        ConstMatrixView a2 = a.sub(0, m/2, (m + 1) / 2, n);
        ConstMatrixView b1 = b.sub(0, 0, p, m / 2);
        ConstMatrixView b2 = b.sub(m / 2, 0, p, (m + 1) / 2);
        binet_recursion(a1, b1, c);
        binet_recursion(a2, b2, c);
    }
}

Matrix Matrix::strassen_multiply(Matrix *a, Matrix *b) {
    if (a->rows != b->rows || a->columns != b->columns || a->rows != a->columns) {
        throw std::invalid_argument("Matrix multiplication matrix is not the same size");
    }
    int n = a->rows;
    Matrix out = Matrix(n, n);

    auto c = MatrixView(out, 0, 0, n, n);
    auto a_view = ConstMatrixView(*a, 0, 0, n, n);
    auto b_view = ConstMatrixView(*b, 0, 0, n, n);
    strassen_recursion(a_view, b_view, c);
    return out;

}

void Matrix::strassen_recursion(ConstMatrixView A, ConstMatrixView B, MatrixView C) {


    unsigned int n = A.n_rows();
    if (n < 2) {
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int j = 0; j < n; ++j) {
                for (unsigned int k = 0; k < n; ++k) {
                    C(i, j) += A(i, k) * B(k, j);
                }
            }
        }
        C.add_flops(n * n * n * 2);
        return;
    }

    if (n % 2 != 0) {
        unsigned int k = n - 1;

        ConstMatrixView A_prime = A.sub(0, 0, k, k);
        ConstMatrixView B_prime = B.sub(0, 0, k, k);
        MatrixView C_prime = C.sub(0, 0, k, k);

        strassen_recursion(A_prime, B_prime, C_prime);


        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                C_prime(i, j) += A(i, k) * B(k, j);
            }
        }
        C.add_flops(k * k * 2);

        double beta = B(k, k);
        for (unsigned int i = 0; i < k; ++i) {
            double sum = 0.0;
            for (unsigned int j = 0; j < k; ++j) {
                sum += A_prime(i, j) * B(j, k); // A' * w
            }
            C(i, k) = sum + beta * A(i, k);
        }
        C.add_flops(k * (k * 2 + 2));

        double alpha = A(k, k);
        for (unsigned int j = 0; j < k; ++j) {
            double sum = 0.0;
            for (unsigned int i = 0; i < k; ++i) {
                sum += A(k, i) * B_prime(i, j); // v^T * B'
            }
            C(k, j) = sum + alpha * B(k, j);
        }

        C.add_flops(k * (k * 2 + 2));

        double gamma = 0.0;
        for (unsigned int i = 0; i < k; ++i) {
            gamma += A(k, i) * B(i, k); // v^T * w
        }
        C(k, k) = gamma + alpha * beta;

        C.add_flops(k * 2 + 2);
        return;
    }

    ConstMatrixView A11 = A.sub(0, 0, n / 2, n/2);
    ConstMatrixView A12 = A.sub(0, n / 2, n / 2, n/2);
    ConstMatrixView A21 = A.sub(n / 2, 0, n/2, n/2);
    ConstMatrixView A22 = A.sub(n / 2, n / 2, n/2, n/2);

    ConstMatrixView B11 = B.sub(0, 0, n/2, n/2);
    ConstMatrixView B12 = B.sub(0, n/2, n/2, n/2);
    ConstMatrixView B21 = B.sub(n / 2, 0, n/2, n/2);
    ConstMatrixView B22 = B.sub(n / 2, n / 2, n/2, n/2);

    MatrixView C11 = C.sub(0, 0, n / 2, n/2);
    MatrixView C12 = C.sub(0, n / 2, n/2, n/2);
    MatrixView C21 = C.sub(n / 2, 0, n/2, n/2);
    MatrixView C22 = C.sub(n / 2, n / 2, n/2, n/2);


    Matrix M1 = compute_helpers(A11 + A22, B11 + B22, n);
    Matrix M2 = compute_helpers(A21 + A22, B11, n);
    Matrix M3 = compute_helpers(A11, B12 - B22, n);
    Matrix M4 = compute_helpers(A22, B21 - B11, n);
    Matrix M5 = compute_helpers(A11 + A12, B22, n);
    Matrix M6 = compute_helpers(A21 - A11, B11 + B12, n);
    Matrix M7 = compute_helpers(A12 - A22, B21 + B22, n);

    C.add_flops(M1.flops + M2.flops + M3.flops + M4.flops + M5.flops + M6.flops + M7.flops);

    C11.write_into(M1 + M4 - M5 + M7);
    C12.write_into(M3 + M5);
    C21.write_into(M2 + M4);
    C22.write_into(M1 - M2 + M3 + M6);

    unsigned int half = n / 2;
    unsigned int additions_flops = 18 * (half * half);
    C.add_flops(additions_flops);

}


Matrix Matrix::compute_helpers(const ConstMatrixView &left, const ConstMatrixView &right, const unsigned n) {
    Matrix M(n/2, n/2);
    strassen_recursion(left, right, MatrixView(M));
    return M;
}

Matrix Matrix::laderman_multiply(Matrix *a, Matrix *b) {
    if (a->rows != b->rows || a->columns != b->columns || a->rows != a->columns) {
        throw std::invalid_argument("Matrix multiplication matrix is not the same size");
    }
    int n = a->rows;
    Matrix out = Matrix(n, n);

    auto c = MatrixView(out, 0, 0, n, n);
    auto a_view = ConstMatrixView(*a, 0, 0, n, n);
    auto b_view = ConstMatrixView(*b, 0, 0, n, n);
    laderman_recursion(a_view, b_view, c);
    return out;
}

void Matrix::laderman_recursion(const ConstMatrixView &A, const ConstMatrixView &B, MatrixView C) {
    unsigned int n = A.n_rows();
    if (n < 3) {
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int j = 0; j < n; ++j) {
                for (unsigned int k = 0; k < n; ++k) {
                    C(i, j) += A(i, k) * B(k, j);
                }
            }
        }
        C.add_flops(n * n * n);
        return;
    }

    unsigned int r = n % 3;


    if (r != 0) {
        unsigned int k = n - r;

        ConstMatrixView A_prime = A.sub(0, 0, k, k);
        ConstMatrixView B_prime = B.sub(0, 0, k, k);
        MatrixView C_prime = C.sub(0, 0, k, k);

        laderman_recursion(A_prime, B_prime, C_prime);

        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                for (unsigned int p = 0; p < r; ++p) {
                    C_prime(i, j) += A(i, k + p) * B(k + p, j);
                }
            }
        }
        C.add_flops(r * k * k * 2);

        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < r; ++j) {
                double sum = 0.0;
                for (unsigned int p = 0; p < n; ++p) {
                    sum += A(i, p) * B(p, k + j);
                }
                C(i, k + j) = sum;
            }
        }
        C.add_flops(k * r * n * 2 + k * r);

        for (unsigned int i = 0; i < r; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                double sum = 0.0;
                for (unsigned int p = 0; p < n; ++p) {
                    sum += A(k + i, p) * B(p, j);
                }
                C(k + i, j) = sum;
            }
        }

        C.add_flops(k * r * n * 2 + k * r);
        for (unsigned int i = 0; i < r; ++i) {
            for (unsigned int j = 0; j < r; ++j) {
                double sum = 0.0;
                for (unsigned int p = 0; p < n; ++p) {
                    sum += A(k + i, p) * B(p, k + j);
                }
                C(k + i, k + j) = sum;
            }
        }

        C.add_flops(r * r * n * 2 + r * r);
        return;
    }

    ConstMatrixView A11 = A.sub(0, 0, n / 3, n / 3);
    ConstMatrixView A12 = A.sub(0, n / 3, n / 3, n / 3);
    ConstMatrixView A13 = A.sub(0, 2 * n / 3, n / 3, n / 3);
    ConstMatrixView A21 = A.sub(n / 3, 0, n / 3, n / 3);
    ConstMatrixView A22 = A.sub(n / 3,  n / 3, n / 3, n / 3);
    ConstMatrixView A23 = A.sub( n / 3, 2 * n / 3, n / 3, n / 3);
    ConstMatrixView A31 = A.sub(2 * n / 3, 0, n / 3, n / 3);
    ConstMatrixView A32 = A.sub(2 * n / 3, n / 3, n / 3, n / 3);
    ConstMatrixView A33 = A.sub(2 * n / 3, 2 * n / 3, n / 3, n / 3);

    ConstMatrixView B11 = B.sub(0, 0, n / 3, n / 3);
    ConstMatrixView B12 = B.sub(0, n / 3, n / 3, n / 3);
    ConstMatrixView B13 = B.sub(0, 2 * n / 3, n / 3, n / 3);
    ConstMatrixView B21 = B.sub(n / 3, 0, n / 3, n / 3);
    ConstMatrixView B22 = B.sub(n / 3,  n / 3, n / 3, n / 3);
    ConstMatrixView B23 = B.sub( n / 3, 2 * n / 3, n / 3, n / 3);
    ConstMatrixView B31 = B.sub(2 * n / 3, 0, n / 3, n / 3);
    ConstMatrixView B32 = B.sub(2 * n / 3, n / 3, n / 3, n / 3);
    ConstMatrixView B33 = B.sub(2 * n / 3, 2 * n / 3, n / 3, n / 3);

    MatrixView C11 = C.sub(0, 0, n / 3, n / 3);
    MatrixView C12 = C.sub(0, n / 3, n / 3, n / 3);
    MatrixView C13 = C.sub(0, 2 * n / 3, n / 3, n / 3);
    MatrixView C21 = C.sub(n / 3, 0, n / 3, n / 3);
    MatrixView C22 = C.sub(n / 3,  n / 3, n / 3, n / 3);
    MatrixView C23 = C.sub( n / 3, 2 * n / 3, n / 3, n / 3);
    MatrixView C31 = C.sub(2 * n / 3, 0, n / 3, n / 3);
    MatrixView C32 = C.sub(2 * n / 3, n / 3, n / 3, n / 3);
    MatrixView C33 = C.sub(2 * n / 3, 2 * n / 3, n / 3, n / 3);

    Matrix M1 = compute_helpers_laderman(A11 + A12 + A13 - A21 - A22 - A32 - A33, B22, n);
    Matrix M2 = compute_helpers_laderman(A11 - A21, B22 - B12, n);
    Matrix M3 = compute_helpers_laderman(A22, B12 - B11 + B21 - B22 - B23 - B31 + B33, n);
    Matrix M4 = compute_helpers_laderman(A21 - A11 + A22, B11 - B12 + B22, n);
    Matrix M5 = compute_helpers_laderman(A21 + A22, B12 - B11,n);
    Matrix M6 = compute_helpers_laderman(A11, B11, n);
    Matrix M7 = compute_helpers_laderman(A31 - A11 + A32, B11 - B13 + B23, n);
    Matrix M8 = compute_helpers_laderman(A31 - A11, B13 - B23,n);
    Matrix M9 = compute_helpers_laderman(A31 + A32, B13 - B11, n);
    Matrix M10 = compute_helpers_laderman(A11 + A12 + A13 - A22 - A23 - A31 - A32, B23, n);
    Matrix M11 = compute_helpers_laderman(A32, B13 - B11+ B21 - B22 - B23 - B31 + B32, n);
    Matrix M12 = compute_helpers_laderman(A32 - A13 + A33, B22 + B31 - B32, n);
    Matrix M13 = compute_helpers_laderman(A13 - A33, B22 - B32, n);
    Matrix M14 = compute_helpers_laderman(A13, B31,n);
    Matrix M15 = compute_helpers_laderman(A32 + A33, B32 - B31, n);
    Matrix M16 = compute_helpers_laderman(A22 - A13 + A23, B23 + B31 - B33, n);
    Matrix M17 = compute_helpers_laderman(A13- A23, B23 - B33, n);
    Matrix M18 = compute_helpers_laderman(A22 + A23, B33 - B31, n);
    Matrix M19 = compute_helpers_laderman(A12, B21, n);
    Matrix M20 = compute_helpers_laderman(A23, B32, n);
    Matrix M21 = compute_helpers_laderman(A21, B13, n);
    Matrix M22 = compute_helpers_laderman(A31, B12, n);
    Matrix M23 = compute_helpers_laderman(A33, B33, n);

    C11.write_into(M6 + M14 + M19);
    C12.write_into(M1 + M4 + M5 + M6 + M12 + M14 + M15);
    C13.write_into(M6 + M7 + M9 + M10 + M14 + M16 + M18);
    C21.write_into(M2 + M3 + M4 + M6 + M14 + M16 + M17);
    C22.write_into(M2 + M4 + M5 + M6 + M20);
    C23.write_into(M14 + M16 + M17 + M18 + M21);
    C31.write_into(M6 + M7 + M8 + M11 + M12 + M13 + M14);
    C32.write_into(M12 + M13 + M14 + M15 + M22);
    C33.write_into(M6 + M7 + M8 + M9 + M23);


    C.add_flops(M1.flops + M2.flops + M3.flops + M4.flops + M5.flops +
                M6.flops + M7.flops + M8.flops + M9.flops + M10.flops +
                M11.flops + M12.flops + M13.flops + M14.flops + M15.flops +
                M16.flops + M17.flops + M18.flops + M19.flops + M20.flops +
                M21.flops + M22.flops + M23.flops);

    unsigned int k = n / 3;
    unsigned int additions_flops = 98 * (k * k);

    C.add_flops(additions_flops);

}

Matrix Matrix::compute_helpers_laderman(const ConstMatrixView &left, const ConstMatrixView &right, unsigned n) {
    Matrix M(n/3, n/3);
    laderman_recursion(left, right, MatrixView(M));
    return M;
}








void Matrix::print() {
    std::cout << "\n";
    for (unsigned int r = 0; r < rows; ++r) {
        std::cout << "| ";
        for (unsigned int c = 0; c < columns; ++c) {
            std::cout << std::setw(11) << std::setprecision(8) << data[r * rows + c] << " ";
        }
        std::cout << "| \n";
    }
    std::cout << "\n";
}


