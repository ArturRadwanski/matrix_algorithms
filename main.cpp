#include <iostream>

#include <vector>

#include <vector>
#include <stdexcept>
#include <cmath>
#include <iomanip>

class Matrix {
public:
    unsigned get_n() const { return n; }
    unsigned get_m() const { return m; }

    Matrix(const unsigned n, const unsigned m)
        : n(n), m(m), mat(n, std::vector<double>(m, 0.0)), operation_count(0) {}

    void set(const unsigned i, const unsigned j, const double val) {
        if (i < n && j < m) {
            mat[i][j] = val;
        } else {
            throw std::out_of_range("Matrix index out of range");
        }
    }

    double get(unsigned i, unsigned j) const {
        if (i < n && j < m) {
            return mat[i][j];
        } else {
            throw std::out_of_range("Matrix index out of range");
        }
    }

    unsigned get_operation_count() const { return operation_count; }

    void print_matrix() {
        std::cout << std::endl;
        for (unsigned i = 0; i < n; ++i) {
            std::cout << "|";
            for (unsigned j = 0; j < m; ++j) {
                std::cout << std::setw(6) << mat[i][j] << " ";
            }
            std::cout << "|\n";
        }
        std::cout << std::endl;
    }

    friend Matrix bine_multiply(Matrix& a, Matrix& b);
    friend Matrix strassen_multiply(Matrix& a, Matrix& b);
    friend Matrix laderman_multiply(Matrix& a, Matrix& b);

private:
    const unsigned n;
    const unsigned m;
    std::vector<std::vector<double>> mat;

protected:
    unsigned operation_count;
    static Matrix add_sub(const Matrix& A, const unsigned rA, const unsigned cA,
                          const Matrix& B, const unsigned rB, const unsigned cB,
                          const unsigned size, unsigned& ops) {
        Matrix res(size, size);
        for (unsigned i = 0; i < size; ++i) {
            for (unsigned j = 0; j < size; ++j) {
                res.mat.at(i).at(j) = A.mat.at(rA + i).at(cA + j) + B.mat.at(rB + i).at(cB + j);
                ops++;
            }
        }
        return res;
    }

    static Matrix bine_sub(const Matrix& A, unsigned rA, unsigned cA,
                           const Matrix& B, unsigned rB, unsigned cB,
                           unsigned size, unsigned& ops) {
        // Baza rekurencji
        if (size <= 2) {
            Matrix res(size, size);
            for (unsigned i = 0; i < size; ++i) {
                for (unsigned j = 0; j < size; ++j) {
                    double sum = 0.0;
                    for (unsigned k = 0; k < size; ++k) {
                        sum += A.mat[rA + i][cA + k] * B.mat[rB + k][cB + j];
                        ops += 2;
                    }
                    res.mat[i][j] = sum;
                }
            }
            return res;
        }

        unsigned k = size / 2;

        Matrix P1 = bine_sub(A, rA, cA,           B, rB, cB,           k, ops);
        Matrix P2 = bine_sub(A, rA, cA + k,      B, rB + k, cB,      k, ops);
        Matrix P3 = bine_sub(A, rA, cA,           B, rB, cB + k,      k, ops);
        Matrix P4 = bine_sub(A, rA, cA + k,      B, rB + k, cB + k,  k, ops);
        Matrix P5 = bine_sub(A, rA + k, cA,      B, rB, cB,           k, ops);
        Matrix P6 = bine_sub(A, rA + k, cA + k,  B, rB + k, cB,      k, ops);
        Matrix P7 = bine_sub(A, rA + k, cA,      B, rB, cB + k,      k, ops);
        Matrix P8 = bine_sub(A, rA + k, cA + k,  B, rB + k, cB + k,  k, ops);

        Matrix C11 = add_sub(P1, 0, 0, P2, 0, 0, k, ops);
        Matrix C12 = add_sub(P3, 0, 0, P4, 0, 0, k, ops);
        Matrix C21 = add_sub(P5, 0, 0, P6, 0, 0, k, ops);
        Matrix C22 = add_sub(P7, 0, 0, P8, 0, 0, k, ops);

        Matrix C(size, size);
        for (unsigned i = 0; i < k; ++i) {
            for (unsigned j = 0; j < k; ++j) {
                C.mat[i][j]         = C11.mat[i][j];
                C.mat[i][j + k]     = C12.mat[i][j];
                C.mat[i + k][j]     = C21.mat[i][j];
                C.mat[i + k][j + k] = C22.mat[i][j];
            }
        }
        return C;
    }
    static Matrix sub_sub(const Matrix& A, const unsigned rA, const unsigned cA,
                          const Matrix& B, const unsigned rB, const unsigned cB,
                          const unsigned size, unsigned& ops) {
        Matrix res(size, size);
        for (unsigned i = 0; i < size; ++i) {
            for (unsigned j = 0; j < size; ++j) {
                res.mat.at(i).at(j) = A.mat.at(rA + i).at(cA + j) - B.mat.at(rB + i).at(cB + j);
                ops++;
            }
        }
        return res;
    }

    static Matrix strassen_sub(const Matrix& A, unsigned rA, unsigned cA,
                               const Matrix& B, unsigned rB, unsigned cB,
                               unsigned size, unsigned& ops) {
        // Baza rekurencji
        if (size <= 2) {
            Matrix res(size, size);
            for (unsigned i = 0; i < size; ++i) {
                for (unsigned j = 0; j < size; ++j) {
                    double sum = 0.0;
                    for (unsigned k = 0; k < size; ++k) {
                        sum += A.mat[rA + i][cA + k] * B.mat[rB + k][cB + j];
                        ops += 2;
                    }
                    res.mat[i][j] = sum;
                }
            }
            return res;
        }

        unsigned k = size / 2;


        // M1 = (A11 + A22) * (B11 + B22)
        Matrix A_11_22 = add_sub(A, rA, cA, A, rA + k, cA + k, k, ops);
        Matrix B_11_22 = add_sub(B, rB, cB, B, rB + k, cB + k, k, ops);
        Matrix M1 = strassen_sub(A_11_22, 0, 0, B_11_22, 0, 0, k, ops);

        // M2 = (A21 + A22) * B11
        Matrix A_21_22 = add_sub(A, rA + k, cA, A, rA + k, cA + k, k, ops);
        Matrix M2 = strassen_sub(A_21_22, 0, 0, B, rB, cB, k, ops);

        // M3 = A11 * (B12 - B22)
        Matrix B_12_22 = sub_sub(B, rB, cB + k, B, rB + k, cB + k, k, ops);
        Matrix M3 = strassen_sub(A, rA, cA, B_12_22, 0, 0, k, ops);

        // M4 = A22 * (B21 - B11)
        Matrix B_21_11 = sub_sub(B, rB + k, cB, B, rB, cB, k, ops);
        Matrix M4 = strassen_sub(A, rA + k, cA + k, B_21_11, 0, 0, k, ops);

        // M5 = (A11 + A12) * B22
        Matrix A_11_12 = add_sub(A, rA, cA, A, rA, cA + k, k, ops);
        // Poprawka: B22 zaczyna się w B na (rB + k, cB + k)
        Matrix M5 = strassen_sub(A_11_12, 0, 0, B, rB + k, cB + k, k, ops);

        // M6 = (A21 - A11) * (B11 + B12)
        Matrix A_21_11 = sub_sub(A, rA + k, cA, A, rA, cA, k, ops);
        Matrix B_11_12 = add_sub(B, rB, cB, B, rB, cB + k, k, ops);
        Matrix M6 = strassen_sub(A_21_11, 0, 0, B_11_12, 0, 0, k, ops);

        // M7 = (A12 - A22) * (B21 + B22)
        Matrix A_12_22 = sub_sub(A, rA, cA + k, A, rA + k, cA + k, k, ops);
        Matrix B_21_22 = add_sub(B, rB + k, cB, B, rB + k, cB + k, k, ops);
        Matrix M7 = strassen_sub(A_12_22, 0, 0, B_21_22, 0, 0, k, ops);


        Matrix C(size, size);
        for (unsigned i = 0; i < k; ++i) {
            for (unsigned j = 0; j < k; ++j) {
                C.mat[i][j]         = M1.mat[i][j] + M4.mat[i][j] - M5.mat[i][j] + M7.mat[i][j];
                C.mat[i][j + k]     = M3.mat[i][j] + M5.mat[i][j];
                C.mat[i + k][j]     = M2.mat[i][j] + M4.mat[i][j];
                C.mat[i + k][j + k] = M1.mat[i][j] - M2.mat[i][j] + M3.mat[i][j] + M6.mat[i][j];
                ops += 12;
            }
        }
        return C;
    }
    static Matrix add(const Matrix& A, const Matrix& B, unsigned& ops) {
        unsigned n = A.get_n();
        Matrix res(n, n);
        for (unsigned i = 0; i < n; ++i) {
            for (unsigned j = 0; j < n; ++j) {
                res.mat[i][j] = A.mat[i][j] + B.mat[i][j];
                ops++;
            }
        }
        return res;
    }
    static Matrix subtract(const Matrix& A, const Matrix& B, unsigned& ops) {
        unsigned n = A.get_n();
        Matrix res(n, n);
        for (unsigned i = 0; i < n; ++i) {
            for (unsigned j = 0; j < n; ++j) {
                res.mat[i][j] = A.mat[i][j] - B.mat[i][j];
                ops++;
            }
        }
        return res;
    }
    static Matrix get_block(const Matrix& src, unsigned row_block, unsigned col_block, unsigned k) {
        Matrix res(k, k);
        unsigned r_start = row_block * k;
        unsigned c_start = col_block * k;
        for (unsigned i = 0; i < k; ++i) {
            for (unsigned j = 0; j < k; ++j) {
                res.mat[i][j] = src.mat[r_start + i][c_start + j];
            }
        }
        return res;
    }

    static void set_block(Matrix& dest, const Matrix& src, unsigned row_block, unsigned col_block, unsigned k) {
        unsigned r_start = row_block * k;
        unsigned c_start = col_block * k;
        for (unsigned i = 0; i < k; ++i) {
            for (unsigned j = 0; j < k; ++j) {
                dest.mat[r_start + i][c_start + j] = src.mat[i][j];
            }
        }
    }
    static Matrix laderman_rec(Matrix& A, Matrix& B, unsigned& adds, unsigned& subs, unsigned& muls) {
    unsigned n = A.get_n();

    if (n == 1) {
        Matrix res(1, 1);
        res.set(0, 0, A.get(0, 0) * B.get(0, 0));
        muls++;
        return res;
    }

    unsigned k = n / 3;

    Matrix A00 = Matrix::get_block(A, 0, 0, k), A01 = Matrix::get_block(A, 0, 1, k), A02 = Matrix::get_block(A, 0, 2, k);
    Matrix A10 = Matrix::get_block(A, 1, 0, k), A11 = Matrix::get_block(A, 1, 1, k), A12 = Matrix::get_block(A, 1, 2, k);
    Matrix A20 = Matrix::get_block(A, 2, 0, k), A21 = Matrix::get_block(A, 2, 1, k), A22 = Matrix::get_block(A, 2, 2, k);

    Matrix B00 = Matrix::get_block(B, 0, 0, k), B01 = Matrix::get_block(B, 0, 1, k), B02 = Matrix::get_block(B, 0, 2, k);
    Matrix B10 = Matrix::get_block(B, 1, 0, k), B11 = Matrix::get_block(B, 1, 1, k), B12 = Matrix::get_block(B, 1, 2, k);
    Matrix B20 = Matrix::get_block(B, 2, 0, k), B21 = Matrix::get_block(B, 2, 1, k), B22 = Matrix::get_block(B, 2, 2, k);

    Matrix U1 = Matrix::add(A20, A22, adds);
    Matrix U2 = Matrix::add(A10, A11, adds);
    Matrix U3 = Matrix::add(A02, U1, adds);
    Matrix U4 = Matrix::subtract(A21, U2, subs);

    Matrix V1 = Matrix::add(B11, B21, adds);
    Matrix V2 = Matrix::subtract(B20, V1, subs);
    Matrix V3 = Matrix::add(B01, V2, adds);
    Matrix V4 = Matrix::subtract(B00, V3, subs);
    Matrix V5 = Matrix::add(B22, V1, adds);
    Matrix V6 = Matrix::subtract(B10, B12, subs);
    Matrix V7 = Matrix::subtract(B01, V5, subs);
    Matrix V8 = Matrix::subtract(V4, V6, subs);

    Matrix M1  = laderman_rec(U1, V5, adds, subs, muls);

    Matrix T2  = Matrix::add(V2, V6, adds);
    Matrix M2  = laderman_rec(U2, T2, adds, subs, muls);

    Matrix M3  = laderman_rec(A21, B12, adds, subs, muls);

    Matrix T4  = Matrix::add(B02, V7, adds);
    Matrix M4  = laderman_rec(A20, T4, adds, subs, muls);

    Matrix M5  = laderman_rec(U3, V7, adds, subs, muls);

    Matrix T6  = Matrix::subtract(A21, A22, subs);
    Matrix M6  = laderman_rec(T6, B11, adds, subs, muls);

    Matrix M7  = laderman_rec(A12, B22, adds, subs, muls);

    Matrix T8  = Matrix::subtract(U1, U2, subs);
    Matrix M8  = laderman_rec(T8, V2, adds, subs, muls);

    Matrix T9  = Matrix::subtract(A01, A02, subs);
    Matrix M9  = laderman_rec(T9, B11, adds, subs, muls);

    Matrix T10 = Matrix::subtract(U3, A10, subs);
    Matrix M10 = laderman_rec(T10, V3, adds, subs, muls);

    Matrix T11 = Matrix::add(A02, A22, adds);
    Matrix T12 = Matrix::subtract(V1, B01, subs);
    Matrix M11 = laderman_rec(T11, T12, adds, subs, muls);

    Matrix M12 = laderman_rec(A02, B22, adds, subs, muls);

    Matrix T13 = Matrix::add(A20, U4, adds);
    Matrix M13 = laderman_rec(T13, V4, adds, subs, muls);

    Matrix M14 = laderman_rec(A00, B02, adds, subs, muls);

    Matrix T15 = Matrix::add(A00, U3, adds);
    Matrix M15 = laderman_rec(T15, B01, adds, subs, muls);

    Matrix T16 = Matrix::add(A02, A12, adds);
    Matrix M16 = laderman_rec(T16, B20, adds, subs, muls);

    Matrix T17 = Matrix::subtract(A11, A21, subs);
    Matrix T18 = Matrix::subtract(V4, B10, subs);
    Matrix M17 = laderman_rec(T17, T18, adds, subs, muls);

    Matrix T19 = Matrix::add(A00, A10, adds);
    Matrix M18 = laderman_rec(T19, B00, adds, subs, muls);

    Matrix M19 = laderman_rec(A01, B12, adds, subs, muls);

    Matrix T20 = Matrix::add(A01, A11, adds);
    Matrix M20 = laderman_rec(T20, B10, adds, subs, muls);

    Matrix T21 = Matrix::subtract(B02, V8, subs);
    Matrix M21 = laderman_rec(A10, T21, adds, subs, muls);

    Matrix T22 = Matrix::subtract(A11, A12, subs);
    Matrix T23 = Matrix::subtract(B20, B21, subs);
    Matrix M22 = laderman_rec(T22, T23, adds, subs, muls);

    Matrix M23 = laderman_rec(U4, V8, adds, subs, muls);

    Matrix W1 = Matrix::add(M5, M12, adds);
    Matrix W2 = Matrix::add(M1, W1, adds);
    Matrix W3 = Matrix::add(M8, W2, adds);
    Matrix W4 = Matrix::subtract(M3, M23, subs);
    Matrix W5 = Matrix::add(M2, W4, adds);
    Matrix W6 = Matrix::add(W3, W5, adds);
    Matrix W7 = Matrix::subtract(M10, W6, subs);
    Matrix W8 = Matrix::add(M17, W7, adds);

    Matrix C00 = Matrix::add(Matrix::add(M18, M20, adds), W8, adds);

    Matrix U_C01 = Matrix::add(M9, M15, adds);
    Matrix C01 = Matrix::subtract(U_C01, W2, subs);

    Matrix C02 = Matrix::add(Matrix::add(M12, M14, adds), M19, adds);

    Matrix C10 = Matrix::subtract(M16, W8, subs);

    Matrix U_C2 = Matrix::add(M16, M22, adds);
    Matrix U_C3 = Matrix::add(U_C2, W3, adds);
    Matrix C11 = Matrix::subtract(U_C3, M10, subs);

    Matrix U_C4 = Matrix::add(Matrix::add(M7, M21, adds), W4, adds);
    Matrix C12 = Matrix::subtract(U_C4, M17, subs);

    Matrix C20 = Matrix::add(Matrix::add(M11, M13, adds), W6, adds);
    Matrix C21 = Matrix::add(Matrix::add(M6, M11, adds), W2, adds);

    Matrix U_C5 = Matrix::add(Matrix::add(M3, M4, adds), M11, adds);
    Matrix C22 = Matrix::subtract(U_C5, W1, subs);

    Matrix C(n, n);
    Matrix::set_block(C, C00, 0, 0, k);
    Matrix::set_block(C, C01, 0, 1, k);
    Matrix::set_block(C, C02, 0, 2, k);
    Matrix::set_block(C, C10, 1, 0, k);
    Matrix::set_block(C, C11, 1, 1, k);
    Matrix::set_block(C, C12, 1, 2, k);
    Matrix::set_block(C, C20, 2, 0, k);
    Matrix::set_block(C, C21, 2, 1, k);
    Matrix::set_block(C, C22, 2, 2, k);

    return C;
}
};

Matrix bine_multiply(Matrix& a, Matrix& b) {
    if (a.get_n() != a.get_m() || b.get_n() != b.get_m() || a.get_m() != b.get_n()) {
        throw std::invalid_argument("It works only for square matrices.");
    }

    unsigned n = a.get_n();

    unsigned p2 = 1;
    while (p2 < n) {
        p2 *= 2;
    }

    Matrix A_padded(p2, p2);
    Matrix B_padded(p2, p2);

    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            A_padded.mat[i][j] = a.mat[i][j];
            B_padded.mat[i][j] = b.mat[i][j];
        }
    }

    unsigned total_ops = 0;
    Matrix res_padded = Matrix::bine_sub(A_padded, 0, 0, B_padded, 0, 0, p2, total_ops);

    Matrix result(n, n);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            result.mat[i][j] = res_padded.mat[i][j];
        }
    }

    result.operation_count = total_ops;
    return result;
}

Matrix strassen_multiply(Matrix& a, Matrix& b) {
    if (a.get_m() != b.get_n() || a.get_m() != a.get_n() || b.get_m() != b.get_n()) {
        throw std::invalid_argument("Strassen works only for square matrices");
    }

    unsigned n = a.get_n();
    unsigned p2 = 1;
    while (p2 < n) {
        p2 *= 2;
    }


    Matrix A_padded(p2, p2);
    Matrix B_padded(p2, p2);

    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            A_padded.set(i, j, a.get(i, j));
            B_padded.set(i, j, b.get(i, j));
        }
    }

    unsigned total_ops = 0;

    Matrix res_padded = Matrix::strassen_sub(A_padded, 0, 0, B_padded, 0, 0, p2, total_ops);

    Matrix result(n, n);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            result.set(i, j, res_padded.get(i, j));
        }
    }

    result.operation_count = total_ops;
    return result;

}

Matrix laderman_multiply(Matrix& a, Matrix& b) {
    unsigned n = a.get_n();
    if (n != a.get_m() || b.get_n() != b.get_m() || n != b.get_n()) {
        throw std::invalid_argument("Matrices must be square and of the same dimensions.");
    }

    unsigned p = 1;
    while (p < n) p *= 3;

    Matrix A_pad(p, p);
    Matrix B_pad(p, p);

    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            A_pad.set(i, j, a.get(i, j));
            B_pad.set(i, j, b.get(i, j));
        }
    }

    unsigned adds = 0, subs = 0, muls = 0;
    Matrix C_pad = Matrix::laderman_rec(A_pad, B_pad, adds, subs, muls);

    Matrix result(n, n);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            result.set(i, j, C_pad.get(i, j));
        }
    }

    result.operation_count = adds + subs + muls; 
    return result;
}

int main() {
    std::cout << "Hello, World!" << std::endl;

    Matrix a(4,4);
    Matrix b(4,4);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            a.set(i,j,i);
            b.set(i,j,j);
        }
    }
    a.print_matrix();
    b.print_matrix();
    Matrix c1 = bine_multiply(a, b), c2 = strassen_multiply(a,b),c3 = laderman_multiply(a, b);
    c2.print_matrix();
    c3.print_matrix();

    return 0;
}