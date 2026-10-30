//
// Created by valleron on 10/29/26.
//
#pragma once
#ifndef ALGORYTMY_MACIERZOWE_MATRIXVIEW_H
#define ALGORYTMY_MACIERZOWE_MATRIXVIEW_H
#include "Matrix.h"

class MatrixView {
    public:
    MatrixView(Matrix& mat)
        : matrix(mat), r(0), c(0), n(mat.n_columns()), m(mat.n_rows()) {}
    MatrixView(Matrix& mat, unsigned int r_start, unsigned int c_start, unsigned int width, unsigned int height)
        : matrix(mat), r(r_start), c(c_start), n(width), m(height) {}

    double& operator()(unsigned int i, unsigned int j) {
        return matrix(r + i, c + j);
    }

    MatrixView sub(unsigned int dr, unsigned int dc, unsigned int new_width, unsigned int new_height) {
        return MatrixView{matrix, r + dr, c + dc, new_width, new_height};
    }

    unsigned int n_cols() const { return n; }
    unsigned int n_rows() const { return m; }

    void write_into(const Matrix& mat) {
        if (mat.n_columns() != n_cols() || mat.n_rows() != n_rows()) {
            throw std::invalid_argument("MatrixView::write_into cant write into matrices with different dimensions");
        }

        for (unsigned int i = 0; i < n_rows(); ++i) {
            for (unsigned int j = 0; j < n_cols(); ++j) {
                (*this)(i, j) = mat(i, j);
            }
        }
    }
    void add_flops(const int count) const {
        matrix.flops+=count;
    }
private:
    Matrix& matrix;
    unsigned int r, c, n, m;
};


#endif //ALGORYTMY_MACIERZOWE_MATRIXVIEW_H