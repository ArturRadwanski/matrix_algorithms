//
// Created by valleron on 10/29/26.
//
#pragma once
#ifndef ALGORYTMY_MACIERZOWE_CONSTMATRIXVIEW_H
#define ALGORYTMY_MACIERZOWE_CONSTMATRIXVIEW_H
#include "Matrix.h"

class ConstMatrixView {
public:
    ConstMatrixView(const Matrix& mat)
        : matrix(mat), r(0), c(0), n(mat.n_columns()), m(mat.n_rows()){}

    ConstMatrixView(const Matrix& mat, unsigned int r_start, unsigned int c_start, unsigned int width, unsigned int height)
        : matrix(mat), r(r_start), c(c_start), n(width), m(height) {}

    double operator()(unsigned int i, unsigned int j) const {
        return matrix(r + i, c + j);
    }

    [[nodiscard]] ConstMatrixView sub(unsigned int dr, unsigned int dc, unsigned int new_width, unsigned int new_height) const {
        return ConstMatrixView{matrix, r + dr, c + dc, new_width, new_height};
    }

    [[nodiscard]] unsigned int n_cols() const { return n; }
    [[nodiscard]] unsigned int n_rows() const { return m; }
    Matrix operator+(const ConstMatrixView& other) const {
        if (this->n_cols() != other.n_cols() || this->n_rows() != other.n_rows()) {
            throw std::invalid_argument("ConstMatrixView::operator+(const ConstMatrixView& other) cant add matrices with different dimensions");
        }

        auto result = Matrix(this -> n_cols(), this -> n_rows());
        for (unsigned int i = 0; i < n_cols(); i++) {
            for (unsigned int j = 0; j < n_rows(); j++) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return result;
    }


    Matrix operator-(const ConstMatrixView& other) const {
        if (this->n_cols() != other.n_cols() || this->n_rows() != other.n_rows()) {
            throw std::invalid_argument("ConstMatrixView::operator+(const ConstMatrixView& other) cant subtract matrices with different dimensions");
        }

        auto result = Matrix(this -> n_cols(), this -> n_rows());
        for (unsigned int i = 0; i < n_cols(); i++) {
            for (unsigned int j = 0; j < n_rows(); j++) {
                result(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return result;
    }

private:
    const Matrix &matrix;
    unsigned int r, c, n, m;
};
#endif //ALGORYTMY_MACIERZOWE_CONSTMATRIXVIEW_H