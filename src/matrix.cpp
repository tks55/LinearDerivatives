#include "matrix.hpp"
#include <stdexcept>
#include <cmath>

Matrix::Matrix(size_t num_rows, size_t num_cols, double starting_val, bool diag): rows_(num_rows), cols_(num_cols) {
    if ((rows_ <= 0) || (cols_ <= 0)) {
        throw std::invalid_argument("INVALID DIMENSION MATRIX!");
    }
    array_ = new double* [rows_];

    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        array_[curr_row] = new double[cols_];
    }

    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < cols_; curr_col++) {
            if (diag == false) {
                array_[curr_row][curr_col] = starting_val;
            } else {
                if (curr_row == curr_col) {
                    array_[curr_row][curr_col] = starting_val;
                } else {
                    array_[curr_row][curr_col] = 0;
                }
            }
        }
    }
}

Matrix::Matrix(size_t num_rows, size_t num_cols, double** init_array): rows_(num_rows), cols_(num_cols) {
    if ((rows_ <= 0) || (cols_ <= 0) || (init_array == nullptr)) {
        throw std::invalid_argument("INVALID DIMENSION MATRIX!");
    }
    array_ = init_array;
}

Matrix::~Matrix() {
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        delete[] array_[curr_row];
    }
    delete[] array_;
    array_ = nullptr;
}

Matrix::Matrix(const Matrix& rhs) {
    rows_ = rhs.rows_;
    cols_ = rhs.cols_;
    array_ = new double*[rows_];
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        array_[curr_row] = new double[cols_];
    }
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < cols_; curr_col++) {
            array_[curr_row][curr_col] = rhs.array_[curr_row][curr_col];
        }
    }
}

Matrix& Matrix::operator=(const Matrix& rhs) {
    if (this == &rhs) {
        return (*this);
    }
    if (array_ != nullptr) {
        for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
            delete[] array_[curr_row];
        }
        delete[] array_;
    }
    rows_ = rhs.rows_;
    cols_ = rhs.cols_;
    array_ = new double*[rows_];
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        array_[curr_row] = new double[cols_];
    }
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < cols_; curr_col++) {
            array_[curr_row][curr_col] = rhs.array_[curr_row][curr_col];
        }
    }
    return (*this);
}

Matrix::Matrix(Matrix&& rhs) {
    rows_ = rhs.rows_;
    cols_ = rhs.cols_;
    array_ = rhs.array_;

    rhs.rows_ = 0;
    rhs.cols_ = 0;
    rhs.array_ = nullptr;
}

Matrix& Matrix::operator=(Matrix&& rhs) {
    if (this == &rhs) {
        return (*this);
    }
    if (array_ != nullptr) {
        for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
            delete[] array_[curr_row];
        }
        delete[] array_;
    }
    rows_ = rhs.rows_;
    cols_ = rhs.cols_;
    array_ = rhs.array_;

    rhs.rows_ = 0;
    rhs.cols_ = 0;
    rhs.array_ = nullptr;
    return (*this);
}

std::string Matrix::ToString() const {
    std::string output = "";
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < cols_; curr_col++) {
            output = output + std::to_string(array_[curr_row][curr_col]);
            if (curr_col != cols_ - 1) {
                output = output + " ";
            }
        }
        if (curr_row != rows_ - 1) {
            output = output + "\n";
        }
    }
    return (output);
}

std::pair<size_t, size_t> Matrix::GetDims() const {
    std::pair<double, double> dims{rows_, cols_};
    return (dims);
}

void Matrix::ChangeEntry(size_t row, size_t col, double val) {
    if (((row < 0) || (row >= rows_)) || ((col < 0) || (col >= cols_))) {
        throw std::invalid_argument("OUT OF BOUNDS ERROR!");
    }
    array_[row][col] = val;
}

Matrix operator*(const Matrix& matrix, const double scalar) {
    Matrix new_matrix = Matrix(matrix.rows_, matrix.cols_);
    for (size_t curr_row = 0; curr_row < matrix.rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < matrix.cols_; curr_col++) {
            new_matrix.array_[curr_row][curr_col] = (matrix.array_[curr_row][curr_col] * scalar);
        }
    }
    return (new_matrix);
}

Matrix operator*(const double scalar, const Matrix& matrix) {
    Matrix new_matrix = matrix * scalar;
    return (new_matrix);
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
    if (lhs.cols_ != rhs.rows_) {
        throw std::invalid_argument("INVALID MATRIX MULTIPLICATION DIMENSIONS!");
    }
    Matrix new_matrix = Matrix(lhs.rows_, rhs.cols_, 0, false);
    for (size_t curr_col = 0; curr_col < lhs.cols_; curr_col++) {
        for (size_t curr_row = 0; curr_row < lhs.rows_; curr_row++) {
            for (size_t int_col = 0; int_col < rhs.cols_; int_col++) {
                new_matrix.array_[curr_row][int_col] += lhs.array_[curr_row][curr_col] * rhs.array_[curr_col][int_col];
            }
        }
    }
    return (new_matrix);
}

Matrix& Matrix::operator*=(const double scalar) {
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < cols_; curr_col++) {
            array_[curr_row][curr_col] = (array_[curr_row][curr_col] * scalar);
        }
    }
    return (*this);
}

Matrix& Matrix::operator*=(const Matrix& rhs) {
    Matrix new_matrix = *this * rhs;
    *this = std::move(new_matrix);
    return (*this);
}

Matrix operator/(const Matrix& matrix, const double scalar) {
    if (scalar == 0) {
        throw new std::invalid_argument("DIVIDE BY ZERO ERROR!");
    }
    Matrix new_matrix = matrix * (1 / scalar);
    return (new_matrix);
}

Matrix& Matrix::operator/=(const double scalar) {
    if (scalar == 0) {
        throw new std::invalid_argument("DIVIDE BY ZERO ERROR!");
    }
    *this = (*this * (1 / scalar));
    return (*this);
}

Matrix Matrix::operator+(const Matrix& rhs) const {
    if ((rows_ != rhs.rows_) || (cols_ != rhs.cols_)) {
        throw std::invalid_argument("INVALID MATRIX ADDITION/SUBTRACTION DIMENSIONS!");
    }
    Matrix new_matrix = Matrix(*this);
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < cols_; curr_col++) {
            new_matrix.array_[curr_row][curr_col] += rhs.array_[curr_row][curr_col];
        }
    }
    return (new_matrix);
}

Matrix& Matrix::operator+=(const Matrix& rhs)  {
    if ((rows_ != rhs.rows_) || (cols_ != rhs.cols_)) {
        throw std::invalid_argument("INVALID MATRIX ADDITION/SUBTRACTION DIMENSIONS!");
    }
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < cols_; curr_col++) {
            array_[curr_row][curr_col] += rhs.array_[curr_row][curr_col];
        }
    }
    return (*this);
}

Matrix Matrix::operator-(const Matrix& rhs) const {
    Matrix new_matrix = Matrix(*this);
    new_matrix += (rhs * -1.0);
    return (new_matrix);
}

Matrix& Matrix::operator-=(const Matrix& rhs) {
    *this += (rhs * - 1.0);
    return (*this);
}

bool Matrix::operator==(const Matrix& rhs) const {
    if ((rows_ != rhs.rows_) || (cols_ != rhs.cols_)) {
        return (false);
    }
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < cols_; curr_col++) {
            if (std::abs(array_[curr_row][curr_col] - rhs.array_[curr_row][curr_col]) > kPrecisionDelta) {
                return (false);
            }
        }
    }
    return (true);
}

Matrix Matrix::T() const {
    Matrix new_matrix = Matrix(cols_, rows_, 0, false);
    for (size_t curr_row = 0; curr_row < rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < cols_; curr_col++) {
            new_matrix.array_[curr_col][curr_row] = array_[curr_row][curr_col];
        }
    }
    return (new_matrix);
}

double Matrix::Det() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("CANNOT DETERMINE DETERMINANT--NONSQUARE MATRIX!");
    } else if (rows_ == 1) {
        return (array_[0][0]);
    } else if (rows_ == 2) {
        return (array_[0][0] * array_[1][1]) - (array_[0][1] * array_[1][0]);
    }
    double determinant = 1.0;
    std::tuple<Matrix, Matrix, Matrix, size_t> matrix_tuple = this->PTREF();
    Matrix permute_matrix = std::get<0>(matrix_tuple);
    Matrix ref_matrix = std::get<2>(matrix_tuple);
    size_t row_swaps = std::get<3>(matrix_tuple);
    for (size_t diag = 0; diag < ref_matrix.rows_; diag++) {
        determinant *= ref_matrix.array_[diag][diag];
        if ((determinant == 0) || (std::isnan(determinant))) {
            return (0);
        }
    }
    if (row_swaps % 2 == 1) {
        determinant *= -1;
    }
    return (determinant);
}

std::tuple<Matrix, Matrix, Matrix, size_t> Matrix::PTREF() const {
    Matrix permutation_matrix = Matrix(rows_, rows_, 1, true);
    Matrix lt_matrix = Matrix(rows_, rows_, 1, true);
    Matrix ut_matrix = Matrix(*this);
    size_t row_exchanges = 0;
    for (size_t root_row = 0; root_row < ut_matrix.rows_ - 1; root_row++) {
        Matrix epermutation_matrix = Matrix(rows_, rows_, 1, true);
        size_t max_row = root_row;
        size_t max_val = std::abs(ut_matrix.array_[root_row][root_row]);
        for (size_t curr_row = root_row + 1; curr_row < ut_matrix.rows_; curr_row++) {
            double curr_val = ut_matrix.array_[curr_row][root_row];
            if (std::abs(curr_val) > max_val) {
                max_row = curr_row;
                max_val = std::abs(curr_val);
            }
        }
        if (max_row != root_row) {
            epermutation_matrix.array_[root_row][root_row] = 0;
            epermutation_matrix.array_[max_row][max_row] = 0;
            epermutation_matrix.array_[root_row][max_row] = 1;
            epermutation_matrix.array_[max_row][root_row] = 1;
            row_exchanges++;
        }
        permutation_matrix = epermutation_matrix * permutation_matrix;
        ut_matrix = epermutation_matrix * ut_matrix;
    }
    for (size_t root_row = 0; root_row < ut_matrix.rows_ - 1; root_row++) {
        for (size_t curr_row = root_row + 1; curr_row < ut_matrix.rows_; curr_row++) {
            double curr_val = ut_matrix.array_[curr_row][root_row];
            double curr_div = ut_matrix.array_[root_row][root_row];
            double curr_scalar = (curr_val / curr_div);
            lt_matrix.array_[curr_row][root_row] += curr_scalar;
            for (size_t curr_col = 0; curr_col < ut_matrix.cols_; curr_col++) {
                ut_matrix.array_[curr_row][curr_col] -= (curr_scalar * ut_matrix.array_[root_row][curr_col]); 
            }
        }
    }
    std::tuple<Matrix, Matrix, Matrix, size_t> matrix_tuple{permutation_matrix, lt_matrix, ut_matrix, row_exchanges};
    return (matrix_tuple);
}

Matrix Matrix::Inv() const {
    std::tuple<Matrix, Matrix, Matrix, double> matrix_tuple = this->PTREF();
    Matrix permutation_matrix = std::get<0>(matrix_tuple);
    Matrix lt_matrix = std::get<1>(matrix_tuple);
    Matrix ut_matrix = std::get<2>(matrix_tuple);
    Matrix ut_matrix_inverse = Matrix(rows_, rows_, 1, true);
    Matrix lt_matrix_inverse = Matrix(rows_, rows_, 1, true);
    for (size_t root_row = 0; root_row < rows_ - 1; root_row++) {
        size_t root_row_inv = ut_matrix.rows_ - root_row - 1;
        for (size_t curr_row = root_row + 1; curr_row < rows_; curr_row++) {
            size_t curr_row_inv = ut_matrix.rows_ - curr_row - 1;
            double curr_val_lt = lt_matrix.array_[curr_row][root_row];
            double curr_div_lt = lt_matrix.array_[root_row][root_row];
            double curr_val_ut = ut_matrix.array_[curr_row_inv][root_row_inv];
            double curr_div_ut = ut_matrix.array_[root_row_inv][root_row_inv];
            if ((curr_div_ut == 0 || std::isnan(curr_div_ut)) || (curr_div_lt == 0 || std::isnan(curr_div_lt))) {
                throw (std::invalid_argument("CANNOT CALCULATE INVERSE--SINGULAR MATRIX!"));
            }
            double curr_scalar_lt = (curr_val_lt / curr_div_lt);
            double curr_scalar_ut = (curr_val_ut / curr_div_ut);
            for (size_t curr_col = 0; curr_col < cols_; curr_col++) {
                lt_matrix.array_[curr_row][curr_col] -= (curr_scalar_lt * lt_matrix.array_[root_row][curr_col]); 
                ut_matrix.array_[curr_row_inv][curr_col] -= (curr_scalar_ut * ut_matrix.array_[root_row_inv][curr_col]); 
                lt_matrix_inverse.array_[curr_row][curr_col] -= (curr_scalar_lt * lt_matrix_inverse.array_[root_row][curr_col]);
                ut_matrix_inverse.array_[curr_row_inv][curr_col] -= (curr_scalar_ut * ut_matrix_inverse.array_[root_row_inv][curr_col]);
            }
        }
    }
    for (size_t diag = 0; diag < ut_matrix.rows_; diag++) {
        double inv_scalar_lt = (1 / lt_matrix.array_[diag][diag]);
        double inv_scalar_ut = (1 / ut_matrix.array_[diag][diag]);
        for (size_t curr_col = 0; curr_col < ut_matrix.cols_; curr_col++) {
            lt_matrix_inverse.array_[diag][curr_col] *= inv_scalar_lt;
            ut_matrix_inverse.array_[diag][curr_col] *= inv_scalar_ut;
        }
    }
    Matrix inverse_matrix = ut_matrix_inverse * lt_matrix_inverse * permutation_matrix;
    return (inverse_matrix);
}

Matrix Matrix::REF() const {
    std::tuple<Matrix, Matrix, Matrix, size_t> matrix_tuple = this->PTREF();
    Matrix ref_matrix = std::get<2>(matrix_tuple);
    for (size_t curr_row = 0; curr_row < ref_matrix.rows_; curr_row++) {
        for (size_t curr_col = 0; curr_col < ref_matrix.cols_; curr_col++) {
            if (std::isnan(ref_matrix.array_[curr_row][curr_col])) {
                throw std::invalid_argument("ZERO COLUMN/ROW VECTOR DETECTED--NOT CURRENTLY ABLE TO DETERMINE REF! (MATRIX IS SINGULAR)");
            }
        }
    }
    return (ref_matrix);
}

double Matrix::GetEntry(size_t row, size_t col) const {
    double entry_value = array_[row][col];
    return (entry_value);
}