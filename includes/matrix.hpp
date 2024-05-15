#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <utility>
#include <string>

class Matrix {
    public:

        //Constructors
        Matrix(size_t num_rows = 2, size_t num_cols = 2, double starting_val = 1, bool diag = true); //default constructs a 2 x 2 identity matrix
        Matrix(size_t num_rows, size_t num_cols, double** init_array); //constructs a matrix from the given initializer array

        //Rule of Five
        Matrix(const Matrix& rhs);
        Matrix& operator=(const Matrix& rhs);
        ~Matrix();
        Matrix(Matrix&& rhs);
        Matrix& operator=(Matrix&& rhs);

        //Interface Functions
        std::string ToString() const;    
        std::pair<size_t, size_t> GetDims() const;
        void ChangeEntry(size_t row, size_t col, double val);

        //Operators
        // Matrix operator*(const double scalar);
        // Matrix operator*(const Matrix& rhs);
        // Matrix& operator*=(const double scalar);
        // Matrix& operator*=(const Matrix& rhs);

        // Matrix operator/(const double scalar);
        // Matrix& operator/=(const double scalar);

        // Matrix operator+(const Matrix& rhs);
        // Matrix& operator+=(const Matrix& rhs);

        // Matrix operator-(const Matrix& rhs);
        // Matrix& operator-=(const Matrix& rhs);

    private:
        size_t rows_;
        size_t cols_;
        double** array_;
};

#endif