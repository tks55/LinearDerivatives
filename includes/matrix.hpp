#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <utility>
#include <string>
#include <tuple>
#include <random>
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
        Matrix T() const;
        double Det() const;
        Matrix REF() const;
        double GetEntry(size_t row, size_t col) const;
        Matrix Inv() const;
        static Matrix RandomMatrixFloat(size_t seed, size_t num_rows = 2, size_t num_cols = 2, double lower_bound = 0, double upper_bound = 10);
        static Matrix RandomMatrixInt(size_t seed, size_t num_rows = 2, size_t num_cols = 2, int lower_bound = 0, int upper_bound = 10);

        //To Implement
        Matrix RREF() const;
        //Operators
        friend Matrix operator*(const Matrix& matrix, const double scalar);
        friend Matrix operator*(const double scalar, const Matrix& matrix);
        friend Matrix operator*(const Matrix& lhs, const Matrix& rhs);
        Matrix& operator*=(const double scalar);
        Matrix& operator*=(const Matrix& rhs);

        friend Matrix operator/(const Matrix& matrix, const double scalar);
        Matrix& operator/=(const double scalar);

        Matrix operator+(const Matrix& rhs) const;
        Matrix& operator+=(const Matrix& rhs);

        Matrix operator-(const Matrix& rhs) const;
        Matrix& operator-=(const Matrix& rhs);

        bool operator==(const Matrix& rhs) const;

    private:
        //Instance Variables
        size_t rows_;
        size_t cols_;
        double** array_;

        //Constants
        const double kPrecisionDelta = 0.00000001;

        //Helper Functions
        std::tuple<Matrix, Matrix, Matrix, size_t> PTREF() const;
};

#endif