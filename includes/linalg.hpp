#ifndef LINALG_HPP
#define LINALG_HPP

#include "matrix.hpp"
#include "vector.hpp"
#include <tuple>
#include <utility>
#include <functional>

class LinAlg {

    public:

        //Class Functions
        static Vector Solve(const Matrix& a, const Vector& b);
        static std::tuple<Matrix, Matrix, Matrix> PLUDecomp(const Matrix& matrix);
        static std::pair<double, Vector> PowerIter(const Matrix& matrix, size_t max_iter = 10000);
        static Vector* ToColumnVectors(const Matrix& matrix);
        static Matrix ToMatrix(Vector* column_vectors, size_t num_vectors);
        static Matrix ApplyFunction(const Matrix& matrix, std::function<double(double)> function);

        //To Implement
        static std::pair<Matrix, Matrix> QRDecomp(const Matrix& matrix);

    private:

};

#endif