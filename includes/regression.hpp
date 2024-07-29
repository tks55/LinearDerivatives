#ifndef REGRESSION_HPP
#define REGRESSION_HPP

#include "matrix.hpp"
#include "vector.hpp"

class Regression {

    public:

        //Class Functions
        static Vector LinearRegression(const Matrix& input_values, const Vector& target_vector);

    private:
};

#endif