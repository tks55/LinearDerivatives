#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "matrix.hpp"

class Vector : public Matrix {
    public:
        //Constructors
        Vector(size_t size = 2, double starting_val = 1);
        Vector(size_t size, double* init_array);

        //Member Functions
        double Dot(const Vector& rhs);
    private:
        //N/A
};

#endif