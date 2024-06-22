#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "matrix.hpp"

class Vector : public Matrix {
    public:
        //Constructors
        Vector(size_t size = 2, double starting_val = 1);
        Vector(size_t size, double* init_array);
        Vector(const Matrix& rhs);

        //Member Functions
        double Dot(const Vector& rhs);
        static double Dot(const Vector& lhs, const Vector& rhs);
        Vector Norm();
        static Vector Norm(const Vector& vect);
        Vector ProjectOnto(const Vector& rhs);
        static Vector ProjectOnto(const Vector& lhs, const Vector& rhs);
        double Magnitude();
        static double Magnitude(const Vector& vect);
    private:
        //N/A
};

#endif