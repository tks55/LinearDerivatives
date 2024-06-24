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
        static Vector RandomVectorFloat(size_t seed, size_t size = 2, double lower_bound = 0, double upper_bound = 10);
        static Vector RandomVectorInt(size_t seed, size_t size = 2, int lower_bound = 0, int upper_bound = 10);
        void ChangeEntry(size_t entry, double val);
        double GetEntry(size_t entry) const;

    private:
        using Matrix::Det;
        using Matrix::REF;
        using Matrix::Inv;
        using Matrix::RandomMatrixFloat;
        using Matrix::RandomMatrixInt;
};

#endif