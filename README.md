# LinearDerivatives

A C++ Program designed by [tks55](https://github.com/tks55) to fit regression lines primarily to stock data using linear algebra and other various numerical methods.

## Current Stage of Development

Pre-Alpha

## Current Features

A dynamic Matrix class that allows for:

- The initalization of a Matrix with various constructors (e.g. default, 2D Array)
- Matrix Addition/Subtraction
- Matrix-Scalar Multiplication
- Matrix-Vector/Matrix-Matrix Multiplication (currently in $\Theta(n^3)$ runtime)
- Calculating the Transpose of a Matrix
- Determining the Determinant of a Matrix
- Determining a Row Echelon Form of a Matrix
- Determining the Inverse of a Matrix
- Generation of a Random Matrix (given a seed value) of either Floats or Integers

A dynamic Vector class allowing for:

- The initalization of a Vector with various constructors (e.g. default, 1D Array)
- Calculating the Vector Dot Product
- Calculating the Projection of a Vector onto another Vector
- Calculating a Normalized/Unitized Vector
- Generation of a Random Vector (given a seed value) of either Floats or Integers

## Future Features

- Calculating the Reduced Row Echelon Form of a Matrix
- Numerically solving a system encoded by a Augmented Matrix
- Calculating various decompsitions (e.g. PLU, SVD)
- Computationally crafting Design Matrices based on User Input
- Solving the Normal Equations of Design Matrices to determine Least-Squared Solutions according to user specification
- Addition of new constructors for the Matrix class to input user data via an input file.
- Addition of new functions to output requested data into an output file
- Requesting stock data from a Stock API to be used in further data analysis

## Improvements/Optimizations

- Matrix Multiplication can further be improved to run in $\Theta(n^{\log_2(7)})$ using Strassen's Algorithm for Matrix Multiplication
- Usage of smart pointers, rather than raw pointers to provide a layer of abstraction and prevent direct access to computer memory, limiting issues faced
- More Comments!

## Known Errors
- The REF Function is currently unable to process certain types of singular matrices (e.g. repeated/zero rows, zero column vectors)
- THE REF, Det, and Inv Functions are unable to process some very rare edge cases (e.g. some lower triangular matrices).
