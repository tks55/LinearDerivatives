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
- The REF Function currently runs into errors if a 0-value is evaluated on the diagonal (can be resolved through the determination of a Permutation Matrix)
- The Determinant Function may run into errors originating as a result of the errors in the REF Function