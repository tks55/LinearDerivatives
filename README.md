# LinearDerivatives

A C++ Program designed by [tks55](https://github.com/tks55) to fit regression lines primarily to stock data using linear algebra and other various numerical methods.

## Current Stage of Development

Alpha

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

A LinAlg class allowing for: 

- Numerically solving a system encoded by a Augmented Matrix
- Calculation of a PLU Decomposition
- Calculation of a dominant eigenvalue and its corresponding eigenvector via power iteration
- Transformation of a Matrix into an array of column Vectors
- Transformation of an array of column Vectors into a Matrix
- The application of an input function (e.g. Quadratic, Sin) on all values of a given matrix

A FunctionList class allowing for:

- The storage of functions of input and output double in a list stored by the object
- The addition of new functions to the predefined list/map
- The deletion of functions from the list/map
- The listing of all functions stored in the object
- Easy generation of functions corresponding to elements to various powers (e.g. quadratic, cubic)

A FileReader class allowing for:

- The reading of a user-defined setup file, set as params.txt
- The reading and processing of .csv files into a matrix
- The ability to apply and modify matrices derived from .csv files according to user input specified in params.txt

A Regression class allowing for:

- Linear Regression to be applied on user-defined and modified data

## Future Features

- Calculating the Reduced Row Echelon Form of a Matrix
- Calculating various other decompsitions (e.g. QR, SVD)
- Addition of new ML Techniques (e.g. K-Means Clustering, PCA) to further analyze input data
- Addition of new functions to output requested data into an output file
- Requesting stock data from a Stock API to be used in further data analysis

## Improvements/Optimizations

- Matrix Multiplication can further be improved to run in $\Theta(n^{\log_2(7)})$ using Strassen's Algorithm for Matrix Multiplication
- Usage of smart pointers, rather than raw pointers to provide a layer of abstraction and prevent direct access to computer memory, limiting issues faced
- More Comments!

## Known Errors
- The REF Function is currently unable to process certain types of singular matrices (e.g. repeated/zero rows, zero column vectors)
- THE REF, Det, and Inv Functions are unable to process some very rare edge cases (e.g. some lower triangular matrices).
- Slight errors with reading input files
