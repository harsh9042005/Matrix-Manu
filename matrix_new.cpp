#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>
#include <cmath>

using namespace std;

// Global matrix size
int M_SIZE;

// Function Prototypes
void inputMatrixFromUser(const string &filename);
void loadMatrixFromFile(const string &filename, vector<vector<double>> &mat);
void saveMatrixToFile(const string &filename, const vector<vector<double>> &mat);
void printMatrix(const vector<vector<double>> &mat);
void addMatrix(const vector<vector<double>> &mat1, const vector<vector<double>> &mat2, vector<vector<double>> &result);
void subtractMatrix(const vector<vector<double>> &mat1, const vector<vector<double>> &mat2, vector<vector<double>> &result);
void multiplyMatrix(const vector<vector<double>> &mat1, const vector<vector<double>> &mat2, vector<vector<double>> &result);
void transposeMatrix(const vector<vector<double>> &mat, vector<vector<double>> &result);
double determinantMatrix(const vector<vector<double>> &mat);
int rankMatrix(vector<vector<double>> mat);
void inverseMatrix(const vector<vector<double>> &mat, vector<vector<double>> &result);
vector<double> eigenvalues(const vector<vector<double>> &mat);

void displayInterface() {
    cout << "====================================\n";
    cout << "       MATRIX CALCULATOR \n";
    cout << "====================================\n";
    cout << "Perform operations on matrices!\n";
    cout << "------------------------------------\n";
}
int main() {
    displayInterface();
    cout << "Enter matrix size (e.g., 2 for 2x2, 3 for 3x3): ";
    cin >> M_SIZE;

    inputMatrixFromUser("matrixA.txt");
    inputMatrixFromUser("matrixB.txt");

    vector<vector<double>> matA(M_SIZE, vector<double>(M_SIZE));
    vector<vector<double>> matB(M_SIZE, vector<double>(M_SIZE));
    vector<vector<double>> result(M_SIZE, vector<double>(M_SIZE, 0));

    loadMatrixFromFile("matrixA.txt", matA);
    loadMatrixFromFile("matrixB.txt", matB);

    int choice;
    do {
        cout << "\n** Matrix Operations Menu **\n"
             << "1 - Addition\n"
             << "2 - Subtraction\n"
             << "3 - Multiplication\n"
             << "4 - Transpose (Matrix A)\n"
             << "5 - Determinant (Matrix A)\n"
             << "6 - Rank of (Matrix A)\n"
             << "7 - Inverse of (Matrix A)\n"
             << "8 - Eigenvalues of (Matrix A)\n"
             << "9 - Transpose (Matrix B)\n"
             << "10 - Determinant (Matrix B)\n"
             << "11 - Rank of (Matrix B)\n"
             << "12 - Inverse of (Matrix B)\n"
             << "13 - Eigenvalues of (Matrix B)\n"
             << "14 - Quit\n"
             << "Enter your choice: ";
        cin >> choice;

        switch (choice) {
            case 1:
                addMatrix(matA, matB, result);
                cout << "Addition Result:\n";
                printMatrix(result);
                saveMatrixToFile("result.txt", result);
                break;
            case 2:
                subtractMatrix(matA, matB, result);
                cout << "Subtraction Result:\n";
                printMatrix(result);
                saveMatrixToFile("result.txt", result);
                break;
            case 3:
                multiplyMatrix(matA, matB, result);
                cout << "Multiplication Result:\n";
                printMatrix(result);
                saveMatrixToFile("result.txt", result);
                break;
            case 4:
                transposeMatrix(matA, result);
                cout << "Transpose of Matrix A:\n";
                printMatrix(result);
                saveMatrixToFile("result.txt", result);
                break;
            case 5:
                cout << "Determinant of Matrix A: " << determinantMatrix(matA) << endl;
                break;
            case 6:
                cout << "Rank of Matrix A: " << rankMatrix(matA) << endl;
                break;
            case 7:
                inverseMatrix(matA, result);
                cout << "Inverse of Matrix A:\n";
                printMatrix(result);
                saveMatrixToFile("result.txt", result);
                break;
            case 8: {
                vector<double> eigVals = eigenvalues(matA);
                cout << "Eigenvalues of Matrix A: ";
                for (double val : eigVals) cout << val << " ";
                cout << endl;
                break;
            }
            case 9:
                transposeMatrix(matB, result);
                cout << "Transpose of Matrix B:\n";
                printMatrix(result);
                saveMatrixToFile("result.txt", result);
                break;
            case 10:
                cout << "Determinant of Matrix B: " << determinantMatrix(matB) << endl;
                break;
            case 11:
                cout << "Rank of Matrix B: " << rankMatrix(matB) << endl;
                break;
            case 12:
                inverseMatrix(matB, result);
                cout << "Inverse of Matrix B:\n";
                printMatrix(result);
                saveMatrixToFile("result.txt", result);
                break;
            case 13: {
                vector<double> eigVals = eigenvalues(matB);
                cout << "Eigenvalues of Matrix B: ";
                for (double val : eigVals) cout << val << " ";
                cout << endl;
                break;
            }
            case 14:
                cout << "Exiting program...\n"
                     << "\n********************************************\n"
                     << "\n  Thank You For Using My Matrix Calculator  \n"
                     << "\n********************************************\n";
                break;
            default:
                cout << "Invalid choice. Try again.\n";
        }
    } while (choice != 14);

    return 0;
}


// Function to take user input and store in file
void inputMatrixFromUser(const string &filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    cout << "Enter elements for " << filename << " (" << M_SIZE * M_SIZE << " elements):\n";
    for (int i = 0; i < M_SIZE; i++) {
        for (int j = 0; j < M_SIZE; j++) {
            double value;
            cin >> value;
            file << value << " ";
        }
        file << "\n";
    }
    file.close();
}

// Function to read matrix from file
void loadMatrixFromFile(const string &filename, vector<vector<double>> &mat) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error reading file: " << filename << endl;
        exit(1);
    }

    for (int i = 0; i < M_SIZE; i++) {
        for (int j = 0; j < M_SIZE; j++) {
            file >> mat[i][j];
        }
    }
    file.close();
}

// Function to save matrix to file
void saveMatrixToFile(const string &filename, const vector<vector<double>> &mat) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error writing to file: " << filename << endl;
        return;
    }

    for (int i = 0; i < M_SIZE; i++) {
        for (int j = 0; j < M_SIZE; j++) {
            file << mat[i][j] << " ";
        }
        file << "\n";
    }
    file.close();
}

// Function to print matrix
void printMatrix(const vector<vector<double>> &mat) {
    cout << "-----------------\n";
    for (int i = 0; i < M_SIZE; i++) {
        cout << "| ";
        for (int j = 0; j < M_SIZE; j++) {
            cout << setw(8) << mat[i][j] << " ";
        }
        cout << "|\n";
    }
    cout << "-----------------\n";
}

// Function for matrix addition
void addMatrix(const vector<vector<double>> &mat1, const vector<vector<double>> &mat2, vector<vector<double>> &result) {
    for (int i = 0; i < M_SIZE; i++) {
        for (int j = 0; j < M_SIZE; j++) {
            result[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
}

// Function for matrix subtraction
void subtractMatrix(const vector<vector<double>> &mat1, const vector<vector<double>> &mat2, vector<vector<double>> &result) {
    for (int i = 0; i < M_SIZE; i++) {
        for (int j = 0; j < M_SIZE; j++) {
            result[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
}

// Function for matrix multiplication
void multiplyMatrix(const vector<vector<double>> &mat1, const vector<vector<double>> &mat2, vector<vector<double>> &result) {
    for (int i = 0; i < M_SIZE; i++) {
        for (int j = 0; j < M_SIZE; j++) {
            result[i][j] = 0;
            for (int k = 0; k < M_SIZE; k++) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}

// Function to compute transpose
void transposeMatrix(const vector<vector<double>> &mat, vector<vector<double>> &result) {
    for (int i = 0; i < M_SIZE; i++) {
        for (int j = 0; j < M_SIZE; j++) {
            result[j][i] = mat[i][j];
        }
    }
}

// Function to compute determinant (only for 3x3 matrices)
double determinantMatrix(const vector<vector<double>> &mat) {
    if (M_SIZE != 3) {
        cout << "Determinant calculation is only supported for 3x3 matrices.\n";
        return 0;
    }

    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
           mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
           mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}
// Function to compute matrix rank using row echelon form
int rankMatrix(vector<vector<double>> mat) {
    int rank = 0;
    for (int row = 0; row < M_SIZE; row++) {
        bool nonZeroRow = false;
        for (int col = 0; col < M_SIZE; col++) {
            if (mat[row][col] != 0) {
                nonZeroRow = true;
                break;
            }
        }
        if (nonZeroRow) rank++;
    }
    return rank;
}

// Function to compute inverse using Gaussian elimination
void inverseMatrix(const vector<vector<double>> &mat, vector<vector<double>> &result) {
    int n = mat.size();
    vector<vector<double>> aug(n, vector<double>(2 * n));

    // Create augmented matrix [A | I]
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = mat[i][j];
        }
        aug[i][n + i] = 1;
    }

    // Gaussian elimination
    for (int i = 0; i < n; i++) {
        if (aug[i][i] == 0) {
            cerr << "Matrix is singular, inverse does not exist.\n";
            return;
        }

        double diagElement = aug[i][i];
        for (int j = 0; j < 2 * n; j++) {
            aug[i][j] /= diagElement;
        }

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = aug[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
    }

    // Extract inverse matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = aug[i][n + j];
        }
    }
}

// Function to compute eigenvalues using a simple iterative method (only works for small matrices)
vector<double> eigenvalues(const vector<vector<double>> &mat) {
    if (M_SIZE != 2) {
        cout << "Eigenvalue calculation is currently supported only for 2x2 matrices.\n";
        return {};
    }

    double a = mat[0][0], b = mat[0][1], c = mat[1][0], d = mat[1][1];
    double trace = a + d;
    double determinant = (a * d) - (b * c);
    double lambda1 = (trace + sqrt(trace * trace - 4 * determinant)) / 2;
    double lambda2 = (trace - sqrt(trace * trace - 4 * determinant)) / 2;

    return {lambda1, lambda2};
}