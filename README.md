This code implements the Simplex method, an algorithm used to solve linear programming problems. 




Header Files

#include <iostream>  // For input and output operations
#include <vector>    // For using dynamic arrays (vectors)
#include <limits>    // For getting infinity and small numerical values

These libraries provide essential functionalities for input/output, dynamic arrays, and handling floating-point precision issues.



---

Constants

const double EPSILON = 1e-10; 
const double INF = numeric_limits<double>::infinity();

EPSILON is a very small number used to avoid precision errors when comparing floating-point values.

INF represents infinity, useful for handling unbounded problems in the Simplex method.



---

Struct for Linear Programming Problem

struct LinearProgrammingProblem {
    vector<vector<double>> A; // Coefficients of constraints
    vector<double> b;         // Right-hand side of constraints
    vector<double> c;         // Coefficients of objective function
};

This structure represents a linear programming problem in standard form:

A: Coefficients of the constraints (m × n matrix).

b: Right-hand side values of constraints (vector of size m).

c: Coefficients of the objective function (vector of size n).




---

Function to Print a Matrix

void printMatrix(const vector<vector<double>>& mat) {
    for (const auto& row : mat) {
        for (double elem : row) {
            cout << elem << "\t";
        }
        cout << endl;
    }
}

This function prints a matrix in tabular form.



---

Function to Read a Matrix from User Input

vector<vector<double>> readMatrix(int rows, int cols) {
    vector<vector<double>> mat(rows, vector<double>(cols));
    cout << "Enter the matrix (" << rows << "x" << cols << "):" << endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cin >> mat[i][j];
        }
    }
    return mat;
}

Prompts the user to input a matrix of size (rows × cols).

Reads values and stores them in a 2D vector mat.



---

Function to Read a Vector from User Input

vector<double> readVector(int size) {
    vector<double> vec(size);
    cout << "Enter the vector (" << size << " elements):" << endl;
    for (int i = 0; i < size; ++i) {
        cin >> vec[i];
    }
    return vec;
}

Reads a vector of a given size from user input.



---

Function to Read an LP Problem from Input

LinearProgrammingProblem readProblem() {
    LinearProgrammingProblem problem;
    int m, n;

    cout << "Enter the number of constraints (m): ";
    cin >> m;
    cout << "Enter the number of variables (n): ";
    cin >> n;

    cout << "Enter the coefficients of the constraints (A matrix):" << endl;
    problem.A = readMatrix(m, n);

    cout << "Enter the right-hand side of the constraints (b vector):" << endl;
    problem.b = readVector(m);

    cout << "Enter the coefficients of the objective function (c vector):" << endl;
    problem.c = readVector(n);

    return problem;
}

Reads a linear programming problem from user input.

Stores it in a LinearProgrammingProblem struct.



---

Simplex Method Implementation

vector<double> simplex(const LinearProgrammingProblem& problem) {
    int m = problem.A.size();
    int n = problem.A[0].size();

Extracts the number of constraints m and number of variables n.


Constructing the Initial Tableau

vector<vector<double>> tableau(m + 1, vector<double>(n + m + 1));

Creates a tableau (augmented matrix) of size (m+1) × (n+m+1).


for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
        tableau[i][j] = problem.A[i][j];
    }
    tableau[i][n + i] = 1;  // Slack variable
    tableau[i][n + m] = problem.b[i];
}

Copies constraint coefficients into the tableau.

Adds slack variables (extra columns for inequalities).

Stores the right-hand side vector b.


for (int j = 0; j < n; ++j) {
    tableau[m][j] = -problem.c[j];
}
tableau[m][n + m] = 0;

Stores the objective function in the last row (negated because Simplex maximizes).



---

Simplex Iterations

while (true) {
    int pivot_col = -1;
    for (int j = 0; j <= n + m; ++j) {
        if (tableau[m][j] < -EPSILON) {
            pivot_col = j;
            break;
        }
    }
    if (pivot_col == -1) break;

Selects the entering variable (most negative coefficient in the last row).

If there is no negative value, the solution is optimal.


int pivot_row = -1;
double min_ratio = INF;
for (int i = 0; i < m; ++i) {
    if (tableau[i][pivot_col] > EPSILON) {
        double ratio = tableau[i][n + m] / tableau[i][pivot_col];
        if (ratio < min_ratio) {
            min_ratio = ratio;
            pivot_row = i;
        }
    }
}
if (pivot_row == -1) {
    return vector<double>(); // Unbounded solution
}

Selects the leaving variable (minimum ratio test).



---

Pivoting Operation

double pivot_elem = tableau[pivot_row][pivot_col];
for (int j = 0; j <= n + m; ++j) {
    tableau[pivot_row][j] /= pivot_elem;
}

Divides the pivot row by the pivot element.


for (int i = 0; i <= m; ++i) {
    if (i != pivot_row) {
        double multiplier = tableau[i][pivot_col];
        for (int j = 0; j <= n + m; ++j) {
            tableau[i][j] -= multiplier * tableau[pivot_row][j];
        }
    }
}

Performs row operations to make all other values in the pivot column zero.



---

Extracting the Solution

vector<double> solution(n, 0);
for (int i = 0; i < m; ++i) {
    int var_index = -1;
    for (int j = 0; j < n; ++j) {
        if (abs(tableau[i][j] - 1) < EPSILON) {
            var_index = j;
            break;
        }
    }
    if (var_index != -1) {
        solution[var_index] = tableau[i][n + m];
    }
}

Extracts values of variables from the final tableau.



---

Main Function

int main() {
    LinearProgrammingProblem problem = readProblem();
    vector<double> solution = simplex(problem);

Reads the LP problem.

Solves it using the Simplex method.


cout << "Optimal solution:" << endl;
for (int i = 0; i < solution.size(); ++i) {
    cout << "x" << i + 1 << " = " << solution[i] << endl;
}

Prints the optimal values of variables.


double objective_value = 0;
for (int i = 0; i < solution.size(); ++i) {
    objective_value += solution[i] * problem.c[i];
}
cout << "Z = " << objective_value << endl;

Computes and prints the optimal objective function value.



---

Summary

This C++ program:

1. Reads a Linear Programming Problem.


2. Builds the Simplex tableau.


3. Performs iterations to find the optimal solution.


4. Extracts and prints the results.



This is a basic implementation and can be improved for efficiency and robustness.

