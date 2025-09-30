//The code in this file, as it currently stands, is not my property, and I stake no claim to it. It is, merely, for the sake of my enrichment.

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

// Define the fixed size for performance
constexpr int N = 10;

/**
 * @brief Solves a 10x10 linear system Ax = b using Gaussian Elimination with Partial Pivoting.
 * * This is optimized for a small, fixed size (10x10) by using a fast array representation
 * and focusing on direct, sequential access which helps the CPU's cache and
 * compiler's optimization (e.g., loop unrolling/vectorization).
 * * @param A The 10x10 matrix (input, will be modified/overwritten).
 * @param b The 10-element right-hand side vector (input/output, becomes the solution x).
 * @return true if a unique solution was found, false otherwise (singular matrix).
 */
bool solveMatrix(double A[N][N], double b[N]) {
    // Forward Elimination with Partial Pivoting
    for (int p = 0; p < N; p++) {
        // 1. Partial Pivoting: Find the row with the largest magnitude element in column p
        int max_row = p;
        for (int i = p + 1; i < N; i++) {
            if (std::abs(A[i][p]) > std::abs(A[max_row][p])) {
                max_row = i;
            }
        }

        // If the pivot is zero, the matrix is singular or near-singular.
        if (std::abs(A[max_row][p]) < 1e-10) { // Check for near-zero pivot
            return false;
        }

        // Swap rows max_row and p in both A and b
        if (max_row != p) {
            std::swap_ranges(A[p], A[p] + N, A[max_row]);
            std::swap(b[p], b[max_row]);
        }

        // 2. Elimination
        for (int i = p + 1; i < N; i++) {
            double factor = A[i][p] / A[p][p];

            // Update b[i]
            b[i] -= factor * b[p];

            // Update A[i] from column p + 1 to N - 1
            // Use pointer arithmetic for potentially faster access/cache utilization
            for (int j = p + 1; j < N; j++) {
                A[i][j] -= factor * A[p][j];
            }
            A[i][p] = 0.0; // The element is now theoretically zero
        }
    }

    // Back Substitution (A is now upper triangular)
    for (int i = N - 1; i >= 0; i--) {
        // Subtract all terms already solved
        for (int j = i + 1; j < N; j++) {
            b[i] -= A[i][j] * b[j];
        }

        // Divide by the diagonal element
        b[i] /= A[i][i];
    }

    // The solution is now stored in b
    return true;
}

void printVector(const double vec[N]) {
    std::cout << "[";
    for (int i = 0; i < N; ++i) {
        std::cout << std::setw(8) << std::fixed << std::setprecision(4) << vec[i] << (i < N - 1 ? ", " : "");
    }
    std::cout << "]" << std::endl;
}

int main() {
    // --- Example 1: Simple Solvable System ---

    // A system where the solution is known to be [1, 2, ..., 10]
    double A[N][N] = {
        {10, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {1, 10, 1, 1, 1, 1, 1, 1, 1, 1},
        // ... (rest of the matrix for a complete test)
        // For a simple test:
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, // A mostly simple matrix for demonstration
        {10, 9, 8, 7, 6, 5, 4, 3, 2, 1}, 
        // Fill the rest with values if needed, for a full 10x10 example:
        // Let's use a simpler structure for the demo:
        {2, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 2, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 2, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 2, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 2, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 2, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 2, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 2, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 2, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 2}
    };

    double b[N] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; // Example right-hand side vector

    std::cout << "--- 10x10 Matrix Solver (Gaussian Elimination) ---" << "\n";

    auto start = std::chrono::high_resolution_clock::now();

    if (solveMatrix(A, b)) {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;

        std::cout << "Matrix solved successfully!" << "\n";
        std::cout << "Solution vector x:" << "\n";
        printVector(b);
        std::cout << "Time taken: " << duration.count() << " ms" << "\n";

    } else {
        std::cout << "Error: Matrix is singular or near-singular. No unique solution found." << "\n";
    }
    
    // Note: The actual runtime for 10x10 is typically microseconds and dominated by I/O.
    // For proper speed testing, run the solver function thousands of times in a loop.

    return 0;
}
