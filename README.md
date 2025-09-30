A simple matrix solver, written with the intention of optimizing for speed.

The code relies on compiler optimizations. I used the Clang/LLVM compiler.
The script to run my code is:
  "clang++ -std=c++20 -O3 -march=native <your_file_name>.cpp -o matrix_solver"
