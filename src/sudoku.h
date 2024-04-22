#ifndef _SUDOKU_H_
#define _SUDOKU_H_

#include <stdio.h>
#include <stdlib.h>

#define N 3
#define NN (N*N)
#define SSIZE (NN*NN)

typedef unsigned int uint;

#define error(cond, ...) {                                                          \
        if (cond) {                                                                 \
                fprintf(stderr, "# ERROR: file %s, line %u\n", __FILE__, __LINE__); \
                fprintf(stderr, __VA_ARGS__);                                       \
                exit(EXIT_FAILURE);                                                 \
        }}

#define malloc_check(size) ({                       \
        void* _ptr = malloc(size);                  \
        error(_ptr == NULL, "Failed to malloc.\n"); \
        (_ptr);                                     \
})

typedef struct MdSudoku {
        uint size;
        int brd[SSIZE];   // Filled with 0-s and original numbers
        int idxs[SSIZE];  // If -1 use brd elements, otherwise use chromosome
        int fill[NN];     // How many allels does chromosome have
        int fill_acc[NN]; // Quick index for getting the right part of chromosome
} MdSudoku;


typedef struct Params {
        uint pop_size;
        uint iters_max;
        uint iters_reset;
        uint restes_max;
        double m_prob;
        double c_prob;
        int weight_const;
} Params;

typedef struct SudokuSolver {
        Params info;
        MdSudoku md;

        int*  pop;
        int* wpop;

        double* weights;
} SudokuSolver;

MdSudoku mdsudoku_load_file(const char* fname);
void mdsudoku_print(MdSudoku* sudoku);
void sudoku_print(MdSudoku* md, int* ch);

int sudoku_get_val(MdSudoku* md, int* sdk, int i, int j);

Params load_params_file(const char* fname);

SudokuSolver solver_alloc(Params* par, MdSudoku* md);
void solver_free(SudokuSolver* ss);

int calc_row_error_val(MdSudoku* md, int* ch, int r, int val);
int calc_col_error_val(MdSudoku* md, int* ch, int r, int val);
int calc_errors(MdSudoku* md, int* ch);

int solver_calc_weights(SudokuSolver* ss, int* min_err);
int solver_next_gen(SudokuSolver* ss, int* min_err);
int solver_next_gen_final(SudokuSolver* ss, int* min_err);

void write_random_pop(int* bpop, int count, SudokuSolver* ss);
void write_best_pop(int* bpop, int count, SudokuSolver* solver);

/* Some utility */
void swap(int* a, int* b);
void permute_random(int n, int* arr);
int sorted_search(double* arr, int len, double val);

#endif /* _SUDOKU_H_ */
