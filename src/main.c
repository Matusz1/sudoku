#include "sudoku.h"
#include <time.h>

void print_usage();
void print_timing(int msec);

int main(int argc, char **argv) {
        if (argc == 2) {
                MdSudoku md = mdsudoku_load_file(argv[1]);
                mdsudoku_print(&md);
                return 0;
        }
        if (argc != 3) {
                print_usage();
                exit(EXIT_FAILURE);
        }


        Params par = load_params_file(argv[1]);
        MdSudoku md = mdsudoku_load_file(argv[2]);
        srand(time(NULL));
        const uint np_size = 20;
        int* full_pop = malloc_check(sizeof(*full_pop)*par.restes_max*np_size*md.size);
        clock_t start = clock();

        int err = 100;
        for (uint k = 0; k < par.restes_max; ++k) {
                SudokuSolver ss = solver_alloc(&par, &md);

                int i = 0;
                int best_idx = 0;
                int reset_counter = 0;
                for (; i < par.iters_max; ++i) {
                        int new_err;
                        best_idx = solver_next_gen(&ss, &new_err);
                        if (new_err == err)
                                ++reset_counter;
                        else
                                reset_counter = 0;
                        err = new_err;
                        if (err == 0)
                                break;
                        if (reset_counter == par.iters_reset) {
                                if (ss.info.weight_const == 0)
                                        break;
                                ss.info.weight_const = 0;
                                reset_counter = 0;
                        }
                }

                if (err != 0)
                        printf("# %2u: Failed to converge after %4d iterations with error %3d.\n", k+1, i, err);
                else
                        printf("# %2u: Converged after %4d iterations.\n", k+1, i);

                if (err == 0) {
                        clock_t diff = clock() - start;
                        int msec = diff * 1000 / CLOCKS_PER_SEC;
                        print_timing(msec);
                        sudoku_print(&md, ss.pop+best_idx*md.size);
                        break;
                }

                write_best_pop(full_pop+k*md.size*np_size, np_size, &ss);
                solver_free(&ss);
        }

        if (err == 0)
                return 0;

        printf("\n# No solutions found.\n");
        printf("# Running combined population ...\n");
        fflush(stdout);

        par.pop_size = np_size*par.restes_max;
        SudokuSolver ss = {
                .wpop = malloc_check(sizeof(*full_pop)*par.restes_max*np_size*md.size),
                .pop = full_pop,
                .md = md,
                .info = par,
                .weights = malloc_check(sizeof(*ss.weights)*np_size*par.restes_max)
        };
        ss.info.weight_const = 0;

        int min_err = 100;
        solver_calc_weights(&ss, &min_err);
        int best_idx = 0;
        for (int k = 0; k < 2*par.iters_max; ++k) {
                best_idx = solver_next_gen_final(&ss, &min_err);
                if (k == par.iters_max)
                        ss.info.weight_const = 0;
                if (min_err == 0)
                        break;

        }
        printf("# Finished with error = %d.\n", min_err);
        clock_t diff = clock() - start;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        print_timing(msec);
        sudoku_print(&md, ss.pop+best_idx*md.size);

        solver_free(&ss);

	return 0;
}

void print_usage() {
        printf("Usage:\n");
        printf(" --- To solve sudoku\n");
        printf("      ./sudoku <params_fname> <sudoku_fname>\n");
        printf(" --- To write sudoku description\n");
        printf("      ./sudoku <sudoku_fname>\n");
}

void print_timing(int msec) {
        int sec = msec / 1000;
        int min = sec / 60;
        printf("# In %2d min %2d sec\n", min, sec%60);
}
