#include "sudoku.h"
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include <gsl/gsl_sort.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

static inline double rand_double() {
        return (double) rand() / RAND_MAX;
}

int __board_id_from_idxs(int i, int j) {
        return N*(i/N) + (j/N);
}

/* Assumes md->brd is filled, calculates rest from that */
void __md_setup(MdSudoku* md) {
        for (uint i = 0; i < NN; ++i)
                md->fill[i] = 0;
        for (uint i = 0; i < SSIZE; ++i)
                md->idxs[i] = -1;
        md->size = 0;

        /* Count up sublengths of chromosomes */
        for (uint i = 0; i < NN; ++i) {
                for (uint j = 0; j < NN; ++j) {
                        int bid = __board_id_from_idxs(i, j);
                        int v = md->brd[i*NN+j];
                        if (v == 0) {
                                ++md->fill[bid];
                                ++md->size;
                        }
                }
        }

        md->fill_acc[0] = 0;
        for (uint i = 1; i < NN; ++i)
                md->fill_acc[i] = md->fill_acc[i-1] + md->fill[i-1];

        /* Fill in idxing table, reusing fill table, will be refilled in the process */
        for (uint i = 0; i < NN; ++i)
                md->fill[i] = 0;
        
        for (uint i = 0; i < NN; ++i) {
                for (uint j = 0; j < NN; ++j) {
                        int bid = __board_id_from_idxs(i, j);
                        int v = md->brd[i*NN+j];
                        if (v == 0) {
                                md->idxs[i*NN+j] = md->fill_acc[bid] + md->fill[bid];
                                ++md->fill[bid];
                        }
                }
        }
}

MdSudoku mdsudoku_load_file(const char* fname) {
        FILE* f = fopen(fname, "r");
        error(f == NULL, "Failed to open file %s\n", fname);

        MdSudoku md;
        int i = 0, j = 0;
        char num;
        while (fscanf(f, "%c", &num) == 1) {
                if ('0' <= num && num <= '9')
                        md.brd[i++] = num - '0';
        }
        error(i != 9*9, "Failed to load full sudoku from %s\n", fname);
        fclose(f);

        __md_setup(&md);
        return md;
}

double __factorial(int n) {
        double r = 1.0;
        for (int i = 2; i <= n; ++i)
                r *= i;
        return r;
}

void mdsudoku_print(MdSudoku* sudoku) {
        uint count = 0;
        printf("+-------+-------+-------+\n");
        for (uint i = 0; i < 9; ++i) {
                putchar('|');
                for (uint j = 0; j < 9; ++j) {
                        int val = sudoku->brd[i*9+j];
                        if (val) {
                                printf(ANSI_COLOR_MAGENTA);
                                printf(" %d", val);
                                printf(ANSI_COLOR_RESET);
                        }
                        else {
                                printf("  ");
                                ++count;
                        }
                        if (j % 3 == 2)
                                printf(" |");
                }
                putchar('\n');
                if (i % 3 == 2)
                        printf("+-------+-------+-------+\n");
        }
        printf("# EMPTY SLOTS  = %u\n", count);
        double v = 1.0;
        for (int i = 0; i < NN; ++i)
                v *= __factorial(sudoku->fill[i]);
        printf("# PERMUTATIONS = %e\n", v);
}

void sudoku_print(MdSudoku* md, int* ch) {
        printf("+-------+-------+-------+\n");
        for (uint i = 0; i < NN; ++i) {
                putchar('|');
                for (uint j = 0; j < NN; ++j) {
                        int val = sudoku_get_val(md, ch, i, j);
                        if (md->brd[i*9+j]) {
                                printf(ANSI_COLOR_MAGENTA);
                        } else {
                                int err = calc_row_error_val(md, ch, i, val);
                                err += calc_col_error_val(md, ch, j, val);
                                if (err > 0)
                                        printf(ANSI_COLOR_RED);
                        }
                        printf(" %d", val);
                        printf(ANSI_COLOR_RESET);
                        if (j % 3 == 2)
                                printf(" |");
                }
                putchar('\n');
                if (i % 3 == 2)
                        printf("+-------+-------+-------+\n");
        }
}

Params load_params_file(const char* fname) {
        FILE* f = fopen(fname, "r");
        error(f == NULL, "Failed to open file %s\n", fname);

        Params par = {
                .c_prob = 0.1,
                .m_prob = 0.1,
                .pop_size = 50,
        };
        char buff[256];
        while (fscanf(f, "%s", buff) == 1) {
                int e = 0;
                if (strcmp(buff, "pop_size") == 0)
                        e = fscanf(f, "%u", &par.pop_size);
                else if (strcmp(buff, "m_prob") == 0)
                        e = fscanf(f, "%lf", &par.m_prob);
                else if (strcmp(buff, "c_prob") == 0)
                        e = fscanf(f, "%lf", &par.c_prob);
                else if (strcmp(buff, "iters_max") == 0)
                        e = fscanf(f, "%u", &par.iters_max);
                else if (strcmp(buff, "iters_reset") == 0)
                        e = fscanf(f, "%u", &par.iters_reset);
                else if (strcmp(buff, "resets_max") == 0)
                        e = fscanf(f, "%u", &par.restes_max);
                else if (strcmp(buff, "weight_const") == 0)
                        e = fscanf(f, "%d", &par.weight_const);
                error(e != 1, "Unknown parameter %s from %s\n", buff, fname);
        }
        error(par.pop_size % 2 == 1, "Odd population size!\n");

        fclose(f);
        return par;
}

int sudoku_get_val(MdSudoku* md, int* sdk, int i, int j) {
        int idx = md->idxs[i*NN+j];
        if (idx == -1)
                return md->brd[i*NN+j];
        return sdk[idx];
}

int __find(int* arr, uint len, int val) {
        for (int i = 0; i < len; ++i)
                if (arr[i] == val)
                        return i;
        return -1;
}

int __is_already_in_box(MdSudoku* md, int val, int i, int j) {
        for (int ii = 0; ii < N; ++ii) {
                for (int jj = 0; jj < N; ++jj) {
                        int iii = i*N + ii;
                        int jjj = j*N + jj;
                        int idx = iii*NN + jjj;
                        if (md->brd[idx] == val)
                                return 0;
                }
        }
        return -1;
}

int solver_calc_weights(SudokuSolver* ss, int* min_err) {
        int min = 100000000;
        int imin = 0;
        const int weight = ss->info.weight_const;
        for (uint i = 0; i < ss->info.pop_size; ++i) {
                int v = calc_errors(&ss->md, ss->pop+i*ss->md.size);
                if (min > v) {
                        min = v;
                        imin = i;
                }
                if (v == 0)
                        break;
                ss->weights[i] = 1.0 / (v + weight);
        }
        *min_err = min;
        return imin;
}

SudokuSolver solver_alloc(Params* par, MdSudoku* md) {
        SudokuSolver ss = {
                .info = *par,
                .md = *md
        };

        ss.pop = malloc_check(sizeof(*ss.pop)*par->pop_size*md->size);
        ss.wpop = malloc_check(sizeof(*ss.wpop)*par->pop_size*md->size);
        ss.weights = malloc_check(sizeof(*ss.weights)*par->pop_size);
        
        uint pi = 0;
        /* Absolute abomination, but hey, it does the job */
        for (uint i = 0; i < N; ++i)
                for (uint j = 0; j < N; ++j)
                        for (int v = 1; v <= NN; ++v)
                                if (__is_already_in_box(md, v, i, j) == -1)
                                        ss.pop[pi++] = v;

        /*
        printf("\n# UNPERMUTED #\n");
        for (uint i = 0; i < NN; ++i) {
                for (uint j = 0; j < md->fill[i]; ++j) {
                        printf("%d ", ss.pop[j+md->fill_acc[i]]);
                }
                putchar('\n');
        }
        */

        /* Copy to everyone */
        for (uint i = 1; i < par->pop_size; ++i)
                memcpy(ss.pop+i*md->size, ss.pop, sizeof(*ss.pop)*md->size);

        /* Permute everyone */
        for (uint i = 0; i < par->pop_size; ++i)
                for (uint j = 0; j < NN; ++j)
                        permute_random(md->fill[j], ss.pop+i*md->size+md->fill_acc[j]);

        int min_err;
        solver_calc_weights(&ss, &min_err);
        error(min_err == 0, "Found solution while allocating.\n");
        /*
        printf("\n# PERMUTED #\n");
        for (uint i = 0; i < 10; ++i) {
                for (uint j = 0; j < md->fill[3]; ++j)
                        printf("%d ", ss.pop[i*md->size + md->fill_acc[3]+j]);
                putchar('\n');
        }

        printf("\n# RECHECKED #\n");
        pi = 0;
        for (uint i = 0; i < 10; ++i) {
                for (uint j = 0; j < md->size; ++j) {
                        printf("%d ", ss.pop[pi]);
                        ++pi;
                }
                putchar('\n');
        }
        */


        return ss;
}

void solver_calc_roulette(SudokuSolver* ss) {
        const uint size = ss->info.pop_size;
        for (uint i = 1; i < size; ++i)
                ss->weights[i+1] += ss->weights[i];

        const double fact = 1.0 / ss->weights[size-1];
        for (uint i = 0; i < size; ++i)
                ss->weights[i] *= fact;

        for (uint i = 1; i < size; ++i)
                ss->weights[size-i] = ss->weights[size-i-1];
        ss->weights[0] = 0.0;
}

void solver_free(SudokuSolver* ss) {
        free(ss->pop);
        free(ss->wpop);
}

void shift(int* arr, uint len) {
        int v = arr[0];
        for (uint i = 1; i < len; ++i)
                arr[i-1] = arr[i];
        if (len > 0)
                arr[len-1] = v;
}

static bool is_in(int x, const int* v, uint len) {
	for (int i = 0; i != len; ++i)
		if (v[i] == x)
			return true;
	return false;
}

#define IS_BETWEEN(x, a, b) ((a) <= (x) && (x) <= (b))

static int where_is(int v, const int* A, uint len) {
	for (int i = 0; i != len; ++i)
		if (A[i] == v)
			return i;
	return -1;
}

void cross_pmx(uint len, int* p1, int* p2, int* c1, int* c2) {
        if (len == 0)
                return;
        int a, b;
        a = rand() % len+1;
        b = rand() % len+1;
        if (b < a) swap(&a, &b);

	// Copy [a, b)
	for (int i = 0; i != len; ++i) {
		c1[i] = IS_BETWEEN(i, a, b-1) ? p1[i] : -1;
		c2[i] = IS_BETWEEN(i, a, b-1) ? p2[i] : -1;
	}

	for (int i = a; i != b; ++i) {
		int input = c2[i];
		if (is_in(input, c1, len))
			continue;
		int next = c1[i];
		while (true) {
			int idx = where_is(next, p2, len);
			if (c1[idx] == -1) {
				c1[idx] = input;
				break;
			}
			next = c1[idx];
		}
	}
	for (int i = 0; i != len; ++i)
		if (c1[i] == -1)
			c1[i] = p2[i];

	for (int i = a; i != b; ++i) {
		int input = c1[i];
		if (is_in(input, c2, len))
			continue;
		int next = c2[i];
		while (true) {
			int idx = where_is(next, p1, len);
			if (c2[idx] == -1) {
				c2[idx] = input;
				break;
			}
			next = c2[idx];
		}
	}
	for (int i = 0; i != len; ++i)
		if (c2[i] == -1)
			c2[i] = p1[i];
}

void cross_cx(uint len, int* p1, int* p2, int* c1, int* c2) {
        if (len == 0)
                return;
        for (int i = 0; i != len; ++i) {
                        c1[i] = -1;
                        c2[i] = -1;
                }
                c1[0] = p1[0];
                c2[0] = p2[0];
                int i = 0;
                while (!is_in(p2[i], c1, len)) {
                        int j = 0;
                        for (int k = 0; k != len; ++k) {
                                if (p1[k] == p2[i]) {
                                        j = k;
                                }
                        }
                        c1[j] = p1[j];
                        c2[j] = p2[j];
                        i = j;
                }

                for (int k = 0; k != len; ++k) {
                        if (c2[k] == -1) {
                                c1[k] = p2[k];
                                c2[k] = p1[k];
                        }
                }
}

void reverse(uint len, int* arr) {
        int beg = 0;
        int end = len-1;
        while ((end-beg) > 0) {
                swap(arr+beg, arr+end);
                ++beg;
                --end;
        }
}

void mutate_swap(uint len, int* arr) {
        if (len == 0)
                return;
        int p1 = rand() % len;
        int p2 = rand() % len;
        while (p1 == p2)
                p2 = rand() % len;
        swap(arr+p1, arr+p2);
}

void mutate_inverse(uint len, int* arr) {
        if (len == 0)
                return;
        int r1 = rand() % len;
        int r2 = rand() % len;
        if (r1 > r2)
                swap(&r1, &r2);
        reverse(r2-r1, arr+r1);
}

int calc_row_error_val(MdSudoku* md, int* ch, int r, int val) {
        int count = 0;
        int is_orig = 0;
        for (uint i = 0; i < NN; ++i) {
                if (sudoku_get_val(md, ch, r, i) == val) {
                        ++count;
                        if (md->idxs[r*NN+i] == -1)
                                is_orig = 1;
                }
        }
        if (count < 2)
                return 0;
        --count;
        return is_orig ? count*3 : count;
}

int calc_col_error_val(MdSudoku* md, int* ch, int c, int val) {
        int count = 0;
        int is_orig = 0;
        for (uint i = 0; i < NN; ++i) {
                if (sudoku_get_val(md, ch, i, c) == val) {
                        ++count;
                        if (md->idxs[i*NN+c] == -1)
                                is_orig = 1;
                }
        }
        if (count < 2)
                return 0;
        --count;
        return is_orig ? count*3 : count;
}

int calc_row_error(MdSudoku* md, int* ch, int r) {
        int err = 0;
        for (int i = 1; i <= NN; ++i)
                err += calc_row_error_val(md, ch, r, i);
        return err;
}

int calc_col_error(MdSudoku* md, int* ch, int c) {
        int err = 0;
        for (int i = 1; i <= NN; ++i)
                err += calc_col_error_val(md, ch, c, i);
        return err;
}

int calc_errors(MdSudoku* md, int* ch) {
        int err = 0;
        for (uint i = 0; i < NN; ++i) {
                err += calc_row_error(md, ch, i);
                err += calc_col_error(md, ch, i);
        }
        return err;
}

int solver_next_gen(SudokuSolver* ss, int* min_err) {
        solver_calc_roulette(ss);
        const uint psize = ss->info.pop_size;
        const uint csize = ss->md.size;
        double mp = ss->info.m_prob;
        double cp = ss->info.c_prob;
        for (uint i = 0; i < psize/2; ++i) {
                int p1 = sorted_search(ss->weights, psize, rand_double());
                int p2 = sorted_search(ss->weights, psize, rand_double());
                int* c1 = ss->wpop+(2*i)*csize;
                int* c2 = ss->wpop+(2*i+1)*csize;
                memcpy(c1, ss->pop+p1*csize, csize*sizeof(*ss->pop));
                memcpy(c2, ss->pop+p2*csize, csize*sizeof(*ss->pop));
                int n = rand() % NN;
                if (rand_double() < cp) {
                        cross_pmx(ss->md.fill[n],
                                 ss->pop+p1*csize+ss->md.fill_acc[n],
                                 ss->pop+p2*csize+ss->md.fill_acc[n],
                                 c1+ss->md.fill_acc[n],
                                 c2+ss->md.fill_acc[n]);
                }
                if (rand_double() < mp)
                        mutate_inverse(ss->md.fill[n], c1+ss->md.fill_acc[n]);
                if (rand_double() < mp)
                        mutate_inverse(ss->md.fill[n], c2+ss->md.fill_acc[n]);
        }
        memcpy(ss->pop, ss->wpop, sizeof(*ss->pop)*ss->info.pop_size*ss->md.size);
        int idx = solver_calc_weights(ss, min_err);
        return idx;
}

int solver_next_gen_final(SudokuSolver* ss, int* min_err) {
        solver_calc_roulette(ss);
        const uint psize = ss->info.pop_size;
        const uint csize = ss->md.size;
        double mp = ss->info.m_prob;
        double cp = ss->info.c_prob;
        for (uint i = 0; i < psize/2; ++i) {
                int p1 = sorted_search(ss->weights, psize, rand_double());
                int p2 = sorted_search(ss->weights, psize, rand_double());
                int* c1 = ss->wpop+(2*i)*csize;
                int* c2 = ss->wpop+(2*i+1)*csize;
                memcpy(c1, ss->pop+p1*csize, csize*sizeof(*ss->pop));
                memcpy(c2, ss->pop+p2*csize, csize*sizeof(*ss->pop));
                int n = rand() % NN;
                if (rand_double() < cp) {
                        for (uint bid = 1+rand() % NN; bid < NN; ++bid) {
                                memcpy(c1+ss->md.fill_acc[bid], ss->pop+p2*csize+ss->md.fill_acc[bid], sizeof(*c1)*ss->md.fill[bid]);
                                memcpy(c2+ss->md.fill_acc[bid], ss->pop+p1*csize+ss->md.fill_acc[bid], sizeof(*c2)*ss->md.fill[bid]);
                        }
                }
                if (rand_double() < mp)
                        mutate_swap(ss->md.fill[n], c1+ss->md.fill_acc[n]);
                if (rand_double() < mp)
                        mutate_swap(ss->md.fill[n], c2+ss->md.fill_acc[n]);
        }
        memcpy(ss->pop, ss->wpop, sizeof(*ss->pop)*ss->info.pop_size*ss->md.size);
        int idx = solver_calc_weights(ss, min_err);
        return idx;
}

void write_random_pop(int* bpop, int count, SudokuSolver* ss) {
        for (int i = 0; i < count; ++i)
                memcpy(bpop+i*ss->md.size, ss->pop+(rand() % ss->info.pop_size)*ss->md.size, sizeof(*bpop)*ss->md.size);
}

void write_best_pop(int* bpop, int count, SudokuSolver* ss) {
        int tmp;
        solver_calc_weights(ss, &tmp);
        size_t p[ss->info.pop_size];
        for (int i = 0; i < ss->info.pop_size; ++i) {
                p[i] = i;
        }

        static_assert(sizeof(double) != sizeof(long double), "Long double and double have different sizes\n");
        gsl_sort_index(p, ss->weights, 1, ss->info.pop_size);

        for (int i = 0; i < count; ++i)
                memcpy(bpop+i*ss->md.size, ss->pop+(p[ss->info.pop_size-1-i])*ss->md.size, sizeof(*bpop)*ss->md.size);
}

/* Some utility */

void swap(int* a, int* b) {
        int tmp = *a;
        *a = *b;
        *b = tmp;
}

void permute_random(int n, int* arr) {
        for (int i = 0; i < n-1; ++i) {
                uint idx = rand() % (n-i);
                swap(arr+i, arr+i+idx);
        }
}

int sorted_search(double* arr, int len, double val) {
        int beg = 0;
        int end = len;
        while (end - beg > 1) {
                int mid = beg + (end - beg) / 2;
                if (val >= arr[mid])
                        beg = mid;
                else
                        end = mid;
        }
        return beg;
}
