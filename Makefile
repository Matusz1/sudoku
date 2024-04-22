CC = gcc
CFLAGS = -std=c11 -O3
LDLIBS = -lm -lgsl -lgslcblas

SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:$(SRC_DIR)/%=$(OBJ_DIR)/%.o)

$(BIN_DIR)/sudoku: obj $(OBJS)
	$(CC) $(OBJS) $(LDLIBS) -o $@

$(OBJ_DIR)/%.c.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

obj:
	mkdir -p obj

clean:
	rm obj/*.o
