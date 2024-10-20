CC = g++ 
FLAGS = -Wall -Werror -Wextra -std=c++17
MATRIX = *.cc
OBJECTS = *.o 
TEST_LIB = ./tests/*.cc
GCOV_FLAGS=-fprofile-arcs -ftest-coverage
LINK_FLAGS = -lgtest -lgtest_main -pthread 
UNAME_S := $(shell uname -s)

all: test 
	rm -rf ./test
	s21_matrix_oop.a

s21_matrix_oop.a:
	$(CC) $(FLAGS) -c $(MATRIX)
	ar rcs s21_matrix_oop.a $(OBJECTS)
	rm -rf *.o

test: clean
	$(CC) $(FLAGS) $(MATRIX) ${TEST_LIB} ${LINK_FLAGS} -o test
	./test

leaks: test
ifeq ($(UNAME_S), Linux) 
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./test
endif
ifeq ($(UNAME_S), Darwin) 
	leaks -quiet -atExit -- ./test
endif

clean: 
	rm -rf *.o *.a *.gch *.gcno *.gcda gcov_test s21_test.info ./report ./test

rebuild: clean s21_matrix_oop.a 

gcov_report: 
	$(CC) $(FLAGS) $(GCOV_FLAGS) ${TEST_LIB} ${MATRIX} ${LINK_FLAGS} -o test
	./test
	lcov -t "stest" -o s21_test.info -c -d . --rc lcov_branch_coverage=1 --ignore-errors mismatch
	genhtml -o report s21_test.info
	open report/index.html	

style:
	clang-format -style=google -n tests/* ./*.cc ./*.h

format:
	clang-format -style=google -i tests/* ./*.cc ./*.h
