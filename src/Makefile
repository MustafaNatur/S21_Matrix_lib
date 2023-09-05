#!/bin/bash
FLAGS=-Wall -Wextra -Werror
CC=gcc
CHECKFLAGS=$(shell pkg-config --cflags --libs check)
LIBNAME=s21_matrix.a
GLFLAGS=--coverage
GCOVFLAGS=-fprofile-arcs -ftest-coverage
sourceTestsFileName=tests/genTests.check
resultTestsFileName=tests/s21_matrix_test.c
buildFolderName=matrix/s21_*.c

ifeq ($(shell uname -s), Linux)
	LFLAG += -lpthread -lrt -lm -lsubunit
endif

all: s21_matrix.a

s21_matrix.a:
	$(CC) $(FLAGS) -c -g $(buildFolderName)
	ar rc s21_matrix.a s21_*.o
	ranlib s21_matrix.a

test: clean
	$(CC) $(FLAGS) $(GLFLAGS) -c -g $(buildFolderName)
	ar rc s21_matrix.a s21_*.o
	ranlib s21_matrix.a
	$(CC) $(CHECKFLAGS) $(GLFLAGS) -g tests/*.c -o test s21_matrix.a -lcheck $(LFLAG)
	./test

leak: rebuild
	$(CC) $(CHECKFLAGS) -g tests/*.c -o test s21_matrix.a -lcheck $(LFLAG)
	CK_FORK=no leaks -quiet -atExit -- ./test
	make clean

gcov_report: test
	lcov -o tests.info -c -d .
	genhtml -o report tests.info
	open report/index.html

style: 
	clang-format -style=google -n $(buildFolderName) $(resultTestsFileName)

clean:
	rm -f *.gcno *.gcda *.info *.gcov 
	rm -f *.o *.a
	rm -rf report *.dSYM
	rm -f *.out
	rm -rf test
	rm -f gcov_report

rebuild: clean all