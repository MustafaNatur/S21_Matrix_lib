#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SUCCESS 1
#define FAILURE 0

enum {
  OK = 0,
  INCORRECT_MATRIX = 1,
  CALC_ERROR = 2,
};

#define S21_EQ_SUCCESS 1
#define S21_EQ_FAILURE 0

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

// Core
void s21_remove_matrix(matrix_t *A);
int s21_create_matrix(int rows, int columns, matrix_t *result);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);

// Helpers
void killMatrix(matrix_t *A);
int s21_check_matrix(matrix_t *A);
int size_match(matrix_t *A, matrix_t *B);
int matrixIsMultipliable(matrix_t *A, matrix_t *B);
matrix_t minor_creator(int i, int j, matrix_t *A);
