#include "../s21_matrix.h"

// Create matrix
int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = OK;
  if (rows > 0 && columns > 0 && result) {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (!result->matrix) {
      status = INCORRECT_MATRIX;
    }

    for (int i = 0; i < rows; i++) {
      if (!(result->matrix[i] = (double *)calloc(columns, sizeof(double)))) {
        s21_remove_matrix(result);
        status = INCORRECT_MATRIX;
        break;
      }
    }
  } else {
    status = INCORRECT_MATRIX;
  }

  return status;
}

// Remove matrix
void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
  }
  killMatrix(A);
  A = NULL;
}

// Eq
int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = SUCCESS;
  if (s21_check_matrix(A) == 0 && s21_check_matrix(B) == 0 &&
      A->columns == B->columns && A->rows == B->rows) {
    for (int i = 0; i < A->rows && status; i++) {
      for (int j = 0; j < A->columns && status; j++) {
        if (round(A->matrix[i][j] * pow(10, 7)) !=
            round(B->matrix[i][j] * pow(10, 7))) {
          status = FAILURE;
        }
      }
    }
  } else {
    status = FAILURE;
  }
  return status;
}

// Sub
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  killMatrix(result);
  if (!s21_check_matrix(A) && !s21_check_matrix(B) && size_match(A, B)) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - (B->matrix[i][j]);
      }
    }
  } else {
    status = size_match(A, B) ? INCORRECT_MATRIX : CALC_ERROR;
  }
  return status;
}

// Sum
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  killMatrix(result);
  if (!s21_check_matrix(A) && !s21_check_matrix(B) && size_match(A, B)) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + (B->matrix[i][j]);
      }
    }
  } else {
    status = size_match(A, B) ? INCORRECT_MATRIX : CALC_ERROR;
  }
  return status;
}

// Mult num
int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int status = OK;
  killMatrix(result);
  if (!s21_check_matrix(A)) {
    if (!s21_create_matrix(A->rows, A->columns, result)) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    }
  } else {
    status = INCORRECT_MATRIX;
  }
  return status;
}

// Mult matrix
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  if (!s21_check_matrix(A) && !s21_check_matrix(B) &&
      matrixIsMultipliable(A, B) &&
      !s21_create_matrix(A->rows, B->columns, result)) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        result->matrix[i][j] = 0;
        for (int k = 0; k < B->rows; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  } else {
    status = matrixIsMultipliable(A, B) ? INCORRECT_MATRIX : CALC_ERROR;
  }
  return status;
}

// Transpose
int s21_transpose(matrix_t *A, matrix_t *result) {
  int status = OK;
  if (!s21_check_matrix(A) && !s21_create_matrix(A->columns, A->rows, result)) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  } else {
    status = INCORRECT_MATRIX;
  }
  return status;
}

// Determination
int s21_determinant(matrix_t *A, double *result) {
  int res = OK;
  if (!result) {
    res = CALC_ERROR;
  } else if (s21_check_matrix(A) != 0) {
    res = INCORRECT_MATRIX;
  } else if (A->columns != A->rows) {
    res = CALC_ERROR;
  } else {
    *result = 0;
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else if (A->rows == 2) {
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    } else {
      for (int j = 0; j < A->columns; j++) {
        matrix_t minor = minor_creator(0, j, A);
        double temp = 0;
        s21_determinant(&minor, &temp);
        *result += pow(-1.0, j) * temp * A->matrix[0][j];
        s21_remove_matrix(&minor);
      }
    }
  }
  return res;
}

// Complements
int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int status = OK;
  if (!s21_check_matrix(A) && (A->rows > 1 && A->columns > 1) &&
      (A->rows == A->columns) &&
      !s21_create_matrix(A->rows, A->columns, result)) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        matrix_t minor = minor_creator(i, j, A);
        double temp = 0;
        s21_determinant(&minor, &temp);
        result->matrix[i][j] = pow(-1, (j + i)) * temp;
        s21_remove_matrix(&minor);
      }
    }
  } else {
    if (s21_check_matrix(A) == 1 ||
        s21_create_matrix(A->rows, A->columns, result)) {
      status = INCORRECT_MATRIX;
    } else if (!((A->rows > 1 && A->columns > 1) && (A->rows == A->columns))) {
      status = CALC_ERROR;
    }
  }
  return status;
}

// Inverse
int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = OK;
  double determinant = 0;
  matrix_t calc_complements;
  matrix_t transpont;
  if (!result) {
    res = CALC_ERROR;
  } else if (s21_check_matrix(A) != 0) {
    res = INCORRECT_MATRIX;
  } else if (A->columns != A->rows) {
    res = CALC_ERROR;
  } else {
    killMatrix(result);
    res = s21_determinant(A, &determinant);
    if (!determinant) {
      res = CALC_ERROR;
    } else if (!res) {
      res = s21_calc_complements(A, &calc_complements);
      if (!res) {
        res = s21_transpose(&calc_complements, &transpont);
        if (!res) {
          res = s21_mult_number(&transpont, 1 / determinant, result);
        } else {
          res = CALC_ERROR;
        }
        s21_remove_matrix(&calc_complements);
      } else {
        res = CALC_ERROR;
      }
      s21_remove_matrix(&transpont);
    }
  }
  return res;
}

// Helpers
void killMatrix(matrix_t *A) {
  A->matrix = NULL;
  A->columns = 0;
  A->rows = 0;
}

int s21_check_matrix(matrix_t *A) {
  return (A && A->rows > 0 && A->columns > 0 && A->matrix) ? 0 : 1;
}

int size_match(matrix_t *A, matrix_t *B) {
  return A->columns == B->columns && A->rows == B->rows ? 1 : 0;
}

int matrixIsMultipliable(matrix_t *A, matrix_t *B) {
  return A->columns == B->rows;
}

matrix_t minor_creator(int i, int j, matrix_t *A) {
  matrix_t minor;
  s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
  int n = 0;
  int k = 0;
  int m = 0;
  for (int l = 0; l < A->rows; l++) {
    for (n = 0; n < A->columns; n++) {
      if (n != j && l != i) {
        minor.matrix[k][m] = A->matrix[l][n];
        m++;
        if (m == A->columns - 1) {
          k++;
          m = 0;
        }
      }
    }
  }
  return minor;
}
