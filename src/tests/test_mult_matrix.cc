#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(function_class, mult_matrix_1) {
  double arr1[] = {12.312, 424.244, 234.324, 24.5443, 12.34, 34.424};
  double arr2[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  double res1[] = {436.556,  436.556, 436.556, 258.8683, 258.8683,
                   258.8683, 46.764,  46.764,  46.764};
  int n1 = 2, m1 = 3;
  int n2 = 3, m2 = 2;
  S21Matrix A(m1, n1), B(m2, n2), res(m1, n2);
  for (int i = 0, count = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++, count++) {
      A[i][j] = arr1[count];
      B[j][i] = arr2[count];
    }
  }
  for (int i = 0, count = 0; i < m1; i++) {
    for (int j = 0; j < n2; j++, count++) {
      res[i][j] = res1[count];
    }
  }
  A.MulMatrix(B);
  for (int i = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, mult_matrix_2) {
  S21Matrix A(4, 4), B(3, 4);
  EXPECT_THROW(A.MulMatrix(B), logic_error);
}

TEST(function_class, mult_matrix_12) {
  S21Matrix A(4, 3), B(4, 4);
  EXPECT_THROW(A.MulMatrix(B), logic_error);
}

TEST(function_class, mult_matrix_4) {
  int n1 = 2, m1 = 2;
  int n2 = 2, m2 = 2;
  S21Matrix A(m1, n1), B(m2, n2), res(m1, n2);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double arr2[] = {24.455, 0.424, 244.522, 13.3321};
  double res1[] = {104038.081328, 5661.2837204, 11732.0147446, 426.58043803};
  for (int i = 0, count = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr2[count];
    }
  }
  for (int i = 0, count = 0; i < m1; i++) {
    for (int j = 0; j < n2; j++, count++) {
      res[i][j] = res1[count];
    }
  }
  A.MulMatrix(B);
  for (int i = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, mult_matrix_5) {
  S21Matrix A(4, 4), B(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      A[i][j] = 1;
      B[i][j] = 1;
    }
  }
  A[0][0] = INFINITY;
  EXPECT_THROW(A.MulMatrix(B), invalid_argument);
}

TEST(function_class, mult_matrix_6) {
  S21Matrix A(4, 4), B(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      A[i][j] = 1;
      B[i][j] = 1;
    }
  }
  A[1][0] = NAN;
  EXPECT_THROW(A.MulMatrix(B), invalid_argument);
}
