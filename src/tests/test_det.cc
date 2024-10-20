#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(function_class, det_matrix_1) {
  int n = 3, m = 3;
  S21Matrix A(n, m);
  double res = 0.0;
  double arr[] = {12.0, 3.0, 3.0, 6.0, 5.0, 4.0, 9.0, 32.0, 5.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  res = A.Determinant();
  ASSERT_TRUE(fabs(round(res * 1000000) - round(-777.0 * 1000000)) < 1);
}

TEST(function_class, det_matrix_2) {
  int n = 2, m = 3;
  S21Matrix A(m, n);
  double arr[] = {12.0, 3.0, 3.0, 6.0, 5.0, 4.0, 9.0, 32.0, 5.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.Determinant(), logic_error);
}

TEST(function_class, det_matrix_3) {
  int n = 2, m = 2;
  S21Matrix A(m, n);
  double res = 0.0;
  double arr[] = {12.0, 3.0, 3.0, 6.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  res = A.Determinant();
  ASSERT_TRUE(fabs(round(res * 1000000) - round(63.0 * 1000000)) < 1);
}

TEST(function_class, det_matrix_4) {
  int n = 1, m = 1;
  S21Matrix A(m, n);
  double res = 0.0;
  double arr[] = {12.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  res = A.Determinant();
  ASSERT_TRUE(fabs(round(res * 1000000) - round(12.0 * 1000000)) < 1);
}

TEST(function_class, det_matrix_5) {
  int n = 1, m = 1;
  S21Matrix A(m, n);
  double arr[] = {12.0 / 0.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.Determinant(), invalid_argument);
}

TEST(function_class, det_matrix_6) {
  int n = 4, m = 4;
  S21Matrix A(m, n);
  double res = 0.0;
  double arr[] = {12.0, 3.0, 14.0, 2.0, 6.0, 5.0, 5.0,  6.0,
                  4.0,  6.0, 8.0,  9.0, 4.0, 6.0, -2.0, -13.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  res = A.Determinant();
  ASSERT_TRUE(fabs(round(res * 1000000) - round(-5220.0 * 1000000)) < 1);
}

TEST(function_class, det_matrix_7) {
  int n = 5, m = 5;
  S21Matrix A(m, n);
  double res = 0.0;
  double arr[] = {12.0,  3.0, 14.0, 2.0, 1.0, 6.0, 5.0, 5.0, 6.0,
                  2.0,   4.0, 6.0,  8.0, 9.0, 3.0, 4.0, 6.0, -2.0,
                  -13.0, 4.0, 5.0,  6.0, 4.0, 3.0, 4.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  res = A.Determinant();
  ASSERT_TRUE(fabs(round(res * 1000000) - round(-5391.0 * 1000000)) < 1);
}

TEST(function_class, det_matrix_8) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  double res = 0.0;
  double arr[] = {12.12, 3.3, 3.4, 6.0, 5.0, 4.0, 9.01, 32.10, 5.4};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  res = A.Determinant();
  ASSERT_TRUE(fabs(round(res * 1000000) - round(-715.286 * 1000000)) < 1);
}

TEST(function_class, det_matrix_9) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  double arr[] = {0.0 / 0.0, 3.3, 3.4, 6.0, 5.0, 4.0, 9.01, 32.10, 5.4};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.Determinant(), invalid_argument);
}

TEST(function_class, det_matrix_10) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  double arr[] = {1.0 / 0.0, 3.3, 3.4, 6.0, 5.0, 4.0, 9.01, 32.10, 5.4};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.Determinant(), invalid_argument);
}