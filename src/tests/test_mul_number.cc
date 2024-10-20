#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(function_class, mult_number_1) {
  int n = 2, m = 2;
  S21Matrix A(m, n), res(m, n);
  double mul = 1.0;
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double res1[] = {12.312, 424.244, 234.324, 24.5443};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      res[i][j] = res1[count];
    }
  }
  A.MulNumber(mul);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, mult_number_2) {
  int n = 2, m = 2;
  S21Matrix A(m, n), res(m, n);
  double mul = 0.0;
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double res1[] = {0.0, 0.0, 0.0, 0.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      res[i][j] = res1[count];
    }
  }
  A.MulNumber(mul);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, mult_number_3) {
  int n = 2, m = 2;
  S21Matrix A(m, n), res(m, n);
  double mul = 2.212;
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double res1[4];
  for (int i = 0; i < 4; i++) {
    res1[i] = arr1[i] * mul;
  }
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      res[i][j] = res1[count];
    }
  }
  A.MulNumber(mul);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, mult_number_4) {
  int n = 2, m = 2;
  S21Matrix A(m, n), res(m, n);
  double mul = INFINITY;
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
    }
  }
  EXPECT_THROW(A.MulNumber(mul), invalid_argument);
}

TEST(function_class, mult_number_5) {
  int n = 2, m = 2;
  S21Matrix A(m, n), res(m, n);
  double mul = NAN;
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
    }
  }
  EXPECT_THROW(A.MulNumber(mul), invalid_argument);
}

TEST(function_class, mult_number_6) {
  int n = 2, m = 2;
  S21Matrix A(m, n), res(m, n);
  double mul = 1.2;
  double arr1[] = {12.312, NAN, 234.324, 24.5443};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
    }
  }
  EXPECT_THROW(A.MulNumber(mul), invalid_argument);
}

TEST(function_class, mult_number_7) {
  int n = 2, m = 2;
  S21Matrix A(m, n), res(m, n);
  double mul = 1.2;
  double arr1[] = {12.312, INFINITY, 234.324, 24.5443};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
    }
  }
  EXPECT_THROW(A.MulNumber(mul), invalid_argument);
}
