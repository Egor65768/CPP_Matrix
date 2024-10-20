#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(function_class, sum_test1) {
  int n = 2, m = 2;
  S21Matrix A(m, n), B(m, n), res(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double arr2[] = {0.123, 12.3121, 9.9999, -123.123};
  double res1[] = {12.435, 436.5561, 244.3239, -98.5787};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr2[count];
      res[i][j] = res1[count];
    }
  }
  A.SumMatrix(B);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, sum_test2) {
  int n = 1, m = 3;
  S21Matrix A(m, n), B(m, n), res(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double arr2[] = {0.124, 12.3121, 9.9999, -123.123};
  double res1[] = {12.436, 436.5561, 244.3239, -98.5787};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr2[count];
      res[i][j] = res1[count];
    }
  }
  A.SumMatrix(B);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, sum_test3) {
  int n = 1, m = 2;
  S21Matrix A(m, n), B(m, m), res(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double arr2[] = {0.124, 12.3121, 9.9999, -123.123};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
    }
  }
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < m; j++, count++) {
      B[i][j] = arr2[count];
    }
  }
  EXPECT_THROW(A.SumMatrix(B), logic_error);
}

TEST(function_class, sum_test4) {
  int n = 2, m = 2;
  S21Matrix A(m, n), B(m, n), res(m, n);
  double arr1[] = {INFINITY, 424.244, 234.324, 24.5443};
  double arr2[] = {1, 12.3121, 9.9999, -123.123};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr2[count];
    }
  }
  EXPECT_THROW(A.SumMatrix(B), invalid_argument);
}

TEST(function_class, sum_test5) {
  int n = 2, m = 2;
  S21Matrix A(m, n), B(m, n), res(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double arr2[] = {12, 12.3121, NAN, -123.123};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr2[count];
    }
  }
  EXPECT_THROW(A.SumMatrix(B), invalid_argument);
}

TEST(function_class, sum_test6) {
  int n = 2, m = 1;
  S21Matrix A(m, n), B(m, m), res(m, n);
  EXPECT_THROW(A.SumMatrix(B), logic_error);
}

TEST(function_class, error_create) {
  int n = 2, m = -1;
  EXPECT_THROW(S21Matrix A(n, m), length_error);
}

TEST(function_class, sub_test1) {
  int n = 2, m = 2;
  S21Matrix A(m, n), B(m, n), res(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double arr2[] = {0.123, 12.3121, 9.9999, -123.123};
  double res1[] = {12.189, 411.9319, 224.3241, 147.6673};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr2[count];
      res[i][j] = res1[count];
    }
  }
  A.SubMatrix(B);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, sub_test2) {
  int n = 1, m = 3;
  S21Matrix A(m, n), B(m, n), res(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double arr2[] = {0.124, 12.3121, 9.9999, -123.123};
  double res1[] = {12.188, 411.9319, 224.3241, 147.6673};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr2[count];
      res[i][j] = res1[count];
    }
  }
  A.SubMatrix(B);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 100000) - round(A[i][j] * 100000)) <
                  1);
    }
  }
}

TEST(function_class, sub_test3) {
  int n = 1, m = 2;
  S21Matrix A(m, n), B(m, m), res(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double arr2[] = {0.124, 12.3121, 9.9999, -123.123};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
    }
  }
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < m; j++, count++) {
      B[i][j] = arr2[count];
    }
  }
  EXPECT_THROW(A.SubMatrix(B), logic_error);
}

TEST(function_class, sub_test4) {
  int n = 2, m = 2;
  S21Matrix A(m, n), B(m, n);
  double arr1[] = {INFINITY, 424.244, 234.324, 24.5443};
  double arr2[] = {1, 12.3121, 9.9999, -123.123};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr2[count];
    }
  }
  EXPECT_THROW(A.SubMatrix(B), invalid_argument);
}

TEST(function_class, sub_test5) {
  int n = 2, m = 2;
  S21Matrix A(m, n), B(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  double arr2[] = {12, 12.3121, NAN, -123.123};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr2[count];
    }
  }
  EXPECT_THROW(A.SubMatrix(B), invalid_argument);
}

TEST(function_class, sub_test6) {
  int n = 2, m = 1;
  S21Matrix A(m, n), B(m, m);
  EXPECT_THROW(A.SubMatrix(B), logic_error);
}