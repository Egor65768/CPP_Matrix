#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(function_class, eq_test1) {
  int n = 3, m = 3;
  S21Matrix A(n, m), B(n, m);
  double arr[] = {64.234455, 1234.432,  213.221, 332.2344, 0.12343,
                  123.333,   88109.331, 1223.13, 0.21233};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
      B[i][j] = arr[count];
    }
  }
  ASSERT_TRUE(A.EqMatrix(B));
}

TEST(function_class, eq_test2) {
  int n = 3, m = 3;
  S21Matrix A(n, m), B(n, m);
  double arr[] = {64.234455, 1234.432,  213.221, 332.2344, 0.12343,
                  123.333,   88109.331, 1223.13, 0.21233};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
      B[i][j] = arr[count];
    }
  }
  B[0][0] = 1.212334;
  ASSERT_FALSE(A.EqMatrix(B));
}

TEST(function_class, eq_test3) {
  int n = 3, m = 3, k = 2;
  S21Matrix A(m, n), B(m, k);
  double arr[] = {64.234455, 1234.432,  213.221, 332.2344, 0.12343,
                  123.333,   88109.331, 1223.13, 0.21233};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < k; j++, count++) {
      B[i][j] = arr[count];
    }
  }
  ASSERT_FALSE(A.EqMatrix(B));
}

TEST(function_class, eq_test4) {
  int n = 3, m = 3, k = 2;
  S21Matrix A(m, n), B(k, n);
  double arr[] = {64.234455, 1234.432,  213.221, 332.2344, 0.12343,
                  123.333,   88109.331, 1223.13, 0.21233};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  for (int i = 0, count = 0; i < k; i++) {
    for (int j = 0; j < n; j++, count++) {
      B[i][j] = arr[count];
    }
  }
  ASSERT_FALSE(A.EqMatrix(B));
}

TEST(function_class, eq_test5) {
  int n = 3, m = 3;
  S21Matrix A(n, m), B(n, m);
  double arr[] = {64.234455, 1234.432,  213.221, 332.2344, 0.12343,
                  123.333,   88109.331, 1223.13, 0.21233};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
      B[i][j] = arr[count] + 0.000001;
    }
  }
  B[0][0] = 1.212334;
  ASSERT_FALSE(A.EqMatrix(B));
}

TEST(function_class, eq_test6) {
  int n = 3, m = 3;
  S21Matrix A(n, m), B(n, m);
  double arr[] = {64.234455, 1234.432,  213.221, 332.2344, 0.12343,
                  123.333,   88109.331, 1223.13, 0.21233};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
      B[i][j] = arr[count];
    }
  }
  B[2][2] += 0.0000001;
  ASSERT_FALSE(A.EqMatrix(B));
}

TEST(function_class, eq_test7) {
  int n = 3, m = 3;
  S21Matrix A(n, m), B(n, m);
  double arr[] = {64.234455, 1234.432,  213.221, 332.2344, 0.12343,
                  123.333,   88109.331, 1223.13, 0.21233};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
      B[i][j] = arr[count];
    }
  }
  B[2][2] += 1e-8;
  ASSERT_TRUE(A.EqMatrix(B));
}

TEST(function_class, eq_test8) {
  int n = 3, m = 3;
  S21Matrix A(n, m), B(n, m);
  double arr[] = {64.234455, 1234.432,  213.221, 332.2344, 0.12343,
                  123.333,   88109.331, 1223.13, 0.21233};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
      B[i][j] = arr[count];
    }
  }
  B[0][0] += 0.000001;
  ASSERT_FALSE(A.EqMatrix(B));
}
