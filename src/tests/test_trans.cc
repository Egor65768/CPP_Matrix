#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(function_class, transpose_1) {
  double arr[] = {13.243, 123.342, 0.1233, 123.34, 1242.34,
                  0.242,  345.4,   34.5,   456.2};
  int n = 3, m = 3;
  S21Matrix A(m, n), res(m, n);

  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
      res[j][i] = arr[count];
    }
  }

  A = A.Transpose();

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, transpose_2) {
  int n = 3, m = 2;
  S21Matrix A(m, n), res(n, m);
  double arr[] = {13.43, 23.342, 0.13, 9823.34, 1442.34, 110.242, 32.2};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
      res[j][i] = arr[count];
    }
  }
  A = A.Transpose();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, transpose_3) {
  int m = 3, n = 4;
  S21Matrix A(m, n), res(n, m);
  double arr[] = {1.23, 32.3,    53.1,  0.01,      244.323, 132.3, 3.1,
                  3.4,  356.632, 0.423, 3204.2423, 432.44,  32.4};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
      res[j][i] = arr[count];
    }
  }
  A = A.Transpose();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}
