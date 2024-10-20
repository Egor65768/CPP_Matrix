#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(function_class, inverse_matrix_1) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  double arr[] = {2.0, 5.0, 7.0, 6.0, 3.0, 4.0, 5.0, -2.0, -3.0};
  double res[] = {1.0, -1.0, 1.0, -38.0, 41.0, -34.0, 27.0, -29.0, 24.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  A = A.InverseMatrix();
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(res[count] * 1000000)) <
                  1);
    }
  }
}

// свойство обратной матрицы произведение обратной на обычную дает единичную
TEST(function_class, inverse_matrix_2) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  S21Matrix B(m, n);
  double arr[] = {2.0, 5.0, 7.0, 6.0, 3.0, 4.0, 5.0, -2.0, -3.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  B = A.InverseMatrix();
  A.MulMatrix(B);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) {
        ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(1.0 * 1000000)) < 1);
      } else {
        ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(0.0 * 1000000)) < 1);
      }
    }
  }
}

TEST(function_class, inverse_matrix_3) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  double arr[] = {2.0, 0.0, 1.0, 0.0, -3.0, -1.0, -2.0, 4.0, 0.0};
  double res[] = {2.0, 2.0, 1.5, 1.0, 1.0, 1.0, -3.0, -4.0, -3.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  A = A.InverseMatrix();
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(res[count] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, inverse_matrix_4) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  double arr[] = {2.0, NAN, 1.0, 0.0, -3.0, -1.0, -2.0, 4.0, 0.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.InverseMatrix(), invalid_argument);
}

TEST(function_class, inverse_matrix_5) {
  int n = 1, m = 3;
  S21Matrix A(m, n);
  double arr[] = {2.0, 12.2, 1.0, 0.0, -3.0, -1.0, -2.0, 4.0, 0.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.InverseMatrix(), logic_error);
}

TEST(function_class, inverse_matrix_6) {
  int n = 1, m = 3;
  S21Matrix A(m, n);
  double arr[] = {1.0, 1.0, 3.0, 1.0, 1.0, 3.0, 2.0, -1.0, 4.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.InverseMatrix(), logic_error);
}

TEST(function_class, inverse_matrix_7) {
  int n = 2, m = 2;
  S21Matrix A(m, n);
  double arr[] = {13.2, 13.0, 3.0, 1.0};
  double res[] = {-5.0 / 129.0, 65.0 / 129.0, 5.0 / 43.0, -22.0 / 43.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  A = A.InverseMatrix();
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(res[count] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, inverse_matrix_8) {
  int n = 1, m = 1;
  S21Matrix A(m, n);
  double arr[] = {2.0};
  A[0][0] = arr[0];
  A = A.InverseMatrix();
  ASSERT_TRUE(fabs(round(A[0][0] * 1000000) - round(0.5 * 1000000)) < 1);
}

TEST(function_class, inverse_matrix_9) {
  int n = 4, m = 4;
  S21Matrix A(m, n);
  double arr[] = {2.0, 1.0, 0.0, 0.0, 3.0, 2.0,  0.0, 0.0,
                  1.0, 1.0, 3.0, 4.0, 2.0, -1.0, 2.0, 3.0};
  double res[] = {2.0,  -1.0,  0.0, 0.0,  -3.0,  2.0,  0.0,  0.0,
                  31.0, -19.0, 3.0, -4.0, -23.0, 14.0, -2.0, 3.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  A = A.InverseMatrix();
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(res[count] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, inverse_matrix_10) {
  int n = 2, m = 2;
  S21Matrix A(m, n);
  double arr[] = {0.0, 0.0, 1.0, 12.3};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.InverseMatrix(), logic_error);
}