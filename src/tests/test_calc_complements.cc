#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(function_class, create_calc_complements_1) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  double arr[] = {1.0, 2.0, 3.0, 0.0, 4.0, 2.0, 5.0, 2.0, 1.0};
  double res[] = {0.0, 10.0, -20.0, 4.0, -14.0, 8.0, -8.0, -2.0, 4.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  A = A.CalcComplements();
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(res[count] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, create_calc_complements_2) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  double arr[] = {5.0, 7.0, 1.0, -4.0, 1.0, 0.0, 2.0, 0.0, 3.0};
  double res[] = {3.0, 12.0, -2.0, -21.0, 13.0, 14.0, -1.0, -4.0, 33.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  A = A.CalcComplements();
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(res[count] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, create_calc_complements_3) {
  int n = 2, m = 3;
  S21Matrix A(m, n);
  double arr[] = {1.0, 2.0, 3.0, 0.0, 4.0, 2.0, 5.0, 2.0, 1.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.CalcComplements(), logic_error);
}

TEST(function_class, create_calc_complements_4) {
  int n = 3, m = 2;
  S21Matrix A(m, n);
  double arr[] = {1.0, 2.0, 3.0, 0.0, 4.0, 2.0, 5.0, 2.0, 1.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.CalcComplements(), logic_error);
}

TEST(function_class, create_calc_complements_5) {
  int n = 1, m = 1;
  S21Matrix A(m, n);
  A[0][0] = 12.4;
  A = A.CalcComplements();
  ASSERT_TRUE(fabs(round(A[0][0] * 1000000) - round(1.0 * 1000000)) < 1);
}

TEST(function_class, create_calc_complements_7) {
  int n = 2, m = 2;
  S21Matrix A(m, n);
  double arr[] = {6.0, 36.0, 96.0, 27.0};
  double res[] = {27.0, -96.0, -36.0, 6.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  A = A.CalcComplements();
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(res[count] * 1000000)) <
                  1);
    }
  }
}

TEST(function_class, create_calc_complements_8) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  double arr[] = {INFINITY, 2.0, 3.0, 0.0, 4.0, 2.0, 5.0, 2.0, 1.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.CalcComplements(), invalid_argument);
}

TEST(function_class, create_calc_complements_9) {
  int n = 3, m = 3;
  S21Matrix A(m, n);
  double arr[] = {23.1, 2.0, 3.0, 0.0, 4.0, NAN, 5.0, 2.0, 1.0};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr[count];
    }
  }
  EXPECT_THROW(A.CalcComplements(), invalid_argument);
}