#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(operator_class, sum_test1) {
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
  S21Matrix otv = A + B;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(
          fabs(round(res[i][j] * 1000000) - round(otv[i][j] * 1000000)) < 1);
    }
  }
}

TEST(operator_class, sum_test2) {
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
  S21Matrix otv = A + B;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(
          fabs(round(res[i][j] * 1000000) - round(otv[i][j] * 1000000)) < 1);
    }
  }
}

TEST(operator_class, sum_test3) {
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
  EXPECT_THROW(A + B, logic_error);
}

TEST(operator_class, sum_test4) {
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
  EXPECT_THROW(A + B, invalid_argument);
}

TEST(operator_class, sum_test5) {
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
  EXPECT_THROW(A + B, invalid_argument);
}

TEST(operator_class, sum_test6) {
  int n = 2, m = 1;
  S21Matrix A(m, n), B(m, m), res(m, n);
  EXPECT_THROW(A + B, logic_error);
}

TEST(operator_class, sub_test1) {
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
  S21Matrix otv = A - B;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(
          fabs(round(res[i][j] * 1000000) - round(otv[i][j] * 1000000)) < 1);
    }
  }
}

TEST(operator_class, sub_test2) {
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
  S21Matrix otv = A - B;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 100000) - round(otv[i][j] * 100000)) <
                  1);
    }
  }
}

TEST(operator_class, sub_test3) {
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
  EXPECT_THROW(A - B, logic_error);
}

TEST(operator_class, sub_test4) {
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
  EXPECT_THROW(A - B, invalid_argument);
}

TEST(operator_class, sub_test5) {
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
  EXPECT_THROW(A - B, invalid_argument);
}

TEST(operator_class, sub_test6) {
  int n = 2, m = 1;
  S21Matrix A(m, n), B(m, m);
  EXPECT_THROW(A - B, logic_error);
}

TEST(operator_class, mult_number_1) {
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
  S21Matrix otv = A * mul;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(
          fabs(round(res[i][j] * 1000000) - round(otv[i][j] * 1000000)) < 1);
    }
  }
}

TEST(operator_class, mult_number_2) {
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
  S21Matrix otv = A * mul;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(
          fabs(round(res[i][j] * 1000000) - round(otv[i][j] * 1000000)) < 1);
    }
  }
}

TEST(operator_class, mult_matrix_1) {
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
  S21Matrix otv = A * B;
  for (int i = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++) {
      ASSERT_TRUE(
          fabs(round(res[i][j] * 1000000) - round(otv[i][j] * 1000000)) < 1);
    }
  }
}

TEST(operator_class, mult_matrix_2) {
  S21Matrix A(4, 4), B(3, 4);
  EXPECT_THROW(A * B, logic_error);
}

TEST(operator_class, mult_matrix_12) {
  S21Matrix A(4, 3), B(4, 4);
  EXPECT_THROW(A * B, logic_error);
}

TEST(operator_class, eq_test1) {
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
  ASSERT_TRUE(A == B);
}

TEST(operator_class, eq_test2) {
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
  ASSERT_FALSE(A == B);
}

TEST(operator_class, eq_test3) {
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
  ASSERT_FALSE(A == B);
}

TEST(operator_class, eq_test4) {
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
  ASSERT_FALSE(A == B);
}

TEST(operator_class, eq_test5) {
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
  ASSERT_FALSE(A == B);
}

TEST(operator_class, eq_test12) {
  int n = 3, m = 3;
  S21Matrix A(n, m), B(n, 2);
  ASSERT_FALSE(A == B);
}

TEST(operator_class, ravno_1) {
  double arr1[] = {12.312, 424.244, 234.324, 24.5443, 12.34, 34.424};
  int n1 = 2, m1 = 3;
  S21Matrix A(m1, n1), B(m1, n1);
  for (int i = 0, count = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr1[count];
    }
  }
  S21Matrix otv;
  otv = A;
  for (int i = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++) {
      ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(B[i][j] * 1000000)) <
                  1);
    }
  }
  for (int i = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++) {
      ASSERT_TRUE(fabs(round(otv[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(operator_class, ravno_2) {
  double arr1[] = {12.312, 424.244, 234.324, 24.5443, 12.34, 34.424};
  int n1 = 2, m1 = 3;
  S21Matrix A(m1, n1), B(m1, n1);
  for (int i = 0, count = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++, count++) {
      A[i][j] = arr1[count];
      B[i][j] = arr1[count] + 2.3;
    }
  }
  S21Matrix otv;
  otv = A;
  B = otv;
  for (int i = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++) {
      ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(B[i][j] * 1000000)) <
                  1);
      ASSERT_TRUE(fabs(round(otv[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
      ASSERT_TRUE(fabs(round(otv[i][j] * 1000000) - round(B[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(operator_class, sum_test1_new) {
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
  A += B;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(operator_class, sum_test2_new) {
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
  A += B;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(operator_class, sum_test3new) {
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
  EXPECT_THROW(A += B, logic_error);
}

TEST(operator_class, sub_test1new) {
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
  A -= B;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(operator_class, sub_test2new) {
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
  A -= B;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 100000) - round(A[i][j] * 100000)) <
                  1);
    }
  }
}

TEST(operator_class, sub_test3new) {
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
  EXPECT_THROW(A -= B, logic_error);
}

TEST(operator_class, mult_number_1new) {
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
  A *= mul;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(operator_class, mult_number_2new) {
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
  A *= mul;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(operator_class, mult_matrix_1new) {
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
  A *= B;
  for (int i = 0; i < m1; i++) {
    for (int j = 0; j < n1; j++) {
      ASSERT_TRUE(fabs(round(res[i][j] * 1000000) - round(A[i][j] * 1000000)) <
                  1);
    }
  }
}

TEST(operator_class, mult_matrix_2new) {
  S21Matrix A(4, 4), B(3, 4);
  EXPECT_THROW(A *= B, logic_error);
}

TEST(operator_class, mult_matrix_12new) {
  S21Matrix A(4, 3), B(4, 4);
  EXPECT_THROW(A *= B, logic_error);
}

TEST(operator_class, set_rows) {
  int n = 2, m = 2;
  S21Matrix A(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
    }
  }
  A.set_rows(1);
  ASSERT_TRUE(1 == A.get_rows());
  for (int i = 0, count = 0; i < A.get_rows(); i++) {
    for (int j = 0; j < A.get_cols(); j++, count++) {
      ASSERT_TRUE(
          fabs(round(A[i][j] * 1000000) - round(arr1[count] * 1000000)) < 1);
    }
  }
}

TEST(operator_class, set_rows1) {
  int n = 2, m = 2;
  S21Matrix A(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
    }
  }
  A.set_rows(4);
  ASSERT_TRUE(4 == A.get_rows());
  for (int i = 0, count = 0; i < A.get_rows(); i++) {
    for (int j = 0; j < A.get_cols(); j++, count++) {
      if (i < 2) {
        ASSERT_TRUE(
            fabs(round(A[i][j] * 1000000) - round(arr1[count] * 1000000)) < 1);
      } else {
        ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(0 * 1000000)) < 1);
      }
    }
  }
}

TEST(operator_class, set_cols) {
  int n = 2, m = 2;
  S21Matrix A(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
    }
  }
  A.set_cols(1);
  ASSERT_TRUE(1 == A.get_cols());
  for (int i = 0, count = 0; i < A.get_rows(); i++) {
    for (int j = 0; j < A.get_cols(); j++, count++) {
      ASSERT_TRUE(fabs(round(A[i][j] * 1000000) -
                       round(arr1[2 * count] * 1000000)) < 1);
    }
  }
}

TEST(operator_class, set_cols1) {
  int n = 2, m = 2;
  S21Matrix A(m, n);
  double arr1[] = {12.312, 424.244, 234.324, 24.5443};
  for (int i = 0, count = 0; i < m; i++) {
    for (int j = 0; j < n; j++, count++) {
      A[i][j] = arr1[count];
    }
  }
  A.set_cols(4);
  ASSERT_TRUE(4 == A.get_cols());
  for (int i = 0, count = 0; i < A.get_rows(); i++) {
    for (int j = 0; j < A.get_cols(); j++) {
      if (j < 2) {
        ASSERT_TRUE(
            fabs(round(A[i][j] * 1000000) - round(arr1[count] * 1000000)) < 1);
        count++;
      } else {
        ASSERT_TRUE(fabs(round(A[i][j] * 1000000) - round(0 * 1000000)) < 1);
      }
    }
  }
}