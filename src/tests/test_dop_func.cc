#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"

TEST(dopfunction_class, set_rows_1) {
  int n = 3, m = 3;
  S21Matrix A(n, m);
  EXPECT_THROW(A.set_rows(-1), logic_error);
}

TEST(dopfunction_class, set_rows_2) {
  int n = 3, m = 3;
  S21Matrix A(n, m);
  A.set_rows(2);
  ASSERT_TRUE(A.get_rows() == 2);
}

TEST(dopfunction_class, set_cols_1) {
  int n = 3, m = 3;
  S21Matrix A(n, m);
  EXPECT_THROW(A.set_cols(-1), logic_error);
}

TEST(dopfunction_class, set_cols_2) {
  int n = 3, m = 3;
  S21Matrix A(n, m);
  A.set_cols(2);
  ASSERT_TRUE(A.get_cols() == 2);
}

TEST(dopfunction_class, bracje1) {
  int n = 3, m = 3;
  S21Matrix A(n, m);
  EXPECT_THROW(A[4][2] = 3, out_of_range);
}

TEST(dopfunction_class, bracje2) {
  int n = 3, m = 3;
  S21Matrix A(n, m);
  A[2][2] = 3.1415;
  ASSERT_TRUE(fabs(round(A[2][2] * 1000000) - round(3.1415 * 1000000)) < 1);
}

TEST(dopfunction_class, bracje3) {
  int n = 3, m = 3;
  S21Matrix A(n, m);
  A[2][2] = 3.1415;
  ASSERT_TRUE(fabs(round(A(2, 2) * 1000000) - round(3.1415 * 1000000)) < 1);
}

TEST(dopfunction_class, bracje4) {
  int n = 3, m = 3;
  S21Matrix A(n, m);
  EXPECT_THROW(A(4, 2), out_of_range);
}

TEST(dopfunction_class, bracje5) {
  int n = 3, m = 3;
  S21Matrix A(n, m);
  EXPECT_THROW(A(-1, 2) = 3, length_error);
}
