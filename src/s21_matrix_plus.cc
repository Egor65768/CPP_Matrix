#include "s21_matrix_oop.h"

/*
 КОНСТРУКТОРЫ И ДЕСТРУКТОР
 */

S21Matrix::S21Matrix() : rows_(DEFAULT_ROWS), cols_(DEFAULT_COLS) {
  Allocate();
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ > 0 && cols_ > 0) {
    Allocate();
  } else {
    throw length_error("Size array error");
  }
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  Allocate();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) {
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() { Deallocate(); }

/*
  ОСНОВНЫЕ ФУНКЦИИ
*/

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  bool res = true;
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    res = false;
  }
  for (int i = 0; i < rows_ && res != false; i++) {
    for (int j = 0; j < cols_ && res != false; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7) {
        res = false;
      }
    }
  }
  return res;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  checkMatrixAddSubValidity(other);
  addOrSubtractMatrix(other, 1);
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  checkMatrixAddSubValidity(other);
  addOrSubtractMatrix(other, -1);
}

void S21Matrix::MulNumber(const double num) {
  if (nan_or_inf() || isnan(num) || isinf(num)) {
    throw invalid_argument("Matrix element cannot be nan or inf");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (matrix_nan_or_inf(other)) {
    throw invalid_argument("Matrix element cannot be nan or inf");
  } else if (cols_ != other.rows_) {
    throw logic_error("Matrix dimensions do not match");
  }
  double res = 0.0;
  S21Matrix buf_matrix(rows_, other.get_cols());
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      res = 0;
      for (int k = 0; k < cols_; k++) {
        res = res + matrix_[i][k] * other.matrix_[k][j];
      }
      buf_matrix[i][j] = res;
    }
  }
  copyMatrix(buf_matrix);
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix trans(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      trans[j][i] = (*this)[i][j];
    }
  }
  return trans;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw logic_error("Matrix dimensions do not match");
  } else if (nan_or_inf()) {
    throw invalid_argument("Matrix element cannot be nan or inf");
  }
  S21Matrix CalcComplementsMatrix(rows_, rows_);
  if (rows_ > 1) {
    S21Matrix CalcComplementsMatrix(rows_, rows_);
    S21Matrix buf_matrix(rows_ - 1, rows_ - 1);
    double buf_determinant = 0.0;
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < rows_; j++) {
        buildMatrixWithoutRowCol(&buf_matrix, i, j);
        buf_determinant = buf_matrix.Determinant();
        CalcComplementsMatrix[i][j] = buf_determinant;
      }
    }
    CalcComplementsMatrix.computeMatrixCofactors();
    return CalcComplementsMatrix;
  }
  CalcComplementsMatrix.set_cols(1);
  CalcComplementsMatrix.set_rows(1);
  CalcComplementsMatrix[0][0] = 1.0;
  return CalcComplementsMatrix;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw logic_error("Invalid matrix");
  } else if (nan_or_inf()) {
    throw invalid_argument("matrix element cannot be nan or inf");
  }
  double res = 0;
  int sign = 1;
  if (rows_ == 3) {
    for (int row = 0; row < rows_; row++) {
      S21Matrix buf(rows_ - 1, cols_ - 1);
      if ((row + 2) % 2 == 0) {
        sign = 1;
      } else {
        sign = -1;
      }
      buildMatrixWithoutRowCol(&buf, 0, row);
      res = res + sign * (*this)[0][row] * compute2x2Minor(buf);
    }
  } else if (rows_ == 2) {
    res = compute2x2Minor(*this);
  } else if (rows_ == 1) {
    res = (*this)[0][0];
  } else if (rows_ > 3) {
    for (int i = 0; i < rows_; i++) {
      S21Matrix buf_matrix(rows_ - 1, cols_ - 1);
      buildMatrixWithoutRowCol(&buf_matrix, 0, i);
      double minor = buf_matrix.Determinant();
      if ((i + 2) % 2 == 0) {
        sign = 1;
      } else {
        sign = -1;
      }
      res = res + sign * (*this)[0][i] * minor;
    }
  }
  return res;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (nan_or_inf()) {
    throw invalid_argument("matrix element cannot be nan or inf");
  } else if (cols_ != rows_) {
    throw logic_error("Invalid matrix");
  }
  double det = Determinant();
  if (det == 0) {
    throw logic_error("The determinant cannot be zero");
  }
  S21Matrix res_matrix(rows_, cols_);
  res_matrix = CalcComplements();
  res_matrix = res_matrix.Transpose();
  det = 1.0 / det;
  res_matrix.MulNumber(det);
  return res_matrix;
}

/*
  ЧИСЛОВЫЕ ОПЕРАТОРЫ
*/

S21Matrix S21Matrix::operator+(const S21Matrix &other) const {
  S21Matrix res = (*this);
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) const {
  S21Matrix res = (*this);
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const double num) const {
  S21Matrix res = (*this);
  res.MulNumber(num);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) const {
  S21Matrix res = (*this);
  res.MulMatrix(other);
  return res;
}

bool S21Matrix::operator==(const S21Matrix &other) const noexcept {
  return this->EqMatrix(other);
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this != &other) {
    copyMatrix(other);
  }
  return *this;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) noexcept {
  if (this != &other) {
    swap(rows_, other.rows_);
    swap(cols_, other.cols_);
    swap(matrix_, other.matrix_);
  }
  return *this;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  this->MulMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double num) {
  this->MulNumber(num);
  return *this;
}

double *S21Matrix::operator[](int row) const {
  if (this->rows_ <= row || this->rows_ < 0) {
    throw out_of_range("There is no element with such an index");
  }
  return matrix_[row];
}

double &S21Matrix::operator()(int row, int col) const {
  if (this->rows_ <= row || this->cols_ <= col) {
    throw out_of_range("There is no element with such an index");
  } else if (row < 0 || col < 0) {
    throw length_error("The dimension of the matrix must be greater than zero");
  }
  return matrix_[row][col];
}

/*
ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
*/

bool S21Matrix::nan_or_inf() {
  bool res = false;
  for (int i = 0; i < this->rows_ && !res; i++) {
    for (int j = 0; j < this->cols_ && !res; j++) {
      if (isnan((*this)[i][j]) || isinf((*this)[i][j])) {
        res = true;
      }
    }
  }
  return res;
}

void S21Matrix::copyMatrix(const S21Matrix &other) {
  Deallocate();
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  Allocate();
  for (int i = 0; i < rows_ && i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      (*this)[i][j] = other[i][j];
    }
  }
}

double S21Matrix::compute2x2Minor(S21Matrix A) {
  double res = 0;
  res = A[0][0] * A[1][1] - A[1][0] * A[0][1];
  return res;
}

void S21Matrix::buildMatrixWithoutRowCol(S21Matrix *matrix, int row, int col) {
  int i1 = 0, j1 = 0;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (i != row && j != col) {
        (*matrix)[j1][i1] = (*this)[i][j];
        i1++;
        if (i1 == matrix->rows_) {
          j1++;
          i1 = 0;
        }
      }
    }
  }
}

void S21Matrix::computeMatrixCofactors() {
  int sign = -1;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if ((i + j + 2) % 2 == 0) {
        sign = 1;
      } else {
        sign = -1;
      }
      (*this)[i][j] *= sign;
    }
  }
}

bool S21Matrix::matrix_nan_or_inf(const S21Matrix &other) {
  bool res = false;
  for (int i = 0; i < rows_ && !res; i++) {
    for (int j = 0; j < cols_ && !res; j++) {
      if (isnan((*this)[i][j]) || isinf((*this)[i][j])) {
        res = true;
      }
    }
  }
  for (int i = 0; i < other.rows_ && !res; i++) {
    for (int j = 0; j < other.cols_ && !res; j++) {
      if (isnan(other[i][j]) || isinf(other[i][j])) {
        res = true;
      }
    }
  }
  return res;
}

void S21Matrix::checkMatrixAddSubValidity(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw logic_error("Matrix dimensions do not match");
  } else if (matrix_nan_or_inf(other)) {
    throw invalid_argument("Matrix element cannot be nan or inf");
  }
}

/*
АКСЕССОРЫ И МУТАТОРЫ
*/

int S21Matrix::get_rows() const { return rows_; }

int S21Matrix::get_cols() const { return cols_; }

void S21Matrix::set_rows(int rows) {
  int old_rows = rows_;
  if (rows <= 0) {
    throw length_error("Size array error");
  }
  if (rows != old_rows) {
    S21Matrix buf_matrix(rows, cols_);
    for (int row = 0; row < rows && row < rows_; row++) {
      for (int col = 0; col < cols_; col++) {
        if (row < old_rows) {
          buf_matrix[row][col] = (*this)[row][col];
        } else {
          buf_matrix[row][col] = 0;
        }
      }
    }
    copyMatrix(buf_matrix);
  }
}

void S21Matrix::set_cols(int cols) {
  if (cols <= 0) {
    throw length_error("Size array error");
  }
  if (cols != cols_) {
    S21Matrix buf_matrix(rows_, cols);
    for (int row = 0; row < rows_; row++) {
      for (int col = 0; col < cols_ && col < cols; col++) {
        if (col < cols_) {
          buf_matrix[row][col] = (*this)[row][col];
        } else {
          buf_matrix[row][col] = 0;
        }
      }
    }
    copyMatrix(buf_matrix);
  }
}

/*
ПРИВАТНЫЕ ФУНКЦИИ
*/

void S21Matrix::addOrSubtractMatrix(const S21Matrix &other, int operation) {
  for (int row = 0; row < rows_; row++) {
    for (int col = 0; col < cols_; col++) {
      this->matrix_[row][col] =
          this->matrix_[row][col] + operation * other.matrix_[row][col];
    }
  }
}

void S21Matrix::Allocate() {
  matrix_ = new double *[rows_];
  for (int row = 0; row < rows_; row++) {
    matrix_[row] = new double[cols_];
  }
  for (int row = 0; row < rows_; row++) {
    for (int col = 0; col < cols_; col++) {
      matrix_[row][col] = 0;
    }
  }
}

void S21Matrix::Deallocate() {
  for (int row = 0; row < rows_; row++) {
    delete[] matrix_[row];
  }
  delete[] matrix_;
  matrix_ = nullptr;
}
