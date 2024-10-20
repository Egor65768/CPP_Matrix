#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#include <algorithm>
#include <cmath>
#include <iostream>

#define DEFAULT_ROWS 2
#define DEFAULT_COLS 2

using namespace std;

class S21Matrix {
 private:
  int rows_, cols_;
  double** matrix_;

  // Выделение памяти под матрицу
  void Allocate();

  // Отчистка памяти
  void Deallocate();

  // Функция осуществляющая сложение и вычитание матриц
  void addOrSubtractMatrix(const S21Matrix& other, int operation);

 public:
  /*
  КОНСТРУКТОРЫ И ДЕСТРУКТОР
  */

  // Конструктор для создание дефолтной матрицы
  S21Matrix();

  // Конструктор для создание матрицы заданной размерности
  S21Matrix(int rows, int cols);

  // Конструктор копирования
  S21Matrix(const S21Matrix& other);

  // Конструктор переноса
  S21Matrix(S21Matrix&& other);

  // деструктор
  ~S21Matrix();

  /*
  ОСНОВНЫЕ ФУНКЦИИ
  */

  // Проверяет матрицы на равенство между собой
  bool EqMatrix(const S21Matrix& other) const;

  // Суммирует две матрицы
  void SumMatrix(const S21Matrix& other);

  // Вычитает из текущей матрицы другую
  void SubMatrix(const S21Matrix& other);

  // Умножает текущую матрицу на число
  void MulNumber(const double num);

  // Умножает текущую матрицу на вторую
  void MulMatrix(const S21Matrix& other);

  // Создает новую транспонированную матрицу из текущей и возвращает ее.
  S21Matrix Transpose();

  // Вычисляет матрицу алгебраических дополнений и возвращает ее.
  S21Matrix CalcComplements();

  // Вычисляет определитель матрницы
  double Determinant();

  // Вычисляет и возвращает обратную матрицу.
  S21Matrix InverseMatrix();

  /*
  ЧИСЛОВЫЕ ОПЕРАТОРЫ
  */

  //+
  S21Matrix operator+(const S21Matrix& other) const;

  //-
  S21Matrix operator-(const S21Matrix& other) const;

  //* матрица на число
  S21Matrix operator*(const double num) const;

  //* на матрицу
  S21Matrix operator*(const S21Matrix& other) const;

  //==
  bool operator==(const S21Matrix& other) const noexcept;

  //=
  S21Matrix& operator=(const S21Matrix& other);

  //=
  S21Matrix& operator=(S21Matrix&& other) noexcept;

  //+=
  S21Matrix& operator+=(const S21Matrix& other);

  //-=
  S21Matrix& operator-=(const S21Matrix& other);

  //*= матрица
  S21Matrix& operator*=(const S21Matrix& other);

  //*= число
  S21Matrix& operator*=(const double num);

  // Оператор для доступа к элементу матрицы
  double* operator[](int row) const;

  // Оператор возвращает элемент матрицы
  double& operator()(int row, int col) const;

  /*
  ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
  */

  // Проверяет текущую матрицу на содержание inf или nan
  bool nan_or_inf();

  // Удаляет матрицу и создает новую в которую записывает передаваемую матрицу
  void copyMatrix(const S21Matrix& other);

  // Вычисляет определитель матрицы 2 на 2
  double compute2x2Minor(S21Matrix A);

  // Создает матрицу, которая получается вычеркивание строки и столбца
  void buildMatrixWithoutRowCol(S21Matrix* matrix, int row, int col);

  // Cтавит нужные знаки для матрицы алгебраических дополнений
  void computeMatrixCofactors();

  // Проверяет текущую матрицу и переданную на содержание inf или nan
  bool matrix_nan_or_inf(const S21Matrix& other);

  // Проверяет валидность матрицы для сложения и вычитания
  void checkMatrixAddSubValidity(const S21Matrix& other);

  /*
  АКСЕССОРЫ И МУТАТОРЫ
  */

  // Возвращает кол-во строк
  int get_rows() const;

  // Возвращает кол-во столбцов
  int get_cols() const;

  // Изменяет кол-во строк
  void set_rows(int rows);

  // Изменяет кол-во столбцов
  void set_cols(int cols);
};

#endif