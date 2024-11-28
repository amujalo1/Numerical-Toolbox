#include <cmath>
#include <exception>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
#include <stdexcept>
#include <typeinfo>
#include <vector>
using namespace std;

class Vector {
  vector<double> v;
  void testIndex(int i) const {
    if (i < 0 || i > (v.size() - 1)) {
      throw range_error("Invalid index");
    }
  }

public:
  explicit Vector(int n);
  Vector(const vector<double> &v_cpy);
  Vector(std::initializer_list<double> l);
  int NElems() const;
  double &operator[](int i);
  double operator[](int i) const;
  double &operator()(int i);
  double operator()(int i) const;
  double Norm() const;
  friend double VectorNorm(const Vector &v);
  double GetEpsilon() const;
  void Print(char separator = '\n', double eps = -1) const;
  friend void PrintVector(const Vector &v, char separator = '\n',
                          double eps = -1) {
    v.Print(separator, eps);
  }
  void swap(Vector &temp) { v.swap(temp.v); }
  friend Vector operator+(const Vector &v1, const Vector &v2);
  Vector &operator+=(const Vector &v);
  friend Vector operator-(const Vector &v1, const Vector &v2);
  Vector &operator-=(const Vector &v);
  friend Vector operator*(double s, const Vector &v);
  friend Vector operator*(const Vector &v, double s);
  Vector &operator*=(double s);
  friend double operator*(const Vector &v1, const Vector &v2);
  friend Vector operator/(const Vector &v, double s);
  Vector &operator/=(double s);

  const std::vector<double> &getV() const { return v; }

  void Chop(double eps = -1) {
    if (eps < 0)
      eps = GetEpsilon();
    for (double &d : v)
      if (fabs(d) < eps)
        d = 0;
  }
  bool EqualTo(const Vector &v, double eps = -1) const {
    if (NElems() != v.NElems())
      return false;
    if (eps < 0)
      eps = GetEpsilon();
    for (int i = 0; i < NElems(); i++) {
      if ((*this)[i] - v[i] > eps) {
        return false;
      }
    }
    return true;
  }
};
Vector::Vector(int n) {
  if (n <= 0)
    throw range_error("Bad dimension");
  v = vector<double>(n, 0.);
}
Vector::Vector(const vector<double> &v_cpy) { v = v_cpy; }

Vector::Vector(std::initializer_list<double> l) {
  if (l.size() == 0)
    throw range_error("Bad dimension");
  v.resize(l.size());
  copy(l.begin(), l.end(), v.begin());
}
int Vector::NElems() const { return v.size(); }
double &Vector::operator[](int i) {
  testIndex(i);
  return v[i];
}
double Vector::operator[](int i) const {
  testIndex(i);
  return v[i];
}
double &Vector::operator()(int i) {
  testIndex(i - 1);
  return v.at(i - 1);
}
double Vector::operator()(int i) const {
  testIndex(i - 1);
  return v.at(i - 1);
}
double Vector::Norm() const {
  double sum = 0.0;
  for (const double &d : v)
    sum += d * d;
  return sqrt(sum);
  ;
}
double VectorNorm(const Vector &v) { return v.Norm(); }
double Vector::GetEpsilon() const {
  return 10 * numeric_limits<double>::epsilon() * Norm();
}
void Vector::Print(char separator, double eps) const {
  if (eps == -1) {
    eps = GetEpsilon();
  }

  for (size_t i = 0; i < v.size(); i++) {
    if (std::abs(v[i]) < eps) {
      std::cout << "0";
    } else {
      std::cout << v[i];
    }
    if (i < v.size() - 1) {
      std::cout << separator;
    }
  }
}
Vector operator+(const Vector &v1, const Vector &v2) {
  Vector rez = v1;
  rez += v2;
  return rez;
}
Vector &Vector::operator+=(const Vector &v) {
  if (this->v.size() != v.v.size())
    throw domain_error("Incompatible formats");
  for (int i = 0; i < this->v.size(); i++)
    this->v[i] += v.v[i];
  return *this;
}
Vector operator-(const Vector &v1, const Vector &v2) {
  Vector rez = v1;
  rez -= v2;
  return rez;
}
Vector &Vector::operator-=(const Vector &v) {
  if (this->v.size() != v.v.size())
    throw domain_error("Incompatible formats");
  for (int i = 0; i < this->v.size(); i++)
    this->v[i] -= v.v[i];
  return *this;
}
Vector operator*(double s, const Vector &v) { return v * s; }
Vector operator*(const Vector &v, double s) {
  Vector rez = v;
  rez *= s;
  return rez;
}
Vector &Vector::operator*=(double s) {
  for (int i = 0; i < this->v.size(); i++)
    this->v[i] *= s;
  return *this;
}
double operator*(const Vector &v1, const Vector &v2) {
  if (v1.v.size() != v2.v.size())
    throw domain_error("Incompatible formats");
  double rez = 0;
  for (int i = 0; i < v1.v.size(); i++)
    rez += v1[i] * v2[i];
  return rez;
}
Vector operator/(const Vector &v, double s) {
  Vector rez = v;
  rez /= s;
  return rez;
}
Vector &Vector::operator/=(double s) {
  if (s < std::numeric_limits<double>::epsilon())
    throw domain_error("Division by zero");
  for (int i = 0; i < NElems(); i++)
    this->v[i] /= s;
  return *this;
}
class Matrix {
  vector<Vector> m;
  void testIndex(int i) const {
    if (i < 0 || i > (m.size() - 1)) {
      throw range_error("Invalid index");
    }
  }
  Matrix VektorKol(Vector v) {
    Matrix mat(v.NElems(), 1);
    for (int i = 0; i < v.NElems(); i++) {
      mat[i][0] = v[i];
    }
    return mat;
  }

public:
  Matrix(int m, int n);
  Matrix(const Vector &v);
  Matrix(std::initializer_list<std::vector<double>> l);
  int NRows() const;
  int NCols() const;
  double *operator[](int i);
  const double *operator[](int i) const;
  const Vector &dajKlasuVector(int i) const;
  double &operator()(int i, int j);
  double operator()(int i, int j) const;
  double Norm() const;
  friend double MatrixNorm(const Matrix &m);
  double GetEpsilon() const;
  void Print(int width = 10, double eps = -1) const;
  friend void PrintMatrix(const Matrix &m, int width = 10, double eps = -1) {
    m.Print(width, eps);
  }
  friend Matrix operator+(const Matrix &m1, const Matrix &m2);
  Matrix &operator+=(const Matrix &m);
  friend Matrix operator-(const Matrix &m1, const Matrix &m2);
  Matrix &operator-=(const Matrix &m);
  friend Matrix operator*(double s, const Matrix &m);
  friend Matrix operator*(const Matrix &m, double s);
  Matrix &operator*=(double s);
  friend Matrix operator*(const Matrix &m1, const Matrix &m2);
  Matrix &operator*=(const Matrix &m);
  friend Vector operator*(const Matrix &m, const Vector &v);
  friend Matrix Transpose(const Matrix &m);
  void Transpose();

  void Chop(double eps = -1) {
    if (eps < 0)
      eps = GetEpsilon();
    for (Vector &red : m)
      red.Chop();
  }
  bool EqualTo(const Matrix &m, double eps = -1) const {

    if (this->NRows() != m.NRows() || this->NCols() != m.NCols()) {
      return false;
    }
    if (eps < 0) {
      eps = GetEpsilon();
    }

    for (int i = 0; i < this->NRows(); i++) {
      if (!this->m[i].EqualTo(m.dajKlasuVector(i), eps)) {
        return false;
      }
    }
    return true;
  }
  friend Matrix LeftDiv(Matrix m1, Matrix m2);
  friend Vector LeftDiv(Matrix m, Vector v);
  friend Matrix operator/(const Matrix &m, double s);
  Matrix &operator/=(double s);
  friend Matrix operator/(Matrix m1, Matrix m2);
  Matrix &operator/=(Matrix m);
  double Det() const;
  friend double Det(Matrix m);
  void Invert();
  friend Matrix Inverse(Matrix m);
  void ReduceToRREF();
  friend Matrix RREF(Matrix m);
  int Rank() const;
  friend int Rank(Matrix m);
};
const Vector &Matrix::dajKlasuVector(int i) const {
  testIndex(i);
  return m[i];
}
Matrix::Matrix(int m, int n) {
  if (m <= 0 || n <= 0)
    throw range_error("Bad dimension");
  this->m = vector<Vector>(m, Vector(n));
}
Matrix::Matrix(const Vector &v) {
  for (int i = 0; i < v.NElems(); i++)
    m.push_back(Vector({v[i]}));
}
Matrix::Matrix(std::initializer_list<std::vector<double>> l) {
  if (l.size() == 0)
    throw range_error("Bad dimension");
  int d = l.begin()->size();
  for (auto a : l) {
    if (a.size() == 0)
      throw range_error("Bad dimension");
    else if (d != a.size())
      throw logic_error("Bad matrix");
    m.push_back(a);
  }
}
int Matrix::NRows() const { return m.size(); }
int Matrix::NCols() const { return m[0].NElems(); }
double *Matrix::operator[](int i) {
  testIndex(i);
  return &m.at(i)[0];
}
const double *Matrix::operator[](int i) const {
  testIndex(i);
  return m[i].getV().data();
}
double &Matrix::operator()(int i, int j) {
  testIndex(i - 1);
  return m.at(i - 1)(j);
}
double Matrix::operator()(int i, int j) const {
  testIndex(i - 1);
  return m.at(i - 1)(j);
}
double Matrix::Norm() const {
  double sum = 0.0;
  for (const Vector &v : m)
    for (const double &val : v.getV())
      sum += val * val;

  return sqrt(sum);
}
double MatrixNorm(const Matrix &m) { return m.Norm(); }
double Matrix::GetEpsilon() const {
  return 10 * Norm() * numeric_limits<double>::epsilon();
}
void Matrix::Print(int width, double eps) const {
  if (eps < 0.) {
    eps = GetEpsilon();
  }

  for (const Vector &v : m) {
    for (const double &val : v.getV()) {
      if (std::fabs(val) < eps) {
        cout << setw(width) << '0';
      } else if (val < 0)
        std::cout << std::setw(width + 1) << val;
      else {
        cout << setw(width) << val;
      }
    }
    std::cout << std::endl;
  }
}

Matrix operator+(const Matrix &m1, const Matrix &m2) {
  Matrix rez = m1;
  rez += m2;
  return rez;
}
Matrix &Matrix::operator+=(const Matrix &m) {
  if (this->m.size() != m.NRows() || this->m.at(0).NElems() != m.NCols()) {
    throw domain_error("Incompatible formats");
    return *this;
  }
  for (int i = 0; i < m.NRows(); i++) {
    this->m[i] += m.m[i];
  }
  return *this;
}

Matrix operator-(const Matrix &m1, const Matrix &m2) {
  Matrix rez = m1;
  rez -= m2;
  return rez;
}
Matrix &Matrix::operator-=(const Matrix &m) {
  if (this->m.size() != m.NRows() || this->m.at(0).NElems() != m.NCols()) {
    throw domain_error("Incompatible formats");
    return *this;
  }
  for (int i = 0; i < m.NRows(); i++) {
    this->m[i] -= m.m[i];
  }
  return *this;
}

Matrix operator*(double s, const Matrix &m) {
  Matrix rez = m;
  rez *= s;
  return rez;
}
Matrix operator*(const Matrix &m, double s) { return s * m; }
Matrix &Matrix::operator*=(double s) {
  for (int i = 0; i < this->NRows(); i++) {
    this->m[i] *= s;
  }
  return *this;
}

Matrix operator*(const Matrix &m1, const Matrix &m2) {
  Matrix rez = m1;
  rez *= m2;
  return rez;
}
Matrix &Matrix::operator*=(const Matrix &m) {
  if (this->NCols() != m.NRows()) {
    throw domain_error("Incompatible formats");
  }
  Matrix rez(this->NRows(), m.NCols());
  for (int i = 0; i < this->NRows(); i++) {
    for (int j = 0; j < m.NCols(); j++) {
      for (int k = 0; k < this->NCols(); k++) {
        rez[i][j] += (*this)[i][k] * m[k][j];
      }
    }
  }
  (*this) = rez;
  return *this;
}
Vector operator*(const Matrix &m, const Vector &v) {
  if (m.NCols() != v.NElems()) {
    throw std::domain_error("Incompatible formats");
  }

  Vector rez(m.NRows());
  for (int i = 0; i < m.NRows(); i++) {
    rez[i] = 0.0;
    for (int j = 0; j < m.NCols(); j++) {
      rez[i] += m[i][j] * v[j];
    }
  }
  return rez;
}

Matrix Transpose(const Matrix &m) {
  Matrix rez = m;
  rez.Transpose();
  return rez;
}
void Matrix::Transpose() {
  if (NCols() == NRows()) {
    for (int i = 0; i < NRows(); i++)
      for (int j = i + 1; j < NCols(); j++)
        std::swap(m[i][j], m[j][i]);
  } else {
    Matrix pomocna(NCols(), NRows());
    for (int i = 0; i < NRows(); i++)
      for (int j = 0; j < NCols(); j++)
        pomocna[j][i] = m[i][j];
    *this = pomocna;
  }
}
Matrix LeftDiv(Matrix m1, Matrix m2) {
  if (m1.NRows() != m1.NCols()) {
    throw domain_error("Divisor matrix is not square");
  }

  if (m1.NRows() != m2.NRows()) {
    throw domain_error("Incompatible formats");
  }
  double eps = m1.GetEpsilon();
  int n = m1.NRows();
  int m = m2.NCols();
  for (int k = 0; k < n; ++k) {
    int p = k;
    for (int i = k + 1; i < n; ++i) {
      if (fabs(m1[i][k]) > fabs(m1[p][k])) {
        p = i;
      }
    }
    if (fabs(m1[p][k]) < eps) {
      throw domain_error("Divisor matrix is singular");
    }
    if (p != k) {
      for (int i = 0; i < m1.NCols(); ++i) {
        swap(m1[k][i], m1[p][i]);
      }
      for (int i = 0; i < m2.NCols(); ++i) {
        swap(m2[k][i], m2[p][i]);
      }
    }
    for (int i = k + 1; i < n; ++i) {
      double u = m1[i][k] / m1[k][k];
      for (int j = k + 1; j < n; ++j) {
        m1[i][j] -= u * m1[k][j];
      }
      for (int j = 0; j < m; ++j) {
        m2[i][j] -= u * m2[k][j];
      }
    }
  }
  Matrix x(n, m);
  for (int k = 0; k < m; ++k) {
    for (int i = n - 1; i >= 0; --i) {
      double s = m2[i][k];
      for (int j = i + 1; j < n; ++j) {
        s -= m1[i][j] * x[j][k];
      }
      x[i][k] = s / m1[i][i];
    }
  }
  return x;
}
Vector LeftDiv(Matrix m, Vector v) {
  Matrix temp(v.NElems(), 1);
  temp = temp.VektorKol(v);
  temp = LeftDiv(m, temp);
  Vector rez(v.NElems());
  for (int i = 0; i < v.NElems(); i++)
    rez[i] = temp[i][0];
  return rez;
}
Matrix operator/(const Matrix &m, double s) {
  Matrix rez = m;
  rez /= s;
  return rez;
}
Matrix &Matrix::operator/=(double s) {
  if (s < std::numeric_limits<double>::epsilon()) {
    throw domain_error("Division by zero");
  }
  for (int i = 0; i < NRows(); ++i) {
    for (int j = 0; j < NCols(); ++j) {
      this->operator[](i)[j] /= s;
    }
  }
  return *this;
}
Matrix &Matrix::operator/=(Matrix m) {
  if (m.NRows() != m.NCols()) {
    throw domain_error("Divisor matrix is not square");
  } else if ((*this).NCols() != m.NCols()) {
    throw domain_error("Incompatible formats");
  } else if (std::fabs(m.Det()) < std::numeric_limits<double>::epsilon())
    throw std::domain_error("Divisor matrix is singular");

  /*Matrix m_ = (*this);
  m_.Transpose();
  m.Transpose();
  (*this) = LeftDiv(m,m_);
  (*this).Transpose();*/
  double eps = (*this).GetEpsilon();
  int n = (*this).NRows();
  Matrix x((*this).NRows(), (*this).NCols());
  for (int k = 0; k < n; ++k) {
    int p = k;
    for (int i = k + 1; i < n; ++i) {
      if (fabs(m[k][i]) > fabs(m[k][p])) {
        p = i;
      }
    }
    if (p != k) {
      for (int i = 0; i < n; ++i) {
        swap(m[i][k], m[i][p]);
      }
      for (int i = 0; i < n; ++i) {
        swap((*this)[i][k], (*this)[i][p]);
      }
    }
    for (int i = k + 1; i < n; ++i) {
      double u = m[k][i] / m[k][k];
      for (int j = k + 1; j < n; ++j) {
        m[j][i] -= u * m[j][k];
      }
      for (int j = 0; j < n; ++j) {
        (*this)[j][i] -= u * (*this)[j][k];
      }
    }
  }
  for (int k = 0; k < n; ++k) {
    for (int i = n - 1; i >= 0; --i) {
      double s = (*this)[k][i];
      for (int j = i + 1; j < n; ++j) {
        s -= m[j][i] * x[k][j];
      }
      x[k][i] = s / m[i][i];
    }
  }
  (*this) = x;
  return (*this);
}
Matrix operator/(Matrix m1, Matrix m2) {
  m1 /= m2;
  return m1;
}

double Det(Matrix m) {
  if (m.NCols() != m.NRows())
    throw domain_error("Matrix is not square");
  double det = 1;
  int n = m.NRows();
  double eps = m.GetEpsilon();
  for (int k = 0; k < n; k++) {
    int p = k;
    for (int i = k + 1; i < n; i++) {
      if (fabs(m[i][k]) > fabs(m[p][k])) {
        p = i;
      }
    }
    if (fabs(m[p][k]) < eps) {
      return 0;
    }
    if (p != k) {
      for (int i = 0; i < n; ++i) {
        swap(m[k][i], m[p][i]);
      }
      det *= -1;
    }
    det = det * m[k][k];
    for (int i = k + 1; i < n; i++) {
      double u = m[i][k] / m[k][k];
      for (int j = k + 1; j < n; j++) {
        m[i][j] -= u * m[k][j];
      }
    }
  }
  return det;
}
double Matrix::Det() const { return ::Det(*this); }
void Matrix::Invert() {
  if (NCols() != NRows())
    throw domain_error("Matrix is not square");
  if (fabs(Det()) < numeric_limits<double>::epsilon())
    throw domain_error("Matrix is singular");
  int n = NRows();
  double mi;

  for (int k = 0; k < n; k++) {

    mi = (*this)[k][k];
    (*this)[k][k] = 1;
    for (int j = 0; j < n; j++)
      (*this)[k][j] /= mi;
    for (int i = 0; i < n; i++) {
      if (i != k) {
        mi = (*this)[i][k];
        (*this)[i][k] = 0;
        for (int j = 0; j < n; j++)
          (*this)[i][j] -= mi * (*this)[k][j];
      }
    }
  }
}
Matrix Inverse(Matrix m) {
  m.Invert();
  return m;
}
void Matrix::ReduceToRREF() {
  int m = NRows();
  int n = NCols();
  const double eps = GetEpsilon();
  vector<bool> w(n);
  int k = -1;
  int l = -1;
  int p = 0;
  for (int j = 0; j < n; j++) {
    w.at(j) = false;
  }
  while (k < m && l < n) {
    double v = 0;
    l++;
    k++;
    while (v < eps && l < n) {
      p = k;
      for (int i = k; i < m; ++i) {
        if (fabs((*this)[i][l]) > v) {
          v = fabs((*this)[i][l]);
          p = i;
        }
      }
      if (v < eps) {
        l++;
      }
    }
    if (l < n) {
      w[l] = true;
      if (p != k) {
        for (int j = 0; j < n; j++) {
          swap((*this)[k][j], (*this)[p][j]);
        }
      }
      double u = (*this)[k][l];
      for (int j = l; j < n; j++) {
        (*this)[k][j] /= u;
      }
      for (int i = 0; i < m; i++) {
        if (i != k) {
          u = (*this)[i][l];
          for (int j = l; j < n; j++) {
            (*this)[i][j] -= u * (*this)[k][j];
          }
        }
      }
    }
  }
}
Matrix RREF(Matrix m) {
  m.ReduceToRREF();
  return m;
}
int Matrix::Rank() const {
  Matrix temp(NRows(), NCols());
  temp += *this;
  int m = temp.NRows();
  int n = temp.NCols();
  const double eps = GetEpsilon();
  vector<bool> w(n);
  int k = -1;
  int l = -1;
  int p = 0;
  for (int j = 0; j < n; j++) {
    w.at(j) = false;
  }
  while (k < m && l < n) {
    l++;
    k++;
    double v = 0;
    while (v < eps && l < n) {
      p = k;
      for (int i = k; i < m; ++i) {
        if (fabs(temp[i][l]) > v) {
          v = fabs(temp[i][l]);
          p = i;
        }
      }
      if (v < eps) {
        l++;
      }
    }
    if (l < n) {
      w[l] = true;
      if (p != k) {
        temp.m[k].swap(temp.m[p]);
      }
      double mi = temp[k][l];
      for (int j = l; j < n; j++)
        temp[k][j] /= mi;
      for (int i = 0; i < m; i++) {
        if (i != k) {
          mi = temp[i][l];
          for (int j = l; j < n; j++)
            temp[i][j] -= mi * temp[k][j];
        }
      }
    }
  }
  return k;
}
int Rank(Matrix m) { return m.Rank(); }

class LUDecomposer {
  Matrix mat;
  Vector w;

public:
  LUDecomposer(Matrix m) : mat(m), w(m.NCols()) {
    if (mat.NCols() != mat.NRows())
      throw std::domain_error("Matrix is not square");

    int n = mat.NRows();
    for (int j = 0; j < n; j++) {

      for (int i = 0; i <= j; i++) {
        double s = mat[i][j];
        for (int k = 0; k < i; k++)
          s -= mat[i][k] * mat[k][j];

        mat[i][j] = s;
      }

      int p = j;
      for (int i = j + 1; i < n; i++) {
        double s = mat[i][j];
        for (int k = 0; k <= j - 1; k++)
          s -= mat[i][k] * mat[k][j];

        mat[i][j] = s;
        if (fabs(s) > fabs(mat[p][j]))
          p = i;
      }
      if (fabs(mat[p][j]) < mat.GetEpsilon())
        throw std::domain_error("Matrix is singular");
      if (p != j) {
        for (int r = 0; r < n; r++) {
          swap(mat[p - 1][r], mat[j][r]);
        }
      }

      w[j] = p;
      double mi = mat[j][j];
      for (int i = j + 1; i < n; i++) {
        mat[i][j] /= mi;
      }
    }
  }

  void Solve(const Vector &b, Vector &x) const;
  Vector Solve(Vector b) const;
  void Solve(Matrix &b, Matrix &x) const;
  Matrix Solve(Matrix b) const;

  Matrix GetCompactLU() const { return mat; }
  Matrix GetL() const;

  Matrix GetU() const;

  Vector GetPermuation() const {
    Vector temp(w.NElems());
    for (int i = 0; i < w.NElems(); i++)
      temp[i] = w[i];
    return temp;
  }
};

Matrix LUDecomposer::Solve(Matrix b) const {
  if (b.NRows() != mat.NRows() || b.NCols() != mat.NCols())
    throw std::domain_error("Incompatible formats");
  Matrix x(b.NRows(), b.NCols());
  Solve(b, x);
  return x;
}

void LUDecomposer::Solve(Matrix &b, Matrix &x) const {
  if (b.NRows() != x.NRows() || b.NCols() != x.NCols())
    throw std::domain_error("Incompatible formats");

  for (int m = 0; m < b.NCols(); m++) {

    for (int i = 0; i < b.NRows(); i++) {
      int p = w[i];
      long double s = b[p][m];
      b[p][m] = b[i][m];

      for (int j = 0; j < i; j++)
        s -= mat[i][j] * b[j][m];
      b[i][m] = s;
    }

    for (int i = b.NRows() - 1; i >= 0; i--) {
      double s = b[i][m];
      for (int j = i + 1; j < b.NRows(); j++)
        s -= mat[i][j] * x[j][m];
      x[i][m] = s / mat[i][i];
    }
  }
}

Vector LUDecomposer::Solve(Vector b) const {
  if (b.NElems() != mat.NRows())
    throw std::domain_error("Incompatible formats");

  for (int i = 0; i < b.NElems(); i++) {
    int p = w[i];
    long double s = b[p];
    b[p] = b[i];

    for (int j = 0; j < i; j++)
      s -= mat[i][j] * b[j];
    b[i] = s;
  }

  for (int i = b.NElems() - 1; i >= 0; i--) {
    double s = b[i];
    for (int j = i + 1; j < b.NElems(); j++)
      s -= mat[i][j] * b[j];
    b[i] = s / mat[i][i];
  }
  return b;
}

void LUDecomposer::Solve(const Vector &b, Vector &x) const {
  if (b.NElems() != mat.NRows())
    throw std::domain_error("Incompatible formats");
  else if (b.NElems() != x.NElems())
    throw std::domain_error("Incompatible formats");
  x = std::move(Solve(b));
}

Matrix LUDecomposer::GetU() const {
  Matrix U(mat.NRows(), mat.NCols());
  for (int i = 0; i < U.NRows(); i++)
    for (int j = 0; j < U.NRows(); j++)
      if (j >= i)
        U[i][j] = mat[i][j];
  return U;
}

Matrix LUDecomposer::GetL() const {
  Matrix L(mat.NRows(), mat.NCols());
  for (int i = 0; i < L.NRows(); i++)
    for (int j = 0; j <= i; j++) {
      if (i == j)
        L[i][j] = 1;
      else
        L[i][j] = mat[i][j];
    }
  return L;
}

class QRDecomposer {
  Matrix mat;
  Matrix V;
  Vector d;

public:
  QRDecomposer(Matrix m);
  void Solve(const Vector &b, Vector &x) const;
  Vector Solve(Vector b) const;
  void Solve(Matrix &b, Matrix &x) const;
  Matrix Solve(Matrix b) const;
  Vector MulQWith(Vector v) const;
  Matrix MulQWith(Matrix m) const;
  Vector MulQTWith(Vector v) const;
  Matrix MulQTWith(Matrix m) const;
  Matrix GetQ() const;
  Matrix GetR() const;
};

QRDecomposer::QRDecomposer(Matrix m)
    : mat(m), V(m.NRows(), m.NCols()), d(mat.NRows()) {

  if (mat.NRows() < mat.NCols())
    throw std::domain_error("Invalid matrix format");
  for (int k = 0; k < mat.NCols(); k++) {
    long double s = 0;
    for (int i = k; i < mat.NRows(); i++)
      s += mat[i][k] * mat[i][k];

    s = std::sqrt(s);
    long double mi = std::sqrt(s * (s + std::fabs(mat[k][k])));

    if (std::fabs(mi) < mat.GetEpsilon())
      throw std::domain_error("Matrix is singular");

    if (mat[k][k] < 0)
      s = -s;
    mat[k][k] = (mat[k][k] + s) / mi;

    for (int i = k + 1; i < mat.NCols(); i++)
      mat[i][k] = mat[i][k] / mi;

    d[k] = -s;

    for (int j = k + 1; j < mat.NCols(); j++) {
      s = 0;
      for (int i = k; i < mat.NRows(); i++)
        s += mat[i][k] * mat[i][j];

      for (int i = k; i < mat.NCols(); i++)
        mat[i][j] -= s * mat[i][k];
    }
  }
}

Matrix QRDecomposer::Solve(Matrix b) const {

  if (mat.NCols() != mat.NRows())
    throw std::domain_error("Matrix is not square");
  if (b.NRows() != mat.NRows() || b.NCols() != mat.NCols())
    throw std::domain_error("Incompatible formats");

  Matrix x(b.NRows(), b.NCols());

  for (int i = 0; i < b.NCols(); i++) {
    Vector v(b.NRows());
    for (int j = 0; j < b.NRows(); j++) {
      v[j] = b[j][i];
    }

    v = Solve(v);
    for (int j = 0; j < b.NRows(); j++) {
      x[j][i] = v[j];
    }
  }
  return x;
}

void QRDecomposer::Solve(Matrix &b, Matrix &x) const {
  if (mat.NCols() != mat.NRows())
    throw std::domain_error("Matrix is not square");
  if (b.NRows() != mat.NRows() || b.NCols() != mat.NCols())
    throw std::domain_error("Incompatible formats");
  if (b.NRows() != x.NRows() || b.NCols() != x.NCols())
    throw std::domain_error("Incompatible formats");
  x = Solve(b);
}

Vector QRDecomposer::Solve(Vector b) const {
  if (mat.NCols() != mat.NRows())
    throw std::domain_error("Matrix is not square");

  else if (b.NElems() != mat.NRows())
    throw std::domain_error("Incompatible formats");

  for (int k = 0; k < b.NElems(); k++) {
    double s = 0;
    for (int i = k; i < mat.NRows(); i++)
      s += mat[i][k] * b[i];

    for (int i = k; i < mat.NRows(); i++)
      b[i] -= s * mat[i][k];
  }

  for (int i = b.NElems() - 1; i >= 0; i--) {
    double s = b[i];
    for (int j = i + 1; j < b.NElems(); j++)
      s -= mat[i][j] * b[j];

    b[i] = s / d[i];
  }
  return b;
}

void QRDecomposer::Solve(const Vector &b, Vector &x) const {
  if (b.NElems() != mat.NRows())
    throw std::domain_error("Incompatible formats");
  else if (b.NElems() != x.NElems())
    throw std::domain_error("Incompatible formats");
  x = std::move(Solve(b));
}

Matrix QRDecomposer::MulQTWith(Matrix m) const {
  if (mat.NCols() != m.NRows())
    throw std::domain_error("Incompatible formats");

  for (int j = 0; j < m.NRows(); j++) {
    for (int k = 0; k < m.NRows(); k++) {
      double s = 0;
      for (int i = k; i < mat.NRows(); i++)
        s += mat[i][k] * m[i][j];

      for (int i = k; i < mat.NRows(); i++)
        m[i][j] -= s * mat[i][k];
    }
  }
  return m;
}

Matrix QRDecomposer::MulQWith(Matrix m) const {
  if (mat.NCols() != m.NRows())
    throw std::domain_error("Incompatible formats");

  for (int j = 0; j < m.NRows(); j++) {
    for (int k = mat.NCols() - 1; k >= 0; k--) {
      double s = 0;
      for (int i = k; i < mat.NRows(); i++)
        s += mat[i][k] * m[i][j];

      for (int i = k; i < mat.NRows(); i++)
        m[i][j] -= s * mat[i][k];
    }
  }
  return m;
}

Matrix QRDecomposer::GetQ() const {
  Matrix Q(mat.NRows(), mat.NRows());

  for (int j = 0; j < mat.NRows(); j++) {
    for (int i = 0; i < mat.NRows(); i++)
      Q[i][j] = 0;
    Q[j][j] = 1;

    for (int k = mat.NCols() - 1; k >= 0; k--) {
      double s = 0;
      for (int i = k; i < mat.NRows(); i++)
        s += mat[i][k] * Q[i][j];

      for (int i = k; i < mat.NRows(); i++)
        Q[i][j] -= s * mat[i][k];
    }
  }
  return Q;
}

Matrix QRDecomposer::GetR() const {
  Matrix R(mat.NRows(), mat.NCols());
  for (int i = 0; i < mat.NRows(); i++)
    for (int j = 0; j < mat.NCols(); j++) {
      if (i == j)
        R[i][j] = d[i];
      else if (i > j)
        R[i][j] = 0;
      else
        R[i][j] = mat[i][j];
    }
  return R;
}

Vector QRDecomposer::MulQWith(Vector v) const {
  if (mat.NCols() != mat.NRows())
    throw std::domain_error("Matrix is not square");

  else if (v.NElems() != mat.NRows())
    throw std::domain_error("Incompatible formats");

  for (int k = v.NElems() - 1; k >= 0; k--) {
    double s = 0;
    for (int i = k; i < mat.NRows(); i++)
      s += mat[i][k] * v[i];

    for (int i = k; i < mat.NRows(); i++)
      v[i] -= s * mat[i][k];
  }
  return v;
}

Vector QRDecomposer::MulQTWith(Vector v) const {
  if (mat.NCols() != mat.NRows())
    throw std::domain_error("Matrix is not square");

  else if (v.NElems() != mat.NRows())
    throw std::domain_error("Incompatible formats");

  for (int k = 0; k < v.NElems(); k++) {
    double s = 0;
    for (int i = k; i < mat.NRows(); i++)
      s += mat[i][k] * v[i];

    for (int i = k; i < mat.NRows(); i++)
      v[i] -= s * mat[i][k];
  }
  return v;
}

int main() {
  try {
    cout << "<<Testiranje lijevo i desno matricno djeljenje, metode i funkciju "
            "inverse...>> "
         << endl;
    Matrix m{{2, 4, 6}, {1, 1, 1}};
    m = m / 2;
    Matrix r{{1, 2, 3}, {0.5, 0.5, 0.5}};
    if (m.EqualTo(r))
      cout << "Tesitranje Matrix /=(skalar) je USPJESNO" << endl;
    else
      cout << "Tesitranje Matrix /=(skalar) je USPJESNO";

    try {
      m /= 0;
      return 0;
    } catch (const domain_error &e) {
      cout << "Tesitranje Matrix /=(skalar) je USPJESNO bacio izuzetak-"
           << e.what() << endl;
    }
    Matrix m2{{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};
    Matrix m1{{1, 1, 1}, {3, 6, 9}, {2, 1, 2}};
    try {
      m1 /= m2;
    } catch (const domain_error &e) {
      cout << "Tesitranje Matrix /=(Matrix) je USPJESNO bacio izuzetak-"
           << e.what() << endl;
    }
    m2 = Matrix({{1, 2, 3}, {5, 5, 5}});
    m1 = Matrix({{1, 1, 11}, {123, 234, 9}, {2, 1, 2}});
    try {
      m1 /= m2;
    } catch (const domain_error &e) {
      cout << "Tesitranje Matrix /=(Matrix) je USPJESNO bacio izuzetak-"
           << e.what() << endl;
    }

    try {
      Matrix m{{1, 1, 1}, {2, 2, 2}};
      Vector v{1, 2, 3};
      Matrix r = LeftDiv(m, v);
    } catch (const domain_error &e) {
      cout << "Tesitranje Matrix::LeftDiv(Matrix,Vector) je USPJESNO bacio "
              "izuzetak-"
           << e.what() << endl;
    }
    try {
      Matrix m{{1, 2, 3}, {1, 2, 3}, {1, 2, 3}};
      Vector v{1, 2, 3, 4};
      Matrix r = LeftDiv(m, v);
    } catch (const domain_error &e) {
      cout << "Tesitranje Matrix::LeftDiv(Matrix,Vector) je USPJESNO bacio "
              "izuzetak-"
           << e.what() << endl;
    }

    m = Matrix({{1, 2, 3}, {3, 2, 1}, {2, 1, 3}});
    Vector v{2, 1, 2};
    r = LeftDiv(m, v);
    Vector cmp{0.0833333, 0.0833333, 0.583333};
    if (r.EqualTo(cmp, 1))
      cout << "Tesitranje Matrix::LeftDiv(Matrix,Vector) je USPJESAN" << endl;
    else
      cout << "Tesitranje Matrix::LeftDiv(Matrix,Vector) nije USPJESAN" << endl;

    try {
      Matrix m1{{1, -2}, {-3, 6}};
      Matrix m2{{1, 2}, {2, 3}};
      Matrix r = LeftDiv(m2, m1);
    } catch (const domain_error &e) {
      cout << "Tesitranje Matrix::LeftDiv(Matrix,Matrix) je USPJESNO bacio "
              "izuzetak-"
           << e.what() << endl;
    }

    {
      Matrix m1{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
      Matrix m2{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
      Matrix r = LeftDiv(m2, m1);
      if (r.EqualTo(m1))
        cout << "Tesitranje Matrix::LeftDiv(Matrix,Matrix) je USPJESAN" << endl;
      else
        cout << "Tesitranje Matrix::LeftDiv(Matrix,Matrix) nije USPJESAN"
             << endl;
    }

    try {
      Matrix m1{{1, 2, 3}, {6, 5, 4}, {7, 6, 7}};
      Matrix m2{{14, 10}, {3, 0}};
      Matrix r = LeftDiv(m1, m2);
    } catch (const domain_error &e) {
      cout << "Tesitranje Matrix::LeftDiv(Matrix,Matrix) je USPJESNO bacio "
              "izuzetak-"
           << e.what() << endl;
    }

    try {
      Matrix testDet{{1, 2, 3}, {4, 5, 6}};
      testDet.Det();
    } catch (const domain_error &e) {
      cout << "Tesitranje Matrix::Det(Matrix) je USPJESNO bacio izuzetak-"
           << e.what() << endl;
    }

    Matrix s{{1, 2, 3}, {4, 25, 6}, {7, 8, 9}};
    if (fabs(s.Det()) - 240 < 1e-10)
      cout << "Tesitranje Matrix::Det(Matrix) je USPJESAN" << endl;
    else
      cout << "Tesitranje Matrix::Det(Matrix) nije USPJESAN" << endl;

    {
      Matrix m{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
      m.ReduceToRREF();
      Matrix r{{1, 0, -1}, {0, 1, 2}, {0, 0, 0}};
      if (m.EqualTo(r))
        cout << "Tesitranje Matrix::ReduceToRREF() je USPJESAN" << endl;
      else
        cout << "Tesitranje Matrix::ReduceToRREF() nije USPJESAN" << endl;

      if (RREF(m).EqualTo(r))
        cout << "Tesitranje firend ReduceToRREF() je USPJESAN" << endl;
      else
        cout << "Tesitranje firend ReduceToRREF() nije USPJESAN" << endl;
    }

    {
      Matrix m{{3, 0, 2}, {0, 0, 1}, {2, -2, 1}};
      Matrix inv{{0.2, -0.2, 0.2}, {0.2, 0.3, -0.3}, {0, 1, 0}};
      m.Invert();
      if (inv.EqualTo(m))
        cout << "Tesitranje Matrix::ReduceToRREF() je USPJESAN" << endl;
      else
        cout << "Tesitranje Matrix::ReduceToRREF() nije USPJESAN" << endl;
    }

    try {
      Matrix m{{3, 2, 1}, {3, 1, 1}};
      m.Invert();
    } catch (const domain_error &e) {
      cout << "Tesitranje Matrix::Invert() je USPJESNO bacio izuzetak-"
           << e.what() << endl;
    }

    try {
      Matrix m{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
      m.Invert();
    } catch (const domain_error &e) {
      cout << "Tesitranje Matrix::Invert() je USPJESNO bacio izuzetak-"
           << e.what() << endl;
    }

    {
      Matrix m{{3, 0, 2}, {0, 0, 1}, {2, -2, 1}};
      Matrix inv{{0.2, -0.2, 0.2}, {0.2, 0.3, -0.3}, {0, 1, 0}};

      if (inv.EqualTo(Inverse(m)))
        cout << "Tesitranje friend ReduceToRREF() je USPJESAN" << endl;
      else
        cout << "Tesitranje friend ReduceToRREF() nije USPJESAN" << endl;
    }
    {
      Matrix m{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
      try {
        Inverse(m);
      } catch (const domain_error &e) {
        cout << "Tesitranje friend Inverse() je USPJESNO bacio izuzetak-"
             << e.what() << endl;
      }
    }

  } catch (...) {
    cerr << "Nesto ne valja" << endl;
  }
  cout << "<<Testiranje Vector::EqualTo>> ";
  Vector equal{1.5, 3.5, 2.75, 0.5};
  if (!equal.EqualTo({1.65, 3.7, 2.75, 0.5}, 0.11)) {
    return 0;
  }
  if (equal.EqualTo({1.4, 3.5, 2.75, 0.5})) {
    return 0;
  }
  cout << "USPJESNO" << endl;

  cout << "<<Testiranje Vecto::Chop>> ";
  Vector chop{0.2, 66, 0.24, 0.3, -0.24, -0.4, -0.2};
  chop.Chop(0.25);
  if (!chop.EqualTo({0, 66, 0, 0.3, 0, -0.4, 0})) {
    return 0;
  }
  if (chop.EqualTo({0, 66, 0, 0.3, 0, -0.4, -0.2})) {
    return 0;
  }
  cout << "USPJESNO" << endl;

  std::cout << "<<Testiranje LU-faktorizacije>>" << endl;
  try {
    Matrix mat1{{0, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    LUDecomposer s(mat1);

    try {
      LUDecomposer(Matrix{{0, 0}});
      return 0;
    } catch (const domain_error &e) {
      cout << "LU-konstruktor je bacio ispravan izuzetak: " << e.what() << endl;
    }

    try {
      LUDecomposer(Matrix{{3, 12}, {2, 8}});
      return 0;
    } catch (const domain_error &e) {
      cout << "LU-konstruktor je bacio ispravan izuzetak: " << e.what() << endl;
    }

    try {
      // los format
      LUDecomposer(Matrix{{1, 2}, {3, 4}})
          .Solve(Vector({1, 2, 3, 4, 5, 6, 7, 8}));
      return 0;
    } catch (const domain_error &e) {
      cout << "LU-konstruktor je bacio ispravan izuzetak: " << e.what() << endl;
    }

    try {
      LUDecomposer(Matrix{{1, 2}, {3, 4}})
          .Solve(Matrix({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}));
      return 0;
    } catch (const domain_error &e) {
      cout << "LU-konstruktor je bacio ispravan izuzetak: " << e.what() << endl;
    }

    Matrix rj{{4123, 12334, 1234}, {4321, 1345, 12345}, {12451, 2321, 123451}};
    LUDecomposer solver(rj);

    if ((solver.GetL() + solver.GetU() -
         Matrix({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}))
            .EqualTo(solver.GetCompactLU()))
      cout << "Test LUDecomposer::CompactLU() je USPJESAN" << endl;
    else {
      cout << "Test LUDecomposer::CompactLU() nije USPJESAN" << endl;
      return 0;
    }

    Vector w = solver.GetPermuation();
    if (std::abs(w[0] - 1) < 5 * std::numeric_limits<double>::epsilon() ||
        std::abs(rj[0][0]) > rj.GetEpsilon())
      std::cout << "Test GetPermutation je USPJESAN" << endl;
    else {
      std::cout << "Test GetPermutation nije USPJESAN" << endl;
      return 0;
    }

    Matrix A{{0, 3, 2}, {4, 6, 1}, {3, 1, 7}};
    Vector x{1, 2, 4};
    Vector rez = A * x;
    LUDecomposer lu(A);
    lu.Solve(rez, rez);
    if (x.EqualTo(rez)) {
      cout << "Test LUDecomposer::Solve() je USPJESAN" << endl;
    } else {
      cout << "Test LUDecomposer::Solve() nije USPJESAN" << endl;
    }
    A = Matrix({{0, 3, 2}, {4, 6, 1}, {3, 1, 7}});
    Matrix X{{4, 1, 5}, {1, 2, 1}, {8, 7, 9}};
    rj = A * x;
    lu = A;
    lu.Solve(rez, rez);
    if (x.EqualTo(rez)) {
      cout << "Test LUDecomposer::Solve() je USPJESAN" << endl;
    } else {
      cout << "Test LUDecomposer::Solve() nije USPJESAN" << endl;
    }
    Matrix B{{2, 3, 5}, {2, 3, 7}, {4, 1, 8}},
        jedinicna{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    LUDecomposer lu2(B);
    if (lu2.GetCompactLU().EqualTo(lu2.GetL() + lu2.GetU() - jedinicna)) {
      cout << "Test LUDecomposer::GetCompactLU() je USPJESAN" << endl;
    } else {
      cout << "Test LUDecomposer::GetCompactLU() ni je USPJESAN" << endl;
    }
    A = Matrix({{0, 3, 2}, {4, 6, 1}, {3, 1, 7}});
    lu = A;
    if (lu.GetPermuation().NElems() == 3) {
      cout << "Test LUDecomposer::GetPermuation() je USPJESAN" << endl;
    } else {
      cout << "Test LUDecomposer::GetPermuation() ni je USPJESAN" << endl;
    }

  } catch (...) {
    cerr << "Nesto ne valja" << endl;
  }
  cout << endl;
  cout << "<Testiranje QR-faktorizacije>" << endl;
  try {
    try {
      Matrix losFormat{{1, 1, 1}, {1, 1, 1}};
      QRDecomposer q(losFormat);
      return 0;
    } catch (const domain_error &e) {
      cout << "LU-konstruktor je bacio ispravan izuzetak: " << e.what() << endl;
    }

    try {
      Matrix singularna{{1, 1, 1}, {0, 1, 0}, {1, 0, 1}};
      QRDecomposer q(singularna);
      return 0;
    } catch (const domain_error &e) {
      cout << "LU-konstruktor je bacio ispravan izuzetak: " << e.what() << endl;
    }

    try {
      Matrix Test{{12, -51}, {6, 167}, {7, 8}};
      QRDecomposer q(Test);
      Vector b{1, 2, 3, 4};
      Vector x{1, 2, 3, 4};
      q.Solve(b, x);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test Solve(b, x) je bacio ispravan izuzetak: " << e.what()
           << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}};
      QRDecomposer q(Test);
      Vector b{1, 2, 3};
      Vector x{1, 2, 3, 4};
      q.Solve(b, x);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test Solve(b, x) je bacio ispravan izuzetak: " << e.what()
           << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}, {14, 58, -33}};
      QRDecomposer q(Test);
      Vector b{1, 2, 3, 4};
      Vector x = q.Solve(b);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test Solve(b bez x) je bacio ispravan izuzetak: " << e.what()
           << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}};
      QRDecomposer q(Test);
      Vector b{1, 2, 3, 4};
      Vector x = q.Solve(b);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test Solve(b bez x) je bacio ispravan izuzetak: " << e.what()
           << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}, {14, 58, -33}};
      QRDecomposer q(Test);
      Matrix b{{-3, 5, 7}, {43, 68, -21}, {13, 36, 23}};
      Matrix x{{6, 3, 5}, {-8, 21, 45}, {31, 4, -12}};
      q.Solve(b, x);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test Solve(matrica b, matrica x) je bacio ispravan izuzetak: "
           << e.what() << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}};
      QRDecomposer q(Test);
      Matrix b{{-3, 5, 7}, {43, 68, -21}};
      Matrix x{{6, 3, 5}, {-8, 21, 45}, {31, 4, -12}};
      q.Solve(b, x);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test Solve(matrica b, matrica x) je bacio ispravan izuzetak: "
           << e.what() << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}, {14, 58, -33}};
      QRDecomposer q(Test);
      Vector b{-3, 5, 7, 21};
      Vector x = q.MulQWith(b);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test MulQWith(v) je bacio ispravan izuzetak: " << e.what()
           << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}};
      QRDecomposer q(Test);
      Vector b{-3, 5, 7, 21};
      Vector x = q.MulQWith(b);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test MulQWith(v) je bacio ispravan izuzetak: " << e.what()
           << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}, {14, 58, -33}};
      QRDecomposer q(Test);
      Matrix b{{-3, 5, 7}, {43, 68, -21}};
      Matrix x = q.MulQWith(b);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test MulQWith(matrica) je bacio ispravan izuzetak: " << e.what()
           << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}, {14, 58, -33}};
      QRDecomposer q(Test);
      Vector b{-3, 5, 7, 21};
      Vector x = q.MulQTWith(b);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test MulQTWith(v) je bacio ispravan izuzetak: " << e.what()
           << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}};
      QRDecomposer q(Test);
      Vector b{-3, 5, 7, 21};
      Vector x = q.MulQTWith(b);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test MulQTWith(v) je bacio ispravan izuzetak: " << e.what()
           << "\n";
    }

    try {
      Matrix Test{{12, -51, 4}, {6, 167, -68}, {7, 8, 9}, {14, 58, -33}};
      QRDecomposer q(Test);
      Matrix b{{-3, 5, 7}, {43, 68, -21}};
      Matrix x = q.MulQTWith(b);
      return 0;
    } catch (const domain_error &e) {
      cout << "Test MulQTWith(matrica) je bacio ispravan izuzetak: " << e.what()
           << "\n";
    }

    cout << endl;

    Matrix A{{5.54, 2.25, 1.56}, {8.437, 1.458, 5.89}, {8.30, 6.12, 1.6}};
    QRDecomposer qr(A);
    if (A.EqualTo(qr.GetQ() * qr.GetR()))
      cout << "Testiranje QR-faktorizacije je USPJESAN" << endl;
    else
      cout << "Testiranje QR-faktorizacije nije USPJESAN" << endl;

    Matrix A1{{3, 3, 4}, {5, 5, 5}, {2, 7, 5}};
    Matrix B1{{2, 2, 3}, {5, 3, 1}, {5, 1, 2}};
    QRDecomposer qr1(A1);
    Matrix x1 = qr1.Solve(B1);
    if (x1.EqualTo({{0.8, 0.56, -1.08}, {1.2, -0.16, -1.12}, {-1, 0.2, 2.4}}))
      cout << "Testiranje QRDecomposer::Solve(Matrix) je USPJESAN" << endl;
    else
      cout << "Testiranje QRDecomposer::Solve(Matrix) nije USPJESAN" << endl;

    A = Matrix({{0, 3, 2}, {4, 6, 1}, {3, 1, 7}});
    Vector x{1, 2, 4};
    Vector rez = A * x;
    qr = A;
    qr.Solve(rez, rez);
    if (x.EqualTo(rez)) {
      cout << "Testiranje QRDecomposer::Solve(Vector) je USPJESAN" << endl;
    } else {
      cout << "Testiranje QRDecomposer::Solve(Vector) nije USPJESAN" << endl;
    }
    A = Matrix({{0, 3, 2}, {4, 6, 1}, {3, 1, 7}});
    x = Vector({1, 2, 4});
    qr = A;
    Matrix Q = qr.MulQTWith(x);
    if (Q.EqualTo(Transpose(qr.GetQ()) * x)) {
      cout << "Testiranje QRDecomposer::MulQTWith(Vector) je USPJESAN" << endl;
    } else {
      cout << "Testiranje QRDecomposer::MulQTWith(Vector) je USPJESAN" << endl;
    }
    A = Matrix({{0, 3, 2}, {4, 6, 1}, {3, 1, 7}});
    Matrix X1{{4, 1, 5}, {1, 2, 1}, {8, 7, 9}};
    qr = A;
    Q = qr.MulQWith(x);
    if (Q.EqualTo(qr.GetQ() * x)) {
      cout << "Testiranje QRDecomposer::MulQTWith(Matrix) je USPJESAN" << endl;
    } else {
      cout << "Testiranje QRDecomposer::MulQTWith(Matrix) je USPJESAN" << endl;
    }
  } catch (...) {
    cerr << "Nesto ne valja" << endl;
  }

  return 0;
}
