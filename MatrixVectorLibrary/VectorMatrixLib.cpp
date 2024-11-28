#include <cmath>
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
  v = vector<double>(l);
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
  if (s == 0)
    throw domain_error("Division by zero");
  for (int i = 0; i < this->v.size(); i++)
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

public:
  Matrix(int m, int n);
  Matrix(const Vector &v);
  Matrix(std::initializer_list<std::vector<double>> l);
  int NRows() const;
  int NCols() const;
  double *operator[](int i);
  const double *operator[](int i) const;
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
};
Matrix::Matrix(int m, int n) {
  if (m <= 0)
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
  testIndex(i-1);
  return m.at(i-1)(j);
}
double Matrix::operator()(int i, int j) const {
  testIndex(i-1);
  return m.at(i-1)(j);
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
      if (std::abs(val) < eps) {
        cout << setw(width) << 0;
      } else {
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
  Matrix rez(NCols(), NRows());
  for (int i = 0; i < NRows(); i++) {
    for (int j = 0; j < NCols(); j++) {
      rez[j][i] = m[i][j];
    }
  }
  *this = rez;
}

int main() {
    try {
        Vector v1(3);
        v1[0] = 1.0;
        v1[1] = 2.0;
        v1[2] = 3.0;

        Vector v2({4.0, 5.0, 6.0});

        cout << "Vektor v1: ";
        v1.Print(' ');
        cout << "\nVektor v2: ";
        v2.Print(' ');

        Vector zbir = v1 + v2;
        cout << "\nVektor v1 + v2: ";
        zbir.Print(' ');

        Vector razlika = v1 - v2;
        cout << "\nVektor v1 - v2: ";
        razlika.Print(' ');

        Vector pomnozen = v1 * 2.0;
        cout << "\nVektor v1 * 2.0: ";
        pomnozen.Print(' ');

        double skalarniProdukt = v1 * v2;
        cout << "\nSkalarni produkt v1 i v2: " << skalarniProdukt << endl;

        double norma_v1 = VectorNorm(v1);
        cout << "Norma v1: " << norma_v1 << endl;

        Matrix m1(2, 3);
        m1[0][0] = 1.0;
        m1[0][1] = 2.0;
        m1[0][2] = 3.0;
        m1[1][0] = 4.0;
        m1[1][1] = 5.0;
        m1[1][2] = 6.0;

        Matrix m2({{7.0, 8.0, 9.0}, {10.0, 11.0, 12.0}});

        cout << "Matrica m1:" << endl;
        m1.Print();

        cout << "Matrica m2:" << endl;
        m2.Print();

        Matrix zbirMatrica = m1 + m2;
        cout << "Zbir matrica m1 i m2:" << endl;
        zbirMatrica.Print();

        Matrix razlikaMatrica = m1 - m2;
        cout << "Razlika matrica m1 i m2:" << endl;
        razlikaMatrica.Print();

        Matrix pomnozenaMatrica = m1 * 2.0;
        cout << "Matrica m1 * 2.0:" << endl;
        pomnozenaMatrica.Print();

        Matrix umnozakMatrica = m1 * m2;
        cout << "Umnožak matrica m1 i m2:" << endl;
        umnozakMatrica.Print();

        Vector v3({1.0, 2.0, 3.0});
        Vector rezultat = m1 * v3;
        cout << "Matrično-vektorsko množenje:" << endl;
        rezultat.Print();
    
        Vector neispravanVektor(0); 
        Matrix neispravnaMatrica({{1.0, 2.0, 3.0}, {4.0, 5.0}});
    } catch (const exception &e) {
        cerr << "Izuzetak: "<< e.what() << endl;
    }

    return 0;
}
