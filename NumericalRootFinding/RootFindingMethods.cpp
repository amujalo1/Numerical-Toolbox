#include <algorithm>
#include <complex>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

enum RegulaFalsiMode { Unmodified, Illinois, Slavic, IllinoisSlavic };
template <typename FunType>
double RegulaFalsiSolve(FunType f, double a, double b,
                        RegulaFalsiMode mode = Slavic, double eps = 1e-10,
                        int maxiter = 100) {
  if (f(a) * f(b) > 0)
    throw range_error("Root must be bracketed");
  if (eps < 0 || maxiter < 0)
    throw domain_error("Invalid parameters");
  double f1, f2, f3, C, Cold;
  int n;
  function<double(double)> help = [](double x) { return x / (1 + abs(x)); };
  switch (mode) {
  case Unmodified:
    f1 = f(a);
    f2 = f(b);
    C = a;
    Cold = b;
    n = 0;
    while (abs(C - Cold) > eps) {
      if (n == maxiter)
        throw logic_error("Given accuracy has not achieved");
      Cold = C;
      C = (a * f2 - b * f1) / (f2 - f1);
      f3 = f(C);
      if (abs(f3) < eps)
        return C;
      if (f1 * f3 < 0) {
        b = a;
        f2 = f1;
      }
      a = C;
      f1 = f3;
      n++;
    }
    return C;
    break;
  case Illinois:
    f1 = f(a);
    f2 = f(b);
    C = a;
    Cold = b;
    n = 0;
    while (abs(C - Cold) > eps) {
      if (n == maxiter)
        throw logic_error("Given accuracy has not achieved");
      Cold = C;
      C = (a * f2 - b * f1) / (f2 - f1);
      f3 = f(C);
      if (abs(f3) < eps)
        return C;
      if (f1 * f3 < 0) {
        b = a;
        f2 = f1;
      } else {
        f2 /= 2;
      }
      a = C;
      f1 = f3;
      n++;
    }
    return C;
    break;
  case Slavic:
    f1 = help(f(a));
    f2 = help(f(b));
    C = a;
    Cold = b;
    n = 0;
    while (abs(C - Cold) > eps) {
      if (n == maxiter)
        throw logic_error("Given accuracy has not achieved");
      Cold = C;
      C = (a * f2 - b * f1) / (f2 - f1);
      double f3 = help(f(C));
      if (abs(f3) < eps)
        return C;
      if (f1 * f3 < 0) {
        b = C;
        f2 = f3;
      } else {
        a = C;
        f1 = f3;
      }
      n++;
    }
    return C;
    break;
  case IllinoisSlavic:
    f1 = help(f(a));
    f2 = help(f(b));
    C = a;
    Cold = b;
    n = 0;
    while (abs(C - Cold) > eps) {
      if (n == maxiter)
        throw logic_error("Given accuracy has not achieved");
      Cold = C;
      C = (a * f2 - b * f1) / (f2 - f1);
      double f3 = help(f(C));
      if (abs(f3) < eps)
        return C;
      if (f1 * f3 < 0) {
        b = C;
        f2 = f1;
      } else {
        f2 /= 2;
      }
      a = C;
      f1 = f3;
      n++;
    }
    return C;
    break;
  }
}
template <typename FunTip>
bool BracketRoot(FunTip f, double x0, double &a, double &b, double hinit = 1e-5,
                 double hmax = 1e10, double lambda = 1.4) {
  if (hinit <= 0 || hmax <= 0 || lambda <= 0)
    throw domain_error("Invalid parameters");
  double temp1 = x0;
  double f1 = f(temp1);
  double h = hinit;
  while (abs(h) < hmax) {
    double temp2 = temp1 + h;
    double f2 = f(temp2);
    if (f1 * f2 <= 0) {
      if (temp1 > temp2)
        swap(temp1, temp2);
      a = temp1;
      b = temp2;
      return true;
    }
    h *= lambda;
    temp1 = temp2;
    f1 = f2;
  }
  h = hinit;
  temp1 = x0;
  f1 = f(temp1);
  while (abs(h) < hmax) {
    double temp2 = temp1 - h;
    double f2 = f(temp2);
    if (f1 * f2 <= 0) {
      if (temp1 > temp2)
        swap(temp1, temp2);
      a = temp1;
      b = temp2;
      return true;
    }
    h *= lambda;
    temp2 = temp1;
    f2 = f1;
  }
  return false;
}

template <typename FunTip>
double RiddersSolve(FunTip f, double a, double b, double eps = 1e-10,
                    int maxiter = 100) {
  if (f(a) * f(b) > 0)
    throw range_error("Root must be bracketed");
  if (eps < 0 || maxiter < 0)
    throw domain_error("Invalid parameters");

  double f1 = f(a), f2 = f(b);
  int br(0);
  while (abs(b - a) >= eps) {
    if (br == maxiter)
      throw logic_error("Given accuracy has not achieved");

    double c = (a + b) / 2, f3 = f(c);
    if (abs(f3) < eps)
      return c;
    int sign;
    if ((f1 - f2) > 0)
      sign = 1;
    else if ((f1 - f2) < 0)
      sign = -1;
    else
      sign = 0;
    double d = c + (f3 * (c - a) * sign) / sqrt(f3 * f3 - f1 * f2);
    double f4 = f(d);
    if (abs(f4) < eps)
      return d;
    if (f3 * f4 <= 0) {
      a = c;
      b = d;
      f1 = f3;
      f2 = f4;
    } else if (f1 * f4 <= 0) {
      b = d;
      f2 = f4;
    } else {
      a = d;
      f1 = f4;
    }
    br++;
  }
  return (a + b) / 2;
}

template <typename FunTip1, typename FunTip2>
double NewtonRaphsonSolve(FunTip1 f, FunTip2 fprim, double x0,
                          double eps = 1e-10, double damping = 0,
                          int maxiter = 100) {
  if (eps < 0 || maxiter < 0)
    throw domain_error("Invalid parameters");
  else if (damping < 0 || damping >= 1)
    throw domain_error("Invalid parameters");
  double dX(numeric_limits<double>::infinity()), v(f(x0)), d(fprim(x0));
  int br(0);

  while (fabs(dX) > eps) {

    if (abs(fprim(x0)) < eps || br == maxiter || !isfinite(x0))
      throw logic_error("Convergence has not achieved");
    if (fabs(v) <= eps)
      return x0;
    dX = v / d;
    double w = v;
    v = f(x0 - dX);
    d = fprim(x0 - dX);
    while (abs(v) > abs(w) || !isfinite(v) || d == 0) {
      dX *= damping;
      v = f(x0 - dX);
      d = fprim(x0 - dX);
    }
    x0 -= dX;
    br++;
  }
  return x0;
}

bool operator==(complex<double> c1, complex<double> c2) {
  return (abs(c1.real() - c2.real()) < numeric_limits<double>::epsilon() &&
          abs(c1.imag() - c2.imag()) < numeric_limits<double>::epsilon());
}

complex<double> operator*(complex<double> c1, complex<double> c2) {
  return {c1.real() * c2.real() - c1.imag() * c2.imag(),
          c1.real() * c2.imag() + c1.imag() * c2.real()};
}

complex<double> operator*(double x, complex<double> c) {
  complex<double> temp(x, 0);
  return temp * c;
}

complex<double> operator*(complex<double> c, double x) { return x * c; }

pair<complex<double>, bool> Laguerre(vector<complex<double>> p, int n,
                                     complex<double> x, double eps,
                                     int maxiter) {
  complex<double> dX(numeric_limits<double>::infinity());
  int k = 1;
  while (abs(dX) > eps && k < maxiter) {
    complex<double> f = p[n];
    complex<double> d = 0;
    complex<double> s = 0;
    for (int i(n - 1); i >= 0; i--) {
      s = s * x + 2 * d;
      d = d * x + f;
      f = f * x + p[i];
    }
    if (f == 0)
      return {x, true};
    complex<double> r = sqrt((n - 1) * ((n - 1) * d * d - n * f * s));
    if (abs(d + r) > abs(d - r))
      dX = n * f / (d + r);
    else
      dX = n * f / (d - r);
    x -= dX;
    k++;
  }
  if (abs(dX) <= eps)
    return {x, true};
  return {x, false};
}

vector<complex<double>> PolyRoots(vector<complex<double>> coefficients,
                                  double eps = 1e-10, int maxiters = 100,
                                  int maxtrials = 10) {
  if (eps < 0 || maxiters < 0 || maxtrials < 0)
    throw domain_error("Invalid parameters");
  int n = coefficients.size() - 1, i(n), it(0);
  vector<complex<double>> a;
  while (i >= 1) {
    if (it == maxiters)
      throw logic_error("Convergence has not achieved");
    int t = 1;
    bool var(false);
    complex<double> x;
    while (!var && (t < maxtrials)) {
      x = {0.0, 0.0};
      pair<complex<double>, bool> pair =
          Laguerre(coefficients, i, x, eps, maxiters);
      var = pair.second;
      x = pair.first;
      t++;
    }
    if (!var)
      throw logic_error("Convergence has not achieved");
    if (abs(x.imag()) <= eps)
      x = x.real();
    a.push_back(x);

    complex<double> v = coefficients[i];

    for (int j = i - 1; j >= 0; j--) {
      complex<double> w = coefficients[j];
      coefficients[j] = v;
      v = w + x * v;
    }
    i--;
    it++;
  }

  sort(a.begin(), a.end(), [](complex<double> x, complex<double> y) {
    if (x.real() < y.real())
      return true;
    else if (x.real() > y.real())
      return false;
    return x.imag() < y.imag();
  });
  return a;
}

vector<complex<double>> PolyRoots(vector<double> coefficients,
                                  double eps = 1e-10, int maxiters = 100,
                                  int maxtrials = 10) {
  if (eps < 0 || maxiters < 0 || maxtrials < 0)
    throw domain_error("Invalid parameters");
  int n = coefficients.size() - 1, i(n), it(0);
  vector<complex<double>> a(n + 1);
  vector<complex<double>> complexCoefficients;
  for (const double &realPart : coefficients) {
    complexCoefficients.push_back(complex<double>(realPart, 0.0));
  }

  while (i >= 1) {
    if (it == maxiters)
      throw logic_error("Convergence has not achieved");
    int t = 1;
    bool var(false);
    complex<double> x;
    while (!var && (t < maxtrials)) {
      x = {1, 1};
      pair<complex<double>, bool> pair =
          Laguerre(complexCoefficients, i, x, eps, maxiters);
      var = pair.second;
      x = pair.first;
      t++;
    }
    if (!var)
      throw logic_error("Convergence has not achieved");
    if (abs(x.imag()) <= eps) {
      x = x.real();
      a[i] = x;
      double v = complexCoefficients[i].real();
      for (int j = i - 1; j >= 0; j--) {
        double w = complexCoefficients[j].real();
        complexCoefficients[j] = v;
        v = w + x.real() * v;
      }
      i--;
    } else {
      a[i] = x;
      a[i - 1] = conj(x);
      double alfa = 2 * x.real(), beta = abs(x);
      beta *= beta;
      double u = complexCoefficients[i].real();
      double v = complexCoefficients[i - 1].real() + alfa * u;
      for (int j = i - 2; j >= 0; j--) {
        double w = complexCoefficients[j].real();
        complexCoefficients[j] = u;
        u = v;
        v = w + alfa * v - beta * complexCoefficients[j].real();
      }
      i -= 2;
    }
    it++;
  }
  a.erase(a.begin());
  sort(a.begin(), a.end(), [](complex<double> c1, complex<double> c2) {
    if (abs(c1.real() - c2.real()) > numeric_limits<double>::epsilon())
      return c1.real() < c2.real();
    return c1.imag() < c2.imag();
  });
  return a;
}

// Funkcija za pronalaženje korijena jednadžbe
double kvadratnaFunkcija(double x) { return x * x * x - 4; }

// Test funkcija za RegulaFalsiSolve sa podrazumijevanim načinom rada (Slavic)
void testirajRegulaFalsiPodrazumijevano() {
  try {
    double korijen = RegulaFalsiSolve(kvadratnaFunkcija, -3, 3);
    cout << "Pronađen korijen sa podrazumijevanim načinom rada: " << korijen
         << endl;
  } catch (const exception &e) {
    cerr << "Izuzetak: " << e.what() << endl;
  }
}

// Test funkcija za RegulaFalsiSolve sa Illinois načinom rada
void testirajRegulaFalsiIllinois() {
  try {
    double korijen = RegulaFalsiSolve(kvadratnaFunkcija, -3, 3, Illinois);
    cout << "Pronađen korijen sa Illinois načinom rada: " << korijen << endl;
  } catch (const exception &e) {
    cerr << "Izuzetak: " << e.what() << endl;
  }
}

// Test funkcija za RegulaFalsiSolve sa Slavic načinom rada
void testirajRegulaFalsiSlavic() {
  try {
    double korijen = RegulaFalsiSolve(kvadratnaFunkcija, -3, 3, Slavic);
    cout << "Pronađen korijen sa Slavic načinom rada: " << korijen << endl;
  } catch (const exception &e) {
    cerr << "Izuzetak: " << e.what() << endl;
  }
}

// Test funkcija za RegulaFalsiSolve sa IllinoisSlavic načinom rada
void testirajRegulaFalsiIllinoisSlavic() {
  try {
    double korijen = RegulaFalsiSolve(kvadratnaFunkcija, -3, 3, IllinoisSlavic);
    cout << "Pronađen korijen sa IllinoisSlavic načinom rada: " << korijen
         << endl;
  } catch (const exception &e) {
    cerr << "Izuzetak: " << e.what() << endl;
  }
}

// Implementacija funkcije za testiranje BracketRoot
void testirajBracketRoot() {
  try {
    double a, b;
    bool success = BracketRoot(kvadratnaFunkcija, -3, a, b);
    if (success) {
      cout << "BracketRoot succeeded. Bracket: [" << a << ", " << b << "]"
           << endl;
    } else {
      cerr << "BracketRoot failed." << endl;
    }
  } catch (const exception &e) {
    cerr << "Izuzetak pri pozivu BracketRoot: " << e.what() << endl;
  }
}

// Implementacija funkcije za testiranje RiddersSolve
void testirajRiddersSolve() {
  try {
    double korijen = RiddersSolve(kvadratnaFunkcija, -3, 3);
    cout << "Pronađen korijen RiddersSolve: " << korijen << endl;
  } catch (const exception &e) {
    cerr << "Izuzetak pri pozivu RiddersSolve: " << e.what() << endl;
  }
}

// Implementacija funkcije za testiranje NewtonRaphsonSolve
void testirajNewtonRaphsonSolve() {
  try {
    double korijen = NewtonRaphsonSolve(
        kvadratnaFunkcija, [](double x) { return 3 * x * x; }, 2);
    cout << "Pronađen korijen NewtonRaphsonSolve: " << korijen << endl;
  } catch (const exception &e) {
    cerr << "Izuzetak pri pozivu NewtonRaphsonSolve: " << e.what() << endl;
  }
}

// Implementacija funkcije za testiranje PolyRoots sa kompleksnim koeficijentima
void testirajPolyRootsComplex() {
  try {
    vector<complex<double>> coefficients = {1, 0, -4};
    vector<complex<double>> roots = PolyRoots(coefficients);
    cout << "Polinomski korijeni (sa kompleksnim koeficijentima): ";
    for (const auto &root : roots) {
      cout << root << " ";
    }
    cout << endl;
  } catch (const exception &e) {
    cerr << "Izuzetak pri pozivu PolyRoots sa kompleksnim koeficijentima: "
         << e.what() << endl;
  }
}

// Implementacija funkcije za testiranje PolyRoots sa realnim koeficijentima
void testirajPolyRootsReal() {
  try {
    vector<double> realCoefficients = {1, 0, -4};
    vector<complex<double>> roots = PolyRoots(realCoefficients);
    cout << "Polinomski korijeni (sa realnim koeficijentima): ";
    for (const auto &root : roots) {
      cout << root << " ";
    }
    cout << endl;
  } catch (const exception &e) {
    cerr << "Izuzetak pri pozivu PolyRoots sa realnim koeficijentima: "
         << e.what() << endl;
  }
}
int main() {
  // Pokrećemo test funkcije
  testirajRegulaFalsiPodrazumijevano();
  testirajRegulaFalsiIllinois();
  testirajRegulaFalsiSlavic();
  testirajRegulaFalsiIllinoisSlavic();
  testirajBracketRoot();
  testirajRiddersSolve();
  testirajNewtonRaphsonSolve();
  testirajPolyRootsComplex();
  testirajPolyRootsReal();
  return 0;
}
