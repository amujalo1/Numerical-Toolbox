#include <algorithm>
#include <cmath>
#include <complex>
#include <exception>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <math.h>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

class AbstractInterpolator {
  mutable int kresirani_index = -1;

protected:
  std::vector<std::pair<double, double>> dots;

  int Locate(double x) const {
    if (x <= dots.at(0).first) {
      kresirani_index = 0;
      return 0;
    } else if (x > dots.at(dots.size() - 1).first) {
      kresirani_index = dots.size() - 2;
      return dots.size();
    }
    if (kresirani_index >= 0 && kresirani_index < dots.size() - 1) {
      if (x >= dots[kresirani_index].first &&
          x < dots[kresirani_index + 1].first) {
        return kresirani_index + 1;
      } else if (x < dots[kresirani_index].first) {
        --kresirani_index;

      } else {
        ++kresirani_index;
      }
    }
    std::pair<double, double> Z = std::make_pair(x, 0);

    auto it = lower_bound(
        dots.begin(), dots.end(), Z,
        [](const pair<double, double> &a, const pair<double, double> &b) {
          return a.first < b.first;
        });
    kresirani_index = distance(dots.begin(), it) - 1;
    return kresirani_index + 1;
  }

public:
  AbstractInterpolator(const std::vector<std::pair<double, double>> &data)
      : dots(data) {
    std::sort(dots.begin(), dots.end(),
              [](const pair<double, double> &a, const pair<double, double> &b) {
                if (a.first < b.first)
                  return true;
                else if (a.first == b.first)
                  throw domain_error("Invalid data set");
                return false;
              });
  }
  virtual double operator()(double x) const = 0;
};

class LinearInterpolator : public AbstractInterpolator {
  double lin(double x, double x1, double y1, double x2, double y2) const {
    return (((x2 - x) / (x2 - x1)) * y1 + ((x - x1) / (x2 - x1)) * y2);
  }

public:
  using AbstractInterpolator::AbstractInterpolator;
  double operator()(double x) const override {
    int index = Locate(x);
    if (index >= 1 && index < dots.size()) {
      return lin(x, dots.at(index - 1).first, dots.at(index - 1).second,
                 dots.at(index).first, dots.at(index).second);
    } else if (index == 0) {
      return lin(x, dots.at(0).first, dots.at(0).second, dots.at(1).first,
                 dots.at(1).second);
    } else {
      return lin(
          x, dots.at(dots.size() - 2).first, dots.at(dots.size() - 2).second,
          dots.at(dots.size() - 1).first, dots.at(dots.size() - 1).second);
    }
  }
};
class PolynomialInterpolator : public AbstractInterpolator {
  std::vector<double> q;

public:
  PolynomialInterpolator(const std::vector<std::pair<double, double>> &data)
      : AbstractInterpolator(data) {
    for (int i = 0; i < dots.size(); i++)
      q.push_back(dots.at(i).second);

    int n = dots.size();
    for (int j = 0; j <= n - 2; j++)
      for (int i = n - 1; i >= j + 1; i--)
        q.at(i) = (q[i] - q[i - 1]) / (dots[i].first - dots[i - j - 1].first);
  }

  double operator()(double x) const override {
    double f = q[dots.size() - 1];
    for (int i = q.size() - 2; i >= 0; i--)
      f = f * (x - dots[i].first) + q[i];
    return f;
  }

  void AddPoint(const std::pair<double, double> &p) {
    for (int i = 0; i < dots.size(); i++)
      if (p.first == dots.at(i).first)
        throw std::domain_error("Invalid point");
    dots.push_back(p);

    std::vector<double> temp;
    for (int i = 0; i < dots.size(); i++)
      temp.push_back(dots.at(i).second);

    int n = dots.size();
    for (int j = 0; j <= n - 2; j++)
      for (int i = n - 1; i >= j + 1; i--)
        temp.at(i) =
            (temp[i] - temp[i - 1]) / (dots[i].first - dots[i - j - 1].first);

    q.push_back(temp[n - 1]);
  }

  std::vector<double> GetCoefficients() const {
    std::vector<double> p(dots.size(), 0);
    std::vector<double> W(dots.size() + 1, 0);

    W.at(0) = 1;

    for (int i = 1; i <= dots.size(); i++) {
      W.at(i) = W.at(i - 1);
      for (int j = i - 1; j >= 1; j--)
        W[j] = W[j - 1] - dots[i - 1].first * W[j];
      W[0] = -dots[i - 1].first * W[0];
    }

    for (int i = 1; i <= dots.size(); i++) {
      double a = 1;
      for (int j = 1; j <= dots.size(); j++)
        if (i != j)
          a = a * (dots[i - 1].first - dots[j - 1].first);

      a = dots[i - 1].second / a;

      std::vector<double> v(W);
      for (int j = 0; j < dots.size(); j++)
        v[j] = W[j];

      for (int j = dots.size() - 1; j >= 0; j--) {
        v[j] += dots[i - 1].first * v[j + 1];
        p[j] += a * v[j + 1];
      }
    }
    return p;
  }
};

class PiecewisePolynomialInterpolator : public AbstractInterpolator {
  int stepen;
  double Lagrange(double x, int poc, int kraj) const {
    double s = 0;
    for (int i = poc - 1; i < kraj; i++) {
      double p = dots.at(i).second;
      for (int j = poc - 1; j < kraj; j++)
        if (j != i)
          p *= (x - dots.at(j).first) / (dots.at(i).first - dots.at(j).first);
      s += p;
    }
    return s;
  }

public:
  PiecewisePolynomialInterpolator(
      const std::vector<std::pair<double, double>> &data, int order)
      : AbstractInterpolator(data) {
    if (order >= dots.size() || order < 1)
      throw domain_error("Invalid order");
    stepen = order;
  }
  double operator()(double x) const override {
    int index = Locate(x);
    int pocetak, kraj;
    if (stepen % 2 != 0) {
      pocetak = index - (stepen - 1) / 2;
      kraj = index + (stepen + 1) / 2;

      if (pocetak < 1) {
        return Lagrange(x, 1, stepen + 1);
      } else if (kraj > dots.size()) {

        return Lagrange(x, dots.size() - stepen, dots.size());
      }
      return Lagrange(x, pocetak, kraj);
    }

    pocetak = index - stepen / 2;
    kraj = index + stepen / 2;
    if (kraj > dots.size()) {
      return Lagrange(x, dots.size() - stepen, dots.size());
    } else if (pocetak < 1) {
      return Lagrange(x, 1, stepen + 1);
    }
    return Lagrange(x, pocetak, kraj);
  }
};

class SplineInterpolator : public AbstractInterpolator {
  vector<double> q;
  vector<double> r;
  vector<double> s;

  double fun(double x, int i) const {
    double t = x - dots.at(i).first;
    return dots.at(i).second + t * (q.at(i) + t * (r.at(i) + t * s.at(i)));
  }

public:
  SplineInterpolator(const std::vector<std::pair<double, double>> &data)
      : AbstractInterpolator(data) {
    int n = dots.size();
    r.resize(n);
    r[0] = 0;
    r[n - 1] = 0;
    for (int i = 1; i < n - 1; i++) {
      double add = 2 * (dots[i + 1].first - dots[i - 1].first);
      s.push_back(add);
      r[i] = 3 * ((dots[i + 1].second - dots[i].second) /
                      (dots[i + 1].first - dots[i].first) -
                  (dots[i].second - dots[i - 1].second) /
                      (dots[i].first - dots[i - 1].first));
    }

    for (int i = 1; i < n - 2; i++) {
      double mi = (dots[i + 1].first - dots[i].first) / s[i - 1];
      s[i] -= mi * (dots[i + 1].first - dots[i].first);
      r[i + 1] -= mi * r[i];
    }

    r[n - 2] /= s[n - 3];
    for (int i = n - 3; i > 0; i--)
      r[i] = (r[i] - (dots[i + 1].first - dots[i].first) * r[i + 1]) / s[i - 1];

    q.resize(n);
    s.resize(n);

    for (int i = 0; i < n - 1; i++) {
      double dX = dots[i + 1].first - dots[i].first;
      s[i] = (r[i + 1] - r[i]) / (3 * dX);
      q[i] = (dots[i + 1].second - dots[i].second) / dX -
             dX * (r[i + 1] + 2 * r[i]) / 3;
    }
  }
  double operator()(double x) const override {
    int index = Locate(x);
    if (index == 0) {
      return fun(x, 0);
    } else if (index == dots.size()) {
      return fun(x, dots.size() - 2);
    }
    return fun(x, index - 1);
  }
};
class BarycentricInterpolator : public AbstractInterpolator {
  int stepen;
  vector<double> tezinskiKoef;

public:
  BarycentricInterpolator(const std::vector<std::pair<double, double>> &data,
                          int order)
      : AbstractInterpolator(data) {
    if (!(order >= 0 && order <= dots.size()))
      throw std::domain_error("Invalid order");
    this->stepen = order;
    double p = 1;
    tezinskiKoef.resize(dots.size());
    for (int i = 1; i <= dots.size(); i++) {
      tezinskiKoef[i - 1] = 0;
      for (int k = (1 < (i - order) ? (i - order) : 1);
           k <= (i < (dots.size() - order) ? i : (dots.size() - order)); k++) {
        p = 1;
        for (int j = k; j <= k + order; j++)
          if (j != i)
            p /= (dots.at(i - 1).first - dots.at(j - 1).first);
        if (k % 2 == 0)
          p *= (-1);
      }
      tezinskiKoef[i - 1] += p;
    }
  }
  double operator()(double x) const override {
    double p = 0;
    double q = 0;
    for (int i = 0; i < dots.size(); i++) {
      if (fabs(x - dots.at(i).first) < numeric_limits<double>::epsilon())
        return dots.at(i).second;
      double u = tezinskiKoef[i] / (x - dots.at(i).first);
      p += u * dots.at(i).second;
      q += u;
    }
    return p / q;
  }
  std::vector<double> GetWeights() const { return tezinskiKoef; }
};

class TrigonometricInterpolator : public AbstractInterpolator {
public:
  TrigonometricInterpolator(const std::vector<std::pair<double, double>> &data)
      : AbstractInterpolator(data) {
    if (dots.front().second != dots.back().second) {
      throw std::domain_error("Function is not periodic");
    }
  }

  double operator()(double x) const override {
    int n = dots.size();
    /*double y = dots[0].second;
    for(int i= 1; i<n; i++)
        y+=dots[i].second*(sin(i*x)+cos(i*x));
    return y;*/
    double sum = 0.0;

    for (int k = 0; k < n; ++k) {
      double term = dots[k].second;
      for (int j = 1; j < n; ++j) {
        if (j != k) {  
          double w = 4*atan(1.) * j / (dots[n - 1].first - dots[0].first);
          double factor = w * (x - dots[k].first);
          term *= ((j % 2 == 1) ? std::sin(factor) : std::cos(factor))
                  * std::sin(w * (x - dots[j].first)) / std::sin(w * (dots[k].first - dots[j].first));
        }
      }
      sum += term;
    }
    return sum;
  }
};
void testLinearInterpolator() {
  cout << "<<Testiranje linearnog interpolatora>>" << endl;
  cout << "f(x) = e^((3 sin(x))/(log(2,x+4))) * x^(3/4)" << endl;
  cout << endl;
  auto f = [](double x) {
    return exp((3 * sin(x)) / log2(x + 4)) * pow(x, 3. / 4);
  };
  try {
    vector<pair<double, double>> data;
    for (double i = 15; i > 0; i -= 0.15) {
      data.push_back(make_pair<double, double>(double(i), f(i)));
    }
    LinearInterpolator lin(data);
    for (double x = 0; x < 20; x += 2.23435) {
      cout << "f(" << x << ") = " << f(x) << ", LinearInterpolator(" << x
           << ") = " << lin(x) << ", Greska = " << fabs(f(x) - lin(x)) << endl;
    }

  } catch (const domain_error &err) {
    cerr << err.what() << endl;
    ;
  }
}
void testPolynomialInterpolator() {
  cout << "<<Testiranje Polinomialnog interpolatora>>" << endl;
  // Lagrangeov polinom
  cout << "f(X) = ((x-1)(x-4.28)(x-10.28))/((7.58-1)(7.58-4.28)(10.28-7.58))"
       << endl;
  cout << endl;
  auto f = [](double x) {
    return ((x - 1) * (x - 4.28) * (x - 10.28)) /
           ((7.58 - 1) * (7.58 - 4.28) * (10.28 - 7.58));
  };
  try {
    vector<pair<double, double>> data;
    for (double i = 15; i >= -15; i -= 1.5) {
      data.push_back(make_pair<double, double>(double(i), f(i)));
    }
    PolynomialInterpolator poly(data);
    for (double x = -20; x < 20; x += 4.23435) {
      cout << "f(" << x << ") = " << f(x) << ", PolynomialInterpolator(" << x
           << ") = " << poly(x) << ", Greska = " << fabs(f(x) - poly(x))
           << endl;
    }
  } catch (const domain_error &err) {
    cerr << err.what() << endl;
    ;
  }
}
void testSplineInterpolator() {
  cout << "<<Testiranje Spline interpolatora>>" << endl;
  cout << "f(x) = sin(x)" << endl;
  cout << endl;
  auto f = [](double x) { return sin(x); };
  try {
    vector<pair<double, double>> data;
    for (double i = 15; i >= -15; i -= 1.5) {
      data.push_back(make_pair<double, double>(double(i), f(i)));
    }
    SplineInterpolator spline(data);
    for (double x = -4 * atan(1.); x <= 4 * atan(1.); x += 4 * atan(1.) / 5) {
      cout << "sin(" << x << ") = " << f(x) << ", SplineInterpolator(" << x
           << ") = " << spline(x) << ", Greska = " << fabs(f(x) - spline(x))
           << endl;
    }
  } catch (const domain_error &err) {
    cerr << err.what() << endl;
    ;
  }
}
void testPiecewisePolynomialInterpolator() {
  cout << "<<Testiranje Piecewise Polynomialnog interpolatora>>" << endl;
  cout << "f(x) = sin(x^(5/3))cos^2(x^(3/2)) / tan(e^(x^3))" << endl;
  cout << endl;
  auto f = [](double x) {
    return atan(exp(pow(x, 3))) * sin(pow(x, 5.0 / 3.0)) *
           pow(cos(pow(x, 3.0 / 2.0)), 2);
  };
  try {
    vector<pair<double, double>> data;
    for (double i = 0; i < 10; i += 0.25) {
      data.push_back(make_pair<double, double>(double(i), f(i)));
    }
    PiecewisePolynomialInterpolator Piecewise(data, 1);
    for (double x = 0; x < 15; x += 15. / 10) {
      cout << "f(" << x << ") = " << f(x)
           << ", PiecewisePolynomialInterpolator(" << x
           << ") = " << Piecewise(x)
           << ", Greska = " << fabs(f(x) - Piecewise(x)) << endl;
    }

  } catch (const domain_error &err) {
    cerr << err.what() << endl;
    ;
  }
}
void testBarycentricInterpolator() {
  cout << "<<Testiranje Baricentricnog interpolatora>>" << endl;
  cout << "f(x) = Riemann_zeta_imag()" << endl;
  cout << endl;
  auto zeta_imaginary_part = [](double t) {
    double imag_part = 0.0;
    for (int n = 1; n <= 100; ++n) {
      imag_part += 1.0 / pow(n, 0.5 + t);
    }
    return imag_part;
  };
  try {
    vector<pair<double, double>> data;
    for (double i = -10; i <= 10; i += 0.5) {
      data.push_back(
          make_pair<double, double>(double(i), zeta_imaginary_part(i)));
    }
    BarycentricInterpolator Barycentric(data, 1);
    for (double x = -15; x <= 15; x += 2.5) {
      cout << "f(" << x << ") = " << zeta_imaginary_part(x)
           << ", BarycentricInterpolator(" << x << ") = " << Barycentric(x)
           << ", Greska = " << fabs(zeta_imaginary_part(x) - Barycentric(x))
           << endl;
    }
  } catch (const domain_error &err) {
    cerr << err.what() << endl;
    ;
  }
}
void testTrigonometricInterpolator() {
  cout << "<<Testiranje Trigonometrijskog interpolatora>>" << endl;
  cout << "sin(x) = x - (x^3)/6 + (x^5)/120 - (x^7)/5040" << endl;
  cout << endl;
  auto f = [](double x) {
    return x - (pow(x, 3)) / 6 + pow(x, 5) / 120 - pow(x, 7) / 5040;
  };
  auto g = [](double x) { return sin(x); };
  try {
    vector<pair<double, double>> data;
    for (double i = 0; i < 4 * atan(1.); i += 0.01) {
      data.push_back(make_pair<double, double>(double(i), f(i)));
    } data.push_back(make_pair<double, double>(4 * atan(1.), 0.));
    TrigonometricInterpolator trig(data);
    for (double x = -4 * atan(1.); x <= 4 * atan(1.); x += 4 * atan(1.) / 5) {
      cout << "sin(" << x << ") = " << f(x) << ", TrigonometricInterpolator("
           << x << ") = " << trig(x) << ", Greska = " << fabs(f(x) - trig(x))
           << endl;
    }
  } catch (const domain_error &err) {
    cerr << err.what() << endl;
    ;
  }
}

int main() {
  testLinearInterpolator();
  cout << endl;
  testPolynomialInterpolator();
  cout << endl;
  testSplineInterpolator();
  cout << endl;
  testPiecewisePolynomialInterpolator();
  cout << endl;
  testBarycentricInterpolator();
  cout << endl;
  testTrigonometricInterpolator();
  return 0;
}