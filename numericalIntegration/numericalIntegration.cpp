#include <iostream>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;
auto operator+ (decltype(make_pair(0.0, false)) a, decltype(make_pair(0.0, false)) b) -> decltype(make_pair(0.0, false)) {
    return make_pair(a.first + b.first, a.second && b.second);
}
template <typename FunType>
pair<double, double> RombergIntegration(FunType f, double a, double b, double eps = 1e-10, int max = 1000000, int min = 50){
    if(eps < 0 || max < 0 || min < 0 || max < min)
        throw domain_error("Bad parameter");
    int N = 2;
    double h = (b - a) / N;
    double s = (f(a) + f(b)) / 2;
    double Iold = s;
    vector<double> I;
    for (int i = 1; N <= max; i++) {
        for (int j = 1; j <= N / 2; j++)
            s += f(a + (2 * j - 1) * h);

        I.push_back(h*s);
        double p = 4;
        for (int k = I.size() - 2; k >= 0; k--) {
            I[k] = (p * I[k+1] - I[k]) / (p - 1);
            p *= 4;
        }
        if (N >= min && fabs(I[0] - Iold) <= eps)
            return make_pair(I[0], true);

        Iold = I[0];
        h /= 2;
        N *= 2;
    }

    return make_pair(Iold, false);
}

template <typename FunType>
pair<double, bool> TanhSinhIntegration(FunType f, double a, double b, double eps = 1e-10, int nmax = 1000000, int nmin = 20, double range = 3.5){
    if (eps < 0 || nmin < 0 || range < 0 || nmax < nmin)
        throw domain_error("Bad parameter");

    int znak = 1;
    if (a > b) {
        znak *= -1;
        swap(a, b);
    }


    double Iold = 0;
    for (double N = 2, h = 2 * range / N, p = (b + a) / 2, q = (b - a) / 2, s =0; N < nmax; N *= 2, h /= 2) {
        for (int i = 1; i <= N / 2; i++) {
            double t = -range + (2 * i - 1) * h;
            double u = M_PI * sinh(t) / 2;
            double v = f(p + q * tanh(u));

            if (isfinite(v))
                s += q * M_PI * cosh(t) * v / (2 * cosh(u) * cosh(u));
        }

        double I = h * s;

        if (N >= nmin && fabs(I - Iold) <= eps)
            return make_pair(znak * I, true);

        Iold = I;
    }

    return make_pair(Iold * znak, false);

}


template <typename FunType>
auto AdaptiveStep(FunType f, double a, double b, double eps, double f1, double f2, double f3, double R) -> decltype(make_pair(0.0, false)) {
    if (!isfinite(f1)) f1 = 0;
    if (!isfinite(f2)) f2 = 0;
    if (!isfinite(f3)) f3 = 0;
    double c((a + b) / 2);
    double I1 = (b - a) * (f1 + 4 * f3 + f2) / 6;
    double f4 = f((a + c) / 2);
    double f5 = f((c + b) / 2);
    if (!isfinite(f4)) f4 = 0;
    if (!isfinite(f5)) f5 = 0;
    double I2 = (b - a) * (f1 + 4 * f4 + 2 * f3 + 4 * f5 + f2) / 12;
    if (fabs(I1 - I2) <= eps) return make_pair(I2, true);
    if (R <= 0) return make_pair(I2, false);
    return AdaptiveStep(f, a, c, eps, f1, f3, f4, R - 1) + AdaptiveStep(f, c, b, eps, f3, f2, f5, R - 1);
}



template <typename FunType>
auto AdaptiveIntegration(FunType f, double a, double b, double eps = 1e-10, int maxdepth = 30, int nmin = 1) -> decltype(make_pair(0.0, false)) {
    if (eps < 0 || maxdepth < 0 || nmin < 0) throw domain_error("Bad parameter");
    double h((b - a) / nmin);
    auto s = make_pair(0.0, true);
    for (int i = 1; i <= nmin; i++) {
        s = s + AdaptiveStep(f, a, a + h, eps, f(a), f(a + h), f((a + h) / 2), maxdepth);
        a += h;
    }
    return s;
}

int main() {
    try {
        // Testiranje RombergIntegration
        auto rombergResult = RombergIntegration([](double x) { return sin(x); }, 0, M_PI);
        cout << "Rezultat Romberg integracije: " << rombergResult.first << " | Uspješno: " << rombergResult.second << endl;
    } catch (const exception &e) {
        cerr << "Greška prilikom Romberg integracije: " << e.what() << endl;
    }

    try {
        // Testiranje TanhSinhIntegration
        auto tanhSinhResult = TanhSinhIntegration([](double x) { return exp(-x * x); }, -2, 2);
        cout << "Rezultat TanhSinh integracije: " << tanhSinhResult.first << " | Uspješno: " << tanhSinhResult.second << endl;
    } catch (const exception &e) {
        cerr << "Greška prilikom TanhSinh integracije: " << e.what() << endl;
    }

    try {
        // Testiranje AdaptiveIntegration
        auto adaptiveResult = AdaptiveIntegration([](double x) { return sqrt(x); }, 0, 4);
        cout << "Rezultat Adaptivne integracije: " << adaptiveResult.first << " | Uspješno: " << adaptiveResult.second << endl;
    } catch (const exception &e) {
        cerr << "Greška prilikom Adaptivne integracije: " << e.what() << endl;
    }

    return 0;
}





