#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>
using namespace std;

class ChebyshevApproximation {
    vector<double> c;
    double xmin, xmax;
    int n, m;

    ChebyshevApproximation(vector<double> c, double xmin, double xmax)
        : c(c), xmin(xmin), xmax(xmax), m(c.size() - 1) {}

public:
    template <typename FunType>
    ChebyshevApproximation(FunType f, double xmin, double xmax, int n)
        : xmax(xmax), xmin(xmin), m(n), n(n) {
        if (xmin >= xmax || n < 1)
            throw std::domain_error("Bad parameters");

        c.resize(n + 1);
        vector<double> w;
        for (int i = 0; i <= 4 * n + 3; i++) {
            double k = cos(M_PI * i / (2 * n + 2));
            w.push_back(k);
        }

        vector<double> v;
        for (int i = 0; i <= n; i++) {
            double k = (xmin + xmax + (xmax - xmin) * w[2 * i + 1]) / 2;
            v.push_back(f(k));
        }

        for (int k = 0; k <= n; k++) {
            double s = 0;
            for (int i = 0; i <= n; i++) {
                int j = (k * (2 * i + 1)) % (4 * n + 4);
                s = s + v[i] * w[j];
            }
            c[k] = 2 * s / (n + 1);
        }
    }

    void set_m(int k) {
        if (!(k > 1 && k < n))
            throw domain_error("Bad order");
        this->m = k;
    }

    void trunc(double eps) {
        if( eps < 0 ) throw domain_error("Bad tolerance");
        for(int i = m; i >= 0; i--) {
            if(i < 1 ) 
                throw domain_error("Bad tolerance");
            else if(fabs(c[i]) > eps) {
                m = i;
                break;
            }
        }
    }

    double operator()(double x) const {
        if (x > xmax || x < xmin)
            throw domain_error("Bad argument");
        double t = (2 * x - xmin - xmax) / (xmax - xmin);
        double p = 1, q = t, s = c[0] / 2 + c[1] * t;
        for (int k = 2; k <= m; k++) {
            double r = 2 * t * q - p;
            p = q;
            q = r;
            s += c[k] * r;
        }
        return s;
    }

    double derivative(double x) const {
        if (x > xmax || x < xmin)
            throw domain_error("Bad argument");
        double t = (2 * x - xmin - xmax) / (xmax - xmin);
        double p = 1, q = 4 * t, s = c[1] + q * c[2], r;
        for (int k = 3; k <= m; k++) {
            r = k * (2 * t * q / (k - 1) - p / (k - 2));
            p = q;
            q = r;
            s += c[k] * r;
        }
        return 2 * s / (xmax - xmin);
    }

    ChebyshevApproximation derivative() const {
        double u = 4. / (xmax - xmin);
        vector<double> K(c.size());
        K[m - 1] = u * m * c[m];
        K[m - 2] = u * (m - 1) * c[m - 1];
        for (int k = m - 3; k >= 0; k--)
            K[k] = K[k + 2] + u * (k + 1) * c[k + 1];
        return ChebyshevApproximation(K, xmin, xmax);
    }

    ChebyshevApproximation antiderivative() const {
        vector<double> K(m + 2);
        K[0] = 0.0;
        for (int k = 1; k <= m + 1; k++)
            K[k] = (xmax - xmin) / (4 * k) * (c[k - 1] - (k < m ? c[k + 1] : 0));
        return ChebyshevApproximation(K, xmin, xmax);
    }

    double integrate(double a, double b) const {
        if (a < xmin || a > xmax || b < xmin || b > xmax)
            throw domain_error("Bad interval");
        ChebyshevApproximation F(antiderivative());
        return F(b) - F(a);
    }

    double integrate() const {
        double s = 0.0;
        for (int k = 1; k <= (m + 1) / 2; k++)
            s += 2 * c[2 * k] / (1 - 4 * k * k);
        s *= (xmax - xmin) / 2;
        s += c[0] * (xmax - xmin) / 2;
        return s;
    }
};

int main() {
    try {
        // Test 1: Kreiranje ChebyshevApproximation koristeći funkciju sin(x)
        auto f1 = [](double x) { return sin(x); };
        ChebyshevApproximation approx1(f1, 0, M_PI, 5);

        // Test 2: Evaluacija aproksimacije u određenoj tački
        double rezultat1 = approx1(M_PI / 2);
        cout << "Test 2 Rezultat: " << rezultat1 << endl;

        // Test 3: Evaluacija derivacije u određenoj tački
        double rezultatDerivacije1 = approx1.derivative(M_PI / 2);
        cout << "Test 3 Rezultat: " << rezultatDerivacije1 << endl;

        // Test 4: Kreiranje ChebyshevApproximation sa lošim parametrima
        try {
            ChebyshevApproximation approx2(f1, M_PI, 0, 5); // Trebalo bi izazvati izuzetak
        } catch (const exception& e) {
            cout << "Test 4 Izuzetak: " << e.what() << endl;
        }

        // Test 5: Postavljanje lošeg reda
        try {
            approx1.set_m(10); // Trebalo bi izazvati izuzetak
        } catch (const exception& e) {
            cout << "Test 5 Izuzetak: " << e.what() << endl;
        }

        // Test 6: Skraćivanje s lošom tolerancijom
        try {
            approx1.trunc(-0.01); // Trebalo bi izazvati izuzetak
        } catch (const exception& e) {
            cout << "Test 6 Izuzetak: " << e.what() << endl;
        }

        // Test 7: Integracija preko nevažećeg intervala
        try {
            double integralniRezultat = approx1.integrate(2, 3); // Trebalo bi izazvati izuzetak
        } catch (const exception& e) {
            cout << "Test 7 Izuzetak: " << e.what() << endl;
        }

        // Test 8: Integracija preko važećeg intervala
        double integralniRezultat2 = approx1.integrate(0, M_PI / 2);
        cout << "Test 8 Rezultat: " << integralniRezultat2 << endl;

        // Test 9: Integracija preko cijelog opsega
        double integralniRezultat3 = approx1.integrate();
        cout << "Test 9 Rezultat: " << integralniRezultat3 << endl;

    } catch (const exception& e) {
        cout << "Neočekivan Izuzetak: " << e.what() << endl;
    }

    return 0;
}