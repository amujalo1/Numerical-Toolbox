#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
template<typename FunTip>
double RK4Step(FunTip f, double x, double y, double h) {
    double k1 = f(x,y);
    double k2 = f(x + h/2, y + h * k1 / 2);
    double k3 = f(x + h/2, y + h * k2 / 2);
    double k4 = f(x + h, y + h * k3);
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

template<typename FunTip>
vector<pair<double,double>> RK4Integrator(FunTip f, double x0, double y0, double xmax, double h, double eps = 1e-8, bool adaptive = false)
{
    if((h > 0 && xmax < x0) || (h <= 0 && xmax > x0)) return {{x0,y0}};
    vector<pair<double,double>>koef;
    if(!adaptive) {
        double x = x0, y = y0;
        if(h > 0) {
            while(x <= xmax + eps) {
                koef.push_back({x,y});
                y = RK4Step(f,x,y,h);
                x += h;
            }
        }
        else {
            while(x >= xmax - eps) {
                koef.push_back({x,y});
                y = RK4Step(f,x,y,h);
                x += h;
            }
        }
    }
    else {
        if(h > 0) {
            double x = x0, y = y0;
            koef.push_back({x,y});
            while(x <= xmax + eps) {
                double u = RK4Step(f,x,y,h/2);
                double v = RK4Step(f, x + h/2,u, h/2);
                double w = RK4Step(f,x,y,h);
                double delta = abs(w - v) / h;
                if(delta <= eps) {
                    x += h;
                    y = v;
                    koef.push_back({x,y});
                }
                h = h * min(5.0,0.9 * pow(eps/delta,1/4.));
            }
            if(abs(koef[koef.size() - 1].first - xmax) > eps) {
                koef.erase(koef.begin());
                h = xmax - koef[koef.size() - 1].first;
                double u = RK4Step(f,x,y,h/2);
                double v = RK4Step(f, x + h / 2, u , h/2);
                double w = RK4Step(f,x,y,h);
                koef[koef.size() - 1] = {xmax,v};
            }
        }
        else {
            double x = x0, y = y0;
            koef.push_back({x,y});
            while(x >= xmax - eps) {
                double u = RK4Step(f,x,y,h/2);
                double v = RK4Step(f, x + h / 2, u , h/2);
                double w = RK4Step(f,x,y,h);
                double delta = abs(w-v) / ((-1) * h);
                if(delta <= eps) {
                    x += h;
                    y = v;
                    koef.push_back({x,y});
                }
                h *= min(5.0,0.9 * pow(eps/delta,1/4.));
            }
            if(abs(koef[koef.size()-1].first - xmax) > eps) {
                koef.erase(koef.begin());
                h = xmax - koef[koef.size() - 1].first;
                double u = RK4Step(f,x,y,h/2);
                double v = RK4Step(f,x + h/2, u, h/2);
                double w = RK4Step(f,x,y,h);
                koef[koef.size() - 1] = {xmax,v};
            }
        }
    }
    return koef;
}

template <typename FunType>
vector<pair<double, vector<double>>> RK4SystemIntegrator(
    FunType f, double x0, vector<double> y0, double xmax, double h) {
    if (y0.size() != f(x0, y0).size()) {
        throw range_error("Incompatible formats");
    }

    vector<pair<double, vector<double>>> result;
    result.push_back(make_pair(x0, y0));

    double x = x0;
    vector<double> y = y0;

    while ((h > 0 && x < xmax) || (h < 0 && x > xmax)) {
        double K1_1 = h * f(x, y)[0];
        double K1_2 = h * f(x, y)[1];

        double K2_1 = h * f(x + h / 2, {y[0] + K1_1 / 2, y[1] + K1_2 / 2})[0];
        double K2_2 = h * f(x + h / 2, {y[0] + K1_1 / 2, y[1] + K1_2 / 2})[1];

        double K3_1 = h * f(x + h / 2, {y[0] + K2_1 / 2, y[1] + K2_2 / 2})[0];
        double K3_2 = h * f(x + h / 2, {y[0] + K2_1 / 2, y[1] + K2_2 / 2})[1];

        double K4_1 = h * f(x + h, {y[0] + K3_1, y[1] + K3_2})[0];
        double K4_2 = h * f(x + h, {y[0] + K3_1, y[1] + K3_2})[1];

        y[0] = y[0] + (K1_1 + 2 * K2_1 + 2 * K3_1 + K4_1) / 6;
        y[1] = y[1] + (K1_2 + 2 * K2_2 + 2 * K3_2 + K4_2) / 6;

        x = x + h;
        result.push_back(make_pair(x, y));
    }

    return result;
}

int main() {
    auto f = [](double x, vector<double> y) -> vector<double> {
        return {2 * y[0] + 3 * y[1] + x, y[0] - 2 * y[1] + 1};
    };
    double x0 = 0.0;
    vector<double> y0 = {1.0, 2.0};
    double xmax = 2.0;
    double h = 0.1;
    auto result = RK4SystemIntegrator(f, x0, y0, xmax, h);

    for (const auto& point : result) {
        cout << "x = " << point.first << ", y = (" << point.second[0] << ", " << point.second[1] << ")\n";
    }
    return 0;
}