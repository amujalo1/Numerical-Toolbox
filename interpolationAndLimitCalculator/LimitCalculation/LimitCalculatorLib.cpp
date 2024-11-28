#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <functional>
#include <iomanip>
using namespace std;
template<typename FunType>
std::pair<double,bool>Limit(FunType f, double x0, double h = 0, double eps = 1e-8, double nmax = 20) {
    if(eps <= 0|| !(nmax > 3 && nmax <= 30)) throw std::domain_error("Invalid parameters");
    if(h < 0.0001 && h >= 0  && x0 >= 0) h = 0.001 * std::max(1.,std::fabs(x0)); 

    double yOld = std::numeric_limits<double>::infinity();
    std::vector<double> y(nmax);
    for(int i = 0; i < nmax; i++) {
        y.at(i) = f(x0 + h);
        double p = 2;
        for(int k = i - 1; k >= 0; k--) {
            y.at(k) = (p * y.at(k + 1) - y.at(k)) / (p - 1);
            p *= 2;
        }
        if(std::fabs(y.at(0) - yOld) < eps) return std::make_pair(y.at(0),true);
        yOld = y.at(0);
        h /= 2.;
    }
    return std::make_pair(y.at(0),false);
}

template <typename FunType>
void testLimitFunction(FunType f, double x0, double ocekivaniLimit) {
    auto rezultat = Limit(f, x0);

    std::cout << "Testiranje limita za funkciju kad x -> " << x0 << std::endl;
    std::cout << "Rezultat: " << std::setprecision(15) << rezultat.first << std::endl;
    std::cout << "Očekivani Limit: " << ocekivaniLimit << std::endl;
    std::cout << "Uspjeh: " << (rezultat.second && std::fabs(rezultat.first - ocekivaniLimit) < 1e-8 ? "Da" : "Ne") << std::endl;
    std::cout << std::endl;
}

int main() {
    // Testiranje Limit funkcije sa funkcijom sin(x)/x kako x prilazi 0
    auto sincPrekoX = [](double x) {
        return (std::fabs(x) < 1e-8) ? 1.0 : std::sin(x) / x;
    };
    // Testiranje za x koji prilazi 0, očekujemo da je limit 1
    testLimitFunction(sincPrekoX, 0.0, 1.0);
    return 0;
}