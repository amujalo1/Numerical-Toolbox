#include <iostream>
#include <cmath>
#include <math.h>
using namespace std;
template <typename FunType>
double FindMinimum(FunType f, double x0, double eps = 1e-8, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4){
    if(eps<=0 || hinit<=0 || hmax<=0 || lambda<=0) 
        throw domain_error("Invalid parameters");
    //"tehnika lovljenja"
    double a = x0 - hinit;
    double b = x0 - hinit;
    double c = x0;
    bool pronadjen = false;
     while(abs(hinit) < hmax) {
        if(f(c + hinit) < f(c)) {
            b = c + hinit;
            a = c - hinit;
        }
        else if(f(c - hinit) < f(c)) {
            b = c - hinit;
            a = b - hinit;
        }
        else {
            a = c - hinit;
            b = c + hinit;
            pronadjen = true;
            break;
        }
        c = b;
        hinit *= lambda;
    }
    if(!pronadjen) 
        throw logic_error("Minimum has not found");
    //metod zlatnog presjeka
    double G = (1+sqrt(5))/2, d;
    if(abs(c-a) < abs(b-c))
        d=b-(b-c)/G;
    else {
        d = c;
        c = a + (c - a)/G;
    }
    double u = f(c);
    double v = f(d);
    while(abs(b-a)>eps){
        if(u < v){
            b = d;
            d = c;
            c = a+ (c - a)/G;
            v = u;
            u = f(c);
        }
        else {
            a = c;
            c = d;
            d = b - (b - d)/G;
            u = v;
            v = f(d);
        }
    }
    return (a + b)/2;
}


// Primer funkcije koju želimo minimizirati
double exampleFunction(double x) {
    return x * x + 3 * sin(x);
}

int main() {
    try {
        // Pozivamo FindMinimum za traženje minimuma funkcije exampleFunction
        double minimum = FindMinimum(exampleFunction, 0.0);

        // Ispisujemo pronađeni minimum
        cout << "Minimum pronadjen u x = " << minimum << endl;
    } catch (const exception& e) {
        // Uhvatimo izuzetak ako nešto pođe po zlu
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
