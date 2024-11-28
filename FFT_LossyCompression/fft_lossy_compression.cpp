#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <complex>
using namespace std;

bool stepen2(int x) {
    double y = log(x) / log(2);
    return fabs(y - (int)y) < 1e-8;
}
void FFT(vector<double>& x, vector<complex<double>>& xk, int N, int s = 0, int d = 0, int t = 1) {
    if (N == 1) {
        xk[d] = x[s];
    } else {
        FFT(x, xk, N / 2, s, d, 2 * t);
        FFT(x, xk, N / 2, s + t, d + N / 2, 2 * t);
        complex<double> mi = 1;
        complex<double> w = {cos(2 * M_PI / N), -sin(2 * M_PI / N)};
        complex<double> u, v;
        for (int k = d; k <= d + N / 2 - 1; k++) {
            u = xk[k];
            v = mi * xk[k + N / 2];
            xk[k] = u + v;
            xk[k + N / 2] = u - v;
            mi = mi * w;
        }
    }
}

void InvFFT(vector<complex<double>>& xk, vector<complex<double>>& x, int N, int s = 0, int d = 0, int t = 1) {
    if (N == 1) {
        x[d] = xk[s];
    } else {
        InvFFT(xk, x, N / 2, s, d, 2 * t);
        InvFFT(xk, x, N / 2, s + t, d + N / 2, 2 * t);
        complex<double> mi = 1.;
        complex<double> w = polar(1., (2 * M_PI) / N);
        complex<double> u, v;
        for (int k = d; k <= d + N / 2 - 1; k++) {
            u = x[k];
            v = mi * x[k + N / 2];
            x[k] = (u + v) / 2.;
            x[k + N / 2] = (u - v) / 2.;
            mi = mi * w;
        }
    }
}
vector<double> LossyCompress(vector<double> data, int new_size) {
    if (!(new_size > 1 && new_size <= data.size())) 
        throw range_error("Bad new size");
    int N = data.size();
    if (!stepen2(N)) 
        throw range_error("Data size must be a power of two");
    vector<double> y(N);
    for (int i = 0; i < N; i++) 
        y[i] = (i < N / 2) ? data[2 * i] : data[2 * (N - i) - 1];
    vector<complex<double>> yk(N);
    FFT(y, yk, N);
    vector<double> xk(new_size);
    for (int k = 0; k < new_size - 1; k++) 
        xk[k] = real(polar(1., -(M_PI * k) / (2 * N)) * yk[k]);
    xk[new_size - 1] = N;
    return xk;
}

vector<double> LossyDecompress(vector<double> compressed) {
    int N = compressed[compressed.size() - 1];
    compressed[compressed.size() - 1] = 0;
    if (!stepen2(N) || N < compressed.size())
        throw logic_error("Bad compressed sequence");
    compressed.resize(N);
    vector<complex<double>> yk(N);
    yk[0] = compressed[0];
    for (int k = 1; k < N; k++)
        yk[k] = 2. * polar(1., (M_PI * k) / (2 * N)) * compressed[k];
    vector<complex<double>> y(N);
    InvFFT(yk, y, N);
    vector<double> x(N);
    for (int n = 0; n < N; n++) 
        x[n] = (n % 2 == 0) ? real(y[n / 2]) : real(y[N - (n + 1) / 2]);
    return x;
}

int main() {
     try {
        LossyCompress({1, 6, 3, 5}, -5);
    }
    catch(range_error& e) {
        cout << "Range error in LossyCompress: '" << e.what() << "'" << endl;
    }

    try {
        LossyCompress({1, 2, 3, 4}, 6);
    }
    catch(range_error& e) {
        cout << "Range error in LossyCompress: '" << e.what() << "'" << endl;
    }

    try {
        LossyCompress({1, 2, 3, 4, 5}, 2);
    }
    catch(range_error& e) {
        cout << "Range error in LossyCompress: '" << e.what() << "'" << endl;
    }

    try {
        LossyDecompress({0, 0, 0, 0, 3});
    }
    catch(logic_error& e) {
        cout << "Logic error in LossyDecompress: '" << e.what() << "'" << endl;
    }
    catch(range_error& e) {
        cout << "Range error in LossyDecompress: '" << e.what() << "'" << endl;
    }

    try {
        LossyDecompress({0, 0, 0, 0, 0, 0, 4});
    }
    catch(logic_error& e) {
        cout << "Logic error in LossyDecompress: '" << e.what() << "'" << endl;
    }
    catch(range_error& e) {
        cout << "Range error in LossyDecompress: '" << e.what() << "'" << endl;
    }

    try {
        LossyDecompress({0, 0, 0, 0, 1.5});
    }
    catch(logic_error& e) {
        cout << "Logic error in LossyDecompress: '" << e.what() << "'" << endl;
    }
    catch(range_error& e) {
        cout << "Range error in LossyDecompress: '" << e.what() << "'" << endl;
    }
    return 0;
}
