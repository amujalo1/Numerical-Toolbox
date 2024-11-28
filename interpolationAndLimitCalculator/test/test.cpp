#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

const double epsilon = 0.0001;

class AbstractInterpolator
{
    mutable int lociran_interval = 0;

protected:

    std::vector<std::pair<double,double>> tacke;
    int Locate(double x) const
    {
        if(x <= tacke.at(0).first) {
            lociran_interval = 0;
            return 0;
        } else if (x > tacke.at(tacke.size() - 1).first) {
            lociran_interval = tacke.size() -2;
            return tacke.size();

        }

        if(lociran_interval + 1 <= tacke.size() - 1 && x > tacke.at(lociran_interval).first && x <= tacke.at(lociran_interval + 1).first)
            return lociran_interval + 1;

        else if(lociran_interval > 0 && x <= tacke.at(lociran_interval).first && x > tacke.at(lociran_interval - 1 ).first) {
            lociran_interval--;
            return lociran_interval + 1;
        } else if(lociran_interval + 2 <= tacke.size() -1 && x > tacke.at(lociran_interval + 1).first && x <= tacke.at(lociran_interval + 2).first) {
            lociran_interval++;
            return lociran_interval + 1;
        }

        std::pair<double,double>Z = std::make_pair(x,std::numeric_limits<int>::min());

        lociran_interval = std::lower_bound(tacke.begin(), tacke.end(), Z,[](const std::pair<double,double>& X, const std::pair<double,double>&Y) {
            return X.first < Y.first;
        }) - tacke.begin();
        lociran_interval -= 1;
        return lociran_interval + 1;
    }
public:
    AbstractInterpolator(const std::vector<std::pair<double,double>> &data) : tacke(data)
    {
        std::sort(tacke.begin(), tacke.end(),[](std::pair<double,double> A, std::pair<double,double> B) {
            if(A.first < B.first) return true;
            else if(A.first == B.first) throw std::domain_error("Invalid data set");
            return false;
        });
    }
    virtual double operator()(double x) const = 0;
};

class LinearInterpolator : public AbstractInterpolator
{
    double proracun(double x, double x1, double y1, double x2, double y2 )const
    {
        return ( (x2 - x) / (x2 - x1) ) * y1 + ( (x - x1) / (x2 - x1) ) * y2;
    }
public:
    LinearInterpolator(const std::vector<std::pair<double,double>> &data) : AbstractInterpolator(data) {}
    double operator() (double x) const override
    {
        int intervalIndex = Locate(x);
        if(intervalIndex == 0) {
            return proracun(x,tacke.at(0).first,tacke.at(0).second,tacke.at(1).first,tacke.at(1).second);
        } else if(intervalIndex == tacke.size()) {
            return proracun(x, tacke.at(tacke.size() - 2).first,tacke.at(tacke.size() - 2).second,tacke.at(tacke.size() - 1).first,tacke.at(tacke.size() - 1).second);
        }
        return proracun(x,tacke.at(intervalIndex - 1).first,tacke.at(intervalIndex - 1).second, tacke.at(intervalIndex).first, tacke.at(intervalIndex).second);
    }
};

class PolynomialInterpolator : public AbstractInterpolator
{
    std::vector<double>y;                                                       // koeficijenti izracunati za Newtona

public:
    PolynomialInterpolator(const std::vector<std::pair<double,double>>& data) : AbstractInterpolator(data)
    {
        for(int i = 0; i < tacke.size(); i++)
            y.push_back(tacke.at(i).second);

        int n = tacke.size();
        for(int j = 0; j <= n - 2; j++)
            for(int i = n - 1; i >= j + 1; i--)
                y[i] = (double)(y[i] - y[i - 1]) / (double)(tacke.at(i).first - tacke.at(i - j - 1).first);
    }

    double operator()(double x) const override
    {
        double f = y[tacke.size() - 1];                                         //Horner algoritam
        for(int i = y.size() - 2; i >= 0; i--)
            f = f * (x - tacke.at(i).first) + y[i];
        return f;
    }

    void AddPoint(const std::pair<double,double>&p)
    {

        for(int i = 0; i < tacke.size(); i++)
            if(tacke.at(i).first == p.first) throw std::domain_error("Invalid point");
        tacke.push_back(p);

        std::vector<double>y2;
        for(int i = 0; i < tacke.size(); i++)
            y2.push_back(tacke.at(i).second);

        int n = tacke.size();
        for(int j = 0; j <= n - 2; j++)
            for(int i = n - 1; i >= j + 1; i--)
                y2[i] = (double)(y2[i] - y2[i - 1]) / (double)(tacke.at(i).first - tacke.at(i - j - 1).first);

        y.push_back(y2[n-1]);

    }
    std::vector<double> GetCoefficients() const
    {
        std::vector<double>p(tacke.size(), 0);
        std::vector<double>w(tacke.size() + 1,0);

        w.at(0) = 1;

        for(int i = 1; i <= tacke.size(); i++) {
            w.at(i) = w.at(i - 1);
            for(int j = i - 1; j >= 1; j--)
                w[j] = w[j-1] - tacke.at(i -1).first * w[j];
            w[0] = - tacke.at(i - 1).first * w[0];
        }

        for(int i = 1; i <= tacke.size(); i++) {
            double alfa = 1;
            for(int j = 1; j <= tacke.size(); j++)
                if(i != j)
                    alfa =alfa* (tacke.at(i - 1).first - tacke.at(j - 1).first);
            alfa = tacke.at(i - 1).second / alfa;

            std::vector<double> v(w);
            for(int j = 0; j < tacke.size(); j++)
                v[j] = w[j];

            for(int j = tacke.size() - 1; j >= 0; j--) {
                v[j] += tacke.at(i - 1).first * v[j + 1];
                p[j] += alfa * v[j + 1];
            }
        }
        return p;
    }
};

class PiecewisePolynomialInterpolator : public AbstractInterpolator
{
    int order;
    double Lagrange(double x, int s, int e) const
    {
        double res = 0;
        for(int i = s - 1; i < e; i++) {
            double p = tacke.at(i).second;
            for(int j = s - 1; j < e; j++)
                if(j != i)p *= (x - tacke.at(j).first) / (tacke.at(i).first - tacke.at(j).first);
            res += p;
        }
        return res;
    }

public:
    PiecewisePolynomialInterpolator(const std::vector<std::pair<double,double>>&data, int order) : AbstractInterpolator(data)
    {
        if(order < 1 || order >= tacke.size()) throw std::domain_error("Invalid order");
        this->order = order;
    }

    double operator()(double x) const override
    {
        int index = Locate(x);

        if(order % 2 != 0) {
            int start = index - (order - 1) / 2;
            int end = index + (order + 1) / 2;

            if(start < 1) {                                                     // tacka blizu lijevog ruba, uzima se prvih k + 1 tacaka
                return Lagrange(x,1, order + 1);
            } else if(end > tacke.size()) {                                     // tacka blizu desnog ruba, uzimamo zadnjih k + 1 tacaka
                return Lagrange(x, tacke.size() - order,tacke.size());
            }
            return Lagrange(x,start,end);
        }

        int start = index - order / 2;
        int end = index + order / 2;

        if(start < 1) {
            return Lagrange(x,1, order + 1);
        } else if(end > tacke.size()) {
            return Lagrange(x, tacke.size() - order,tacke.size());
        }
        return Lagrange(x,start,end);
    }
};
class SplineInterpolator : public AbstractInterpolator
{
    std::vector<double> r;                                                      //splajn koeficijenti
    std::vector<double> q;
    std::vector<double> a;

    double proracun(int i, double x) const
    {
        double t = x - tacke[i].first;
        return tacke[i].second + t*(q[i] + t*(r[i] + t*a[i]));
    }

public:
    SplineInterpolator(const std::vector<std::pair<double,double>> &data) : AbstractInterpolator(data)
    {
        int n = tacke.size();

        r.resize(n);
        r[0] = 0;
        r[n - 1] = 0;

        for(int i = 1; i < n - 1; i++) {
            double dodaj = 2 * (tacke[i + 1].first - tacke[i - 1].first);
            a.push_back(dodaj);
            r[i] = 3 * ((tacke[i+1].second - tacke[i].second) / (tacke[i+1].first - tacke[i].first) - (tacke[i].second - tacke[i-1].second) / (tacke[i].first - tacke[i-1].first));
        }

        for(int i = 1; i < n - 2; i++) {
            double mi = (tacke[i + 1].first - tacke[i].first) / a[i - 1];
            a[i] -= mi * (tacke[i + 1].first - tacke[i].first);
            r[i + 1] -= mi * r[i];
        }

        r[n - 2] /= a[n - 3];
        for(int i = n - 3; i >0; i--)
            r[i] = (r[i] - (tacke[i + 1].first - tacke[i].first) * r[i + 1]) / a[i-1];

        q.resize(n);
        a.resize(n);

        for(int i = 0; i < n - 1; i++) {
            double deltaX = tacke[i + 1].first - tacke[i].first;
            a[i] = (r[i + 1] - r[i]) / (3*deltaX);
            q[i] = (tacke[i+1].second - tacke[i].second)/deltaX - deltaX*(r[i + 1] + 2*r[i])/3;
        }

    }

    double operator() (double x) const override
    {
        int index = Locate(x);                                                  // vrati mi matematicki indeks

        if(index == 0) {
            return proracun(0, x);
        } else if(index == tacke.size()) {
            return proracun(tacke.size() - 2,x);
        }

        return proracun(index - 1, x);
    }
};

class BarycentricInterpolator  : public AbstractInterpolator

{
    int order;
    std::vector<double> tezinski;
    int max(int x, int y)
    {
        if(x < y) return y;
        return x;
    }

    int min(int x, int y)
    {
        if(x < y) return x;
        return y;
    }

public:
    BarycentricInterpolator (const std::vector<std::pair<double,double>>&data, int order) : AbstractInterpolator(data)
    {
        if(!(order >= 0 && order <= tacke.size() )) throw std::domain_error("Invalid order");
        this->order = order;
        double p = 1;
        tezinski.resize(tacke.size());
        for(int i = 1; i <= tacke.size(); i++) {
            tezinski[i - 1] = 0;
            for(int k = max(1, i - order); k <= min(i, tacke.size() - order); k++) {
                p = 1;
                for(int j = k ; j <= k + order ; j++)
                    if(j != i) p /= (tacke.at(i - 1).first - tacke.at(j - 1).first);
                if(k % 2 == 0) p *= (-1);
            }
            tezinski[i - 1] += p;
        }
    }

    double operator()(double x) const override
    {
        double p = 0;
        double q = 0;
        for(int i = 0; i < tacke.size(); i++) {
            if(x == tacke.at(i).first) return tacke.at(i).second;
            double u = tezinski[i] / (x - tacke.at(i).first);
            p += u * tacke.at(i).second;
            q += u;
        }
        return p/q;
    }

    std::vector<double>GetWeights() const
    {
        return tezinski;
    }
};
class TrigonometricInterpolator : public AbstractInterpolator {
public:
    TrigonometricInterpolator(const std::vector<std::pair<double, double>> &data)
        : AbstractInterpolator(data) {
        if (tacke.front().second != tacke.back().second) {
            throw std::domain_error("Function is not periodic");
        }
    }

    double operator()(double x) const override {
        int n = tacke.size();
        double sum = 0.0;

        for (int k = 0; k < n; ++k) {
            double term = tacke[k].second;

            for (int j = 1; j <= n / 2; ++j) {
                double factor = 2.0 * M_PI * j * (x - tacke[k].first) /
                                (tacke[n - 1].first - tacke[0].first);
                term += 2.0 * tacke[k].second * std::cos(factor);
            }

            sum += term;
        }

        return sum;
    }
};
int main()
{
    {
        //AddPoint test
        PolynomialInterpolator L({{1,2},{3,4},{6,6},{8,5}});
        L.AddPoint({5,7});
        if(L(4.) - 6.28571 < 0.01)
            std::cout <<"OK\n" ;
        else std::cout << "NOT OK";
    }

    {
        //operator() test
        std::cout << "Test operatora() :\n";
        {
            LinearInterpolator L({{1,2},{3,4},{6,6},{8,5}});
            if(L(4.) - 4.94286 < epsilon)
                std::cout <<"OK\n" ;
            else std::cout << "NOT OK";
        }

        {
            PolynomialInterpolator L({{1,2},{3,4},{6,6},{8,5}});
            if(L(4.) - 4.94286 < epsilon)
                std::cout <<"OK\n" ;
            else std::cout << "NOT OK";
        }

        {
            PiecewisePolynomialInterpolator L({{1,2},{3,4},{6,6},{8,5}},3);
            if(L(4.) - 4.94286 < epsilon)
                std::cout <<"OK\n" ;
            else std::cout << "NOT OK";
        }

        {
            SplineInterpolator L({{1,2},{3,4},{6,6},{8,5}});
            if(L(4.) - 4.97314 < epsilon)
                std::cout <<"OK\n" ;
            else std::cout << "NOT OK";
        }

        {
            BarycentricInterpolator L({{1,2},{3,4},{6,6},{8,5}},3);
            if(L(4.) - 4.94286 < epsilon)
                std::cout <<"OK\n" ;
            else std::cout << "NOT OK";
        }
    }


    {
        //testiranje vrijednosti
        //5 x^5 + 4 x^4 + 3 x^3 + 2 x^2 + x
        std::vector<std::pair<double, double>> data{{1,15},{2,258},{3,1641},{4,6372},{5,18555},{6,44790}};
        double x1 = 1.5, x2 = 2.5, x3 = 4.5;
        LinearInterpolator i1(data);
        PolynomialInterpolator i2(data);
        PiecewisePolynomialInterpolator i3(data,3);
        SplineInterpolator i4(data);
        BarycentricInterpolator i5(data,3);

        std::cout << "Actual:" << std::endl;
        std::cout << x1 << ": " << 74.3438 << std::endl;
        std::cout << x2 << ": " << 706.406 << std::endl;
        std::cout << x3 << ": " << 11185. << std::endl;

        std::cout<<"LinearInterpolator:"<<std::endl;
        std::cout << x1 << ": " << i1(x1) << std::endl;
        std::cout << x2 << ": " << i1(x2) << std::endl;
        std::cout << x3 << ": " << i1(x3) << std::endl;

        std::cout<<"PolynomialInterpolator:"<<std::endl;
        std::cout << x1 << ": " << i2(x1) << std::endl;
        std::cout << x2 << ": " << i2(x2) << std::endl;
        std::cout << x3 << ": " << i2(x3) << std::endl;

        std::cout<<"PiecewisePolynomialInterpolator:"<<std::endl;
        std::cout << x1 << ": " << i3(x1) << std::endl;
        std::cout << x2 << ": " << i3(x2) << std::endl;
        std::cout << x3 << ": " << i3(x3) << std::endl;

        std::cout<<"SplineInterpolator:"<<std::endl;
        std::cout << x1 << ": " << i4(x1) << std::endl;
        std::cout << x2 << ": " << i4(x2) << std::endl;
        std::cout << x3 << ": " << i4(x3) << std::endl;

        std::cout<<"BarycentricInterpolator:"<<std::endl;
        std::cout << x1 << ": " << i5(x1) << std::endl;
        std::cout << x2 << ": " << i5(x2) << std::endl;
        std::cout << x3 << ": " << i5(x3) << std::endl;
    }


    {
        //Test izuzetaka
        std::cout<<"Test izuzetaka:" << std::endl;
        try {
            LinearInterpolator({{1,2},{1,2}});
            std::cout<<"NOT OK" << std::endl;
        } catch(std::domain_error e) {
            std::cout<<"OK " <<e.what()<< std::endl;
        }

        try {
            PolynomialInterpolator p({{1,2},{2,2}});
            p.AddPoint({2,3});
            std::cout<<"NOT OK" << std::endl;
        } catch(std::domain_error e) {
            std::cout<<"OK " <<e.what()<< std::endl;
        }

        try {
            PiecewisePolynomialInterpolator({{1,2},{2,2}}, 3);
            std::cout<<"NOT OK" << std::endl;
        } catch(std::domain_error e) {
            std::cout<<"OK " <<e.what()<< std::endl;
        }

        try {
            PiecewisePolynomialInterpolator({{1,2},{2,2}}, 0);
            std::cout<<"NOT OK" << std::endl;
        } catch(std::domain_error e) {
            std::cout<<"OK " <<e.what()<< std::endl;
        }

        try {
            BarycentricInterpolator({{1,2},{2,2}}, 3);
            std::cout<<"NOT OK" << std::endl;
        } catch(std::domain_error e) {
            std::cout<<"OK " <<e.what()<< std::endl;
        }

        try {
            BarycentricInterpolator({{1,2},{2,2}}, -1);
            std::cout<<"NOT OK" << std::endl;
        } catch(std::domain_error e) {
            std::cout<<"OK " << e.what()<<std::endl;
        }
    }

    {
        std::vector<std::pair<double, double>> data{{1,15},{2,258},{3,1641},{4,6372},{5,18555}};
        double x1 = 1.5, x2 = 2.5, x3 = 4.5;
        PolynomialInterpolator p(data);
        BarycentricInterpolator b({{1,15},{2,258},{3,1641},{4,6372},{5,18555},{6,44790}}, 3);
        p.AddPoint({6,44790});

        std::cout<<"PolynomialInterpolator sa AddPoint:"<<std::endl;
        std::cout << x1 << ": " << p(x1) << std::endl;
        std::cout << x2 << ": " << p(x2) << std::endl;
        std::cout << x3 << ": " << p(x3) << std::endl;

        std::vector<double>expected{0,1,2,3,4,5};
        if(expected == p.GetCoefficients()) std::cout <<"OK\n";

        std::cout<<"PolynomialInterpolator GetCoefficients:"<<std::endl;
        for(auto x : p.GetCoefficients()) std::cout << x << " ";
        std::cout << std::endl;

        std::cout<<"BarycentricInterpolator GetWeights:"<<std::endl;
        for(auto x : b.GetWeights()) std::cout << x << " ";
        std::cout << std::endl;
    }
    return 0;
}