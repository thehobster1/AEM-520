#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return pow(x, 2);
}

int main() {
    double x0 = 2;
    double h = 0.0001;
    double derivative = central_difference(f, x0, h);
    cout << "The derivative of f(x) = x^2 at x = 2 is: " << derivative << endl;
    return 0;
}

double central_difference(double (*f)(double), double x0, double h) {
    double derivative = (f(x0 + h) - f(x0 - h)) / (2 * h);
    return derivative;
}