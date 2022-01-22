#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

double ComputeNorm(const vector<double>& x) {
    double result = 0.0;
    
    for (auto a : x) {
        result += a*a;
    }
    return sqrt(result);
}
