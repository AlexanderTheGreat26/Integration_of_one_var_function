#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <utility>
#include <algorithm>


std::random_device rd;  // Will be used to obtain a seed for the random number engine.
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd().


const double pi = std::acos(-1);
const double left_border = 0;
const double right_border = 2*pi;
const int N = 1e6;
const double acc = (right_border - left_border) / N;


double function (double & x) {return std::pow(x,2) * std::cos(x);}

std::vector <double> mesh (double left_border, const double & right_border, const double & step);

std::vector<double> function_nodes (std::vector<double> & xx, double f(double & x));

std::vector<int> zeros_indexes (std::vector<double> & f);

std::vector<double> heights (std::vector<double> & function_nodes, std::vector<int> & zero_indexes);

std::vector<int> points_distribution (std::vector<double> & x_mesh, std::vector<double> & height,
                                      std::vector<int> & zero_indexes, const int & number_of_points);

double MC_integration (double f(double & x), double & x_min, double & x_max, double & height, const double & points_number);


int main () {
    std::vector<double> x_mesh = std::move(mesh(left_border, right_border, acc));
    std::vector<double> f_nodes = std::move(function_nodes(x_mesh, function));
    std::vector<int> zero_ind = std::move(zeros_indexes(f_nodes));
    std::vector<double> height = std::move(heights(f_nodes, zero_ind));

    std::vector<int> points = std::move(points_distribution(x_mesh, height, zero_ind, N));

    double sum = 0;
    for (int i = 0; i < zero_ind.size()-1; ++i)
            sum += MC_integration(function, x_mesh[zero_ind[i]], x_mesh[zero_ind[i + 1]], height[i], points[i]);

    std::cout << "Integrand:\t x^2*cos(x)\n"
                 "Left border:\t" << left_border <<
                 "\nRight border:\t" << right_border <<
                 "\nNumber of points:\t" << N <<
                 "\nResult:\t" << sum << '\n';

    return 0;
}


// Creates mesh from left border to right border with given step.
std::vector <double> mesh (double left_border, const double & right_border, const double & step) {
    std::vector <double> xx ((right_border-left_border) / step + 1);
    xx[0] = left_border;
    std::generate(xx.begin()+1, xx.end(), [&] {left_border += step; return left_border;});
    return xx;
}


// Returns std::vector of function values according to the given argument mesh (xx).
std::vector<double> function_nodes (std::vector<double> & xx, double f(double & x)) {
    std::vector<double> ff (xx.size());
    for (int i = 0; i < xx.size(); ++i)
        ff[i] = f(xx[i]);
    return ff;
}


std::vector<int> zeros_indexes (std::vector<double> & f) {
    std::vector<int> indexes;
    indexes.emplace_back(0);
    for (int i = 1; i < f.size()-1; ++i)
        if ((f[i] > 0 && f[i+1] < 0) || (f[i] < 0 && f[i+1] > 0))
            indexes.emplace_back(i+1);
    indexes.emplace_back(f.size()-1);
    return indexes;
}


std::vector<double> heights (std::vector<double> & function_nodes, std::vector<int> & zero_indexes) {
    std::vector<double> H(zero_indexes.size()-1);
    for (int i = 0; i < zero_indexes.size()-1; ++i) {
            const auto node_min = function_nodes.begin() + zero_indexes[i];
            const auto node_max = function_nodes.begin() + zero_indexes[i+1];
            double min = *std::min_element(node_min, node_max);
            double max = *std::max_element(node_min, node_max);
            H[i] = (std::fabs(max) > std::fabs(min)) ? max : min;
        }
    return H;
}


std::vector<int> points_distribution (std::vector<double> & x_mesh, std::vector<double> & height, std::vector<int> & zero_indexes,
                                      const int & number_of_points) {
    std::vector<int> points_count(zero_indexes.size()-1);
    std::vector<double> areas;
    double sum_area = 0;
    for (int i = 0; i < zero_indexes.size()-1; ++i) {
        areas.emplace_back((x_mesh[zero_indexes[i+1]] - x_mesh[zero_indexes[i]]) * std::fabs(height[i]));
        sum_area += areas[i];
    }
    for (int i = 0; i < points_count.size(); ++i)
        points_count[i] = number_of_points * areas[i]/sum_area;
    return points_count;
}


double MC_integration (double  f(double & x), double & x_min, double & x_max, double & height, const double & points_number) {
    int n = 0;
    std::uniform_real_distribution<double> dis_x (x_min, x_max);
    std::uniform_real_distribution<double> dis_y (0.0, height);
    for (int k = 0; k < points_number; ++k) {
        double x = dis_x(gen);
        double y = dis_y(gen);
         if ((f(x) > 0 && y < f(x) && y > 0) || (f(x) < 0 && y > f(x) && y < 0))
            ++n;
    }
    return (x_max - x_min) * height * n / points_number;
}


double rectangle_method (double f (double & x), double left_border, double right_border, const double & eps){
    int n = 10;
    double I2 = 1.0e20;
    double I1 = 0;
    while (std::fabs(I2 - I1) >= eps) {
        double dx = (left_border - right_border) / n;
        I1 = I2;
        I2 = (f(left_border) + f(right_border)) / 2.0;
        for (int i = 1; i < n; ++i) {
            double x = left_border + i*dx;
            I2 += f(x);
        }
        I2 *= dx;
        ++n;
    }
    return I2;
}


double trapezium_method (double f (double & x), const double & left_border, const double & right_border, const double & eps) {
    int n = 10;
    double I2 = 1.0e20;
    double I1 = 0;
    while (std::fabs(I2 - I1) >= eps) {
        double dx = (left_border - right_border) / n;
        I1 = I2;
        I2 = 0;
        for (int i = 0; i < n-1; ++i) {
            double x = left_border + i*dx;
            I2 += f(x);
        }
        I2 *= dx;
        ++n;
    }
    return I2;
}