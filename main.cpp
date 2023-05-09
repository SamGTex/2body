#include "vector_functions.h"
#include <algorithm>
#include <chrono>
#include <cstdlib> //for exit(1) to catch errors
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

// best decision ever
using namespace std;
using vec = vector<double>;

// calc acceleration vector a
vec get_a(const vec &r1, const vec &r2, const double M) {
    // a = F(r1, r2)/m
    vec r_dif = vec_sub(r1, r2);

    // check if r1-r2 == 0
    if (vec_abs(r_dif) == 0) {
        cerr << "Division by zero: r1 and r2 must be different!";
        exit(1);
    }

    double factor = -M / pow(vec_abs(r_dif), 3);
    return vec_scalar(r_dif, factor);
}

// calc velocity vector v
vec get_v(const vec &v_old, const vec &a, const double h) {
    // v + a*h
    return vec_add(v_old, vec_scalar(a, h));
}

// calc space vector r
vec get_r(const vec &r_old, const vec &v, const double h) {
    // r + v*h
    return vec_add(r_old, vec_scalar(v, h));
}

// Euler algorithm
// results (space vectors r1, r2) in csv file
// return last r1, r2, v1, v2 vector as tuple
tuple<vec, vec, vec, vec> euler(vec r1, vec r2, vec v1, vec v2, double m1, double m2, double h, double T_max) {
    // var m1
    vec r1_old{0.0, 0.0};
    vec v1_old{0.0, 0.0};
    vec a1{0.0, 0.0};
    vec a1_old{0.0, 0.0};

    // var m2
    vec r2_old{0.0, 0.0};
    vec v2_old{0.0, 0.0};
    vec a2{0.0, 0.0};
    vec a2_old{0.0, 0.0};

    // open csv files to save the space vectors r1, r2
    // use stepsize h in filename
    string h_string = to_string(h);
    replace(h_string.begin(), h_string.end(), '.', '_');

    string path_r1 = string("build/euler_r1_h") + h_string + string(".csv");
    string path_r2 = string("build/euler_r2_h") + h_string + string(".csv");

    ofstream file_r1(path_r1);
    ofstream file_r2(path_r2);

    // algorithm
    for (int i = 1; i * h < T_max; i++) { // t = i*h, i=1,2,..
        // m1
        a1_old = a1;
        a1 = get_a(r1, r2, m2); // update acceleration a

        v1_old = v1;
        v1 = get_v(v1_old, a1, h); // calc new velocity

        r1_old = r1;
        r1 = get_r(r1_old, v1, h); // calc new position

        file_r1 << i*h << ", " << r1[0] << ", " << r1[1] << endl; // write to csv

        // m2
        a2_old = a2;
        a2 = get_a(r2, r1_old, m1); // update acceleration a

        v2_old = v2;
        v2 = get_v(v2_old, a2, h); // calc new velocity

        r2_old = r2;
        r2 = get_r(r2_old, v2, h); // calc new position

        file_r2 << i*h << ", " << r2[0] << ", " << r2[1] << endl; // write to csv
    }

    // close csv files
    file_r1.close();
    file_r2.close();

    return make_tuple(r1, r2, v1, v2);
}

// Verlet algorithm
tuple<vec, vec, vec, vec> verlet(vec r1, vec r2, vec v1, vec v2, double m1, double m2, double h, double T_max) {
    // var m1
    vec r1_old{0.0, 0.0};  // r_n
    vec r1_oold{0.0, 0.0}; // r_n-1

    vec a1{0.0, 0.0};
    vec a1_old{0.0, 0.0};

    // var m2
    vec r2_old{0.0, 0.0};
    vec r2_oold{0.0, 0.0};

    vec a2{0.0, 0.0};
    vec a2_old{0.0, 0.0};

    // calc initial value for r_-1
    // r_-1 = r_0 - v_0*h + 0.5*a_0*h^2
    r1_old = vec_add(vec_sub(r1, vec_scalar(v1, h)), vec_scalar(get_a(r1, r2, m2), 0.5 * h * h));
    r2_old = vec_add(vec_sub(r2, vec_scalar(v2, h)), vec_scalar(get_a(r2, r1, m1), 0.5 * h * h));

    // open csv files to save the space vectors r1, r2
    // use stepsize h in filename
    string h_string = to_string(h);
    replace(h_string.begin(), h_string.end(), '.', '_');

    string path_r1 = string("build/verlet_r1_h") + h_string + string(".csv");
    string path_r2 = string("build/verlet_r2_h") + h_string + string(".csv");

    ofstream file_r1(path_r1);
    ofstream file_r2(path_r2);

    // algorithm
    for (int i = 1; i * h < T_max; i++) { // t = i*h, i=1,2,..
        // m1
        a1_old = a1;
        a1 = get_a(r1, r2, m2); // update a

        r1_oold = r1_old;
        r1_old = r1;

        // r_n+1 = 2*r_n - r_n-1 + a_n*h^2
        r1 = vec_add(vec_sub(vec_scalar(r1_old, 2), r1_oold), vec_scalar(a1, h * h));

        file_r1 << i*h << ", " << r1[0] << ", " << r1[1] << endl; // write to csv

        // m2
        a2_old = a2;
        a2 = get_a(r2, r1_old, m1); // update a

        r2_oold = r2_old;
        r2_old = r2;

        // r_n+1 = 2*r_n - r_n-1 + a_n*h^2
        r2 = vec_add(vec_sub(vec_scalar(r2_old, 2), r2_oold), vec_scalar(a2, h * h));

        file_r2 << i*h << ", " << r2[0] << ", " << r2[1] << endl; // write to csv
    }

    // close csv files
    file_r1.close();
    file_r2.close();

    return make_tuple(r1, r2, v1, v2);
}

int main() {
    // Var mass 1: initial conditions
    const double m1 = 1.0;
    const vec r1{0.0, 1.0};
    const vec v1{0.8, 0.0}; //0.8

    // Var mass 2: initial conditions
    const double m2 = 2.0; // =2
    const vec r2{0.0, -0.5};
    const vec v2{-0.4, 0.0};

    // Var time
    // t = i*h = 1*h, 2*h, ..., T_max
    double T_max = 100.0;

    // Euler algorithm >>>
    // a)
    euler(r1, r2, v1, v2, m1, m2, 1.0, T_max);
    euler(r1, r2, v1, v2, m1, m2, 0.1, T_max);

    // b) time measurement for h=0.01
    auto start_euler = chrono::high_resolution_clock::now();
    euler(r1, r2, v1, v2, m1, m2, 0.01, T_max);
    auto end_euler = chrono::high_resolution_clock::now();
    auto duration_euler = chrono::duration_cast<chrono::milliseconds>(end_euler - start_euler);
    cout << "Euler algorithm duration: " << duration_euler.count() << " ms" << endl;

    // c) h=0.05 and back with h=-0.05
    // save last r1, r2, v1, v2 vector
    tuple<vec, vec, vec, vec> result_euler = euler(r1, r2, v1, v2, m1, m2, 0.05, T_max);

    vec r1_new_euler = get<0>(result_euler);
    vec r2_new_euler = get<1>(result_euler);
    vec v1_new_euler = get<2>(result_euler);
    vec v2_new_euler = get<3>(result_euler);

    // use the latest positions/velocities and use now h -> -h
    // doesnt work
    // euler(r1_new_euler, r2_new_euler, v1_new_euler, v2_new_euler, m1, m2, -0.05, T_max);

    // <<<

    // Verlet algorithm >>>
    verlet(r1, r2, v1, v2, m1, m2, 1.0, T_max);
    verlet(r1, r2, v1, v2, m1, m2, 0.1, T_max);

    // time measurement for h=0.01
    auto start_verlet = chrono::high_resolution_clock::now();
    verlet(r1, r2, v1, v2, m1, m2, 0.01, T_max);
    auto end_verlet = chrono::high_resolution_clock::now();
    auto duration_verlet = chrono::duration_cast<chrono::milliseconds>(end_verlet - start_verlet);
    cout << "Verlet algorithm duration: " << duration_verlet.count() << " ms" << endl;

    // c) h=0.05 and back with h=-0.05
    // save last r1, r2, v1, v2 vector
    tuple<vec, vec, vec, vec> result_verlet = verlet(r1, r2, v1, v2, m1, m2, 0.05, T_max);

    vec r1_new_verlet = get<0>(result_verlet);
    vec r2_new_verlet = get<1>(result_verlet);
    vec v1_new_verlet = get<2>(result_verlet);
    vec v2_new_verlet = get<3>(result_verlet);

    // use the latest positions/velocities and use now h -> -h, going back in time T_max -> -T_max
    // doesent work
    // verlet(r1_new_verlet, r2_new_verlet, v1_new_verlet, v2_new_verlet, m1, m2, -0.05, T_max);

    // <<<

    return 0;
}
