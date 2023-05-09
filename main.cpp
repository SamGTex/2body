#include "vector_functions.h"
#include <algorithm>
#include <chrono>
#include <cstdlib> //for exit(1) to catch errors
#include <fstream>
#include <sstream>
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

// read config file, get parameters
tuple<vec, vec, vec, vec, double, double, double, double, int> get_config(string path){
    ifstream configFile(path); // open the config file for reading

    // define variables to store the parameters
    vector<double> r1, r2, v1, v2;
    double m1, m2, h, T_max;
    int algorithm;

    // read the parameters from the config file
    string line;
    char delim = '=';
    while (getline(configFile, line)) {
        // split line into paramName and paramValue using "="
        stringstream ss(line);
        string paramName, paramValue;
        getline(ss, paramName, delim);
        getline(ss, paramValue);

        // split paramValue into two parts using ","
        stringstream ss2(paramValue);
        string part;
        vector<double> values;
        while (getline(ss2, part, ',')) {
            values.push_back(stod(part));
        }

        if (paramName == "r1") {
            r1 = values;
        } else if (paramName == "r2") {
            r2 = values;
        } else if (paramName == "v1") {
            v1 = values;
        } else if (paramName == "v2") {
            v2 = values;
        } else if (paramName == "m1") {
            m1 = stod(paramValue);
        } else if (paramName == "m2") {
            m2 = stod(paramValue);
        } else if (paramName == "h") {
            h = stod(paramValue);
        } else if (paramName == "T_max") {
            T_max = stod(paramValue);
        } else if (paramName == "algorithm") {
            algorithm = stoi(paramValue);
        }
    }

    configFile.close(); // close the config file

    return make_tuple(r1, r2, v1, v2, m1, m2, h, T_max, algorithm);
}


int main() {
    const string path_config = "config.txt";

    // read in parameters from config file
    tuple<vec, vec, vec, vec, double, double, double, double, int> params = get_config(path_config);
  
    vec r1 = get<0>(params);
    vec r2 = get<1>(params);
    vec v1 = get<2>(params);
    vec v2 = get<3>(params);
    double m1 = get<4>(params);
    double m2 = get<5>(params);
    double h = get<6>(params);
    double T_max = get<7>(params);
    int algorithm = get<8>(params);

    // start time measurement
    auto start_euler = chrono::high_resolution_clock::now(); 
    
    // start algorithm depending on config file parameter "algorithm"
    if (algorithm==0){
        cout << "Starting Euler algorithm" << endl;
        euler(r1, r2, v1, v2, m1, m2, h, T_max);

    }else if (algorithm==1){
        cout << "Starting Verlet algorithm" << endl;
        verlet(r1, r2, v1, v2, m1, m2, h, T_max);
    } else{
        cerr << "Algorithm not found!" << endl;
    }

    // stop time and print duration
    auto end_euler = chrono::high_resolution_clock::now();
    auto duration_euler = chrono::duration_cast<chrono::milliseconds>(end_euler - start_euler);
    cout << "Algorithm duration: " << duration_euler.count() << " ms" << endl;

    return 0;
}
