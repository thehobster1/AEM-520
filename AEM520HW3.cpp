#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>

using namespace std;


// double DerX(double *, double *, int *);
// double DerX4(double*, double*,int*);
// double DerX6(double*, double*,int*);
// double DerX_2(double*, double*,int*);
// double DerX4_2(double*, double*,int*);
// double DerX6_2(double*, double*,int*);
// double* fill_ghostNodes(double*, int,int);
// void create(string,double*,double*,int);

// input nodes, ghost cells, and order from text file

//vector<double> calSin

class funcValues {
    public:

        vector<double> xValue; 
        vector<double> yValue;
        int gc;
        double dx;

        funcValues(double n){
            //xValue.resize(n);
            //yValue.resize(n);
            //double x = 0;
            dx = 2.0/n;
            for (int i = 0; i < n; i++){
                xValue.push_back(dx * i + dx*0.5);
            }
        }
};

// double nodes(int n) {

// }

// double sin() {

// }

void println(string fmt, funcValues &func)
{
    cout << '\n' << fmt << '\n';
    //for (int i = 0 + func.gc/2; i < func.xValue.size() - func.gc/2; i++)
    for (int i = 0; i < func.xValue.size(); i++)
        cout << func.xValue.at(i) << ' ' << func.yValue.at(i) << '\n';
    cout << '\n';
};

// void printlny(string fmt, funcValues &func)
// {
//     //for (int i = 0 + func.gc/2; i < func.xValue.size() - func.gc/2; i++)
//     cout << "Y Values: " << '\n';
//     for (int i = 0; i < func.yValue.size(); i++)
//         cout << func.yValue.at(i) << ' ' << '\n';//<< func.yValue.at(i) << '\n';
//     cout << '\n';
// };

void sinFunc(funcValues &func){
    //for (int i; i < func.xValue.size() + 1; i++)
    for (auto & e : func.xValue){
        func.yValue.push_back(sin(e*M_PI));
    }
};

void cosFunc(funcValues &func){
    //for (int i; i < func.xValue.size() + 1; i++)
    for (auto & e : func.xValue){
        func.yValue.push_back(cos(e*M_PI));
    }
};

void negsinFunc(funcValues &func){
    //for (int i; i < func.xValue.size() + 1; i++)
    for (auto & e : func.xValue){
        func.yValue.push_back(-1*sin(e*M_PI));
    }
};

void ghostNodes(int n, funcValues &func){
    func.gc = n;
    double xMin = func.xValue.front();
    double xMax = func.xValue.back();
    //cout <<"xMin: " << xMin <<'\n';
    //cout << "xMax: " << xMax << '\n';
    for (int i = 0; i < n/2; i++){
        xMin -= func.dx; 
        func.xValue.insert(func.xValue.begin(), xMin);
        xMax += func.dx;
        //cout <<"xMax Changing: "<< xMax << ' ';
        func.xValue.push_back(xMax);
        func.yValue.push_back(func.yValue.at(2*i));
        func.yValue.insert(func.yValue.begin(), func.yValue.at(-2*(i)+func.yValue.size()-2));

    }
    //cout << '\n';
};

void d1o2 (funcValues &func, funcValues &funcDer){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        double val = (-1.0/2.0*func.yValue.at(i-1) + 1.0/2.0*func.yValue.at(i+1))/(func.dx*M_PI);
        funcDer.yValue.push_back(val);
    }
}

void d1o4 (funcValues &func, funcValues &funcDer){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        double val = (1.0/12.0*func.yValue.at(i-2) - 2.0/3.0*func.yValue.at(i-1) + 2.0/3.0*func.yValue.at(i+1)-1.0/12.0*func.yValue.at(i+2))/(func.dx*M_PI);
        funcDer.yValue.push_back(val);
    }
}

void d1o6 (funcValues &func, funcValues &funcDer){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        double val = (-1.0/60.0*func.yValue.at(i-3) + 3.0/20.0*func.yValue.at(i-2) - 3.0/4.0*func.yValue.at(i-1) + 3.0/4.0*func.yValue.at(i+1)-3.0/20.0*func.yValue.at(i+2) + 1.0/60.0*func.yValue.at(i+3))/(func.dx*M_PI);
        funcDer.yValue.push_back(val);
    }
}

void d2o2 (funcValues &func, funcValues &funcDer){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        double val = (1.0*func.yValue.at(i-1) -2.0*func.yValue.at(i) + 1.0*func.yValue.at(i+1))/(pow(func.dx*M_PI,2.0));
        funcDer.yValue.push_back(val);
    }
}

void d2o4 (funcValues &func, funcValues &funcDer){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        double val = (-1.0/12.0*func.yValue.at(i-2) + 4.0/3.0*func.yValue.at(i-1) -5.0/2.0*func.yValue.at(i)+ 4.0/3.0*func.yValue.at(i+1)-1.0/12.0*func.yValue.at(i+2))/(pow(func.dx*M_PI,2.0));
        funcDer.yValue.push_back(val);
    }
}

void d2o6 (funcValues &func, funcValues &funcDer){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        double val = (1.0/90.0*func.yValue.at(i-3) - 3.0/20.0*func.yValue.at(i-2) + 3.0/2.0*func.yValue.at(i-1) - 49.0/18.0*func.yValue.at(i) + 3.0/2.0*func.yValue.at(i+1)-3.0/20.0*func.yValue.at(i+2) + 1.0/90.0*func.yValue.at(i+3))/(pow(func.dx*M_PI,2.0));
        funcDer.yValue.push_back(val);
    }
}

void create(funcValues &func, string fmt){
    fstream fout;

    fout.open(fmt + ".csv", ios::out);
    
    for (int i = 0; i < func.xValue.size(); i++){
        fout << func.xValue.at(i) << ", " << func.yValue.at(i) << '\n';
    }
}

void error(funcValues &sample, funcValues &func, funcValues &error){
    for (int i = 0; i < func.xValue.size(); i++){
        double err = abs(sample.yValue.at(i) - func.yValue.at(i));
        error.yValue.push_back(err);
    }
}


int main(){
    double xMin = 0;
    double xHigh = 2;
    // double cfl = 0.01;
    // double n = 10;
    // double dt = 1.0/n*cfl;

    // funcValues sinWaveInitial(n);
    // //sinWave.yValue = generate(sinWave.yValue.begin(), sinWave.yValue.end(), )
    // //println("x: ", sinWaveInitial.xValue);
    // sinFunc(sinWaveInitial);
    // ghostNodes(6, sinWaveInitial);
    // println("Sin Wave Initial", sinWaveInitial);

    // funcValues der1n10o6(n);

    // d1o4(sinWaveInitial, der1n10o6);

    // println("D1O6", der1n10o6);
    // //printlny("Sin Wave Initial", sinWaveInitial);

    funcValues sin_10_2(10);
    funcValues sin_100_2(100);
    funcValues sin_1000_2(1000);

    sinFunc(sin_10_2);
    sinFunc(sin_100_2);
    sinFunc(sin_1000_2);

    create(sin_10_2, "sin_10");
    create(sin_100_2, "sin_100");
    create(sin_1000_2, "sin_1000");

    ghostNodes(2, sin_10_2);
    ghostNodes(2, sin_100_2);
    ghostNodes(2, sin_1000_2);

    funcValues sin_10_4(10);
    funcValues sin_100_4(100);
    funcValues sin_1000_4(1000);

    sinFunc(sin_10_4);
    sinFunc(sin_100_4);
    sinFunc(sin_1000_4);

    ghostNodes(4, sin_10_4);
    ghostNodes(4, sin_100_4);
    ghostNodes(4, sin_1000_4);

    funcValues sin_10_6(10);
    funcValues sin_100_6(100);
    funcValues sin_1000_6(1000);

    sinFunc(sin_10_6);
    sinFunc(sin_100_6);
    sinFunc(sin_1000_6);

    ghostNodes(6, sin_10_6);
    ghostNodes(6, sin_100_6);
    ghostNodes(6, sin_1000_6);



    funcValues cos_10(10);
    funcValues cos_100(100);
    funcValues cos_1000(1000);

    cosFunc(cos_10);
    cosFunc(cos_100);
    cosFunc(cos_1000);


    create(cos_10, "cos_10");
    create(cos_100, "cos_100");
    create(cos_1000, "cos_1000");


    funcValues negsin_10(10);
    funcValues negsin_100(100);
    funcValues negsin_1000(1000);

    negsinFunc(negsin_10);
    negsinFunc(negsin_100);
    negsinFunc(negsin_1000);

//Der 1 All Nodes and Orders
    funcValues der_1_2_10(10);
    funcValues der_1_2_100(100);
    funcValues der_1_2_1000(1000);

    d1o2(sin_10_2, der_1_2_10);
    d1o2(sin_100_2, der_1_2_100);
    d1o2(sin_1000_2, der_1_2_1000);

    create(der_1_2_10, "der_1_2_10");
    create(der_1_2_100, "der_1_2_100");
    create(der_1_2_1000, "der_1_2_1000");

    funcValues der_1_2_err_10(10);
    funcValues der_1_2_err_100(100);
    funcValues der_1_2_err_1000(1000);

    error(cos_10, der_1_2_10, der_1_2_err_10);
    error(cos_100, der_1_2_100, der_1_2_err_100);
    error(cos_1000, der_1_2_1000, der_1_2_err_1000);

    create(der_1_2_err_10, "der_1_2_err_10");
    create(der_1_2_err_100, "der_1_2_err_100");
    create(der_1_2_err_1000, "der_1_2_err_1000");


    funcValues der_1_4_10(10);
    funcValues der_1_4_100(100);
    funcValues der_1_4_1000(1000);

    d1o4(sin_10_4, der_1_4_10);
    d1o4(sin_100_4, der_1_4_100);
    d1o4(sin_1000_4, der_1_4_1000);

    create(der_1_4_10, "der_1_4_10");
    create(der_1_4_100, "der_1_4_100");
    create(der_1_4_1000, "der_1_4_1000");

    funcValues der_1_4_err_10(10);
    funcValues der_1_4_err_100(100);
    funcValues der_1_4_err_1000(1000);

    error(cos_10, der_1_4_10, der_1_4_err_10);
    error(cos_100, der_1_4_100, der_1_4_err_100);
    error(cos_1000, der_1_4_1000, der_1_4_err_1000);

    create(der_1_4_err_10, "der_1_4_err_10");
    create(der_1_4_err_100, "der_1_4_err_100");
    create(der_1_4_err_1000, "der_1_4_err_1000");


    funcValues der_1_6_10(10);
    funcValues der_1_6_100(100);
    funcValues der_1_6_1000(1000);

    d1o6(sin_10_6, der_1_6_10);
    d1o6(sin_100_6, der_1_6_100);
    d1o6(sin_1000_6, der_1_6_1000);

    create(der_1_6_10, "der_1_6_10");
    create(der_1_6_100, "der_1_6_100");
    create(der_1_6_1000, "der_1_6_1000");

    funcValues der_1_6_err_10(10);
    funcValues der_1_6_err_100(100);
    funcValues der_1_6_err_1000(1000);

    error(cos_10, der_1_6_10, der_1_6_err_10);
    error(cos_100, der_1_6_100, der_1_6_err_100);
    error(cos_1000, der_1_6_1000, der_1_6_err_1000);

    create(der_1_6_err_10, "der_1_6_err_10");
    create(der_1_6_err_100, "der_1_6_err_100");
    create(der_1_6_err_1000, "der_1_6_err_1000");


// //Der 2 Non-c

    funcValues der_2_2_10(10);
    funcValues der_2_2_100(100);
    funcValues der_2_2_1000(1000);

    d2o2(sin_10_2, der_2_2_10);
    d2o2(sin_100_2, der_2_2_100);
    d2o2(sin_1000_2, der_2_2_1000);

    create(der_2_2_10, "der_2_2_10");
    create(der_2_2_100, "der_2_2_100");
    create(der_2_2_1000, "der_2_2_1000");

    funcValues der_2_2_err_10(10);
    funcValues der_2_2_err_100(100);
    funcValues der_2_2_err_1000(1000);

    error(negsin_10, der_2_2_10, der_2_2_err_10);
    error(negsin_100, der_2_2_100, der_2_2_err_100);
    error(negsin_1000, der_2_2_1000, der_2_2_err_1000);

    create(der_2_2_err_10, "der_2_2_err_10");
    create(der_2_2_err_100, "der_2_2_err_100");
    create(der_2_2_err_1000, "der_2_2_err_1000");


    funcValues der_2_4_10(10);
    funcValues der_2_4_100(100);
    funcValues der_2_4_1000(1000);

    d2o4(sin_10_4, der_2_4_10);
    d2o4(sin_100_4, der_2_4_100);
    d2o4(sin_1000_4, der_2_4_1000);

    create(der_2_4_10, "der_2_4_10");
    create(der_2_4_100, "der_2_4_100");
    create(der_2_4_1000, "der_2_4_1000");

    funcValues der_2_4_err_10(10);
    funcValues der_2_4_err_100(100);
    funcValues der_2_4_err_1000(1000);

    error(negsin_10, der_2_4_10, der_2_4_err_10);
    error(negsin_100, der_2_4_100, der_2_4_err_100);
    error(negsin_1000, der_2_4_1000, der_2_4_err_1000);

    create(der_2_4_err_10, "der_2_4_err_10");
    create(der_2_4_err_100, "der_2_4_err_100");
    create(der_2_4_err_1000, "der_2_4_err_1000");


    funcValues der_2_6_10(10);
    funcValues der_2_6_100(100);
    funcValues der_2_6_1000(1000);

    d2o6(sin_10_6, der_2_6_10);
    d2o6(sin_100_6, der_2_6_100);
    d2o6(sin_1000_6, der_2_6_1000);

    create(der_2_6_10, "der_2_6_10");
    create(der_2_6_100, "der_2_6_100");
    create(der_2_6_1000, "der_2_6_1000");

    funcValues der_2_6_err_10(10);
    funcValues der_2_6_err_100(100);
    funcValues der_2_6_err_1000(1000);

    error(negsin_10, der_2_6_10, der_2_6_err_10);
    error(negsin_100, der_2_6_100, der_2_6_err_100);
    error(negsin_1000, der_2_6_1000, der_2_6_err_1000);

    create(der_2_6_err_10, "der_2_6_err_10");
    create(der_2_6_err_100, "der_2_6_err_100");
    create(der_2_6_err_1000, "der_2_6_err_1000");


//  Conservative Der 2
    ghostNodes(2, der_1_2_10);
    ghostNodes(2, der_1_2_100);
    ghostNodes(2, der_1_2_1000);

    ghostNodes(4, der_1_4_10);
    ghostNodes(4, der_1_4_100);
    ghostNodes(4, der_1_4_1000);

    ghostNodes(6, der_1_6_10);
    ghostNodes(6, der_1_6_100);
    ghostNodes(6, der_1_6_1000);

    println("Sin: ", sin_10_4);

    println("der_1_4_10: ", der_1_4_10);

    funcValues der_2_2_cons_10(10);
    funcValues der_2_2_cons_100(100);
    funcValues der_2_2_cons_1000(1000);

    d1o2(der_1_2_10, der_2_2_cons_10);
    d1o2(der_1_2_100, der_2_2_cons_100);
    d1o2(der_1_2_1000, der_2_2_cons_1000);

    create(der_2_2_cons_10, "der_2_2_cons_10");
    create(der_2_2_cons_100, "der_2_2_cons_100");
    create(der_2_2_cons_1000, "der_2_2_cons_1000");

    funcValues der_2_2_cons_err_10(10);
    funcValues der_2_2_cons_err_100(100);
    funcValues der_2_2_cons_err_1000(1000);

    error(negsin_10, der_2_2_cons_10, der_2_2_cons_err_10);
    error(negsin_100, der_2_2_cons_100, der_2_2_cons_err_100);
    error(negsin_1000, der_2_2_cons_1000, der_2_2_cons_err_1000);

    create(der_2_2_cons_err_10, "der_2_2_cons_err_10");
    create(der_2_2_cons_err_100, "der_2_2_cons_err_100");
    create(der_2_2_cons_err_1000, "der_2_2_cons_err_1000");


    funcValues der_2_4_cons_10(10);
    funcValues der_2_4_cons_100(100);
    funcValues der_2_4_cons_1000(1000);

    d1o4(der_1_4_10, der_2_4_cons_10);
    d1o4(der_1_4_100, der_2_4_cons_100);
    d1o4(der_1_4_1000, der_2_4_cons_1000);
    println("der_2_4_10: ", der_2_4_10);
    println("der_2_4_cons_10: ", der_2_4_cons_10);

    create(der_2_4_cons_10, "der_2_4_cons_10");
    create(der_2_4_cons_100, "der_2_4_cons_100");
    create(der_2_4_cons_1000, "der_2_4_cons_1000");

    funcValues der_2_4_cons_err_10(10);
    funcValues der_2_4_cons_err_100(100);
    funcValues der_2_4_cons_err_1000(1000);

    error(negsin_10, der_2_4_cons_10, der_2_4_cons_err_10);
    error(negsin_100, der_2_4_cons_100, der_2_4_cons_err_100);
    error(negsin_1000, der_2_4_cons_1000, der_2_4_cons_err_1000);

    create(der_2_4_cons_err_10, "der_2_4_cons_err_10");
    create(der_2_4_cons_err_100, "der_2_4_cons_err_100");
    create(der_2_4_cons_err_1000, "der_2_4_cons_err_1000");


    funcValues der_2_6_cons_10(10);
    funcValues der_2_6_cons_100(100);
    funcValues der_2_6_cons_1000(1000);

    d1o6(der_1_6_10, der_2_6_cons_10);
    d1o6(der_1_6_100, der_2_6_cons_100);
    d1o6(der_1_6_1000, der_2_6_cons_1000);

    create(der_2_6_cons_10, "der_2_6_cons_10");
    create(der_2_6_cons_100, "der_2_6_cons_100");
    create(der_2_6_cons_1000, "der_2_6_cons_1000");

    funcValues der_2_6_cons_err_10(10);
    funcValues der_2_6_cons_err_100(100);
    funcValues der_2_6_cons_err_1000(1000);

    error(negsin_10, der_2_6_cons_10, der_2_6_cons_err_10);
    error(negsin_100, der_2_6_cons_100, der_2_6_cons_err_100);
    error(negsin_1000, der_2_6_cons_1000, der_2_6_cons_err_1000);

    create(der_2_6_cons_err_10, "der_2_6_cons_err_10");
    create(der_2_6_cons_err_100, "der_2_6_cons_err_100");
    create(der_2_6_cons_err_1000, "der_2_6_cons_err_1000");

    return 69;

};