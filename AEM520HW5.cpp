#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <filesystem>

#include <sstream>
#include <iomanip>

#include <ctime>

std::stringstream ss;



using namespace std;

//namespace fs = std::filesystem;

class funcValues {
    public:

        vector<double> xValue; 
        vector<double> yValue;
        vector<double> yDer;
        vector<double> yDerDer;
        int gc = 0;
        double x0;
        double xf;
        double dx;
        double t;
        double dt;

        funcValues(double n, double initial, double final){
            //xValue.resize(n);
            //yValue.resize(n);
            //double x = 0;
            x0 = initial;
            xf = final;
            dx = (xf - x0)/n;
            //t = time;
            for (int i = 0; i < n; i++){
                xValue.push_back(dx * i + dx*0.5);
            }
        }
};

void create(funcValues &func, string fmt){
    fstream fout;

    fout.open(fmt + ".csv", ios::out);

    int n = func.gc/2;
    
    for (int i = n; i < func.xValue.size() - n; i++){

        //cout << i << '\n';
        fout << func.xValue.at(i) << ", " << func.yValue.at(i) << '\n';
    }
}

void createFolder(vector<funcValues> &func, string folder){

    string s;

    //filesystem fs;

    std::filesystem::create_directory(folder);
    std::filesystem::current_path(folder);

    for (int i = 0; i < func.size(); i+=5){

        

        if (i == 0){

            s = folder + "_00000";
            s = s + to_string(i);
        }
        
        else if (i < 100){

            s = folder + "_0000";
            s = s + to_string(i);
        }

        else if (i < 1000){

            s = folder + "_000";
            s = s + to_string(i);
        }

        else if (i < 10000){

            s = folder + "_00";
            s = s + to_string(i);
        }

        else if (i < 100000){

            s = folder + "_0";
            s = s + to_string(i);
        }

        else {
            s = folder + "_";
            s = s + to_string(i);
        }


        create(func[i], s);
    

    }

    std::filesystem::current_path("..");
}

void sinFunc(funcValues &func){
    //for (int i; i < func.xValue.size() + 1; i++)
    for (auto & e : func.xValue){
        func.yValue.push_back(sin(e)+2);
    }
};

void specFunc(funcValues &func){
    //for (int i; i < func.xValue.size() + 1; i++)
    for (auto & e : func.xValue){
        if (e < M_PI) {
            func.yValue.push_back(e);
        }
        else {
            func.yValue.push_back(-e + 2*M_PI);
        }
    }
};

void printY(string fmt, funcValues &func)
{
    cout << '\n' << fmt << '\n';
    //for (int i = 0 + func.gc/2; i < func.xValue.size() - func.gc/2; i++)
    for (int i = 0; i < func.xValue.size(); i++)
        cout << func.xValue.at(i) << ' ' << func.yValue.at(i) << '\n';
    cout << '\n';
};

void printDer(string fmt, funcValues &func)
{
    cout << '\n' << fmt << '\n';
    //for (int i = 0 + func.gc/2; i < func.xValue.size() - func.gc/2; i++)
    for (int i = 0; i < func.yDer.size(); i++)
        cout << func.xValue.at(i+func.gc/2) << ' ' << func.yDer.at(i) << '\n';
    cout << '\n';
};


void d1o2 (funcValues &func){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        double val = (-1.0/2.0*func.yValue.at(i-1) + 1.0/2.0*func.yValue.at(i+1))/(func.dx*M_PI);
        func.yDer.push_back(val);
    }
}

void d1o4 (funcValues &func){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        double val = (1.0/12.0*func.yValue.at(i-2) - 2.0/3.0*func.yValue.at(i-1) + 2.0/3.0*func.yValue.at(i+1)-1.0/12.0*func.yValue.at(i+2))/(func.dx*M_PI);
        func.yDer.push_back(val);
    }
}

void d1o6 (funcValues &func){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        double val = (-1.0/60.0*func.yValue.at(i-3) + 3.0/20.0*func.yValue.at(i-2) - 3.0/4.0*func.yValue.at(i-1) + 3.0/4.0*func.yValue.at(i+1)-3.0/20.0*func.yValue.at(i+2) + 1.0/60.0*func.yValue.at(i+3))/(func.dx*M_PI);
        func.yDer.push_back(val);
    }
}

void d2o2 (funcValues &func){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        double val = (-1.0/2.0*func.yDer.at(i-1) + 1.0/2.0*func.yDer.at(i+1))/(func.dx*M_PI);
        func.yDerDer.push_back(val);
    }
}



void nextVecLin(funcValues &func, funcValues &next, double c, double v, double cfl){
    

    int n = func.gc/2;

    for(int i = n; i < (func.yValue.size()-n); i++){
        next.yValue.push_back(func.yValue.at(i) - c*func.yDer.at(i-n)*func.dt + v*func.yDerDer.at(i-n)*func.dt);
        //cout << '\n' << next.yValue.at(i-n) << '\n';
    }

    //double uMax = * max_element(next.yValue.begin(), next.yValue.end());

    next.dt = next.dx*cfl/c;

    next.t = next.dt + func.t;

}

void nextVecNon(funcValues &func, funcValues &next, double c, double v, double cfl){
    

    int n = func.gc/2;

    for(int i = n; i < (func.yValue.size()-n); i++){
        next.yValue.push_back(func.yValue.at(i) - c*func.yValue.at(i-n)*func.yDer.at(i-n)*func.dt + v*func.yDerDer.at(i-n)*func.dt);
        //cout << '\n' << next.yValue.at(i-n) << '\n';
    }

    double uMax = * max_element(next.yValue.begin(), next.yValue.end());

    next.dt = next.dx*cfl/(c*uMax);

    next.t = next.dt + func.t;

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


void ghostNodes2(int n, funcValues &func){
    func.gc = n;
    double xMin = func.xValue.front();
    double xMax = func.xValue.back();
    //cout <<"xMin: " << xMin <<'\n';
    //cout << "xMax: " << xMax << '\n';
    for (int i = 0; i < n/2; i++){
        // xMin -= func.dx; 
        // func.xValue.insert(func.xValue.begin(), xMin);
        // xMax += func.dx;
        //cout <<"xMax Changing: "<< xMax << ' ';
        //func.xValue.push_back(xMax);
        func.yDer.push_back(func.yDer.at(2*i));
        func.yDer.insert(func.yDer.begin(), func.yDer.at(-2*(i)+func.yDer.size()-2));

    }
    //cout << '\n';
};

void getDer(int n, funcValues &func){
    ghostNodes(n, func);

    if (n == 2){
        d1o2(func);
        ghostNodes2(n, func);
        d2o2(func);
    }
    else if (n == 4){
        d1o4(func);
    }
    else {
        d1o6(func);
    }

    //printDer("getDer", func);
}

// void getDerDer(int n, funcValues &func){
//     ghostNodes2(n, func);

//     if (n == 2){
//         d1o2(func);
//     }
//     else if (n == 4){
//         d1o4(func);
//     }
//     else {
//         d1o6(func);
//     }

//     //printDer("getDer", func);
// }

int main(){


    clock_t start = clock();

    //vector <funcValues> linWave;
    vector <funcValues> nonWave;
    int n;
    int order;
    double initial;
    double final;
    
    double cfl;
    
    double c;
    double v;

	ifstream input;


    input.open("input.txt", ios::in);
	string line;

    // while (getline(input,line)) {
    //     getline(input,line);
    //     if(line.find("#") == 0) {
    //         continue;
    //     }
    //     input >> n;
    //     //cout << n << endl;
    //     input >> order;
    //     //cout << order << endl;
    //     input >> initial;
    //     //cout << initial << endl;
    //     input >> final;
    //     input >> cfl;
    //     //input >> i_max;
    //     input >> c;
        
        
    //     if (input.eof()){
    //         break;
    //     }
        
    // }

    

    input >> n >> order >> initial >> final >> cfl >> c >> v;

    int tNodes = n*40;

    initial = initial*M_PI;
    final = final*M_PI;
    
    //funcValues linWaveZero(n, initial, final);
    //linWaveZero.t = 0;
    
    //specFunc(linWaveZero);

    //double uMax = * max_element(linWaveZero.yValue.begin(), linWaveZero.yValue.end());
    //linWaveZero.dt = linWaveZero.dx*cfl/c;

    funcValues nonWaveZero(n, initial, final);
    nonWaveZero.t = 0;
    specFunc(nonWaveZero);
    double uMax = * max_element(nonWaveZero.yValue.begin(), nonWaveZero.yValue.end());
    nonWaveZero.dt = nonWaveZero.dx*cfl/(c*uMax);

    //linWave.push_back(linWaveZero);
    nonWave.push_back(nonWaveZero);



    for (int i = 1; i < tNodes; i++){

        //printY("t = " + to_string((i-1)*dt) +": ", linWave[i-1]);

        //funcValues linWaveLast = linWave.at(i-1);

        //getDer(order, linWave[i-1]);

        //funcValues linWaveNext(n, initial, final);

        //nextVecLin(linWave[i-1], linWaveNext, c, v, cfl);

        //cout<<"HERE";

        //printY("t = " + to_string((i)*dt) +": ", linWaveNext);

        //linWave.push_back(linWaveNext);


        getDer(order, nonWave[i-1]);

        funcValues nonWaveNext(n, initial, final);

        nextVecNon(nonWave[i-1], nonWaveNext, c, v, cfl);

        //cout<<"HERE";

        //printY("t = " + to_string((i)*dt) +": ", linWaveNext);

        nonWave.push_back(nonWaveNext);
        
    }

    clock_t here = clock() - start;

    cout << ((float)here)/CLOCKS_PER_SEC << '\n';

    //cout<< "TOTAL linWave: "<< linWave.size()<< '\n';

    //createFolder(linWave, "linWave");
    createFolder(nonWave, "nonWave");

    start = clock()-start;

    cout << ((float)start)/CLOCKS_PER_SEC;
    
    return 420;
}