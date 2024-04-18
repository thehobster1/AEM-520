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
        vector<vector<vector<double>>> uValue;
        //vector<double> uDer;
        //vector<double> uDerDer;
        // vector<double> vValue;
        // vector<double> vDer;
        // vector<double> vDerDer;
        int gc = 0;
        double x0;
        double xf;
        double dx;
        double y0;
        double yf;
        double dy;
        double t;
        double dt;
        double rho;
        double err;

        funcValues(double n, double initialx, double finalx, double initialy, double finaly, double rho0, double err0){
            
            x0 = initialx;
            xf = finalx;
            y0 = initialy;
            yf = finaly;
            dx = (xf - x0)/n;
            dy = (yf - y0)/n;
            rho = rho0;
            err = err0;
            uValue = vector<vector<vector<double>>>(n, vector<vector<double>>(n, vector<double>(12)));
            //t = time;
            for (int i = 0; i < n; i++){
                xValue.push_back(dx * i + dx*0.5);
                yValue.push_back(dy * i + dy*0.5);
            }
        }
};

void create(funcValues &func, string fmt){
    char buff[64];
    fstream fout;
    

    fout.open(fmt + ".csv", ios::out);
    

    int n = func.gc/2;
    
    for (int i = n; i < func.xValue.size() - n; i++){
        for (int j = n; j < func.yValue.size() - n; j++){
        //cout << i << '\n';

            //sprintf_s(buff, "%f,%f,%f,%f\n", func.xValue.at(i), func.yValue.at(j), func.uValue[i][j][0], func.uValue[i][j][1]);
            //sprintf_s(buff, "%f,%f,%f\n", i-n, j-n, abs(func.uValue[i][j][0] + func.uValue[i][j][1]));
            fout << abs(func.uValue[i][j][0] + func.uValue[i][j][1]) << ' ';
            //fout << func.xValue.at(i) << ", " << func.yValue.at(j) << ", " << func.uValue[i][j][0] << ", " << func.uValue[i][j][1] <<'\n';
        }
        fout << '\n';
    }
    fout.close();
}

void createPressure(funcValues &func, string fmt){
    char buff[64];
    fstream fout;
    

    fout.open(fmt + ".csv", ios::out);
    

    int n = func.gc/2;
    
    for (int i = n; i < func.xValue.size() - n; i++){
        for (int j = n; j < func.yValue.size() - n; j++){
        
            fout << abs(func.uValue[i][j][11]) << ' ';
            
        }
        fout << '\n';
    }
    fout.close();
}

void createAxes(funcValues &func, string fmt){
    char buff[64];
    fstream fout;

    fout.open(fmt + ".csv", ios::out);

    int n = func.gc/2;
    
    for (int i = n; i < func.xValue.size() - n; i++){
        
        //cout << i << '\n';

            //sprintf_s(buff, "%f,%f,%f,%f\n", func.xValue.at(i), func.yValue.at(j), func.uValue[i][j][0], func.uValue[i][j][1]);
            sprintf_s(buff, "%f,%f\n", func.xValue.at(i), func.yValue.at(i));
            fout << buff;
            //fout << func.xValue.at(i) << ", " << func.yValue.at(j) << ", " << func.uValue[i][j][0] << ", " << func.uValue[i][j][1] <<'\n';
        
    }
    fout.close();
}

void createFolder(vector<funcValues> &func, string folder){

    string s;

    //filesystem fs;

    std::filesystem::create_directory(folder);
    std::filesystem::current_path(folder);

    int count = 1;
    if (func.size()>500){
        count = func.size()/500;
    }

    int x = func.size();

    for (int i = 0; i < func.size(); i+=count){

        

        // if (i == 0){

        //     s = folder + "_000000";
        //     s = s + to_string(i);
        // }

        if (i < 10){
            s = folder + "_000000";
            s = s + to_string(i);
        }
        
        else if (i < 100){

            s = folder + "_00000";
            s = s + to_string(i);
        }

        else if (i < 1000){

            s = folder + "_0000";
            s = s + to_string(i);
        }

        else if (i < 10000){

            s = folder + "_000";
            s = s + to_string(i);
        }

        else if (i < 100000){

            s = folder + "_00";
            s = s + to_string(i);
        }

        else if (i < 1000000){
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


void createPressureFolder(vector<funcValues> &func, string folder){

    string s;

    //filesystem fs;

    std::filesystem::create_directory(folder);
    std::filesystem::current_path(folder);

    int count = 1;
    if (func.size()>500){
        count = func.size()/500;
    }

    int x = func.size();

    for (int i = 0; i < func.size(); i+=count){

        

        // if (i == 0){

        //     s = folder + "_000000";
        //     s = s + to_string(i);
        // }

        if (i < 10){
            s = folder + "_000000";
            s = s + to_string(i);
        }
        
        else if (i < 100){

            s = folder + "_00000";
            s = s + to_string(i);
        }

        else if (i < 1000){

            s = folder + "_0000";
            s = s + to_string(i);
        }

        else if (i < 10000){

            s = folder + "_000";
            s = s + to_string(i);
        }

        else if (i < 100000){

            s = folder + "_00";
            s = s + to_string(i);
        }

        else if (i < 1000000){
            s = folder + "_0";
            s = s + to_string(i);
        }

        else {
            s = folder + "_";
            s = s + to_string(i);
        }


        createPressure(func[i], s);
    

    }

    std::filesystem::current_path("..");
}

void createFolderTime(vector<funcValues> &func, string folder){

    string s;

    //filesystem fs;

    std::filesystem::create_directory(folder);
    std::filesystem::current_path(folder);

    int count = 1;
    // if (func.size()>1000){
    //     count = func.size()/1000;
    // }

    int x = func.size();

    double counter = 0;

    string time;

    for (int i = 0; i < func.size(); i+=count){

        if (func[i].t >= counter){

            time = to_string(counter);

            // time.erase(time.find_last_not_of('0') + 1, time.length()-1);
            // time.erase(time.find_last_not_of('.') + 1, time.length()-1);
            create(func[i], time);
            counter += 0.5;
            
        }

    }

    std::filesystem::current_path("..");
}

// void sinFunc(funcValues &func){
//     //for (int i; i < func.xValue.size() + 1; i++)
//     for (auto & e : func.xValue){
//         for (auto & e : func.yValue){
//             func.uValue.push_back(sin(e));
//         }
//     }
// };

void specFunc(funcValues &func, double cfl, double c, double nu){
    double uMax = 0;
    for (int i = 0; i < func.xValue.size(); i++){
        for (int j = 0; j < func.yValue.size(); j++){
            
            double u = sin(func.xValue.at(i))*sin(func.yValue.at(j))+1;
            double v = 0.1;
            if(uMax < abs(u)){
            uMax = abs(u); 
            }

            if(uMax < abs(v)){
                uMax = abs(v);
            }
            
            func.uValue[i][j][0] = u;
            func.uValue[i][j][1] = v;
        }
    }
    double dta = func.dx*cfl/(c*abs(uMax));

    double dtb = pow(func.dx, 2)/nu*cfl;

    if (dta < dtb){
        func.dt = dta;
    }
    else{
        func.dt = dtb;
    }
        
};

void printU(string fmt, funcValues &func)
{
    cout << '\n' << fmt << '\n';
    //for (int i = 0 + func.gc/2; i < func.xValue.size() - func.gc/2; i++)
    for (int i = 0; i < func.xValue.size(); i++){
        for (int j = 0; j < func.yValue.size(); i++){
            cout << func.xValue.at(i) << ' ' << func.yValue.at(j) << ' ' << func.uValue[i][j][0] << ' ' << func.uValue[i][j][1] << '\n';
        }
    }
    cout << '\n';
};

// void printV(string fmt, funcValues &func)
// {
//     cout << '\n' << fmt << '\n';
//     //for (int i = 0 + func.gc/2; i < func.xValue.size() - func.gc/2; i++)
//     for (int i = 0; i < func.xValue.size(); i++)
//         cout << func.xValue.at(i) << ' ' << func.uValue[1].at(i) << '\n';
//     cout << '\n';
// };

// void printDer(string fmt, funcValues &func)
// {
//     cout << '\n' << fmt << '\n';
//     //for (int i = 0 + func.gc/2; i < func.xValue.size() - func.gc/2; i++)
//     for (int i = 0; i < func.uDer.size(); i++)
//         cout << func.xValue.at(i+func.gc/2) << ' ' << func.uDer.at(i) << '\n';
//     cout << '\n';
// };




void d1o2 (funcValues &func){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        for(int j = n; j< func.yValue.size() - n; j++){
            func.uValue[i][j][2] = (-1.0/2.0*func.uValue[i-1][j][0] + 1.0/2.0*func.uValue[i+1][j][0])/(func.dx);
            func.uValue[i][j][3] = (-1.0/2.0*func.uValue[i][j-1][1] + 1.0/2.0*func.uValue[i][j+1][1])/(func.dy);
            func.uValue[i][j][4] = (-1.0/2.0*func.uValue[i-1][j][0] + 1.0/2.0*func.uValue[i+1][j][0])/(func.dx);
            func.uValue[i][j][5] = (-1.0/2.0*func.uValue[i][j-1][1] + 1.0/2.0*func.uValue[i][j+1][1])/(func.dy);

            func.uValue[i][j][10] = -1.0*func.rho*(func.uValue[i][j][2]*func.uValue[i][j][2] + 2*func.uValue[i][j][2]*func.uValue[i][j][3] +func.uValue[i][j][3]*func.uValue[i][j][3]);
            
            func.uValue[i][j][11] = 0;
            //func.uDer.push_back(val);
        }
    }
}




void d2o2 (funcValues &func){
    int n = func.gc/2;
    for (int i = n; i < func.xValue.size()-n; i++){
        for(int j = n; j< func.yValue.size() - n; j++){
            func.uValue[i][j][6] = (-1.0/2.0*func.uValue[i-1][j][2] + 1.0/2.0*func.uValue[i+1][j][2])/(func.dx);
            func.uValue[i][j][7] = (-1.0/2.0*func.uValue[i][j-1][3] + 1.0/2.0*func.uValue[i][j+1][3])/(func.dy);
            func.uValue[i][j][8] = (-1.0/2.0*func.uValue[i-1][j][4] + 1.0/2.0*func.uValue[i+1][j][4])/(func.dx);
            func.uValue[i][j][9] = (-1.0/2.0*func.uValue[i][j-1][5] + 1.0/2.0*func.uValue[i][j+1][5])/(func.dy);

            //func.uDer.push_back(val);
        }
    }
}

void nextVecBurger(funcValues &func, funcValues &next, double c, double nu, double cfl){
    

    int n = func.gc/2;

    double uMax = 0;

    double u;
    double v;


    for(int i = n; i < (func.xValue.size()-n); i++){
        int help = func.yValue.size() - n - 1;
        for (int j = n; j < (func.yValue.size() - n); j++){
        
            //u = (func.uValue[i-n][j-n].at(0) - c*func.uValue[i-n][j-n].at(0)*func.uValue[i-n][j-n][2]*func.dt - c*func.uValue[i-n][j-n].at(1)*func.uValue[i-n][j-n][3]*func.dt + nu*(func.uValue[i-n][j-n][6]- func.uValue[i-n][j-n][7])*func.dt);
            //v = (func.uValue[i-n][j-n].at(1) - c*func.uValue[i-n][j-n].at(0)*func.uValue[i-n][j-n][4]*func.dt - c*func.uValue[i-n][j-n].at(1)*func.uValue[i-n][j-n][5]*func.dt + nu*(func.uValue[i-n][j-n][8]- func.uValue[i-n][j-n][9])*func.dt);
            
            //u = (func.uValue[i][j].at(0) - c*func.uValue[i][j].at(0)*func.uValue[i][j][2]*func.dt - c*func.uValue[i][j].at(1)*func.uValue[i][j][3]*func.dt +  + nu*(func.uValue[i][j][6] + func.uValue[i][j][7])*func.dt - func.dt*1.0/(2.0*func.rho)*func.uValue[i+1][j][11] - func.uValue[i-1][j][11]);
            //v = (func.uValue[i][j].at(1) - c*func.uValue[i][j].at(0)*func.uValue[i][j][4]*func.dt - c*func.uValue[i][j].at(1)*func.uValue[i][j][5]*func.dt + nu*(func.uValue[i][j][8] + func.uValue[i][j][9])*func.dt - func.dt*1.0/(2.0*func.rho)*(func.uValue[i][j+1][11] - func.uValue[i][j-1][11]));

            u = func.uValue[i][j][0] - func.dt*(func.uValue[i][j][0]*func.uValue[i][j][2] + func.uValue[i][j][1]*func.uValue[i][j][3] + 1.0/func.rho*func.uValue[i][j][11]/func.dx - nu*(func.uValue[i][j][6] + func.uValue[i][j][7]));
            v = func.uValue[i][j][1] - func.dt*(func.uValue[i][j][0]*func.uValue[i][j][4] + func.uValue[i][j][1]*func.uValue[i][j][5] + 1.0/func.rho*func.uValue[i][j][11]/func.dy - nu*(func.uValue[i][j][8] + func.uValue[i][j][9]));
            //u = func.uValue[i][j][0] - func.dt*(func.uValue[i][j][2]*(func.uValue[i][j][2] - func.uValue[i-1][j][2]) + func.uValue[i][j][3]*(func.uValue[i][j][2] - func.uValue[i][j-1][2]) + 1.0/(2.0*func.rho)*(func.uValue[i+1][j][10] - func.uValue[i-1][j][10])) + nu*func.dt*(func.uValue[i+1][j][])
            if(uMax < abs(u)){
                uMax = abs(u);
            }

            if(uMax < abs(v)){
                uMax = abs(v);
            }
            
            next.uValue[i-n][j-n][0] = u;
            next.uValue[i-n][j-n][1] = v;
        }
        //cout << i << '\n';
    }



    // double uMax = * max_element(next.uValue.begin(), next.uValue.end());
    // double uMin = * min_element(next.uValue.begin(), next.uValue.end());

    // if (abs(uMin) > abs(uMax)){
    //     uMax = uMin;
    // }

    double dta = next.dx*cfl/(c*abs(uMax));

    double dtb = pow(next.dx, 2)/nu;



    if (dta < dtb){
        next.dt = dta;
    }
    else{
        next.dt = dtb;
    }

    //cout << next.dt << '\n';
    


    next.t = next.dt + func.t;

};

void ghostNodes(int n, funcValues &func){
    func.gc = n;
    double xMin = func.xValue.front();
    double xMax = func.xValue.back();
    double yMin = func.yValue.front();
    double yMax = func.yValue.back();
    
    
    for (int i = 0; i < n/2; i++){
        xMin -= func.dx; 
        func.xValue.insert(func.xValue.begin(), xMin);
        xMax += func.dx;
        func.xValue.push_back(xMax);

        yMin -= func.dy; 
        func.yValue.insert(func.yValue.begin(), yMin);
        yMax += func.dy;
        func.yValue.push_back(yMax);

        //vector <vector<double>> xGhost0;
        //vector <vector<double>> xGhost1;

        for (int j = i; j < func.uValue.size() - i; j++){

            // int help = func.uValue.size()-1-i;

            // vector <double> test = func.uValue[j][help];

            func.uValue[j].insert(func.uValue[j].begin(), func.uValue[j][func.uValue.size() - 1 - i]);
            
            func.uValue[j].push_back(func.uValue[j][i+1]);

            
        }
        
        vector <vector<double>> xGhost0 = func.uValue[func.uValue.size() - 1 -i];
        xGhost0[0] = {0,0,0,0,0,0,0,0,0,0};
        xGhost0[xGhost0.size()-1] = {0,0,0,0,0,0,0,0,0,0};
        vector <vector<double>> xGhost1 = func.uValue[i];
        xGhost1[0] = {0,0,0,0,0,0,0,0,0,0};
        xGhost1[xGhost1.size()-1] = {0,0,0,0,0,0,0,0,0,0};
        func.uValue.insert(func.uValue.begin(), xGhost0);

        func.uValue.push_back(xGhost1);

    }
    
};


void ghostNodes2(int n, funcValues &func){
    func.gc = n;
    
    
    for (int i = 0; i < n/2; i++){
        //vector <vector<double>> xGhost0;
        //vector <vector<double>> xGhost1;

        for (int j = 1; j < func.xValue.size() + i - 1; j++){

            
            int help = func.uValue.size()-1-(n/2-i);
            

            func.uValue[j][n/2-i-1][2] = func.uValue[j][func.uValue.size()-1-(n/2-i)][2];
            func.uValue[j][func.uValue.size()-(n/2-i)][2] = func.uValue[j][n/2+i][2];

            func.uValue[j][n/2-i-1][3] = func.uValue[j][func.uValue.size()-1-(n/2-i)][3];
            func.uValue[j][func.uValue.size()-(n/2-i)][3] = func.uValue[j][n/2+i][3];

            func.uValue[j][n/2-i-1][4] = func.uValue[j][func.uValue.size()-1-(n/2-i)][4];
            func.uValue[j][func.uValue.size()-(n/2-i)][4] = func.uValue[j][n/2+i][4];

            func.uValue[j][n/2-i-1][5] = func.uValue[j][func.uValue.size()-1-(n/2-i)][5];
            func.uValue[j][func.uValue.size()-(n/2-i)][5] = func.uValue[j][n/2+i][5];

            

            func.uValue[n/2-i-1][j][2] = func.uValue[func.uValue.size()-1-(n/2-i)][j][2];
            func.uValue[func.uValue.size()-(n/2-i)][j][2] = func.uValue[n/2+i][j][2];

            func.uValue[n/2-i-1][j][3] = func.uValue[func.uValue.size()-1-(n/2-i)][j][3];
            func.uValue[func.uValue.size()-(n/2-i)][j][3] = func.uValue[n/2+i][j][3];

            func.uValue[n/2-i-1][j][4] = func.uValue[func.uValue.size()-1-(n/2-i)][j][4];
            func.uValue[func.uValue.size()-(n/2-i)][j][4] = func.uValue[n/2+i][j][4];

            func.uValue[n/2-i-1][j][5] = func.uValue[func.uValue.size()-1-(n/2-i)][j][5];
            func.uValue[func.uValue.size()-(n/2-i)][j][5] = func.uValue[n/2+i][j][5];

            
            
        }
        

    }
    
};

void ghostNodesPressure(funcValues &func){
    int n = func.gc/2;
    
    
    for (int i = 0; i < n/2; i++){

        for (int j = 1; j < func.xValue.size() + i - 1; j++){

            
            int help = func.uValue.size()-1-(n/2-i);
            

            func.uValue[j][n/2-i-1][11] = func.uValue[j][func.uValue.size()-1-(n/2-i)][11];
            func.uValue[j][func.uValue.size()-(n/2-i)][11] = func.uValue[j][n/2+i][11];
            
            func.uValue[n/2-i-1][j][11] = func.uValue[func.uValue.size()-1-(n/2-i)][j][11];
            func.uValue[func.uValue.size()-(n/2-i)][j][11] = func.uValue[n/2+i][j][11];
        }
        

    }
    
};

void solvePoisson (funcValues &func){
    double err = 1;
    int n = func.gc/2;
    

    ghostNodesPressure(func);
    
    double sample;
    while (err > func.err){
        
        double maxErr = 0;
        for (int i = n; i < func.xValue.size()-n; i++){
            for(int j = n; j< func.yValue.size() - n; j++){
                sample = func.uValue[i][j][11];
                func.uValue[i][j][11] = (func.dx*func.dx*(func.uValue[i][j+1][11] + func.uValue[i][j-1][11]) + func.dy*func.dy*(func.uValue[i+1][j][11] + func.uValue[i-1][j][11]) - func.dx*func.dx*func.dy*func.dy*func.uValue[i][j][10])/(2.0*(func.dy*func.dy + func.dx*func.dx));
                
                if (abs (func.uValue[i][j][11] - sample)/func.uValue[i][j][11] > maxErr) {
                    maxErr = abs (func.uValue[i][j][11] - sample)/func.uValue[i][j][11];
                }
            }
        }
        
        err = maxErr;
        //cout << err << '\n';
    }
    //createPressure(func, "Pressure" + to_string(func.t));

}

void getDer(int n, funcValues &func){
    ghostNodes(n, func);

    if (n == 2){
        d1o2(func);
        solvePoisson(func);
        ghostNodes2(n, func);
        d2o2(func);
    }
    

    //printDer("getDer", func);
}



int main(){


    clock_t start = clock();

    //vector <funcValues> linWave;
    vector <funcValues> burger;
    int n;
    int order;
    double initialx;
    double finalx;
    double initialy;
    double finaly;
    double cfl;
    
    double c;
    double nu;
    double rho;
    double err;

	ifstream input;


    input.open("input.txt", ios::in);
	string line;



    

    input >> n >> order >> initialx >> finalx >> initialy >> finaly >> cfl >> c >> nu >> rho >> err;

    int tNodes = n*2;
    //int tNodes = 1;

    initialx = initialx*M_PI;
    finalx = finalx*M_PI;

    initialy = initialy*M_PI;
    finaly = finaly*M_PI;
    

    funcValues burgerZero(n, initialx, finalx, initialy, finaly, rho, err);
    burgerZero.t = 0;
    specFunc(burgerZero, cfl, c, nu);

    createAxes(burgerZero, "xy");
    // double uMax = * max_element(burgerZero.uValue.begin(), burgerZero.uValue.end());
    // double uMin = * min_element(burgerZero.uValue.begin(), burgerZero.uValue.end());

    // if (abs(uMin) > abs(uMax)){
    //     uMax = uMin;
    // }

    // double dta = burgerZero.dx*cfl/(c*abs(uMax));

    // double dtb = pow(burgerZero.dx, 2)/v*cfl;

    // if (dta < dtb){
    //     burgerZero.dt = dta;
    // }
    // else{
    //     burgerZero.dt = dtb;
    // }
    

    //linWave.push_back(linWaveZero);
    burger.push_back(burgerZero);



    for (int i = 1; i < tNodes + 1; i++){



        getDer(order, burger[i-1]);

        funcValues burgerNext(n, initialx, finalx, initialy, finaly, rho, err);

        nextVecBurger(burger[i-1], burgerNext, c, nu, cfl);

        burger.push_back(burgerNext);

        //cout<< i<<'\n';
        
    }
    //cout<< "Made it";

    clock_t here = clock() - start;

    cout << ((float)here)/CLOCKS_PER_SEC << '\n';

    string str = to_string(nu);

    str.erase(str.find_last_not_of('0') + 1, str.length()-1);
    str.erase(str.find_last_not_of('.') + 1, str.length()-1);

    createFolder(burger, "2Dburger" + str);


    createFolderTime(burger, "2DburgerTime" + str);

    createPressureFolder(burger, "Pressure" + str);

    start = clock()-start;

    cout << ((float)start)/CLOCKS_PER_SEC;
    
    return 2319;
}