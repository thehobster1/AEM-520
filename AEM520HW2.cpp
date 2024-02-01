#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main() {
  std::cout << "Hello, World!" << std::endl;
  

  // file pointer 
    fstream fout; 
  
    // opens an existing csv file or creates a new file. 
    fout.open("reportcard.csv", ios::out | ios::app); 
  
    cout << "Enter the details of 5 students:"
         << " roll name maths phy chem bio" 
    << endl; 
  
    int i, roll, phy, chem, math, bio; 
    string name; 
  
    // Read the input 
    for (i = 0; i < 5; i++) { 
  
        cin >> roll 
            >> name 
            >> math 
            >> phy 
            >> chem 
            >> bio; 
  
        // Insert the data to file 
        fout << roll << ", "
             << name << ", "
             << math << ", "
             << phy << ", "
             << chem << ", "
             << bio 
             << "\n"; 
    }

    return 0;
}

void create(string file, double x[], double y[]) 
{ 
    // file pointer 
    fstream fout; 
  
    // opens an existing csv file or creates a new file. 
    fout.open(file, ios::out | ios::app); 
  
    
  
    int i;
  
    // Read the input 
    for (i = 0; i < sizeof(x); i++) { 
  
        fout << x[i] << ", "
            << y[i] << ", "
            <<"\n"; 
  
    } 
}