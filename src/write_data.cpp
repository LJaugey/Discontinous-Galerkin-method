#include <fstream>
#include <iostream>
#include "parameters.hpp"
#include "write_data.hpp"

using namespace std;


void write_data(double data, string filename, string ext)
{
    ofstream file;
    file.open("output/"+filename+"."+ext, ios::app);

    file.precision(10);
    file<<std::fixed;
    
    file<<data<<"\t";

    file.close();
}


void clear_data(string filename, string ext)
{
    ofstream file;
    file.open("output/"+filename+"."+ext, ios::trunc);
    file.close();
}