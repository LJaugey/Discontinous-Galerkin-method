#ifndef WRITE_DATA_HPP
#define WRITE_DATA_HPP

#include <fstream>
#include <iostream>
#include "ND_array.hpp"
#include "parameters.hpp"


template<class T>
void write_data(T M, string filename, string ext="out")
{
    ofstream file;
    file.open("output/"+filename+"."+ext, ios::app);

    file.precision(10);
    file<<std::fixed;
    
    file<<M;

    file.close();
}

void write_data(double data, string filename, string ext="out");


void clear_data(std::string filename, std::string ext="out");



#endif