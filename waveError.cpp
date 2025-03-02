#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <thread>
#include <chrono>
#include <dirent.h>
#include <random>
#include <algorithm>
#include <string>
#include <cstring>

// for I/O
#include "adios2.h"

// for lossy compression and decompression
#include "mgard/compress_x.hpp"



// TODO:
// ADD DEBUG STATEMENTS LOTS OF THEM ALL THE TIME YOU ARE NOT THAT GOOD AT CODING LOTS OF THEM AND THE VAR--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// fidn out how to get file name of bp file
// find out what time step either through the bp file or through the code ??
// read it in u_data
// then compress it  subtract U - U~
// and write it out


std::vector<double> reader(const std::string& bp_file_path, int steps)
{

        std::cout<<"DEBUG: Starting to read in U_data" <<std::endl;

        adios2::ADIOS adios;
        auto inIO = adios.DeclareIO("Reader");
        inIO.SetEngine("BP");
        // maybe
        auto reader = inIO.Open(bp_file_path, adios2::Mode::Read);


        if (!reader) {
        std::cerr << "ERROR: Failed to open BP file: " << bp_file_path << std::endl;
        return std::vector<double>();
        }

        std::cout<<"DEBUG:: Turning the reader engine on" << std::endl;
        std::vector<double> u_data;
        adios2::Variable<double> var_udata;

        int stepCounter = 0;

        // this should work ??????
        for (int i = 0; i < steps; i++) {
        std::cout << "DEBUG: Reading step " << i << std::endl;
        auto status = reader.BeginStep();

        if (status != adios2::StepStatus::OK) {
                break;
        }

        var_udata = inIO.InquireVariable<double>("u_data");

        if (i == 0) {
                std::vector<std::size_t> shape = var_udata.Shape();
                size_t num_data = shape[0] * shape[1];
                u_data.resize(num_data);
        }

        // Get the data
        reader.Get(var_udata, u_data);

        // PerformGets if needed for deferred reading
        reader.PerformGets();

        reader.EndStep();
        stepCounter++;
        }

        reader.Close();
        std::cout<<"DEBUG: done reading all steps. The number of steps read is: "<< stepCounter << std::endl;


         return u_data;

}

void wrighter()
{
 // this is where I will right out the differences

}


void compression()
{
        // I will compress u_data here

}


void uncompression()
{
        // I will uncompress u_data here

}

void calU_data_Error()
{
        double error = 0.0;
        // I will calculate the U-U~ here and right out here
        // error call wrighter here
}



#ifdef POST_PROCESSING_EXEC
int main(int argc, char** argv){

    std::cout << "DEBUG: Starting program" << std::endl;
        // Check command line arguments
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_bp_file>  <number of steps>  input file will be the same as the output file " << std::endl;
        return 1;
    }

    std::string input_bp_file = argv[1];
    int steps = std::stoi(argv[2]);


    std::cout << "DEBUG: Input BP file: " << input_bp_file << std::endl;
    std::cout << "DEBUG: Output BP file: " << input_bp_file << std::endl;
    std::cout << "DEBUG: number of steps being read is: " << steps<< std::endl;
    // Define error bounds to test
    std::vector<double> error_bounds = {0.1, 0.05, 0.01, 0.005, 0.001, 0.0005};
    std::cout << "DEBUG: Will test error bounds: " << steps <<std::endl;

    for (size_t i = 0; i < error_bounds.size(); i++) {
        std::cout << error_bounds[i];
        if (i < error_bounds.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;


    // adds tsteps and convert to for loop
    std::vector<double> u_data = reader(input_bp_file, steps);
    std::cout<<"DEBUG: u_data read in successfuly"<< std::endl;
    // call reader
        // get data U_data all time steps
        // copy data
        // lossy compress data
        // subtract data
        //
        // uncompress data
        // write out diff
        // write out uncompress_U_data


    std::cout << "DEBUG: All processing completed" << std::endl;
    std::cout <<"That is very good yes" << std::endl;

        return 0;
}
#endif
