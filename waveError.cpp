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


void reader()
{


    // I will read u_data in here


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
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_bp_file> <output_bp_file> [dh]" << std::endl;
        return 1;
    }

    std::string input_bp_file = argv[1];


    std::cout << "DEBUG: Input BP file: " << input_bp_file << std::endl;
    std::cout << "DEBUG: Output BP file: " << input_bp_file << std::endl;

    // Define error bounds to test
    std::vector<double> error_bounds = {0.1, 0.05, 0.01, 0.005, 0.001, 0.0005};
    std::cout << "DEBUG: Will test error bounds: ";
    for (size_t i = 0; i < error_bounds.size(); i++) {
        std::cout << error_bounds[i];
        if (i < error_bounds.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
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


        return 0;
}
#endif
