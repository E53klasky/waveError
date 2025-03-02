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
#include <cstdlib>
#include <tuple>
// for I/O
#include "adios2.h"

// for lossy compression and decompression
#include "mgard/compress_x.hpp"
//#include "/home/eklasky/local/mgard_install/include/mgard/compress_x.hpp"


// TODO:
// ADD DEBUG STATEMENTS LOTS OF THEM ALL THE TIME YOU ARE NOT THAT GOOD AT CODING LOTS OF THEM AND THE VAR--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// fidn out how to get file name of bp file
// find out what time step either through the bp file or through the code ??
// read it in u_data
// then compress it  subtract U - U~
// and write it out


std::tuple<std::vector<double>, int> reader(const std::string& bp_file_path, int steps)
{

        std::cout<<"DEBUG: Starting to read in U_data" <<std::endl;
        int dims = 0;
        adios2::ADIOS adios;
        auto inIO = adios.DeclareIO("Reader");
        inIO.SetEngine("BP");
        // maybe
        auto reader = inIO.Open(bp_file_path, adios2::Mode::Read);


        if (!reader) {
        std::cerr << "ERROR: Failed to open BP file: " << bp_file_path << std::endl;
          return std::make_tuple(std::vector<double>(), -1);
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
                dims = shape[0];
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

         return std::make_tuple(u_data, dims);
         //return u_data;

}



void compression_experiment(const std::vector<double>& u_data,
                            const std::vector<std::size_t>& shape)
{
    // Define the set of absolute error bounds
    std::vector<double> error_bounds = {0.1, 0.05, 0.01, 0.005, 0.001, 0.0005};

    // Storage for decompressed and diff data
    std::vector<double> decompressed_data(u_data.size());
    std::vector<double> diff_data(u_data.size());

    // Loop over each error bound
    for (double error_bound : error_bounds) {
        std::cout << "DEBUG: Running compression/decompression with error bound: " << error_bound << std::endl;

        size_t compressed_size = 0;
        void* compressed_array_cpu = NULL;

        // MGARD config
        mgard_x::Config config;

        // OG size
         std::cout << "DEBUG: Original size before compression: "
              << (u_data.size() * sizeof(double))  << " Bytes" << std::endl;


        // Compress
        //
        //  WHAT IS THE VAR BEING CALLED  HERE ------------------------------------------------------------------------------------
        mgard_x::compress(2, mgard_x::data_type::Double, shape, error_bound, 1.0,
                          mgard_x::error_bound_type::ABS, u_data.data(),
                          compressed_array_cpu, compressed_size, config, false);

        std::cout << "DEBUG: Compression completed. Compressed size: " << compressed_size/(1024.0*1024.0) << " MiB" << std::endl;

        // Decompress
        void* decompressed_array_cpu = NULL;

        // fix print and store in a vector ------------------------------------------------------------------------------------
        mgard_x::decompress(compressed_array_cpu, compressed_size,
                            decompressed_array_cpu, config, false);
        // something

        std::cout << "DEBUG: Decompression completed " << std::endl;

//      std::cout << "DEBUG: Decompressed size: "
//              << (u_data.size() * sizeof(double)) / (1024.0 * 1024.0) << " MiB" << std::endl;
     std::cout << "DEBUG: Decompressed size: " << (u_data.size() * sizeof(double)) << " Bytes " << std::endl;

      //  here!!!!!!!!!! -----------------------------------------------------------------------------------------------


        std::cout << "" << std::endl;
        // Copy decompressed data
        std::memcpy(decompressed_data.data(), decompressed_array_cpu,
                    u_data.size() * sizeof(double));

        // Compute the difference (u - u')
        for (size_t i = 0; i < u_data.size(); i++) {
            diff_data[i] = u_data[i] - decompressed_data[i];
        }

        // Free memory
        if (compressed_array_cpu) free(compressed_array_cpu);
        if (decompressed_array_cpu) free(decompressed_array_cpu);

        // Compute and print potential energy
        //double potential_energy = calc_PE(diff_data, shape);
        //std::cout << "DEBUG: PE for error bound " << error_bound << " = " << potential_energy << std::endl;
    }
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
    // Assuming reader returns a tuple<std::vector<double>, std::size_t>
    auto [u_data, num_data] = reader(input_bp_file, steps); // Unpacking the tuple
    std::cout << "DEBUG: u_data read in successfully" << std::endl;

    // Since num_data represents the total size, and we assume it's a square shape
    std::size_t num_rows = num_data;
    std::size_t num_cols = num_data;

    // Create shape vector
    std::vector<std::size_t> shape = {num_rows, num_cols};
    std::cout << "u_data shape: " << shape[0] << " x " << shape[1] << " Y " << std::endl;

    // Call compression_experiment
    std::cout<<""<<std::endl;

    compression_experiment(u_data, shape);

    std::cout<<"DEBUG: finshed getting the difference"<<std::endl;

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
