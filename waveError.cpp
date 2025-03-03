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


// potential energy
double calc_PE(double *u, double dh, size_t Nx, size_t Ny)
{
    double PE = 0.0;
    double ux, uy;
    size_t r_curr, k;
    for (size_t r=1; r<Nx; r++) {
        r_curr = r * Ny;
        for (size_t c=1; c<Ny; c++) {
            k  = r_curr + c;
            ux = (u[k] - u[k-Ny]) / dh;
            uy = (u[k] - u[k-1])/dh;
            PE += ux*ux + uy*uy;
        }
    }
    return PE * 0.5;
}


std::tuple<std::vector<double>, int> reader(const std::string& bp_file_path, int steps)
{

        //std::cout<<"DEBUG: Starting to read in U_data" <<std::endl;
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

        //std::cout<<"DEBUG:: Turning the reader engine on" << std::endl;
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
                            const std::vector<std::size_t>& shape, double dh)
{
    // Define the set of absolute error bounds
    std::vector<double> error_bounds = {0.1, 0.05, 0.01, 0.005, 0.001, 0.0005};

    // Storage for decompressed and diff data
    std::vector<double> decompressed_data(u_data.size());
    double* diff_data = new double[u_data.size()];
std::vector<std::tuple<double, double, double>> results;
    // Loop over each error bound
    for (double error_bound : error_bounds) {
        std::cout << "DEBUG: Running compression/decompression with error bound: " << error_bound << std::endl;

        size_t compressed_size = 0;
        void* compressed_array_cpu = NULL;

        // MGARD config
        mgard_x::Config config;

        // OG size
         //std::cout << "DEBUG: Original size before compression: "
          //    << (u_data.size() * sizeof(double))  << " Bytes" << std::endl;





        // Compress
        mgard_x::compress_status_type status = mgard_x::compress(2, mgard_x::data_type::Double, shape, error_bound, 1.0,
                          mgard_x::error_bound_type::ABS, u_data.data(),
                          compressed_array_cpu, compressed_size, config, false);

        if (status != mgard_x::compress_status_type::Success) {
                std::cerr << "ERROR: Compression failed with status code " << static_cast<int>(status) << std::endl;
                return;
        }

         //std::cout << "DEBUG: u_data size: " << u_data.size() << " elements" << std::endl;
        //std::cout << "DEBUG: u_data total bytes: " << compressed_size << " bytes" << std::endl;
        //std::cout << "DEBUG: Memory address of u_data: " << static_cast<const void*>(u_data.data()) << std::endl;


        // Decompress
        void* decompressed_array_cpu = NULL;





        mgard_x::decompress(compressed_array_cpu, compressed_size,
                            decompressed_array_cpu, config, false);


       // std::cout << "DEBUG: Decompression completed " << std::endl;
        //std::cout << "DEBUG: Decompressed data memory address: " << static_cast<void*>(decompressed_array_cpu) << std::endl;

        //std::cout << "DEBUG: Decompressed size: "
       //       << (u_data.size() * sizeof(double)) / (1024.0 * 1024.0) << " MiB" << std::endl;
     //std::cout << "DEBUG: Decompressed size: " << (u_data.size() * sizeof(double)) << " Bytes " << std::endl;

      //  here!!!!!!!!!! -----------------------------------------------------------------------------------------------


        //std::cout << "" << std::endl;
        // Copy decompressed data


        std::memcpy(decompressed_data.data(), decompressed_array_cpu,
                    u_data.size() * sizeof(double));


        // Compute the difference (u - u')

        // Compute the difference (u - u')
       // std::cout << "DEBUG: Computing difference between original and decompressed data" << std::endl;
        for (size_t i = 0; i < u_data.size(); i++) {
            diff_data[i] = u_data[i] - decompressed_data[i];
        }

        double PE = calc_PE(diff_data, dh, shape[0], shape[1]);
        double sqrt_PE = std::sqrt(PE);
        double ratio = error_bound / sqrt_PE;


        std::cout << "DEBUG: PE = " << PE << ", sqrt(PE) = " << sqrt_PE << ", eb/sqrt(PE) = " << ratio << std::endl;
        std::cout << "\n" << std::endl;
        // Store results


        results.emplace_back(error_bound, sqrt_PE, ratio);

        //delete[] diff_data;
        // Free memory
        if (compressed_array_cpu) free(compressed_array_cpu);
        if (decompressed_array_cpu) free(decompressed_array_cpu);

        // Compute and print potential energy
        //double potential_energy = calc_PE(diff_data, shape);
        //std::cout << "DEBUG: PE for error bound " << error_bound << " = " << potential_energy << std::endl;

    }
    delete[] diff_data;
}





#ifdef POST_PROCESSING_EXEC
int main(int argc, char** argv){

    std::cout << "DEBUG: Starting program" << std::endl;
        // Check command line arguments
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_bp_file>  <number of steps>  <dh>  input file will be the same as the output file " << std::endl;
        return 1;
    }

    std::string input_bp_file = argv[1];
    int steps = std::stoi(argv[2]);
    double dh = std::stod(argv[3]);


    std::cout << "DEBUG: Input BP file: " << input_bp_file << std::endl;
    std::cout << "DEBUG: Output BP file: " << input_bp_file << std::endl;
    std::cout << "DEBUG: number of steps being read is: " << steps<< std::endl;
    std::cout << "DEBUG: dh is " <<dh <<std::endl;
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

    compression_experiment(u_data, shape, dh);

    std::cout<<"DEBUG: finshed getting the difference"<<std::endl;


        // write out diff
        // write out uncompress_U_data


    std::cout << "All processing completed" << std::endl;
    //std::cout <<"That is very good yes" << std::endl;

        return 0;
}
