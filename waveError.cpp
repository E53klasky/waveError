#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <dirent.h>
#include <random>
#include <algorithm>  // For std::min_element and std::max_element
#include <string.h>
#include <iomanip>
#include "adios2.h"
#include "mgard/compress_x.hpp"

// potential energy
double calc_PE(double *u, double *pe, double dh, size_t Nx, size_t Ny)
{
    double PE = 0.0;
    double ux, uy;
    size_t r_curr, k;
    for (size_t r=1; r<Nx; r++) {
        r_curr = r * Ny;
        for (size_t c=1; c<Ny; c++) {
            k     = r_curr + c;
            ux    = (u[k] - u[k-Ny]) / dh;
            uy    = (u[k] - u[k-1])/dh;
            pe[k] = ux*ux + uy*uy;
            PE   += pe[k];
        }
    }
    return PE * 0.5;
}

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

// root-of-mean-square error
double calc_rmse(double *data_f, double *data_g, double *diff, size_t num_data)
{
    double rmse = 0.0;
    for (size_t i=0; i<num_data; i++) {
        diff[i] = data_f[i] - data_g[i];
        rmse += (diff[i]*diff[i]) / (double) num_data;
    }
    return std::sqrt(rmse);
}

int main(int argc, char **argv) {
    int cnt_argv = 1;
    std::string fname_f(argv[cnt_argv++]);
    // simulation spatial resolution
    double dh      = std::stof(argv[cnt_argv++]);
    // from the init_ts step in frame_f to compare against frame_g
    size_t init_ts = std::stoi(argv[cnt_argv++]);
    std::string fname_err(argv[cnt_argv++]);
    
    // Define error bounds to test
    std::vector<double> error_bounds = {0.1, 0.05, 0.01, 0.005, 0.001, 0.0005};
    
    // Keep tol_s0 always 0
    double tol_s0 = 0.0;
    
    std::cout << "original file: " << fname_f.c_str() << "\n";
    std::cout << "output file: " << fname_err.c_str() << "\n";
    std::cout << "dh = " << dh << ", init_ts = " << init_ts << ", tol(s0) = " << tol_s0 << "\n";
    std::cout << "Testing multiple error bounds for tol(s1): ";
    for (auto& eb : error_bounds) {
        std::cout << eb << " ";
    }
    std::cout << "\n\n";

    adios2::ADIOS ad;
    adios2::IO reader_io_f = ad.DeclareIO("Original");
    reader_io_f.SetEngine("BP");
    adios2::Engine reader_f = reader_io_f.Open(fname_f, adios2::Mode::ReadRandomAccess);
    adios2::Variable<double> variable_f = reader_io_f.InquireVariable<double>("u_data");

    size_t Nx = variable_f.Shape()[0];
    size_t Ny = variable_f.Shape()[1];
    size_t num_data = Nx * Ny;
    std::vector<double> var_f(num_data);
    // difference data
    std::vector<double> err_s1(num_data);
    std::vector<double> err_s0(num_data);
    std::vector<double> PE_s0(num_data);
    std::vector<double> PE_s1(num_data);
    std::vector<double> diff(num_data);

    // compression parameters
    mgard_x::Config config;
    config.lossless = mgard_x::lossless_type::Huffman_Zstd;
    config.dev_type = mgard_x::device_type::SERIAL;
    //config.dev_type = mgard_x::device_type::CUDA;

    // Set the shape of the data
    std::vector<mgard_x::SIZE> shape{Nx, Ny};

    variable_f.SetStepSelection({init_ts, 1});
    reader_f.Get(variable_f, var_f);
    reader_f.PerformGets();

    double v_max = *std::max_element(var_f.begin(), var_f.end());
    double v_min = *std::min_element(var_f.begin(), var_f.end());
    std::cout << "data value ranges: {" << v_max << ", " << v_min << "}\n\n";

    // s=0 compression (only done once as tol_s0 is fixed at 0)
    void *compressed_s0 = NULL;
    void *decompressed_s0 = NULL;
    size_t compressed_size_s0;
    
    mgard_x::compress(2, mgard_x::data_type::Double, shape, tol_s0, (double)0.0,
                      mgard_x::error_bound_type::ABS, var_f.data(),
                      compressed_s0, compressed_size_s0, config, false);

    mgard_x::decompress(compressed_s0, compressed_size_s0, decompressed_s0, config, false);

    double rmse_s0 = calc_rmse(var_f.data(), (double *)decompressed_s0, err_s0.data(), num_data);
    double PE_err_s0 = calc_PE(err_s0.data(), PE_s0.data(), dh, shape[0], shape[1]);
    
    double *data_s0 = static_cast<double*>(decompressed_s0);
    
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "s=0 (fixed): L2 = " << rmse_s0
              << " (rel " << rmse_s0 / (v_max - v_min) << "), "
              << "PE_e = " << PE_err_s0 << ", "
              << "sqrt(PE_e) = " << std::sqrt(PE_err_s0) << ", "
              << "eb/sqr = " << tol_s0/std::sqrt(PE_err_s0) << ", "
              << "compression ratio = " << static_cast<double>(num_data * sizeof(double)) / compressed_size_s0
              << "\n\n";
    
    // Loop through different error bounds for s=1
    std::cout << "=== Results for different error bounds (tol_s1) ===\n";
    for (auto& tol_s1 : error_bounds) {
        std::cout << "\n--- Testing tol_s1 = " << tol_s1 << " ---\n";
        
        void *compressed_s1 = NULL;
        void *decompressed_s1 = NULL;
        size_t compressed_size_s1;
        
        // s=1 compression with current error bound
        mgard_x::compress(2, mgard_x::data_type::Double, shape, tol_s1, (double)1.0,
                          mgard_x::error_bound_type::ABS, var_f.data(),
                          compressed_s1, compressed_size_s1, config, false);

        mgard_x::decompress(compressed_s1, compressed_size_s1, decompressed_s1, config, false);
        
        double *data_s1 = static_cast<double*>(decompressed_s1);
        
        // Calculate metrics
        for (size_t i = 0; i < num_data; i++) {
            diff[i] = var_f[i] - data_s1[i];
        }
        
        double PE_diff = calc_PE(diff.data(), 0.5, shape[0], shape[1]);
        double rmse_diff = calc_rmse(diff.data(), (double *)decompressed_s1, err_s1.data(), num_data);
        double rmse_s1 = calc_rmse(var_f.data(), (double *)decompressed_s1, err_s1.data(), num_data);
        double PE_err_s1 = calc_PE(err_s1.data(), PE_s1.data(), dh, shape[0], shape[1]);
        
        // Print results for current error bound
        std::cout << "s=1: L2 = " << rmse_s1
                  << " (rel " << rmse_s1 / (v_max - v_min) << "), "
                  << "PE_e = " << PE_err_s1 << ", "
                  << "sqrt(PE_e) = " << std::sqrt(PE_err_s1) << ", "
                  << "eb/sqr = " << tol_s1/std::sqrt(PE_err_s1) << ", "
                  << "compression ratio = " << static_cast<double>(num_data * sizeof(double)) / compressed_size_s1
                  << "\n";

        std::cout << "diff: L2 = " << rmse_diff
                  << " (rel " << rmse_diff / (v_max - v_min) << "), "
                  << "PE_e = " << PE_diff << ", "
                  << "sqrt(PE_e) = " << std::sqrt(PE_diff) << ", "
                  << "eb/sqr = " << tol_s1 / std::sqrt(PE_diff) << ", "
                  << "compression ratio = " << static_cast<double>(num_data * sizeof(double)) / (compressed_size_s0 + compressed_size_s1)
                  << "\n";
        
        // Free memory for current iteration
        delete[] data_s1;
    }
    
    delete[] data_s0;
    reader_f.Close();
    
    return 0;
}
