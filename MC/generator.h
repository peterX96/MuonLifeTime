#pragma once

#include "constants.h"

TH1D* generator(int B=2, int seed=0, float min_edge_hist=min_edge, float max_edge_hist=max_edge, int used_ch = 1880);
TH1D* generator_exp(float tau, int integral, int seed=0, float min_edge_hist=min_edge, float max_edge_hist=max_edge, int used_ch = 1880);
template <typename T>
T lin_interpol(T x, T* x_arr, T* y_arr, size_t size);
double generator_theta(TRandom3 *rg, double *theta_arr, double *cdf_arr, size_t arr_size);
