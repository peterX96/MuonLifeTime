#pragma once

const int Nchannels = 2048;
const double min_edge = 0., max_edge = 2048.;
const double cal_par1 = 0.01182, cal_par0 = 0.02727; //time[\mu s] = cal_par1*time[ch] + cal_par0

const double tauplus_values[4] = {2., 2.1, 2.2, 2.3}; //0
const double tauminus_values[4] = {0.7, 0.8, 0.9, 1.}; //1
const double R_values[4] = {1, 1.3, 1.8, 2.3}; //2
const double Nweeks_ratio[4] = {1, 3, 4.5, 6}; //3
const double B_values[3] = {.7, 1, 1.3}; //4, to be multiplied by B_data

//slab geometry
const double slab_length = 182., slab_width = 21., slab_height = 2.59; //in cm
