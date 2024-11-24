#pragma once

const int Nchannels = 2048;
const double min_edge = 0., max_edge = 2048.;
const double cal_par1 = 0.01182, cal_par0 = 0.02727; //time[\mu s] = cal_par1*time[ch] + cal_par0
const double tauplus = 2.2, tauminus = 0.9; // in \mu s
