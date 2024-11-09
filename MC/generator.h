#pragma once

#include <iostream>

const int Nchannels = 2048;
const double min_edge = 0., max_edge = 2048.;

TH1D* generator(int B=2);
