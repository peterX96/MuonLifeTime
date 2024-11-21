#pragma once

#include <iostream>

TH1D* generator(int B=2, int seed=0);
TH1D* generator_exp(float tau, int integral, int seed=0);
