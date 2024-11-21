The new code work like this:
./main.exe N create a Tree inside PartialFit.root call Data_N, inside each Data_N 
there are 4 leafes, each one contain for that run the distribution of fit_bin,fit_err,avg_bin,dev.
(dev) deviation is the distribution of (MEANfit_bin-mu)/SIGMAfit_bin for the histograms.
Has expected the last one is symmetric respect to 0.

PartialFit.root has already inside Data_N trees from N=300,350,400,...,800.

PlotComp is a macro for plot a single histogram extracting it from a specific branch.

For sure, there is room for improvement, that is an early example.
