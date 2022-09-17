README

Source code for manuscript "On the use of cross-validation for the calibration of the adaptive lasso"
by N. Ballout, L. Etievant and V. Viallon

For questions, comments or remarks about the code please contact (ALL or just one) (@xxx.xx)

sessionInfo()
The code has been excuted using
R version 4.0.3 (2020-10-10)
(Platform: x86_64-w64-mingw32/x64 (64-bit))
Running under: Windows 10 x64 (build 19043)
witpackage versions ****

To reproduce the results presented in the manuscript,
run the xx for the simulation and the xx for the application.
All figures will be stored in the figures subfolder.

Note that:
-the code can run in parallel on linux, to do, change the value of the ncores in the parameters in the simulation.R script (same for application.R)
-to reproduce only the figures, skip the Run section in simulations.R, load xxx.Rdata into the intermediate_results subfolder and run only the Figures section (same for applications.R)

Information on the data set can be found in xxxx.