INTRODUCTION
------------

** Caution:

 * Our provided code depends on a crucial MATLAB-based modeling system called 'CVX', developed by Professor Stephen Boyd. More details, including guidelines, can be found at https://cvxr.com/cvx/.


SIMULATION STUDIES
------------------
 * Matlab code for simulation studies are under the folder “simulation”.

 Main code:
 * “mainpart_debias_spa_lala.m” gives the Matlab code for running our method with Lasso-Lasso regularizor, and “mainpart_debias_spa_scsc.m” gives the Matlab code for running our method with SCAD-SCAD regularizor.
 
 Other code: 
 * “Samples_generation_despa.m” gives the Matlab code for generating the Toeplitz design samples, and “Samples_generation_despa_scalefree.m”gives the Matlab code for generating the  scalefree-graph design samples.

 * "only_lasso.m", "nodewise_lasso.m" gives the Matlab code for solving the estimators under the regularized framework with Lasso penalty.

 * "lla_spa.m", "nodewise_lla_spa.m" gives the Matlab code for solving the estimators under the regularized framework with SCAD penalty.

 * "debiased_lasso_spa.m", "debiased_lasso_spa_scad.m" gives the core-part Matlab code for solving the de-biasing estimators (for single variable test) along with other crucial resuls under the regularized framework with Lasso and SCAD, respectively.

 * "debiased_lasso_groupG.m", "debiased_scad_groupG.m" gives the core-part Matlab code for solving the de-biasing estimators (for group variables test) along with other crucial resuls under the regularized framework with Lasso and SCAD, respectively.



REAL DATA STUDIES
-----------------
1. Finance data

* The data file used here is 'clean_2023-10-monthly.csv'. The original macro-data and the related data cleaning code were provided by Michael W. McCracken. For more information, please visit the following URL: https://research.stlouisfed.org/econ/mccracken/fred-databases/

* The file 'main_md.m' provides the Matlab code for running our method with the aforementioned finance macro data.


2. Gene data

# Gene data used in this part is sourced from Fan et al. 2017, where the raw data has been organized into "x_real.csv" and "y_real.csv".
# For more details, refer to https://doi.org/10.1111/rssb.12166

# The file 'main_gene.m' provides the Matlab code for running our method with the aforementioned gene data.

