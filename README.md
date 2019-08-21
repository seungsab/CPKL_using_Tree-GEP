# CPKL_using_Tree-GEP

Compositional Kernel Learning using Tree-based Genetic Programming (CPKL_Tree-GEP) for Gaussian Proecess Regression

The source codes in this repository are uploaded for the replication of results section.
These source codes are only working for the 2-dimensional Mathematical test function ("Branin function")

Three methods can be implemented by running following m-files
- Run_Singel_kernel_GPML.m: Gaussian Process using single kernels
- Run_AUTO_GPR_GPML_HCPKL.m: Gaussian Process using Hierarchical CPK
- Run_AUTO_GPR_GPML_GPlab.m: Gaussian Process using CPKL_Tree-GEP

All results given in the manuscript are in the folders as "1_Test_functions" and "2_Physical_models"
- To see the results from Matehmatical test functions, run "Plot_Results_Test_function.m"
- To see the results from Physcial models, run "Plot_Results_Physcial_models.m"

If you have any questions, please feel free to contact me: seungsab@gmail.com (Seung-Seop Jin).



Before you run the codes, you need the following MATLAB toolboxes. 

[1] GPML toolbox (for Gaussian process): http://www.gaussianprocess.org/gpml/code/matlab/doc/
- C.E. Rasmussen, C.K.I. Williams, Gaussian processes for machine learning, MIT Press, Cambridge, Mass., 2006.

[2] GPlab toolbox (for Tree-based Genetic programing): http://gplab.sourceforge.net/
- William E, Northern J Genetic Programming Lab (GPLab) Tool Set Version 3.0. In: 2008 IEEE Region 5 Conference, 17-20 April 2008 2008. pp 1-6. doi:10.1109/TPSD.2008.4562729
