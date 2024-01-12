Matlab files for the experiments of the paper 'Seeded graph matching for the correlated Gaussian Wigner model via the projected power method' by E. Araya, G.Braun and H.Tyagi. 

Contain the following folders: 
-Synthetic data: corresponds to experiments in Sections 5.1 and 5.2. 
                -Fig.2 can be reproduced with Boosting_seedless_ppmgm_rev.m and
		different_initializations_ppmgm_rev.m
	        -Fig.3 can be reproduced with Boosting_grampa_fullqp_ppmgm_rev.m
	        -Fig.4 can be reproduced with different_iterations_ppmgm_rev.m
	        -Fig.5 can be reproduced with thresholding_experiments_rev.m
	        -Fig.6 can be reproduced with plot_heatmap.m (uses Data_fig6.mat)

-Real data: corresponds to experiments in Sections 5.3 (uses Shrec'16 database:)
	   
            -Fig.7 can be reproduced with plot_shape.m 
	    -Fig.8 can be reproduced with CDF_shrec16.m 

-Helpers: contain a variety of auxiliary methods used in the rest of the code. For exemple, the greedy matching GMWM algorithm (Alg.1 in the paper). 

-Graph matching algorithms: contain the following graph matching methods. 
                            -grampa algorithm --> matching_robust_spectral.m 
			   -umeyama algorithme --> matching_umeyama.m
			   -convex. relaxation to the Birkhoff polytope->matching_full_qp.m
			   -ppmgm-->matching_ppmgm.m

Code written by Ernesto ARAYA VALDIVIA, based on the code from different authors. Each author is credited in the individual .m files. 
	