# Overview of files 

-----------

* `cww_compute_j0.m` - Computes the minimum wavelet scale j_0, used for computing DWT.
* `cww_compute_wave_kernels.m` - Extracts the internal wavelet kernel structure used for computation of DWT.
* `cww_extract_is_per_from_bd_mode.m` - converts a string to a bool value, based on whether or not to use periodic boundary extension.  
* `cww_extract_vm_from_wname.m` - Returns the number of vanishing moments based on a text string with wavelet names.
* `cww_fastwht_2d.m` - Two dimensional (unitary) fast Walsh-Hadamard transform.
* `cww_findim_handle_1d_cs.m` - Computes the one dimensional matrix-vector multiplication with the matrix P_{Omega}*U_{fwht}*U_{idwt}.
* `cww_get_phi_pieces.m` - Computes approximations to the scaling and wavelet functions using the cascade algorithm.
* `cww_get_phi_walsh_pieces.m` - Computes the Walsh transform of the wavelet and scaling function.
* `cww_get_scaling_matrix.m` - Mapping from wavelet coeff. to pointwise samples of the wavelet.
* `cww_handle_1d_cs.m` - Function handle for matrix-vector multiplication with a subsampled Walsh sampling and wavelet reconstruction matrix in one dimension. 
* `cww_handle_1d.m` - Function handle for matrix-vector multiplication with a Walsh sampling and wavelet reconstruction matrix in one dimension. 
* `cww_handle_2d_cs.m` - Function handle for matrix-vector multiplication with a subsampled Walsh sampling and wavelet reconstruction matrix in two dimensions. 
* `cww_handle_2d.m` - Function handle for matrix-vector multiplication with a Walsh sampling and wavelet reconstruction matrix in two dimensions. 
* `cww_kernel_CWW_bd.m` - Kernel function used for computation of Walsh sampling and boundary corrected wavelet reconstruction.  
* `cww_kernel_CWW_per.m` - Kernel function used for computation of Walsh sampling and periodic wavelet reconstruction.  
* `cww_map_wcoeff_to_func_vals_1d.m` - Maps wavelet coeff. to pointwise samples of the wavelet in one dimension.
* `cww_map_wcoeff_to_func_vals_2d.m` - Maps wavelet coeff. to pointwise samples of the wavelet in two dimensions.
* `cww_sample_walsh_1d.m` - Samples a function f:[0,1] -> R using walsh functions.
* `cww_sample_walsh_2d.m` - Samples a function f:[0,1]^2 -> R using walsh functions.
* `cww_sph1_power_law.m` - One dimensional sampling pattern based on power law sampling
* `cww_visualize_1d_samp_patt.m` - Converts a one dimensional sampling pattern to an image


-----------


