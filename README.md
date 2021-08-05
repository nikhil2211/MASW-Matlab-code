# MASW-Matlab-code
3. MASWaves Inversion 
3.1 Conducting the inversion analysis (MASWaves_inversion_MC) 
The function MASWaves_inversion_MC carries out the inversion analysis using a Monte Carlo based 
search procedure.
1) The program computes the theoretical (fundamental mode) dispersion curve for the stratified 
layer model defined by n, alpha_initial, beta_initial, rho and h_initial. The initial value of the 
dispersion misfit function is evaluated.
2) A Monte Carlo based search process is used in search of the shear wave velocity profile (i.e., 
value of beta and h for each layer) that provides the closest fit to the experimental data.
3) Inversion results are provided in the form of:
a. All sampled profiles
b. All sampled profiles whose theoretical dispersion curves fall within the upper/lower 
boundaries of the experimental curve ('accepted profiles') [optional]
Subfunctions:
MASWaves_theoretical_dispersion_curve_FDMA
MASWaves_misfit_MC
Input arguments
c_min Minimum testing Rayleigh wave phase velocity vector [m/s]
c_max Maximum testing Rayleigh wave phase velocity vector [m/s]
c_step Testing Rayleigh wave phase velocity increment [m/s]
delta_c Zero search initiation parameter [m/s] 
At wave number k_i the zero search is initiated at a phase velocity of 
max{c_(i-1)-delta_c , c_min}, where c_(i-1) is the theoretical Rayleigh wave 
phase velocity value at wave number k_(i-1)
n Number of finite thickness layers
n_unsat Number of unsaturated soil layers
(n_unsat = 0 for a fully saturated soil profile)
(n_unsat = n+1 for a fully unsaturated soil profile)
alpha_initial Initial estimate of compressional wave velocity [m/s] (array of length n+1)
nu_unsat Poisson's ratio
beta_initial Initial estimate of shear wave velocity [m/s] (array of length n+1)
rho Mass density vector [kg/m3
] (array of length n+1)
h_initial Initial estimate of layer thicknesses [m] (array of length n)
N_reversals Number of soil layers where velocity reversals are permitted 
(N_reversals = 0 for a normally dispersive profile).
9
Experimental fundamental mode dispersion curve
c_OBS Phase velocity [m/s]
lambda_OBS Wavelength [m]
up_low_boundary - 'yes': Upper/lower boundaries for the experimental dispersion curve are 
available.
- 'no': Upper/lower boundaries for the experimental dispersion curve are 
not available.
c_OBS_up Phase velocity, upper bound curve [m/s]
(can be assigned as 'nan' or [] if up_low_boundary = 'no').
c_OBS_low Phase velocity, lower bound curve [m/s]
(can be assigned as 'nan' or [] if up_low_boundary = 'no').
Search-control parameters
b_S Shear wave velocity search-control parameter
b_h Layer thickness search-control parameter
N_max Maximum number of iterations
e_max Maximum misfit (optional stopping criterion for MC search)
Output arguments
store_all All sampled profiles (cell array)
 For iteration no. i
 store_all{1,i}: Shear wave velocity vector [m/s]
 store_all{2,i}: Layer thickness vector [m]
 store_all{3,i}: Compressional wave velocity vector [m/s] 
 store_all{4,i}: Rayleigh wave velocity vector [m/s]
 store_all{5,i}: Wavelength vector [m]
 store_all{6,i}: Dispersion misfit value [%]
elapsedTime Elapsed time (MC simulations)
store_accepted If up_low_boundary = 'yes'
- Sampled profiles whose DC are within the upper/lower boundaries of the 
 experimental data (cell array)
 For profile no. j
 store_accepted{1,j}: Shear wave velocity vector [m/s]
 store_accepted{2,j}: Layer thickness vector [m]
 store_accepted{3,j}: Compressional wave velocity vector [m/s]
 store_accepted{4,j}: Rayleigh wave velocity vector [m/s]
 store_accepted{5,j}: Wavelength vector [m]
 store_accepted{6,j}: Dispersion misfit value [%]
If up_low_boundary = 'no'
- store_accepted is returned as NaN
10
3.2 Computing theoretical dispersion curves 
(MASWaves_theoretical_dispersion_curve_FDMA) 
The function MASWaves_theoretical_dispersion_curve computes the theoretical fundamental mode 
dispersion curve for the stratified soil model defined by n, alpha, beta, rho and h at wavelengths lambda.
Subfunctions:
MASWaves_FDMA
Input arguments
c_min Minimum testing Rayleigh wave phase velocity vector [m/s]
c_max Maximum testing Rayleigh wave phase velocity vector [m/s]
c_step Testing Rayleigh wave phase velocity increment [m/s]
delta_c Zero search initiation parameter [m/s] 
At wave number k_i the zero search is initiated at a phase velocity of 
max{c_(i-1)-delta_c , c_min}, where c_(i-1) is the theoretical Rayleigh wave 
phase velocity value at wave number k_(i-1)
lambda Wavelength vector [m]
n Number of finite thickness layers
alpha Compressional wave velocity vector [m/s] (array of length n+1)
beta Shear wave velocity vector [m/s] (array of length n+1)
rho Mass density vector [kg/m3
] (array of length n+1)
h Layer thickness vector [m] (array of length n)
Output arguments
c_t Rayleigh wave phase velocity vector (theoretical dispersion curve) [m/s]
lambda_t Rayleigh wave wavelength vector (theoretical dispersion curve) [m]
3.2.1 MASWaves_FDMA
The function MASWaves_FDMA computes the value of the Rayleigh wave dispersion function F for the 
ordered couple (c,k). Computations are conducted using the fast delta matrix algorithm (FDMA).
3 The 
stratified soil model is described in terms of shear wave velocity, compressional wave velocity, mass 
density and layer thicknesses. 
Subfunctions:
MASWaves_FDMA_recursion
3 Buchen, P.W. & Ben-Hador, R. (1996). Free-mode surface-wave computations. Geophysical Journal International,
124(3), 869–887. doi:10.1111/j.1365-246X.1996.tb05642.x.
11
Input arguments
c Rayleigh wave phase velocity [m/s]
k Wave number
n Number of finite thickness layers
alpha Compressional wave velocity vector [m/s] (array of length n+1)
beta Shear wave velocity vector [m/s] (array of length n+1)
rho Mass density vector [kg/m3
] (array of length n+1)
h Layer thickness vector [m] (array of length n)
Output arguments
F F(c,k). Dispersion function value for the ordered couple (c,k).
3.2.2 MASWaves_FDMA_recursion
The function MASWaves_FDMA_recursion conductes the layer recursion of the fast delta matrix 
algorithm.
Subfunctions:
(none)
Input arguments
c Rayleigh wave phase velocity [m/s]
k Wave number
alpha Compressional wave velocity of layer i [m/s]
beta1 Shear wave velocity of layer i [m/s]
beta2 Shear wave velocity of layer (i+1) [m/s]
rho1 Mass density of layer i [kg/m3
]
rho2 Mass density of layer (i+1) [kg/m3
]
h Thickness of layer i [m]
X Recursion vector
Output arguments
X Recursion vector
12
3.3 Evaluating the misfit between theoretical/experimental curves 
(MASWaves_misfit_MC) 
The function MASWaves_misfit_MC is used to evaluate the misfit between the theoretical and the 
experimental fundamental mode dispersion curves. The theoretical and experimental curves must be 
evaluated at the same wavelengths lambda. 
Subfunctions:
(none)
Input arguments
c_t Rayleigh wave phase velocity vector (theoretical dispersion curve) [m/s]
c_OBS Rayleigh wave phase velocity vector (experimental/observed dispersion 
curve) [m/s]
Output arguments
e Dispersion misfit [%]
3.4 Visualizing the inversion results (MASWaves_inversion_MC_plot) 
The function MASWaves_inversion_MC_plot displays the inversion results. The simulated profiles are 
visualized as
1) Sampled Vs profiles/dispersion curves presented using a color scale defined based on dispersion 
misfit values
2) [Optional] Accepted Vs profiles/dispersion curves presented using a color scale defined based on 
dispersion misfit values. Remaining trial profiles shown in gray.
Subfunctions:
(None)
Input arguments
Experimental fundamental mode dispersion curve
c_OBS Phase velocity [m/s]
lambda_OBS Wavelength [m]
up_low_boundary - 'yes': Upper/lower boundaries for the experimental dispersion curve are 
available.
- 'no': Upper/lower boundaries for the experimental dispersion curve are 
not available.
c_OBS_up Phase velocity, upper bound curve [m/s]
(can be assigned as 'nan' or [] if up_low_boundary = 'no').
13
c_OBS_low Phase velocity, lower bound curve [m/s]
(can be assigned as 'nan' or [] if up_low_boundary = 'no').
Inversion results
n Number of finite thickness layers
store_all All sampled profiles (cell array)
 For iteration no. i
 store_all{1,i}: Shear wave velocity vector [m/s]
 store_all{2,i}: Layer thickness vector [m]
 store_all{3,i}: Compressional wave velocity vector [m/s]
 store_all{4,i}: Rayleigh wave velocity vector [m/s]
 store_all{5,i}: Wavelength vector [m]
 store_all{6,i}: Dispersion misfit value [%]
store_accepted Required if up_low_boundary = 'yes'
- Sampled profiles whose DC are within the upper/lower boundaries of the
 experimental data (cell array)
 For profile no. j
 store_accepted{1,j}: Shear wave velocity vector [m/s]
 store_accepted{2,j}: Layer thickness vector [m]
 store_accepted{3,j}: Compressional wave velocity vector [m/s]
 store_accepted{4,j}: Rayleigh wave velocity vector [m/s]
 store_accepted{5,j}: Wavelength vector [m]
 store_accepted{6,j}: Dispersion misfit value [%]
(can be assigned as 'nan' or [] if up_low_boundary = 'no').
MaxDepth Maximum depth for simulated shear wave velocity profiles [m]
