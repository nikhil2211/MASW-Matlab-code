
%%
%  [store_all,elapsedTime,store_accepted] = MASWaves_inversion_MC...
%    (c_min,c_max,c_step,delta_c,n,n_unsat,alpha_initial,nu_unsat,...
%    beta_initial,rho,h_initial,N_reversals,c_OBS,lambda_OBS,...
%    up_low_boundary,c_OBS_up,c_OBS_low,b_S,b_h,N_max,e_max)
%
%%
%  The function MASWaves_inversion_MC carries out the inversion analysis
%  using a Monte Carlo based search procedure.
%  (1) The program computes the theoretical (fundamental mode) dispersion
%      curve for the stratified layer model defined by n, alpha_initial,
%      beta_initial, rho and h_initial. The initial value of the dispersion
%      misfit function is evaluated.
%  (2) A Monte Carlo based search process is used in search of the shear
%      wave velocity profile (i.e., value of beta and h for each layer)
%      that provides the closest fit to the experimental data.
%  (3) Inversion results are provided in the form of
%      - All sampled profiles
%      - All sampled profiles whose theoretical dispersion curves fall
%        within the upper/lower boundaries of the experimental curve
%        ('accepted profiles') [optional]
%
%% Input
%  c_min             Minimum testing Rayleigh wave phase velocity [m/s]
%  c_max             Maximum testing Rayleigh wave phase velocity [m/s]
%  c_step            Testing Rayleigh wave phase velocity increment [m/s]
%  delta_c           Zero search initiation parameter [m/s] 
%                    At wave number k_i the zero search is initiated at  
%                    a phase velocity of max{c_(i-1)-delta_c , c_min}, where 
%                    c_(i-1) is the theoretical Rayleigh wave phase velocity 
%                    value at wave number k_(i-1)
%  n                 Number of finite thickness layers
%  n_unsat           Number of unsaturated soil layers
%                    (n_unsat = 0 for a fully saturated soil profile)
%                    (n_unsat = n+1 for a fully unsaturated soil profile)
%  alpha_initial     Initial estimate of compressional wave velocity [m/s]
%                    (array of length n+1)
%  nu_unsat          Poisson's ratio
%  beta_initial      Initial estimate of shear wave velocity [m/s]
%                    (array of length n+1)
%  rho               Mass density vector [kg/m^3] (array of length n+1)
%  h_initial         Initial estimate of layer thicknesses [m]
%                    (array of length n)
%  N_reversals       Number of soil layers where velocity reversals are
%                    permitted (N_reversals = 0 for a normally dispersive
%                    profile).
%
%  Experimental fundamental mode dispersion curve
%  c_OBS             Phase velocity [m/s]
%  lambda_OBS        Wavelength [m]
%  up_low_boundary   - 'yes'   Upper/lower boundaries for the experimental
%                              dispersion curve are available.
%                    - 'no'    Upper/lower boundaries for the experimental
%                              dispersion curve are not available.
%  c_OBS_up          Phase velocity, upper bound curve [m/s]
%                    (can be assigned as 'nan' or [] if up_low_boundary = 'no').
%  c_OBS_low         Phase velocity, lower bound curve [m/s]
%                    (can be assigned as 'nan' or [] if up_low_boundary = 'no').
%
%  Search-control parameters
%  b_S               Shear wave velocity search-control parameter
%  b_h               Layer thickness search-control parameter
%  N_max             Maximum number of iterations
%  e_max             Maximum misfit (optional stopping criterion for MC search)
%
%% Output
%  store_all         All sampled profiles (cell array)
%                      For iteration no. i
%                      store_all{1,i}: Shear wave velocity vector [m/s]
%                      store_all{2,i}: Layer thickness vector [m]
%                      store_all{3,i}: Compressional wave velocity vector [m/s]                   
%                      store_all{4,i}: Rayleigh wave velocity vector [m/s]
%                      store_all{5,i}: Wavelength vector [m]
%                      store_all{6,i}: Dispersion misfit value [%]
%  elapsedTime       Elapsed time (MC simulations)
%  store_accepted    If up_low_boundary = 'yes'
%                    - Sampled profiles whose DC are within the upper/lower
%                      boundaries of the experimental data (cell array)
%                        For profile no. j
%                        store_accepted{1,j}: Shear wave velocity vector [m/s]
%                        store_accepted{2,j}: Layer thickness vector [m]
%                        store_accepted{3,j}: Compressional wave velocity vector [m/s]
%                        store_accepted{4,j}: Rayleigh wave velocity vector [m/s]
%                        store_accepted{5,j}: Wavelength vector [m]
%                        store_accepted{6,j}: Dispersion misfit value [%]
%                    If up_low_boundary = 'no'
%                    - store_accepted is returned as NaN
%
%% Subfunctions
%  MASWaves_theoretical_dispersion_curve_FDMA
%  MASWaves_misfit_MC
%

%%
function [store_all,elapsedTime,store_accepted] = MASWaves_inversion_MC...
    (c_min,c_max,c_step,delta_c,n,n_unsat,alpha_initial,nu_unsat,beta_initial,rho,h_initial,...
    N_reversals,c_OBS,lambda_OBS,up_low_boundary,c_OBS_up,c_OBS_low,...
    b_S,b_h,N_max,e_max)

% Theoretical dispersion curve and dispersion misfit valye for the initial
% set of model parameters
[c_t,lambda_t] = MASWaves_theoretical_dispersion_curve_FDMA...
    (c_min,c_max,c_step,lambda_OBS,n,alpha_initial,beta_initial,rho,h_initial,delta_c);
e_initial = MASWaves_misfit_MC(c_t,c_OBS);

timerVal = tic; % Start stopwatch timer

% Initialization
store_all = cell(6,N_max);
beta_opt = beta_initial;
h_opt = h_initial;
e_opt = e_initial;
Nend = N_max;

% Iteration
for w = 1:N_max
    [w] % Number of iteration
    
    % Testing shear wave velocity vector
    beta_test = beta_opt + unifrnd(-(b_S/100).*beta_opt,(b_S/100).*beta_opt);
    while issorted(beta_test((N_reversals+1):end)) == 0
        beta_test = beta_opt + unifrnd(-(b_S/100).*beta_opt,(b_S/100).*beta_opt);
    end
    
    % Reevaluate testing compressional wave velocity vector
    alpha_unsat = sqrt((2*(1-nu_unsat))/(1-2*nu_unsat))*beta_test;
    alpha_test = 1440*ones(size(beta_test));
    if n_unsat ~= 0
        alpha_test(1:n_unsat) = alpha_unsat(1:n_unsat);
    end
    
    % Testing layer thickness vector
    h_test = h_opt + unifrnd(-(b_h/100).*h_opt,(b_h/100).*h_opt);
    [c_t,lambda_t] = MASWaves_theoretical_dispersion_curve_FDMA...
        (c_min,c_max,c_step,lambda_OBS,n,alpha_test,beta_test,rho,h_test,delta_c);
    e_test = MASWaves_misfit_MC(c_t,c_OBS);
    
    store_all{1,w} = beta_test;
    store_all{2,w} = h_test;
    store_all{3,w} = alpha_test;
    store_all{4,w} = c_t;
    store_all{5,w} = lambda_t;
    store_all{6,w} = e_test;
    
    if e_test <= e_opt
        e_opt = e_test;
        beta_opt = beta_test;
        h_opt = h_test;
    end
    if e_opt < e_max
        Nend = w;
        break
    end
end
elapsedTime = toc(timerVal);  % Stop stopwatch timer

% Find sampled profiles whose theoretical dispersion curves fall within the
% upper and lower boundaries of the experimental data
if strcmp(up_low_boundary,'yes') == 1
    Vec = zeros(Nend,1);
    count = 0;
    for i = 1:Nend
        temp = zeros(length(c_OBS_low),1);
        for j = 1:length(c_OBS_low)
            temp(j) = (c_OBS_low(j) < store_all{4,i}(j) && store_all{4,i}(j) < c_OBS_up(j));
        end
        if sum(temp) == length(c_OBS_low)
            Vec(i) = 1;
            count = count + 1;
        end
    end
    
    store_accepted = cell(6,count);
    j = 0;
    for i = 1:Nend
        if Vec(i) == 1
            j = j+1;
            store_accepted{1,j} = store_all{1,i};
            store_accepted{2,j} = store_all{2,i};
            store_accepted{3,j} = store_all{3,i};
            store_accepted{4,j} = store_all{4,i};
            store_accepted{5,j} = store_all{5,i};
            store_accepted{6,j} = store_all{6,i};
        end
    end
else
    store_accepted = NaN;
end

end