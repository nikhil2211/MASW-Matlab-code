
%%
%  [c_t,lambda_t] = MASWaves_theoretical_dispersion_curve_FDMA...
%    (c_min,c_max,c_step,lambda,n,alpha,beta,rho,h,delta_c)
%%
%  The function MASWaves_theoretical_dispersion_curve computes the
%  theoretical fundamental mode dispersion curve for the stratified
%  soil model defined by n, alpha, beta, rho and h at wavelengths lambda.
%
%% Input
%  c_min         Minimum testing Rayleigh wave phase velocity [m/s]
%  c_max         Maximum testing Rayleigh wave phase velocity [m/s]
%  c_step        Testing Rayleigh wave phase velocity increment [m/s]
%  lambda        Wavelength vector [m]
%  n             Number of finite thickness layers
%  alpha         Compressional wave velocity vector [m/s] (array of length n+1)
%  beta          Shear wave velocity vector [m/s] (array of length n+1)
%  rho           Mass density vector [kg/m^3] (array of length n+1)
%  h             Layer thickness vector [m] (array of length n)
%  delta_c       Zero search initiation parameter [m/s] 
%                At wave number k_i the zero search is initiated at  
%                a phase velocity of max{c_(i-1)-delta_c , c_min}, where 
%                c_(i-1) is the theoretical Rayleigh wave phase velocity 
%                value at wave number k_(i-1)
%
%% Output
%  c_t           Rayleigh wave phase velocity vector (theoretical
%                dispersion curve) [m/s] 
%  lambda_t      Rayleigh wave wavelength vector (theoretical dispersion
%                curve) [m]
%
%% Subfunctions
%  MASWaves_FDMA
%

%%
function [c_t,lambda_t] = MASWaves_theoretical_dispersion_curve_FDMA...
    (c_min,c_max,c_step,lambda,n,alpha,beta,rho,h,delta_c)

% Determine testing phase velocity values
c_test = c_min:c_step:c_max;

% Wave numbers that correspond to wavelengths lambda
k = (2*pi)./lambda;

% Number of modes (modes = 1, fundamental mode)
modes = 1;

% Initialization
D = zeros(length(c_test),length(k));
c_t = zeros(length(k),modes);
lambda_t = zeros(length(k),modes);

% For each wave number k, compute the dispersion function using 
% increasing values of c_test until its value has a sign change.
m_loc = 1;
delta_m = round(delta_c/c_step);
for j = 1:length(k)
    for m = m_loc:length(c_test)
        D(j,m) = MASWaves_FDMA(c_test(m),k(j),n,alpha,beta,rho,h);
        if m==m_loc
            sign_old = sign(MASWaves_FDMA(c_test(1),k(j),n,alpha,beta,rho,h));
        else
            sign_old = signD;
        end
        signD = sign(D(j,m));
        if sign_old*signD == -1
            c_t(j) = c_test(m);
            lambda_t(j) = 2*pi/k(j);
            m_loc = m - delta_m;
            if m_loc <= 0
                m_loc = 1;
            end
            break
        end
    end
end
end