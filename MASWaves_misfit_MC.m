
%%
%  e = MASWaves_misfit_MC(c_t,c_OBS)
%
%%
%  The function MASWaves_misfit_MC is used to evaluate the misfit between  
%  the theoretical and the experimental fundamental mode dispersion curves. 
%  The theoretical and experimental curves must be evaluated at the same 
%  wavelengths lambda. 
%
%% Input
%
%  Theoretical fundamental mode dispersion curve
%  c_t               Rayleigh wave phase velocity vector [m/s]
%
%  Experimental (observed) fundamental mode dispersion curve
%  c_OBS             Rayleigh wave phase velocity vector [m/s]
%
%% Output
%  e                 Dispersion misfit [%]
%
%% Subfunctions
%  (None)
%

%%
function  e = MASWaves_misfit_MC(c_t,c_OBS)

% Number of data points in theoretical and experimental dispersion curves.
Q = length(c_t); 

% Dispersion misfit
e = (1/Q)*sum(sqrt((c_OBS-c_t).^2)./c_OBS)*100;

end
