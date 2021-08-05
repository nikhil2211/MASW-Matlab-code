
%%
%  X = MASWaves_FDMA_recursion(c,k,alpha,beta1,beta2,rho1,rho2,h,X)
%
%%
%  The function MASWaves_FDMA_recursion conductes the layer recursion
%  of the fast delta matrix algorithm.
%
%% Input
%  c             Rayleigh wave phase velocity [m/s]
%  k             Wave number
%  alpha         Compressional wave velocity of layer i [m/s]
%  beta1         Shear wave velocity of layer i [m/s]
%  beta2         Shear wave velocity of layer (i+1) [m/s]
%  rho1          Mass density of layer i [kg/m^3]
%  rho2          Mass density of layer (i+1) [kg/m^3]
%  h             Thickness of layer i [m]
%  X             Recursion array
%
%% Output
%  X             Recursion array
%
%% Subfunctions
%  (None)
%

%%
function X = MASWaves_FDMA_recursion(c,k,alpha,beta1,beta2,rho1,rho2,h,X)

% Layer parameters
gamma1 = beta1^2/c^2;
gamma2 = beta2^2/c^2;
epsilon = rho2/rho1;
eta = 2*(gamma1-epsilon*gamma2);

a = epsilon+eta;
ak = a-1;
b = 1-eta;
bk = b-1;

% Layer eigenfunctions
r = sqrt(1-(c^2/alpha^2));
C_alpha = cosh(k*r*h);
S_alpha = sinh(k*r*h);

s = sqrt(1-(c^2/beta1^2));
C_beta = cosh(k*s*h);
S_beta = sinh(k*s*h);

% Layer recursion
p1 = C_beta*X(2)+s*S_beta*X(3);
p2 = C_beta*X(4)+s*S_beta*X(5);
p3 = (1/s)*S_beta*X(2)+C_beta*X(3);
p4 = (1/s)*S_beta*X(4)+C_beta*X(5);

q1 = C_alpha*p1-r*S_alpha*p2;
q2 = -(1/r)*S_alpha*p3+C_alpha*p4;
q3 = C_alpha*p3-r*S_alpha*p4;
q4 = -(1/r)*S_alpha*p1+C_alpha*p2;

y1 = ak*X(1)+a*q1;
y2 = a*X(1)+ak*q2;

z1 = b*X(1)+bk*q1;
z2 = bk*X(1)+b*q2;

X(1) = bk*y1+b*y2;
X(2) = a*y1+ak*y2;
X(3) = epsilon*q3;
X(4) = epsilon*q4;
X(5) = bk*z1+b*z2;

end

