load finaldata.mat
load muted_data.mat


u1 = Data ;
u1 = fliplr(u1);%fliping
u2 = Muteddata;
u2= flipud(Muteddata);

fs = 7500; % Hz
N = 12;
x1 = 5; % m
dx = 0.5; % m
T = 0.000125:0.000125:1;
Tmax = max(T);
L = (N-1)*dx;
x = x1:dx:(L+x1);


du1 = 1/50;
du2 = 1/50;
FigWidth = 6; % cm
FigHeight = 8; % cm
FigFontSize = 8; % pt

figure
MASWaves_plot_data(u1,N,dx,x1,L,T,Tmax,du1,FigWidth,FigHeight,FigFontSize)
figure
MASWaves_plot_data(u2,N,dx,x1,L,T,Tmax,du2,FigWidth,FigHeight,FigFontSize)
%% 
cT_min = 1; % m/s
cT_max = 1500; % m/s
delta_cT = 1; % m/s

[f1,c1,A1] = MASWaves_dispersion_imaging(u1,N,x,fs,cT_min,cT_max,delta_cT);
[f2,c2,A2] = MASWaves_dispersion_imaging(u2,N,x,fs,cT_min,cT_max,delta_cT);

%% 
resolution = 100;
fmin = 1; % Hz
fmax = 80; % Hz
FigWidth = 7; % cm
FigHeight = 7; % cm
FigFontSize = 8; % pt
figure
[fplot1,cplot1,Aplot1] = MASWaves_plot_dispersion_image_2D(f1,c1,A1,fmin,fmax,...
    resolution,FigWidth,FigHeight,FigFontSize);
figure
[fplot2,cplot2,Aplot2] = MASWaves_plot_dispersion_image_2D(f2,c2,A2,fmin,fmax,...
    resolution,FigWidth,FigHeight,FigFontSize);
%% 
fmin = 1; % Hz
fmax = 80; % Hz
FigWidth = 10; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
[fplot1,cplot1,Aplot1] = MASWaves_plot_dispersion_image_3D(f1,c1,A1,fmin,fmax,...
    FigWidth,FigHeight,FigFontSize);
[fplot2,cplot2,Aplot2] = MASWaves_plot_dispersion_image_3D(f2,c2,A2,fmin,fmax,...
    FigWidth,FigHeight,FigFontSize);
%% 
f_receivers = 4.5; % Hz
select = 'numbers';
up_low_boundary = 'yes'; 
p = 99; % Percentage
[f_curve01,c_curve01,lambda_curve01,...
    f_curve01_up,c_curve01_up,lambda_curve01_up,...
    f_curve01_low,c_curve01_low,lambda_curve01_low] = ...
    MASWaves_extract_dispersion_curve(f1,c1,A1,fmin,fmax,f_receivers,...
    select,up_low_boundary,p);

[f_curve02,c_curve02,lambda_curve02,...
    f_curve02_up,c_curve02_up,lambda_curve02_up,...
    f_curve02_low,c_curve02_low,lambda_curve02_low] = ...
    MASWaves_extract_dispersion_curve(f2,c2,A2,fmin,fmax,f_receivers,...
    select,up_low_boundary,p);
%% 
FigWidth = 9; % cm
FigHeight = 6; % cm
FigFontSize = 8; % pt
type = 'f_c';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve01,c_curve01,lambda_curve01,...
     f_curve01_up,c_curve01_up,lambda_curve01_up,f_curve01_low,c_curve01_low,...
     lambda_curve01_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)
figure
MASWaves_plot_dispersion_curve(f_curve02,c_curve02,lambda_curve02,...
     f_curve02_up,c_curve02_up,lambda_curve02_up,f_curve02_low,c_curve02_low,...
     lambda_curve02_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)

FigWidth = 7; % cm
FigHeight = 9; % cm
FigFontSize = 8; % pt
type = 'c_lambda';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve01,c_curve01,lambda_curve01,...
     f_curve01_up,c_curve01_up,lambda_curve01_up,f_curve01_low,c_curve01_low,...
     lambda_curve01_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)
figure
MASWaves_plot_dispersion_curve(f_curve02,c_curve02,lambda_curve02,...
     f_curve02_up,c_curve02_up,lambda_curve02_up,f_curve02_low,c_curve02_low,...
     lambda_curve02_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)



