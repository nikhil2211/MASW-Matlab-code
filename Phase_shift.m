
%%
% Note: 
% This script is written in cell mode, i.e. the code is divided into
% several code sections. Each code section begins with two comment characters (%%).
% Below are instructions on evaluating code sections.
%
% - Run the code in the current section
% (1) Place the cursor in the code section. 
% (2) On the Editor tab, click Run Section OR hit Crtl and Return on the
% keyboard
% 
% - Run the code in the current section, and then move to the next section.
% (1) Place the cursor in the code section. 
% (2) On the Editor tab, click Run and Advance
%
%%
clear all
close all
clc
%% ------------- Dispersion analysis ---------------
%load finaldata.mat
%u = Data;
%u = fliplr(u);%fliping
load muted_data.mat
u= flipud(Muteddata); % my x-t domain signals defined for this script
 
fs = 7500; % Hz
N = 12;
x1 = 5; % m
dx = 0.5; % m
T = 0.000125:0.000125:1;
Tmax = max(T);
L = (N-1)*dx;
x = x1:dx:(L+x1);

%%
du = 1/50; 
FigWidth = 6; % cm
FigHeight = 8; % cm
FigFontSize = 8; % pt

figure
MASWaves_plot_data(u,N,dx,x1,L,T,Tmax,du,FigWidth,FigHeight,FigFontSize)

%%
cT_min = 1; % m/s
cT_max = 1500; % m/s
delta_cT = 1; % m/s

[f,c,A] = MASWaves_dispersion_imaging(u,N,x,fs,cT_min,cT_max,delta_cT);

%%
resolution = 100;
fmin = 1; % Hz
fmax = 80; % Hz
FigWidth = 7; % cm
FigHeight = 7; % cm
FigFontSize = 8; % pt
figure
[fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_2D(f,c,A,fmin,fmax,...
    resolution,FigWidth,FigHeight,FigFontSize);

%%
fmin = 1; % Hz
fmax = 80; % Hz
FigWidth = 10; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
[fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_3D(f,c,A,fmin,fmax,...
    FigWidth,FigHeight,FigFontSize);

%%
f_receivers = 4.5; % Hz
select = 'numbers';
up_low_boundary = 'yes'; 
p = 99; % Percentage
[f_curve0,c_curve0,lambda_curve0,...
    f_curve0_up,c_curve0_up,lambda_curve0_up,...
    f_curve0_low,c_curve0_low,lambda_curve0_low] = ...
    MASWaves_extract_dispersion_curve(f,c,A,fmin,fmax,f_receivers,...
    select,up_low_boundary,p);

%%
FigWidth = 9; % cm
FigHeight = 6; % cm
FigFontSize = 8; % pt
type = 'f_c';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
     f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,...
     lambda_curve0_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)

FigWidth = 7; % cm
FigHeight = 9; % cm
FigFontSize = 8; % pt
type = 'c_lambda';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
     f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,...
     lambda_curve0_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)
 
 %% Inversion analysis%%
 
%%
% Repeated use of MASWaves_theoretical_dispersion_curve.m, MASWaves_misfit.m 
% and MASWaves_plot_theor_exp_dispersion_curves.m
% (For iteration, the layer parameters should be updated and this code section run again). 

c_test_min = 0; % m/s
c_test_max = 800; % m/+s
delta_c_test = 1; % m/s
c_test = c_test_min:delta_c_test:c_test_max; % m/s

% Layer parameters
n = 10;
alpha = [1440 1440 1440 1440 1440 1440 1440 1440 1440 1440 1440]; % m/s
h = [1 1 2 2 4 5 6 7 8 9 10 Inf]; % m
beta = [75 90 150 180 240 290 290 300 320 330 340]; % m/s
rho = [1850 1850 1850 1850 1850 1850 1850 1850 1850 1850 1850]; % kg/m^3

up_low_boundary = 'yes';
[c_t,lambda_t] = MASWaves_theoretical_dispersion_curve...
    (c_test,lambda_curve0,h,alpha,beta,rho,n);

up_low_boundary = 'yes';
FigWidth = 8; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
    c_curve0_low,lambda_curve0_low,up_low_boundary,...
    FigWidth,FigHeight,FigFontSize)

e = MASWaves_misfit(c_t,c_curve0);

%%
% Use of MASWaves_inversion

c_test_min = 0; % m/s
c_test_max = 800; % m/+s
delta_c_test = 1; % m/s
c_test = c_test_min:delta_c_test:c_test_max; % m/s

% Layer parameters
n = 10;
alpha = [1440 1440 1440 1440 1440 1440 1440 1440 1440 1440 1440]; % m/s
h = [1 1 2 2 4 5 6 7 8 9 10 Inf]; % m
beta = [75 90 150 180 240 290 290 300 320 330 340]; % m/s
rho = [1850 1850 1850 1850 1850 1850 1850 1850 1850 1850 1850]; % kg/m^3

up_low_boundary = 'yes';
[c_t,lambda_t,e] = MASWaves_inversion(c_test,h,alpha,beta,rho,n,...
    up_low_boundary,c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
    c_curve0_low,lambda_curve0_low);

% View the results
up_low_boundary = 'yes';
FigWidth = 8; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
    c_curve0_low,lambda_curve0_low,up_low_boundary,...
    FigWidth,FigHeight,FigFontSize)


