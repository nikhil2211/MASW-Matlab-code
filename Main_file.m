%%
clear all
close all
clc
%% ------------- Dispersion analysis ---------------
%load finaldata.mat
%u = Data;
%u = fliplr(u);%fliping
load Modsmoothvel-101.mat
%u= flipud(Modsmoothvel); % my x-t domain signals defined for this script
u = Modsmoothvel; 
fs = 8000; % Hz
N = 24;
x1 = 15; % m
dx = 1; % m
T = 0.000125:0.000125:1;
Tmax = max(T);
L = (N-1)*dx;
x = x1:dx:(L+x1);

%%
du = 1/1; 
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
 
%% Load and plot sample data
Data = table([c_curve0],[c_curve0_low],[c_curve0_up],[lambda_curve0],...
                'VariableNames',{'c_OBS' 'c_OBS_low' 'c_OBS_up' 'lambda_OBS'});  
SampleData = sortrows(Data,4);      
% SampleData = readtable('SampleData_MASWavesInversion.csv');
c_OBS = SampleData{1:end,1};
up_low_boundary = 'yes'; % Upper/lower boundary curves are available
c_OBS_low = SampleData{1:end,2};
c_OBS_up = SampleData{1:end,3};
lambda_OBS = SampleData{1:end,4};

figure, hold on
exp_curve = plot(c_OBS,lambda_OBS,'kx-','MarkerFaceColor',...
    'k','MarkerSize', 2,'LineWidth',1);
exp_low = plot(c_OBS_low,lambda_OBS,'r--','LineWidth',1);
exp_up = plot(c_OBS_up,lambda_OBS,'r--','LineWidth',1);

set(gca,'FontSize',9),axis ij, grid on, box off
xlim([0 max(c_OBS_up)]), ylim([0 max(lambda_OBS)+2])
xlabel('Rayleigh wave velocity [m/s]','FontSize',9,'Fontweight','normal')
ylabel('Wavelength [m]','FontSize',9,'Fontweight','normal')
legend([exp_curve,exp_low],'mean','mean \pm std','Location','SouthWest')

set(gcf,'units','centimeters')
pos = [2, 2, 8, 10];
set(gcf,'Position',pos)

%% Initialization

% Testing phase velocity
c_min = 75; % m/s
c_max = 500; % m/s
c_step = 0.9; % m/s
delta_c = 2; % m/s

% Initial values for model parameters
n = 10;
h_initial = [1 1 2 2 4 5 4 6 5 7]; % m (array(n))
beta_initial = [75 90 150 180 240 290 290 300 320 330 340];%array(n+1) % m/s
n_unsat = n+1;
nu_unsat = 0.3;
alpha_temp = sqrt((2*(1-nu_unsat))/(1-2*nu_unsat))*beta_initial; % m/s
alpha_initial = [alpha_temp(1) 1440 1440 1440 1440 1440 1440 ...
    1440 1440 1440 1440]; % m/s array(n+1)
rho = [1850 1850 1850 1850 1900 1900 1900 1900 1950 1950 1950]; % kg/m^3

%% Run inversion analysis (one initation)

% Number of velocity reversals 
N_reversals = 0; % Normally dispersive analysis

% Search-control parameters
b_S = 5; 
b_h = 10; 
N_max = 1000; 
e_max = 0;
tic
% Run inversion
[store_all,elapsedTime,store_accepted] = MASWaves_inversion_MC...
    (c_min,c_max,c_step,delta_c,n,n_unsat,alpha_initial,nu_unsat,...
    beta_initial,rho,h_initial,N_reversals,c_OBS,lambda_OBS,...
    up_low_boundary,c_OBS_up,c_OBS_low,b_S,b_h,N_max,e_max);
toc
%% Display inversion results 
MaxDepth = sum(h_initial);

MASWaves_inversion_MC_plot(c_OBS,lambda_OBS,up_low_boundary,c_OBS_up,...
    c_OBS_low,n,store_all,store_accepted,MaxDepth);

