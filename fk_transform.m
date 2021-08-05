
clear all
close all

%% Edit the parameters below that represent your dataset:

load finaldata.mat % Load the synthetic dataset
U=Data;% my x-t domain signals defined for this script
U= fliplr(U); % reversing my matrix (specifically for my data ) 
pick='none'; % This sets the function for picking the dispersion curve. 
             % Pick setting options are 'manual' , 'auto' or 'none'
waterlevel=0; % Dispersion Image waterlevel expressed as a percent for display/picking 
              % purposes. Amplitudes in the dispersion image that are less
              % than "waterlevel" percent of the peak amplitudes will be 
              % set to NANs. E.g., a waterlevel=0 will allow all Dispersion
              % Image data to be displayed, while a waterlevel=95 will
              % only allow the top 5 percent of amplitudes to be
              % displayed. High values may be desired, especially when
              % 'manual' picking is set. 
            
ct=80:1:1500; % A vector for the velocity of interest in m/s.
               % This is the velocity "scanning" range
               
freqlimits=[1 300]; % For plotting purposes only. Sets the frequency (Hz)
                    % plotting range for the dispersion images. Enter ['none']
% Input Signal Properties %
dt = 0.000125; %sample rate[s]
n_samples = 8000; %no. of samples per receiver
%fs=15000; % sampling frequency Hz
numk = 5000;%Number of velocities or wavenumbers to consider in calculations
min_x = 5; % min offset [m] of reciever spread
max_x = 10.5; % max offset [m] of reciever spread
dx = 0.5; % distance between receivers [m]
multiple = 1; %Used when zeros are padded such that down-sampling every 
              %"multiple" samples is necessary
min_frequency = freqlimits(1,1);
max_frequency = freqlimits(1,2);

fnyq = 1/(2*dt);%Nyquist frequency
df = 1/(n_samples*dt);%Sampling in frequency domain
freq = 0:df:fnyq-1;
freq = freq; %
kres = 2*pi / dx ;%Maximum resolvable wavenumber
%space
dk = 2*pi/(numk*dx);%2*pi/(numk*dx)
k_vals = dk:dk:kres;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% End Inputs and Run Script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Perform two-dimensional FFT...............................................
fk = fft2(U,n_samples,numk)*dx*dt;
fk = fliplr(fk);
%droping column;
%[M1,m1] = size(fk);
%fk = fk(:,1:m1);
%Remove frequencies above/below specificied max/min
%frequencies and downsample (if required by zero padding)
[m,fminID] = min(abs(freq-min_frequency));
[M,fmaxID] = min(abs(freq-max_frequency));
%freq_id = fminID:multiple:fmaxID;
freq = freq(:,fminID:multiple:fmaxID);
fk = fk(fminID:multiple:fmaxID,: );
%% 
%Identify wavenumber associated with maximum in fk domain..................
    
pnorm = zeros(size(fk));
k_peak = zeros(size(freq));
for k = 1:size(fk,1)
        % Normalize by largest number in fk domain at each frequency
        pnorm(k,:) = abs(fk(k,:)) / max(abs( fk(k,:)));
        %# Find peak
       [Value,pk_id] = max(pnorm(k,:));
        k_peak(1,k) = k_vals(1,pk_id);
end

pnorm = pnorm';%transposing matrix
pnorm = pnorm./max(pnorm);
%% Plotting parameters considered in analysis  
% Grid parameters
[freq_grid,wavenum_grid] = meshgrid(freq,k_vals);
vel_grid = 2*pi*(freq_grid./wavenum_grid);
wavel_grid = 2*pi./wavenum_grid;

% Peaks
wavel_peak = 2*pi ./k_peak;
v_peak = wavel_peak.*freq;

%% plot

%Frequency-wavenumber
xgrid = freq_grid;
ygrid = wavenum_grid;

% subplot(1,2,2)
h=pcolor(xgrid,ygrid,pnorm); % plot the normalized dispersion image
set(h, 'EdgeColor', 'none');
colormap(jet)
cb = colorbar;
ylabel('wavenumber')
xlabel('Freq (Hz)')
title('Dispersion Image-fk method')
set(gca,'FontSize',25,'fontweight','bold')



