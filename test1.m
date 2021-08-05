

clear all
close all

%% Edit the parameters below that represent your dataset:

load finaldata.mat % Load the synthetic dataset
U=Data;% my x-t domain signals defined for this script
U = fliplr(U);
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
            
ct=1:1:1000; % A vector for the velocity of interest in m/s.
               % This is the velocity "scanning" range
               
freqlimits=[0 50]; % For plotting purposes only. Sets the frequency (Hz)
                    % plotting range for the dispersion images. Enter ['none']
                    % if you desire the natural limits based on time sampling

% Input Signal Properties %
fs=1500; % sampling frequency Hz
min_x = 5; % min offset [m] of reciever spread
max_x = 10.5; % max offset [m] of reciever spread
d_x = 0.5; % distance between receivers [m]

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% End Inputs and Run Script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = min_x:d_x:max_x; % spatial sampling vector (receiver position)
%t = 0.000125:0.000125:1;
t=1/fs:1/fs:length(U(:,1))/fs;  % time vector in seconds

%%%%% Find the dispersion image %%%%%%
%% 1. First, find the fft of each trace
% determine the proper frequency vector
L=length(t); % number os samples
T=max(t); % max time
if ~mod(L,2)
    f = (1/T)*(0:L/2-1); % n is even
else
    f = (1/T)*(0:(L-1)/2); % n is odd
end

Ufft=fft(U); % fft is performed on each column (the time vector)
Ufft=Ufft(1:length(f),:); % only keep the fft results for the positve frequencies
w=f*2*pi; % convert the freq vector to rads/sec

%% 2. Normalize the spectrum
Rnorm=Ufft./abs(Ufft); % Normalize the spectrum

%% 3. Perfrom the summation over N in the frequency domain
for ii=1:length(w) % send in frequencies one at a time
    As=zeros(length(x),length(ct)); %initialize As or overwite the previous As
    for n=1:length(x) % send in a positon (x) one at a time
        % Iterate As over position
        As(n,:)=exp(1j*w(ii)*(x(n)./ct))*Rnorm(ii,n);
        % Note: The sign convention of MATLAB FFT is negative, so the minus
        % sign is dropped from the "-j" seen in Ryden 2004
    end
     AsSum(ii,:)=abs(sum(As)); % Sum up As for a given frequency. 
end                              % The final result will be a matrix of 
                                 % amplitudes changing with phase 
                                 % velocity in the x direction and 
                                 % frequency in the y direction

AsSum=AsSum'; % transpose the matrix so velocity is on vertical and freq is on horizontal
normed=AsSum./max(AsSum);% normalize the dispersion image so max column values=1

%%
%%%% This next section will pick the dispersion curve base on the "pick"
%%%%%%% setting ('auto' 'manual' or 'none')

%%%%%%%%%%%% curve autopicking %%%%%%%%%%%%%%%%%%%
% This section will auto pick the dispersion curve based on the peak 
% value in each column. This will need updating when 
% higher modes are introduced
if strcmp('auto',pick) == 1 % if auto picking is set
    remaindat=normed; % create a new variable called remaindat from the normalized image
    remaindat(remaindat < 1) = 0;% set all data amplitudes in the new variable that are less than 1 to zero
    [row,col]=find(remaindat); % find the indices where there is non-zero data (the peak that we set equal to 1)
    autovel=ct(row); % autopicked velocity indices that contain the peak value
    autofreq=f(col); % autopicked frequency indices that contain the peak value
    % get rid of all the extra picks at zero frequency
    nz=numel(autofreq)-nnz(autofreq); % determine the number of repeated picks at zero
    autofreq=autofreq(nz:end); % remove the repeated zero picks
    autovel=autovel(nz:end); % remove the repeated zero picks
    DispersionVelocity=autovel; % need to get the rlowess working again, I no longer have access to the toolbox that allows it
    normed(normed < waterlevel/100)=nan; % throw out data less than waterlevel
%%%%%%%%%%%%%%%% end autopicking %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% manual curve picking %%%%%%%%%%%%%
elseif strcmp('manual',pick) == 1 % if manual picking is set
    normed(normed < waterlevel/100)=nan; % throw out data less than waterlevel
    figure('units','normalized','outerposition',[0 0 1 1])
    h=pcolor(f,ct,normed);
    set(h, 'EdgeColor', 'none');
    colormap([jet;1 1 1]);
    ylabel('Phase Velocity (m/s)')
    xlabel('Freq (Hz)')
        if strcmp('none',freqlimits) == 1

        else
            xlim(freqlimits)
        end
    title({'1) Use your mouse to pick the dispersion curve using as many points as possible.';...
    '2) Press ENTER when finished ';'3) MATLAB will interpolate between points';'DUPLICATE PICKS WILL RESULT IN AN ERROR.'})
    set(gca,'FontSize',15,'fontweight','bold')
    [pick_f,pick_ct]=ginput;
    DispersionVelocity=interp1(pick_f,pick_ct,f,'spline');
    hold on
    plot(pick_f,pick_ct,'ko','markersize',10,'linewidth',5)
    plot(f,DispersionVelocity,'k','linewidth',2)
    title('Press ENTER to continue')
    pause
    close
%%%%%%%%% end manual curve picking %%%%%%%%%%%%

% No curve picking
elseif strcmp('none',pick) == 1 % if no picking is desired i.e., 'none' is set
    DispersionVelocity=nan(length(f));
    normed(normed < waterlevel/100)=nan; % throw out data less than waterlevel
end

DispersionCurve=[f' DispersionVelocity']; % this only has meaning if the curve was picked using 'manual' or 'auto'

%%
%%% End calculations. Plot everything %%%%
% Plot the normalized x-t data
%cqi_plotmatrix(U,'dt',0.000125,'t0',0,'scale',1,'skip',0,'fillco','b','linewidth',0.01);
%wiggle(U);

% figure;
% subplot(1,2,1)
% datanorm=U./max(U);
% imagesc(x,t,datanorm)
% colormap(jet)
% xlabel('Receiver Position (m)')
% set(gca,'xaxisLocation','top')
% ylabel('time (sec)')
% %title('x-t domain')
% %ylim([0 0.05])
% set(gca,'FontSize',25,'fontweight','bold')

% Plot the dispersion results
subplot(1,1,1)
h=pcolor(f,ct,normed); % plot the normalized dispersion image
set(h, 'EdgeColor', 'none');
colormap(jet)
cb = colorbar;
hold on
plot(DispersionCurve(:,1),DispersionCurve(:,2),'k','linewidth',2) % Plot the picked dispersion curve. This is nothing if 'none' picking was set
if strcmp('none',freqlimits) == 1 % check for plotting limits
    % if 'none' then ignore the xlim
else
    xlim(freqlimits) % set the xlim if specified
end
ylabel('Phase Velocity (m/s)')
xlabel('Freq (Hz)')
title('Dispersion Image-phase shift method')
set(gca,'FontSize',25,'fontweight','bold')

