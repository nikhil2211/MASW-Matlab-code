
%%
%  MASWaves_inversion_MC_plot(c_OBS,lambda_OBS,up_low_boundary,c_OBS_up,...
%      c_OBS_low,n,store_all,store_accepted,MaxDepth)
%
%%
%  The function MASWaves_inversion_MC_plot displays the inversion
%  results.
%  The simulated profiles are visualized as
%  (1) Sampled Vs profiles/dispersion curves presented using a color scale
%      defined based on dispersion misfit values
%  (2) [Optional] 'Accepted' Vs profiles/dispersion curves presented using 
%      a color scale defined based on dispersion misfit values. Remaining 
%      trial profiles shown in gray.
%
%% Input
%  Experimental fundamental mode dispersion curve
%  c_OBS             Phase velocity [m/s]
%  lambda_OBS        Wavelength [m]
%  up_low_boundary   - 'yes'   Upper/lower boundaries for the experimental
%                              dispersion curve are avaliable.
%                    - 'no'    Upper/lower boundaries for the experimental
%                              dispersion curve are not avaliable.
%  c_OBS_up          Phase velocity, upper bound curve [m/s]
%                    (can be assigned as 'nan' or [] if up_low_boundary = 'no').
%  c_OBS_low         Phase velocity, lower bound curve [m/s]
%                    (can be assigned as 'nan' or [] if up_low_boundary = 'no').
%
%  Inversion results
%  n                 Number of finite thickness layers
%  store_all         All sampled profiles (cell array)
%                      For iteration no. i
%                      store_all{1,i}: Shear wave velocity vector [m/s]
%                      store_all{2,i}: Layer thickness vector [m]
%                      store_all{3,i}: Compressional wave velocity vector [m/s]
%                      store_all{4,i}: Rayleigh wave velocity vector [m/s]
%                      store_all{5,i}: Wavelength vector [m]
%                      store_all{6,i}: Dispersion misfit value [%]
%  store_accepted    Required if up_low_boundary = 'yes'
%                    - Sampled profiles whose DC are within the upper/lower
%                      boundaries of the experimental data (cell array)
%                        For profile no. j
%                        store_accepted{1,j}: Shear wave velocity vector [m/s]
%                        store_accepted{2,j}: Layer thickness vector [m]
%                        store_accepted{3,j}: Compressional wave velocity vector [m/s]
%                        store_accepted{4,j}: Rayleigh wave velocity vector [m/s]
%                        store_accepted{5,j}: Wavelength vector [m]
%                        store_accepted{6,j}: Dispersion misfit value [%]
%                     (can be assigned as 'nan' or [] if up_low_boundary = 'no').
%  MaxDepth          Maximum depth for simulated shear wave velocity
%                    profiles [m]
%
%%  Output
%  (none)
%
%%
function MASWaves_inversion_MC_plot(c_OBS,lambda_OBS,up_low_boundary,...
    c_OBS_up,c_OBS_low,n,store_all,store_accepted,MaxDepth)

%%%% Figure no. 1 %%%%
% Sort sampled profiles by dispersion misfit values
[~,NoPlot_all] = size(store_all);
store_e_all = zeros(1,NoPlot_all);
for i = 1:NoPlot_all
    store_e_all(i) = store_all{6,i};
end
[sort_e_all_plot,order_all_plot]=sort(store_e_all);
sort_e_all_plot = sort_e_all_plot(NoPlot_all:-1:1);
order_all_plot = order_all_plot(NoPlot_all:-1:1);

% Figure settings & colormap
figure
set(gcf,'units','centimeters')
figwidth = 16;
figheight = 10;
pos = [2, 2, figwidth, figheight];
set(gcf,'Position',pos)
vec = [100; 75; 50; 25; 10; 5; 0];
hex =['#FFC77E';'#FFAC6D';'#E7905C';'#BB754B';'#8F5A39';'#5F3C26';'#000000'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
linecolors = interp1(vec,raw,linspace(100,0,NoPlot_all),'pchip');
mymap = colormap(linecolors);

%%% DISPERSION CURVES %%%
ax(1)=subplot(1,2,1);
hold on

% All sampled dispersion curves
for j = 1:1:NoPlot_all
    c_plot = store_all{4,order_all_plot(j)};
    lambda_plot = store_all{5,order_all_plot(j)};
    plot(c_plot,lambda_plot,'color',linecolors(j,:))
end

% Experimental (observed) dispersion curve
plot(c_OBS,lambda_OBS,'r--','LineWidth',1)

% Visual
FigFontSize = 8;
set(gca, 'FontSize', FigFontSize)
axis ij, grid on, box off
set(gca,'TickDir','out')
xlabel('Phase velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
ylabel('Wavelength [m]','FontSize',FigFontSize,'Fontweight','normal')
set(gca,'TickDir','out')

%%% Vs PROFILES %%%
ax(2)=subplot(1,2,2);
hold on
% All sampled Vs profiles
for j = 1:1:NoPlot_all
    clear plot_layer_depth
    clear plot_beta
    h_plot = [0 store_all{2,order_all_plot(j)} max(MaxDepth-sum(store_all{2,order_all_plot(j)}),0)];
    beta_plot = store_all{1,order_all_plot(j)};
    plot_layer_depth = zeros(2*(n+1),1);
    plot_beta = zeros(2*(n+1),1);
    for i = 1:(n+1)
        plot_layer_depth(2*i-1) = sum(h_plot(1:i));
        plot_layer_depth(2*i) = sum(h_plot(1:i+1));
        plot_beta(2*i-1) = beta_plot(i);
        plot_beta(2*i) = beta_plot(i);
    end
    plot(plot_beta, plot_layer_depth,'color',linecolors(j,:));
end

% Visual
FigFontSize = 8;
set(gca, 'FontSize', FigFontSize)
axis ij, grid on, box off
set(gca,'TickDir','out')
xlabel('Shear wave velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
ylabel('Depth [m]','FontSize',FigFontSize,'Fontweight','normal'), ylim([0 MaxDepth])
set(gca,'TickDir','out')

%%% COLORBAR %%%
c = colorbar;
set(c, 'Position', [0.925*.8914  1.5*.11 .0281 0.95*.8150])
for i = 1:2
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [0.9*pos(1) 1.5*pos(2) 0.85*pos(3) 0.95*pos(4)]);
end
c.FontSize = FigFontSize;
c.Label.FontSize = FigFontSize;
c.Label.String = 'Misfit [%]';
c.Ticks = 0:0.25:1;
Y = prctile(sort_e_all_plot,0:25:100);
Y = round(Y,1);
c.TickLabels = Y(end:-1:1);

if strcmp(up_low_boundary,'yes') == 1
    %%%% Figure no. 2 %%%%
    % Sort 'accepted' profiles by dispersion misfit values
    [~,NoPlot] = size(store_accepted);
    store_e = zeros(1,NoPlot);
    for i = 1:NoPlot
        store_e(i) = store_accepted{6,i};
    end
    [sort_e_plot,order_plot]=sort(store_e);
    sort_e_plot = sort_e_plot(NoPlot:-1:1);
    order_plot = order_plot(NoPlot:-1:1);
    
    % Figure settings & colormap
    figure
    set(gcf,'units','centimeters')
    figwidth = 16;
    figheight = 10;
    pos = [2, 2, figwidth, figheight];
    set(gcf,'Position',pos)
    vec = [100; 75; 50; 25; 10; 5; 0];
    hex =['#FFC77E';'#FFAC6D';'#E7905C';'#BB754B';'#8F5A39';'#5F3C26';'#000000'];
    raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
    linecolors = interp1(vec,raw,linspace(100,0,NoPlot),'pchip');
    mymap = colormap(linecolors);
    
    %%% DISPERSION CURVES %%%
    ax(1)=subplot(1,2,1);
    hold on
    
    % All sampled dispersion curves
    for j = 1:1:NoPlot_all
        c_plot = store_all{4,order_all_plot(j)};
        lambda_plot = store_all{5,order_all_plot(j)};
        plot(c_plot,lambda_plot,'color',[0.8 0.8 0.8])
    end
    
    % Upper/lower boundaries for the experimental (observed) 
    % dispersion curve
    e = errorbar(c_OBS,lambda_OBS,c_OBS-c_OBS_low,c_OBS_up-c_OBS,...
        'horizontal','.','MarkerSize',0.5);
    e.LineWidth = 0.25;
    e.Color = 'r';    
    
    % 'Accepted' dispersion curves
    for j = 1:1:NoPlot
        c_plot = store_accepted{4,order_plot(j)};
        lambda_plot = store_accepted{5,order_plot(j)};
        plot(c_plot,lambda_plot,'color',linecolors(j,:))
    end
    
    % Experimental (observed) dispersion curve
    plot(c_OBS,lambda_OBS,'r--','LineWidth',1)
    
    % Visual
    FigFontSize = 8;
    set(gca, 'FontSize', FigFontSize)
    axis ij, grid on, box off
    set(gca,'TickDir','out')
    xlabel('Phase velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
    ylabel('Wavelength [m]','FontSize',FigFontSize,'Fontweight','normal')
    set(gca,'TickDir','out')
    
    %%% Vs PROFILES %%%
    ax(2)=subplot(1,2,2);
    hold on
    % All sampled Vs profiles
    for j = 1:1:NoPlot_all
        clear plot_layer_depth
        clear plot_beta
        h_plot = [0 store_all{2,order_all_plot(j)} max(MaxDepth-sum(store_all{2,order_all_plot(j)}),0)];
        beta_plot = store_all{1,order_all_plot(j)};
        plot_layer_depth = zeros(2*(n+1),1);
        plot_beta = zeros(2*(n+1),1);
        for i = 1:(n+1)
            plot_layer_depth(2*i-1) = sum(h_plot(1:i));
            plot_layer_depth(2*i) = sum(h_plot(1:i+1));
            plot_beta(2*i-1) = beta_plot(i);
            plot_beta(2*i) = beta_plot(i);
        end
        plot(plot_beta, plot_layer_depth,'color',[0.8 0.8 0.8]);
    end
    
    % 'Accepted' Vs profiles
    for j = 1:1:NoPlot
        clear plot_layer_depth
        clear plot_beta
        h_plot = [0 store_accepted{2,order_plot(j)} max(MaxDepth-sum(store_accepted{2,order_plot(j)}),0)];
        beta_plot = store_accepted{1,order_plot(j)};
        plot_layer_depth = zeros(2*(n+1),1);
        plot_beta = zeros(2*(n+1),1);
        for i = 1:(n+1)
            plot_layer_depth(2*i-1) = sum(h_plot(1:i));
            plot_layer_depth(2*i) = sum(h_plot(1:i+1));
            plot_beta(2*i-1) = beta_plot(i);
            plot_beta(2*i) = beta_plot(i);
        end
        plot(plot_beta, plot_layer_depth,'color',linecolors(j,:));
    end
           
    % Visual
    FigFontSize = 8;
    set(gca, 'FontSize', FigFontSize)
    axis ij, grid on, box off
    set(gca,'TickDir','out')
    xlabel('Shear wave velocity [m/s]','FontSize',FigFontSize,'Fontweight','normal')
    ylabel('Depth [m]','FontSize',FigFontSize,'Fontweight','normal'), ylim([0 MaxDepth])
    set(gca,'TickDir','out')
    
    %%% COLORBAR %%%
    colormap(mymap)
    c = colorbar;
    set(c, 'Position', [0.925*.8914  1.5*.11 .0281 0.95*.8150])
    for i = 1:2
        pos=get(ax(i), 'Position');
        set(ax(i), 'Position', [0.9*pos(1) 1.5*pos(2) 0.85*pos(3) 0.95*pos(4)]);
    end
    c.FontSize = FigFontSize;
    c.Label.FontSize = FigFontSize;
    c.Label.String = 'Misfit [%]';
    c.Ticks = 0:0.25:1;
    Y = prctile(sort_e_plot,0:25:100);
    Y = round(Y,1);
    c.TickLabels = Y(end:-1:1);
end
end