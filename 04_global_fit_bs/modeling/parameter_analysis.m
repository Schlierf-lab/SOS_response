%%%%%%%%%%%INPUT%%%%%%%%%%%%%
dataset = 3; % cip0p5:  '1'  cip3:  '2'  cip20:  '3'  uv10:  '4'
%%%%%%%%%%%INPUT%%%%%%%%%%%%%

x_cip = [0 10 25 45 60 90 120 180 210];% x-values for cip curves (min)
x_uv = [0 10 20 35 55 70 100 130 190]; % x-values for uv curves (min)

boolLen = 0;
boolConc = 1;

numBS = 100;

if dataset == 4
    x = x_uv;
else
    x = x_cip;
end

% import data
files = {'co_cip05';...
         'co_cip3';...
         'co_cip20';...
         'co_uv10';...
         };

file_list = dir([files{dataset} '\*.mat*']);

params_bs = zeros(length(file_list), 8);
BS_all = cell(length(file_list), 4);
for i = 1:length(file_list)
    load([files{dataset} '\' file_list(i).name])    % load cell_params variable for condition
    
    %cell_params variable:
    % (:,1): LexA counts per cell 
    % (:,2): cell length
    % (:,3): cell contour length
    % (:,4): cell area
    % (:,5): cell dividing (0)?
    % (:,6): from which movie

    sel = logical(cell_params(:,5));   % cells that were sorted out before have zeros
    conc = cell_params(sel, 1)./cell_params(sel, 4);  % calculate # of LexA/µm^2
    length_all = cell_params(sel, 2);
    
    if boolConc
        fc = figure('position', [488.0000  633.8000  174.6000  128.2000]); hold on;
    end
    
    if boolLen
        fl = figure('position', [669.8000  633.8000  174.4000  128.0000]); hold on;
    end
    
    binsize_c = 20;
    edges_c = 0:binsize_c:1000;
    hD_c = histcounts(conc, edges_c, 'normalization', 'probability');
    
    binsize_l = 0.3;
    edges_l = 0:binsize_l:10;
    hD_l = histcounts(length_all, edges_l, 'normalization', 'probability');

    
    % bootstrapping concentration
    logfit_BS = zeros(numBS, 4);
    for boot = 1:numBS
        conc_bs = datasample(conc, length(conc), 'replace', true); 
        logfit = lognfit(conc_bs);     % fit log-normal distribution to LexA concentration
        logfit_BS(boot,1:2) = logfit;
        fitline_temp = lognpdf(0:1:1000, logfit(1), logfit(2));
        fac_temp = trapz(edges_c(1:end-1), hD_c)/trapz(0:1:1000, fitline_temp);
        if boolConc
            figure(fc)
            plot(0:1:1000, fac_temp.*fitline_temp, '-r', 'linewidth', 0.5)
        end
        
        length_bs = datasample(length_all, length(length_all), 'replace', true); 
        logfit = lognfit(length_bs);     % fit log-normal distribution to LexA concentration
        logfit_BS(boot,3:4) = logfit;
        fitline_temp = lognpdf(0:0.01:10, logfit(1), logfit(2));
        fac_temp = trapz(edges_l(1:end-1), hD_l)/trapz(0:0.01:10, fitline_temp);
        if boolLen
            figure(fl)
            plot(0:0.01:10, fac_temp.*fitline_temp, '-r', 'linewidth', 0.5)
        end
    end
    
    
    % get parameters and plot
    params_bs(i,1:4) = [mean(logfit_BS(:,1)), std(logfit_BS(:,1)), mean(logfit_BS(:,2)), std(logfit_BS(:,2))];
    params_bs(i,5:8) = [mean(logfit_BS(:,3)), std(logfit_BS(:,3)), mean(logfit_BS(:,4)), std(logfit_BS(:,4))];
    BS_all{i,1} = logfit_BS(:,1);
    BS_all{i,2} = logfit_BS(:,2);
    BS_all{i,3} = logfit_BS(:,3);
    BS_all{i,4} = logfit_BS(:,4);
    
        
    fitline_c = lognpdf(0:1:1000, params_bs(i,1), params_bs(i,3));        
    fac_c = trapz(edges_c(1:end-1), hD_c)/trapz(0:1:1000, fitline_c);
    if boolConc    
        figure(fc)
        stairs(edges_c(1:end-1), hD_c)
        plot(0:1:1000, fac_c.*fitline_c, '-k', 'linewidth', 2)
        title([num2str(x(i)) ' min, n = ' num2str(length(conc))])
        xlim([0 500])
        ylim([0 0.4])
    end
    
    fitline_l = lognpdf(0:0.01:10, params_bs(i,5), params_bs(i,7));
    fac_l = trapz(edges_l(1:end-1), hD_l)/trapz(0:0.01:10, fitline_l);
    if boolLen
        figure(fl)
        stairs(edges_l(1:end-1), hD_l)
        plot(0:0.01:10, fac_l.*fitline_l, '-k', 'linewidth', 2)
        title([num2str(x(i)) ' min, n = ' num2str(length(conc))])
    end
end

NumLexA_ex = zeros(length(BS_all), 4);
CV = zeros(length(BS_all), 4);
for j = 1:length(BS_all)
    NumLexA_ex(j,1) = mean(exp(BS_all{j,1}+(BS_all{j,2}.^2)/2));
    NumLexA_ex(j,2) = std(exp(BS_all{j,1}+(BS_all{j,2}.^2)/2));
    NumLexA_ex(j,3) = mean(exp(BS_all{j,3}+(BS_all{j,4}.^2)/2));
    NumLexA_ex(j,4) = std(exp(BS_all{j,3}+(BS_all{j,4}.^2)/2));
%     CV(j,1) = mean(sqrt(exp(BS_all{j,2})-1));
%     CV(j,2) = std(sqrt(exp(BS_all{j,2})-1));
%     CV(j,3) = mean(sqrt(exp(BS_all{j,4})-1));
%     CV(j,4) = std(sqrt(exp(BS_all{j,4})-1));
    CV(j,1) = mean(sqrt(exp(BS_all{j,2}.^2)-1));
    CV(j,2) = std(sqrt(exp(BS_all{j,2}.^2)-1));
    CV(j,3) = mean(sqrt(exp(BS_all{j,4}.^2)-1));
    CV(j,4) = std(sqrt(exp(BS_all{j,4}.^2)-1));
end
x = x(1:length(file_list));

figure('position', [257.8000  477.8000  368.0000  260.0000]);
yyaxis left
errorbar(x, NumLexA_ex(:,1), NumLexA_ex(:,2), '-x', 'linewidth', 2);
ylabel('LexA/µm^2')
yyaxis right
errorbar(x, NumLexA_ex(:,3), NumLexA_ex(:,4), '-x', 'linewidth', 2);
ylabel('cell length (µm)')
title([files{dataset}(4:end) ': LexA concentration'])
xlabel('incubation time')
% xlim([-5 215])
% ylim([30 500])

figure('position', [813.6667  500.3333  261.3333  263.3333]); hold on;
errorbar(x, CV(:,1), CV(:,2), '-', 'linewidth', 2);
errorbar(x, CV(:,3), CV(:,4), '-', 'linewidth', 2);
ylabel('coefficient of variation')
title([files{dataset}(4:end) ': extrinsic noise'])
xlabel('incubation time')
%legend('LexA concentration', 'cell length', 'location', 'southeast');#
% ylim([0.35 1.2])

figure('position', [257.8000  477.8000  368.0000  260.0000]);
errorbar(NumLexA_ex(:,1), NumLexA_ex(:,3), NumLexA_ex(:,4), NumLexA_ex(:,4), NumLexA_ex(:,2), NumLexA_ex(:,2), 'linestyle', 'none');
text(NumLexA_ex(:,1)-5, NumLexA_ex(:,3)+0.05, num2str(x(:)), 'HorizontalAlignment', 'right');
xlabel('LexA/µm^2')
ylabel('cell length')

