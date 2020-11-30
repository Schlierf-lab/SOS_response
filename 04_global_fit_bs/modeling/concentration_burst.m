%%%%%%%%%%%INPUT%%%%%%%%%%%%%
dataset = 2; % cip0p5:  '1'  cip3:  '2'  cip20:  '3'  uv10:  '4'
%%%%%%%%%%%INPUT%%%%%%%%%%%%%

x_cip = [0 10 25 45 60 90 120 180 210];% x-values for cip curves (min)
x_uv = [0 10 20 35 55 70 100 130 190]; % x-values for uv curves (min)

numBS = 10;
corr_fac = 2.09;

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

NumLexA = zeros(length(file_list), 6);
BS_all = cell(length(file_list), 3);
raw_BS_all = cell(length(file_list), 2);
area_all = zeros(length(file_list), 2);


% Set up fittype and options.
P = @(x, C, a, b) C*x.^(a-1).*exp(-x./b);
ft = fittype( 'C*x^(a-1)*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.Lower = [0 0 0];
opts.StartPoint = [10 15 15];

for i = 1:length(file_list)
    load([files{dataset} '\' file_list(i).name])    % load cell_params variable for condition
    
    %cell_params variable:
    % (:,1): LexA counts per cell 
    % (:,2): cell length
    % (:,3): cell contour length
    % (:,4): cell area
    % (:,5): cell dividing (0)?
    % (:,6): from which movie

    hascount = cell_params(:,1) ~= 0;   % cells that were sorted out before have zeros
    conc = corr_fac.*cell_params(hascount, 1)./cell_params(hascount, 4);  % calculate # of LexA/µm^2
    num = corr_fac.*cell_params(hascount, 1);
    area = cell_params(hascount, 4);
    
    area_all(i, 1) = mean(area);
    area_all(i, 2) = std(area);
    
    binsize = 35;
    edges = 0:binsize:2000;
    hD = histcounts(conc, edges, 'normalization', 'probability');
    
    fit_BS = zeros(numBS, 3);
    
    figure('position', [488.0000  633.8000  174.6000  128.2000]); hold on;
    
    for boot = 1:numBS
        conc_bs = datasample(conc, length(conc), 'replace', true); 
        
        hD_bs = histcounts(conc_bs, edges, 'normalization', 'probability');
        [xData, yData] = prepareCurveData(edges(1:end-1)+0.5*binsize, hD_bs);
        [fitres_bs, ~] = fit(xData, yData, ft, opts);
        coefs_bs = coeffvalues(fitres_bs);
        
        C = coefs_bs(1);
        a = coefs_bs(2);
        b = coefs_bs(3);
        fit_BS(boot,:) = [C, a, b];
        
        fitline_temp = P(0:1:2000, C, a, b);
        
        fac_temp = trapz(edges(1:end-1), hD)/trapz(0:1:2000, fitline_temp);
        plot(0:1:2000, fac_temp.*fitline_temp, '-r', 'linewidth', 0.5)
  
    end
    
    NumLexA(i,:) = [mean(fit_BS(:,1)), std(fit_BS(:,1)), mean(fit_BS(:,2)), std(fit_BS(:,2)), mean(fit_BS(:,3)), std(fit_BS(:,3))];
    BS_all{i,1} = fit_BS(:,1);
    BS_all{i,2} = fit_BS(:,2);
    BS_all{i,3} = fit_BS(:,3);
    
        
    fitline_l = P(0:1:2000, NumLexA(i,1), NumLexA(i,3), NumLexA(i,5));
        
    fac_l = trapz(edges(1:end-1), hD)/trapz(0:1:2000, fitline_l);
        
    stairs(edges(1:end-1), hD)
    plot(0:1:2000, fac_l.*fitline_l, '-k', 'linewidth', 2)
%     plot(1:1:1000, fac_g.*fitline_g, '-g', 'linewidth', 2)
    title([num2str(x(i)) ' min, n = ' num2str(length(conc))])
%     set(gca, 'xscale', 'log')
    xlim([0 1100])
    ylim([0 0.48])
end

figure('position', [257.8000  477.8000  368.0000  260.0000]);
errorbar(x, NumLexA(:,3), NumLexA(:,4), '-x', 'linewidth', 2);
xlabel('incubation time')
ylabel('burst frequency')

figure('position', [257.8000  477.8000  368.0000  260.0000]);
errorbar(x, NumLexA(:,5), NumLexA(:,6), '-x', 'linewidth', 2);
xlabel('incubation time')
ylabel('burst size')

