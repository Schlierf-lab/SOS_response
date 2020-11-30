[files, path] = uigetfile('*out.mat','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

um_px = 0.1067;
numBS = 100;

all_params = [];
for movie = 1:file_count
    in_filename = [path files{movie}];
    
    load(in_filename);
    
    
    params = zeros(cellListN, 4); % col1: not dividing (1) or dividing (0), col2: included or excluded, col3: length, col4: area, col5: contour length, col6: movieNR
    for i = 1:cellListN
        model = cellList.meshData{1, 1}{1, i}.model;
        mesh = cellList.meshData{1, 1}{1, i}.mesh;        
        if mesh(1) ~= 0
            midx = mesh(:,1) + (mesh(:,3)-mesh(:,1)) ./ 2;
            midy = mesh(:,2) + (mesh(:,4)-mesh(:,2)) ./ 2;
            d = diff([midx, midy]);
            e = diff([model(:,2), model(:,2)]);
            celllength = sum(sqrt(sum(d.*d,2)))*um_px;
            contlength = sum(sqrt(sum(e.*e,2)))*um_px;
            area = polyarea(model(:,1).*um_px, model(:,2).*um_px);
            params(i,:) = [celllength contlength area movie];
        end
    end
    all_params = [all_params; params];
end

length_all = all_params(all_params(:,1)>0,1);

figure('position', [488.0000  633.8000  174.6000  128.2000]); hold on;
binsize = 0.25;
edges = 0:binsize:10;
hD = histcounts(length_all, edges, 'normalization', 'probability');

logfit_BS = zeros(numBS, 2);
for boot = 1:numBS
    length_bs = datasample(length_all, length(length_all), 'replace', true); 
    logfit = lognfit(length_bs);     % fit log-normal distribution to LexA concentration
    logfit_BS(boot,:) = logfit;
    fitline_temp = lognpdf(0:0.01:10, logfit(1), logfit(2));
    fac_temp = trapz(edges(1:end-1), hD)/trapz(0:0.01:10, fitline_temp);
    plot(0:0.01:10, fac_temp.*fitline_temp, '-r', 'linewidth', 0.5)
end

fitline_l = lognpdf(0:0.01:10, mean(logfit_BS(:,1)), mean(logfit_BS(:,2)));

fac_l = trapz(edges(1:end-1), hD)/trapz(0:0.01:10, fitline_l);

stairs(edges(1:end-1), hD)
plot(0:0.01:10, fac_l.*fitline_l, '-k', 'linewidth', 2)

NumLexA_mean = mean(exp(logfit_BS(:,1)+(logfit_BS(:,2).^2)/2));
NumLexA_std = std(exp(logfit_BS(:,1)+(logfit_BS(:,2).^2)/2));

title([num2str(NumLexA_mean, '%.2f') ' ± ' num2str(NumLexA_std, '%.2f') ' µm'])
