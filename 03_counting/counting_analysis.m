% counting results analysis
um_px = 0.1067;

[files, path] = uigetfile('*Count.mat','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

cell_params = [];
counts_per_cluster = [];
for movie = 1:file_count

    load([path files{movie}]);
    load([path files{movie}(1:end-9) 'out_corr.mat'])
    
    cell_count = zeros(cellListN, 6);
    for cell = 1:cellListN
        model = cellList.meshData{1, 1}{1, cell}.model;
        mesh = cellList.meshData{1, 1}{1, cell}.mesh;
        
        if cellList.selected2(cell) && mesh(1) ~= 0
            in = inpolygon(countArr(:,1), countArr(:,2), model(:,1).*106.7, model(:,2).*106.7);
            count_single = sum(countArr(in,3));
        
            midx = mesh(:,1) + (mesh(:,3)-mesh(:,1)) ./ 2;
            midy = mesh(:,2) + (mesh(:,4)-mesh(:,2)) ./ 2;
            d = diff([midx, midy]);
            e = diff([model(:,2), model(:,2)]);
            celllength = sum(sqrt(sum(d.*d,2)))*um_px;
            contlength = sum(sqrt(sum(e.*e,2)))*um_px;
            area = polyarea(model(:,1).*um_px, model(:,2).*um_px);
            cell_count(cell,:) = [count_single celllength contlength area cellList.selected1(cell) movie];
        end
    end
    cell_params = [cell_params; cell_count];
    counts_per_cluster = [counts_per_cluster; countArr(:,3)];
end

%% plot
figure('Position', [180 514 1090 213]);
hascount = cell_params(:,1) ~= 0;

subplot(1,4,1)
histogram(cell_params(hascount,1), 'BinWidth', 50);
title(['n = ' num2str(length(cell_params(hascount,1)))]);
xlabel('molecules per cell')
lim = get(gca, 'xlim');

subplot(1,4,2)
histogram(cell_params(hascount, 1)./cell_params(hascount, 4), 'BinWidth', 25);
text(0.6, 0.8, ['mean: ' num2str(mean(cell_params(hascount, 1)./cell_params(hascount, 4)), '%.0f')], 'Units', 'normalized')
xlabel('molecules per µm^2')
lim = get(gca, 'xlim');


subplot(1,4,3)
histogram(counts_per_cluster, 'BinWidth', 1, 'normalization', 'probability');
xlim([0 15]);
ylim([0 0.7]);
xlabel('molecules per cluster');
subplot(1,4,4)
scatter(cell_params(hascount, 2), cell_params(hascount, 1)./cell_params(hascount, 4), [], [0 0 0], 'filled', 'markerfacealpha', 0.5)
lin = lsline; lin.Color = 'r'; lin.LineWidth = 2;
title('concentration vs. cell length')
ylabel('molecules per µm^2')
xlabel('cell length (µm)')

%% plot Result
figure;

minN=1;
maxN = 5;


cmap=brewermap(maxN-1, 'spectral');
scatter(countArr(:,1)./106.7,countArr(:,2)./106.7,120,countArr(:,3),'filled');
colorbar;
caxis([minN maxN]);
colormap(cmap);
hold on;

for iterO=1:cellListN
    
    border=cellList.meshData{1}{iterO}.model;
    
    plot(border(:,1),border(:,2),'-k','LineWidth',0.5);
end

set(gca,'ydir', 'reverse');
axis equal;

xlim([0, 512])
ylim([0, 512])
