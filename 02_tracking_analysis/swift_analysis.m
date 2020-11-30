%% import swift data
%%%%%INPUT%%%%%
minlength_MJD = 3; % for MJD analysis
minlength_JD = 3; % for JD analysis
binsize_MJD = 10;
binsize_JD = 5;
threshold = [55]; % set thresholds for MJD plot

plot_thresh = [105 175]; % set thresholds for coloring the tracks
plot_movienr = 1; % which movie do you want to plot the tracks + outlines from?
minlength_plot = 1; % minimum length for plotting (no. of localizatios)

thresh = 3; % for D* analysis, number of steps (no. of localizations - 1)

boolStat = 0; % wanna get single-cell statistics? 
boolCell = 0; % wanna do MJD analysis for certain cell length? boolStat must be '1'. slows it down a lot.
length_spec = [0 5]; % cell length for MJD analysis
%%%%%INPUT%%%%%

s_frame = 0.02;
um_px = 0.1067;

[files, path] = uigetfile('*filt.csv','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

datall = struct('tracks', [], 'segs', [], 'dynamics', [], 'outlines', []);
for movie = 1:file_count
    in_filename = [path files{movie}];
    
    tracks = importdata([in_filename(1:end-4) '.tracked.loc.txt'], ',', 1);
    segs = readtable([in_filename(1:end-4) '.tracked.seg.txt']);
    outlines = open([in_filename(1:end-9) '_out_corr.mat']);
    datall(movie).outlines = outlines.cellList.meshData{1, 1};
    datall(movie).tracks = tracks.data;
    datall(movie).segs = table2array(segs(2:end,1:5));
    datall(movie).dynamics = table2cell(segs(:,6));
end


%% MJD and JD Analysis 

% JD analysis
JD_all = [];
for rep = 1:length(datall)

    trajectories = cell(max(datall(rep).tracks(:,18)),2);
    for i = 1:max(datall(rep).tracks(:,18))
        if sum(datall(rep).tracks(:,18)==i) >= (minlength_JD+1)
            x_track = datall(rep).tracks(datall(rep).tracks(:,18)==i, 13);
            y_track = datall(rep).tracks(datall(rep).tracks(:,18)==i, 14);
            trajectories(i,1) = num2cell(i);
            trajectories(i,2) = mat2cell([x_track, y_track], length(x_track));
        end
    end
    
    trajectories_filt = trajectories(~cellfun('isempty',trajectories(:,1)),:);
    
    JD_file = [];
    for w = 1:length(trajectories_filt)
        JD = [];
        for d = 2:length(trajectories_filt{w,2})
            dist = pdist([trajectories_filt{w,2}(d,1) trajectories_filt{w,2}(d,2); trajectories_filt{w,2}(d-1,1) trajectories_filt{w,2}(d-1,2)], 'euclidean');
            JD(end+1) = dist;
        end
        JD_file = [JD_file JD];
    end
    JD_all = [JD_all JD_file];
end

% figure;
% edges_JD = 0:binsize_JD:1000;
% hJD = histc(JD_all, edges_JD);
% bar(edges_JD+0.5*binsize_JD,hJD./sum(hJD), 1, 'EdgeColor', 'none', 'LineStyle', 'none');
% xlim([0 1000]);
% xlabel('Jumping Distance (nm)')
% ylabel('frequency')
% title(['n = ' num2str(length(JD_all))])

%rayleigh_JD = rayleigh(JD_all, binsize_JD, displaymodels);

allsegs = [];
for i = 1:length(datall)
    allsegs = [allsegs; datall(i).segs];
end

% MJD analysis
MJDn = allsegs(:,5);
sel = MJDn >= minlength_MJD;
MJDn = MJDn(sel);
MJD = allsegs(:,3);
MJD = MJD(sel);
MJDe = allsegs(:,4);
MJDe = MJDe(sel);


MJDw = MJD;
for i = 1:length(MJD)
    if MJDn(i) > 1
        weight = repmat(MJD(i), MJDn(i)-1, 1);
        MJDw = [MJDw; weight];
    end
end

figure;
edges = 0:binsize_MJD:1000;
hMJDw = histc(MJDw, edges);
bar(edges+0.5*binsize_MJD,hMJDw./sum(hMJDw), 1);
pw = gca;
hold on;
% for jj = 2:2:2*displaymodel_MJD
%     plot([rayleigh_JD(displaymodel_MJD).Parameters(jj) rayleigh_JD(displaymodel_MJD).Parameters(jj)].*sqrt(pi/2), [0 max(hMJDw./sum(hMJDw))], '-m', 'linewidth', 2)
% end
xlim([0 800]);
xlabel('Mean Jumping Distance (nm)')
ylabel('frequency')
title(['weighted, n = ' num2str(length(MJD))])

figure;
hMJD = histc(MJD, edges);
bar(edges+0.5*binsize_MJD,hMJD./sum(hMJD), 1);
puw = gca;
hold on;
% for jj = 2:2:2*displaymodel_MJD
%     plot([rayleigh_JD(displaymodel_MJD).Parameters(jj) rayleigh_JD(displaymodel_MJD).Parameters(jj)].*sqrt(pi/2), [0 max(hMJD./sum(hMJD))], '-m', 'linewidth', 2)
% end
xlim([0 800]);
xlabel('Mean Jumping Distance (nm)')
ylabel('frequency')
title('unweighted')

% threshold analysis
if ~isempty(threshold)
    
    % unweighted
    prob = zeros(length(threshold)+1,1);
    prob(1) = sum(hMJDw(edges<threshold(1))./sum(hMJDw));
    plot(pw, [threshold(1) threshold(1)], [0 max(hMJDw./sum(hMJDw))], '-r', 'linewidth', 2)
    for j = 2:length(threshold)
       prob(j) =  sum(hMJDw(edges<threshold(j)&edges>=threshold(j-1))./sum(hMJDw));
       plot(pw, [threshold(j) threshold(j)], [0 max(hMJDw./sum(hMJDw))], '-r', 'linewidth', 2)
    end
    prob(end) = sum(hMJDw(edges>=threshold(end))./sum(hMJDw));
    text(pw, 0.8, 0.8, num2str(100.*prob,'%.2f'), 'Units', 'normalized')
    
    % weighted
    prob = zeros(length(threshold)+1,1);
    prob(1) = sum(hMJD(edges<threshold(1))./sum(hMJD));
    plot(puw, [threshold(1) threshold(1)], [0 max(hMJD./sum(hMJD))], '-r', 'linewidth', 2)
    for j = 2:length(threshold)
       prob(j) =  sum(hMJD(edges<threshold(j)&edges>=threshold(j-1))./sum(hMJD));
       plot(puw, [threshold(j) threshold(j)], [0 max(hMJD./sum(hMJD))], '-r', 'linewidth', 2)
    end
    prob(end) = sum(hMJD(edges>=threshold(end))./sum(hMJD));
    text(puw, 0.8, 0.8, num2str(100.*prob,'%.2f'), 'Units', 'normalized')
end

%% cell statistics
if boolStat
    cellData = struct('cellID', [], 'length', [], 'MJD_cell', []);
    for rep = 1:length(datall) 
        lengths = zeros(length(datall(rep).outlines),1);
        IDs = zeros(length(datall(rep).outlines),1);
        for i = 1:length(datall(rep).outlines)
            if ~isempty(datall(rep).outlines{1, i}.model) && sum(sum(datall(rep).outlines{1, i}.mesh)) > 0
                model = datall(rep).outlines{1, i}.model;
                mesh = datall(rep).outlines{1, i}.mesh;
                midx = mesh(:,1) + (mesh(:,3)-mesh(:,1)) ./ 2;
                midy = mesh(:,2) + (mesh(:,4)-mesh(:,2)) ./ 2;

                d = diff([midx, midy]); 
                lengths(i) = sum(sqrt(sum(d.*d,2)))*um_px;
                IDs(i) = i;
            end
        end
        cellData(rep).length = lengths;
        cellData(rep).cellID = IDs;

%         % get MJDs per cell (not really useful)
%         cellData(rep).MJD_cell = cell(length(cellData(rep).cellID), 1);
%         cellData(rep).MJD_cell(cellfun('isempty',cellData(rep).MJD_cell)) = {0};
%         for segment = 1:length(datall(rep).segs)
%             if datall(rep).segs(segment,5) >= minlength_MJD
%                 IDsel = unique(datall(rep).tracks(datall(rep).tracks(:,18) == datall(rep).segs(segment,1), 17));
%                 temp = cellData(rep).MJD_cell(cellData(rep).cellID == IDsel);
%                 cellData(rep).MJD_cell(cellData(rep).cellID == IDsel) = {[temp{:} datall(rep).segs(segment,3)]};
%             end
%         end
    end

    % get lengths from all movies
    lengths_all = vertcat(cellData(:).length);

    % mean JD per cell
    MJD_cell_all = [];
    for rep = 1:length(datall)
        MJD_cell = zeros(length(datall(rep).outlines),3);
        for cellNR = 1:length(datall(rep).outlines)
            if ~isempty(datall(rep).outlines{1, cellNR}.model) && sum(sum(datall(rep).outlines{1, cellNR}.mesh)) > 0
                cell_tracks = unique(datall(rep).tracks(datall(rep).tracks(:,17) == cellNR, 18));
                trajectories = cell(length(cell_tracks),2);

                iter = 1;
                for i = cell_tracks'
                    if sum(datall(rep).tracks(:,18)==i) >= (minlength_JD+1)
                        x_track = datall(rep).tracks(datall(rep).tracks(:,18)==i, 13);
                        y_track = datall(rep).tracks(datall(rep).tracks(:,18)==i, 14);
                        trajectories(iter,1) = num2cell(i);
                        trajectories(iter,2) = mat2cell([x_track, y_track], length(x_track));
                        iter = iter + 1;
                    end
                end
                trajectories_filt = trajectories(~cellfun('isempty',trajectories(:,1)),:);

                JD_cell = [];
                kk = size(trajectories_filt);
                for w = 1:kk(1)
                    JD = [];
                    for d = 2:length(trajectories_filt{w,2})
                        dist = pdist([trajectories_filt{w,2}(d,1) trajectories_filt{w,2}(d,2); trajectories_filt{w,2}(d-1,1) trajectories_filt{w,2}(d-1,2)], 'euclidean');
                        JD(end+1) = dist;
                    end
                    JD_cell = [JD_cell JD];
                end
                MJD_cell_single = mean(JD_cell);
                MJDstd_cell_single = std(JD_cell);
                MJDn_cell_single = length(JD_cell);

                MJD_cell(cellNR,:) = [MJD_cell_single MJDstd_cell_single MJDn_cell_single];
            end
        end
        MJD_cell_all = [MJD_cell_all; MJD_cell];
    end

%     % get mean cell MJDs for all movies (not used later)
%     mean_MJD_cell = [];
%     for rep = 1:length(cellData)
%         cellmean = zeros(length(cellData(rep).MJD_cell), 1);
%         for i = 1:length(cellData(rep).MJD_cell)
%             temp = cellData(rep).MJD_cell(i);
%             cellmean(i) = mean(nonzeros(temp{:}));
%         end
%         mean_MJD_cell = [mean_MJD_cell; cellmean];
%     end

    figure;
    subplot(2,2,1);
    histogram(lengths_all, 50);
    xlim([0 12])
    xlabel('cell length (µm)')
    ylabel('# of cells')
    title('single-cell length distribution')
    text(0.7, 0.8, {['n = ' num2str(length(lengths_all))]; ['l = ' num2str(mean(lengths_all), '%.2f') ' ± ' num2str(std(lengths_all), '%.2f') ' µm']}, 'Units', 'normalized')

    subplot(2,2,2);
    histogram(MJD_cell_all(:,1), 50);
    xlim([0 750])
    xlabel('mean MJD/single cell')
    ylabel('# of cells')
    title('single-cell MJD')

    subplot(2,2,[3,4]);
    scatter(lengths_all, MJD_cell_all(:,1), 'filled', 'MarkerFaceAlpha', 0.4)
    lin = lsline; lin.Color = 'r'; 
    ylim([0 750])
    xlim([0 12])
    xlabel('cell length')
    ylabel('single-cell MJD')
    title('cell length vs single-cell MJD')

    %sgtitle(['single-cell statistics: n = ' num2str(length(lengths_all))])
end




% extract MJDs from selected cells only
if boolCell

    allsegs_sel = [];
    for rep = 1:length(datall)
        for segment = 1:length(datall(rep).segs)
            IDsel = unique(datall(rep).tracks(datall(rep).tracks(:,18) == datall(rep).segs(segment,1), 17));
            if ~isempty(cellData(rep).length(cellData(rep).cellID == IDsel)) && cellData(rep).length(cellData(rep).cellID == IDsel) >= length_spec(1) && cellData(rep).length(cellData(rep).cellID == IDsel) <= length_spec(2)
                allsegs_sel(end+1,:) = datall(rep).segs(segment, :);
            end
            cellData(rep).MJD_cell = datall(rep).segs(segment, :);
        end
    end
    % MJD analysis
    MJDn = allsegs_sel(:,5);
    sel = MJDn >= minlength_MJD;
    MJDn = MJDn(sel);
    MJD = allsegs_sel(:,3);
    MJD = MJD(sel);
    MJDe = allsegs_sel(:,4);
    MJDe = MJDe(sel);


    MJDw = MJD;
    for i = 1:length(MJD)
        if MJDn(i) > 1
            weight = repmat(MJD(i), MJDn(i)-1, 1);
            MJDw = [MJDw; weight];
        end
    end

    figure;
    edges = 0:binsize_MJD:1000;
    hMJDw = histc(MJDw, edges);
    bar(edges+0.5*binsize_MJD,hMJDw./sum(hMJDw), 1);
    pw = gca;
    hold on;
    xlim([0 1000]);
    xlabel('Mean Jumping Distance (nm)')
    ylabel('frequency')
    title(['cells between ' num2str(length_spec(1)) ' and ' num2str(length_spec(2)) ' µm, weighted, n = ' num2str(length(MJD))])

    figure;
    hMJD = histc(MJD, edges);
    bar(edges+0.5*binsize_MJD,hMJD./sum(hMJD), 1);
    puw = gca;
    hold on;
    xlim([0 1000]);
    xlabel('Mean Jumping Distance (nm)')
    ylabel('frequency')
    title(['cells between ' num2str(length_spec(1)) ' and ' num2str(length_spec(2)) ' µm, unweighted'])


    % JD analysis
    
    JD_all = [];
    for rep = 1:length(datall)

        trajectories = cell(max(datall(rep).tracks(:,18)),2);
        for i = 1:max(datall(rep).tracks(:,18))
            if sum(datall(rep).tracks(:,18)==i)>=minlength_MJD+1 && ~isempty(cellData(rep).length(cellData(rep).cellID == unique(datall(rep).tracks(datall(rep).tracks(:,18)==i,17)))) && cellData(rep).length(cellData(rep).cellID == unique(datall(rep).tracks(datall(rep).tracks(:,18)==i,17))) >= length_spec(1) && cellData(rep).length(cellData(rep).cellID == unique(datall(rep).tracks(datall(rep).tracks(:,18)==i,17))) <= length_spec(2)
                x_track = datall(rep).tracks(datall(rep).tracks(:,18)==i, 13);
                y_track = datall(rep).tracks(datall(rep).tracks(:,18)==i, 14);
                trajectories(i,1) = num2cell(i);
                trajectories(i,2) = mat2cell([x_track, y_track], length(x_track));
            end
        end

        trajectories_filt = trajectories(~cellfun('isempty',trajectories(:,1)),:);

        JD_file = [];
        for w = 1:length(trajectories_filt)
            JD = [];
            for d = 2:length(trajectories_filt{w,2})
                dist = pdist([trajectories_filt{w,2}(d,1) trajectories_filt{w,2}(d,2); trajectories_filt{w,2}(d-1,1) trajectories_filt{w,2}(d-1,2)], 'euclidean');
                JD(end+1) = dist;
            end
            JD_file = [JD_file JD];
        end
        JD_all = [JD_all JD_file];
    end

    figure;
    edges_JD = 0:binsize_JD:1000;
    hJD = histc(JD_all, edges_JD);
    bar(edges_JD+0.5*binsize_JD,hJD./sum(hJD), 1, 'EdgeColor', 'none', 'LineStyle', 'none');
    xlim([0 1000]);
    xlabel('Jumping Distance (nm)')
    ylabel('frequency')
    title(['cells between ' num2str(length_spec(1)) ' and ' num2str(length_spec(2)) ' µm, n = ' num2str(length(JD_all))])
end

%% plot tracks
color = [0.6392 0.0784 0.1804];

figure('renderer', 'painters'); hold on; axis equal
for i = 1:max(datall(plot_movienr).tracks(:,18))
    if sum(datall(plot_movienr).tracks(:,18)==i)>=minlength_plot && datall(plot_movienr).segs(i,3) ~= 0
        if datall(plot_movienr).segs(i,3) > plot_thresh(1) && datall(plot_movienr).segs(i,3) < plot_thresh(2)
            color = [0.5216 0.8314 0.7216];
        elseif datall(plot_movienr).segs(i,3) <= plot_thresh(1)
            color = [0.6392 0.0784 0.1804];
        else
            color = [1 0 1];
        end
        plot(datall(plot_movienr).tracks(datall(plot_movienr).tracks(:,18)==i, 13),datall(plot_movienr).tracks(datall(plot_movienr).tracks(:,18)==i, 14), '-o', 'color', color, 'markerfacecolor', color, 'LineWidth', 1.2, 'markersize', 2)
    end
end

if ~boolCell
    for i = 1:length(datall(plot_movienr).outlines)
        if ~isempty(datall(plot_movienr).outlines{1, i}.model)
            plot(datall(plot_movienr).outlines{1, i}.model(:,1).*106.7, datall(plot_movienr).outlines{1, i}.model(:,2).*106.7, 'color', 'k', 'LineWidth', 2)
        end
    end
end

if boolCell
    for i = 1:length(datall(plot_movienr).outlines)
        if ~isempty(datall(plot_movienr).outlines{1, i}.model)
            if ~isempty(cellData(plot_movienr).length(cellData(plot_movienr).cellID == i)) && cellData(plot_movienr).length(cellData(plot_movienr).cellID == i) >= length_spec(1) && cellData(plot_movienr).length(cellData(plot_movienr).cellID == i) <= length_spec(2)
                color = 'k';
            else
                color = 'r';
            end
            plot(datall(plot_movienr).outlines{1, i}.model(:,1).*106.7, datall(plot_movienr).outlines{1, i}.model(:,2).*106.7, 'color', color, 'LineWidth', 2)
        end
    end
end
xlim([0 512*106.7])
ylim([0 512*106.7])
set(gca, 'ydir', 'reverse')
scalebar('ScaleLengthRatio', 0.25, 'Colour', [0 0 0], 'Bold', 1, 'Unit', 'nm', 'Location', 'southeast');
