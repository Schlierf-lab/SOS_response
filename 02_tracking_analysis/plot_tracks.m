minlength_plot = 5;

[files, path] = uigetfile('*filt.csv');
% if only one file is selected

in_filename = [path files];
tracks = importdata([in_filename(1:end-4) '.tracked.loc.txt'], ',', 1);
outlines = open([in_filename(1:end-9) '_out_corr.mat']);
outlines = outlines.cellList.meshData{1, 1};
tracks = tracks.data;

%% plot tracks

figure('renderer', 'painters'); hold on; axis equal
for i = 1:max(tracks(:,18))
    if sum(tracks(:,18)==i)>=minlength_plot
        plot(tracks(tracks(:,18)==i, 2),tracks(tracks(:,18)==i, 3), '-o', 'LineWidth', 1.2, 'markersize', 2)
    end
end

for i = 1:length(outlines)
    if ~isempty(outlines{1, i}.model)
       plot(outlines{1, i}.model(:,1), outlines{1, i}.model(:,2), 'color', 'k', 'LineWidth', 2)
    end
end
% xlim([30 70])
% ylim([160 200])
set(gca, 'ydir', 'reverse')
scalebar('ScaleLengthRatio', 0.25, 'Colour', [0 0 0], 'Bold', 1, 'Unit', 'nm', 'Location', 'southeast');
