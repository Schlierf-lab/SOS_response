[files, path] = uigetfile('.csv');
filename = [path files];

raw_output = importdata([filename(1:end-4) '.csv']);
%raw_output = importdata('W:\Leo\LexA_fixed\2019-01-24_Dd2_fermi\Dd2_fermi_30_001.csv');
% outlines = open([path [files(1:end-4) '_out_corr.mat']]);
outlines = open([path [files(1:end-4) '_out_corr.mat']]);
%BF = imread('C:\Users\Leo\Documents\Group_Meeting\Apr9\overlay.jpg');
outlines = outlines.cellList.meshData{1, 1};

frame = raw_output.data(:,1);
x = raw_output.data(:,2);
e = raw_output.data(:,15);
y = raw_output.data(:,3);
%y_e = raw_output.data(:,16);
%z = raw_output.data(:,4);
photons = raw_output.data(:,5);

sel = frame > 500;

figure('Renderer', 'Painters'); hold on
axis equal
% map = colormap(brewermap(256, 'Spectral'));
h = scatter(x(sel), y(sel), 15, frame(sel), 'filled', 'r');
set(gca, 'ydir', 'reverse')
colorbar

for i = 1:length(outlines)
    if ~isempty(outlines{1, i}.model)
        plot(outlines{1, i}.model(:,1), outlines{1, i}.model(:,2), 'color', 'k', 'LineWidth', 1.5)
    end
end
xlim([-10 522])
ylim([-10 522])
scalebar('ScaleLengthRatio', 0.25, 'Colour', [0 0 0], 'Bold', 1, 'Unit', 'µm', 'Location', 'southeast');