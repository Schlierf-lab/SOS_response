%%
[files, path] = uigetfile('.csv');
filename = [path files];

raw_output = importdata([filename(1:end-4) '.csv']);

frame = raw_output.data(:,2);
x = raw_output.data(:,3);
e = raw_output.data(:,16);
y = raw_output.data(:,4);
photons = raw_output.data(:,6);


%%
sel = y < 20;

k = figure('position', [1981.0 -0.094 329.6 136.8], 'renderer', 'painters');
%set(k,'paperunits','inches','paperposition',[0 0 50 50])
hold on; axis equal;
xlim([0 56]);
ylim([0 25]);
% plot([0 0], [0 512], '-k')
% plot([0 512], [512 512], '-k')
% plot([512 512], [0 512], '-k')
% plot([0 512], [0 0], '-k')

map = flipud(brewermap(500, 'Spectral'));
%map = [ones(15,3); map];
colormap(map);
h = dscatter(x(sel), y(sel), 'MARKER', 'o', 'BINS', [100 100], 'plottype', 'surf');
%h = scatter(x, y);
colorbar
%set(gca,'Box','on','Color',[0 0 0]);

%%
%print(k, 'test.tif', '-dtiff', '-r100');



