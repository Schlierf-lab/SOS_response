[files, path] = uigetfile('*filt.tracked.loc.txt','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

mols_per_cell = zeros(file_count, 1);
for rep = 1:file_count
    filename = [path files{rep}];

    tracks = importdata(filename, ',', 1);
    tracks = tracks.data;
    load([filename(1:end-20) 'out_corr.mat'], 'cellListN');
    
    disp([num2str(rep) '/' num2str(file_count)]);

    loc_frame = zeros(max(tracks(:,1)), 1);
    for i = 1:max(tracks(:,1))
        loc_frame(i) = sum(tracks(:,1) == i);
    end

    mols_per_cell(rep) = mean(loc_frame(end-50:end))/cellListN;
end

figure;
errorbar(1, mean(mols_per_cell), std(mols_per_cell), 'o', 'linewidth', 2);

figure;
errorbar(1, 1/mean(mols_per_cell), 1/std(mols_per_cell), 'o', 'linewidth', 2);