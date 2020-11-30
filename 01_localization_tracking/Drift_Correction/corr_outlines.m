% parameters
nm_px = 106.7;
firstFrame = 501;
% import swift data
[files, path] = uigetfile('.csv','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

for movie = 1:file_count
    
    disp(['Working on file ' num2str(movie) '/' num2str(file_count) '...']);
    % import .csv file from SMAP and outline file
    file_location = [path files{movie}];

    SMAP = importdata(file_location);
    SMAP_dat = SMAP.data(SMAP.data(:,1)>=firstFrame,:);
    SMAP_txt = SMAP.textdata;
    
    load([file_location(1:end-4) '_out.mat'], 'cellList', 'cellListN');
    [cellList, shift_xy] = optCellBorder(SMAP_dat, cellList, cellListN, -5, 5);
    save([file_location(1:end-4) '_out_corr.mat'], 'cellList', 'cellListN', 'shift_xy');
    
    disp(['shifted ' num2str(shift_xy(1)) ' in x and ' num2str(shift_xy(2)) ' in y.']);
end