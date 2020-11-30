% Converts localization table from SMAP to a swift input file, including a new column for adding the cell identity.
% Also corrects the offset of the outlines from oufti.
% --Leonard Schärfen

% parameters
nm_px = 106.7;
firstFrame = 101;


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
    SMAP_dat = SMAP.data;
    SMAP_txt = SMAP.textdata;

    
    if ~exist([file_location(1:end-4) '_out_corr.mat'], 'file')
        load([file_location(1:end-4) '_out.mat'], 'cellList', 'cellListN');
        [cellList, shift_xy] = optCellBorder(SMAP_dat, cellList, cellListN, -4, 4);    
        save([file_location(1:end-4) '_out_corr.mat'], 'cellList', 'cellListN', 'shift_xy');
        disp(['shifted ' num2str(shift_xy(1)) ' in x and ' num2str(shift_xy(2)) ' in y.']);
    else
        load([file_location(1:end-4) '_out_corr.mat'], 'cellList', 'cellListN');
    end
    % add cellID for each localization, remove outside localizations

    SMAP_txt{17} = 'cell_id';
    SMAP_dat(:,17) = 0;
    for i = 1:cellListN
        if ~isempty(cellList.meshData{1,1}{1,i}.model)% && cellList.selected2(i)
            [in, on] = inpolygon(SMAP_dat(:,13), SMAP_dat(:,14), nm_px.*cellList.meshData{1,1}{1,i}.model(:,1), nm_px.*cellList.meshData{1,1}{1,i}.model(:,2));
            SMAP_dat(in|on,17) = i;
        end
    end
    SMAP_dat(SMAP_dat(:,17)==0, :) = [];
    SMAP_dat(SMAP_dat(:,1)<firstFrame, :) = [];

    % rename relevant variables for swift and save the new .csv

    SMAP_txt{13} = 'X';
    SMAP_txt{14} = 'Y';
    SMAP_txt{2} = 'x_pix';
    SMAP_txt{3} = 'y_pix';
    SMAP_txt{5} = 'Intensity';

    fid = fopen([file_location(1:end-4) '_filt.csv'], 'w');
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', SMAP_txt{:});
    for k=2:length(SMAP_dat)
        fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',SMAP_dat(k,:));
    end
    fclose(fid);
end
disp('Done.');