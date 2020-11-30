% press enter to go to next localization file

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

figure('units','normalized','outerposition',[0 0 1 1]);
colormap(brewermap(256, 'Spectral'));
xlim([0 512])
ylim([0 512])
for movie = 1:file_count
    
    % import .csv file from SMAP and outline file
    file_location = [path files{movie}];

    SMAP = importdata(file_location);
    SMAP_dat = SMAP.data(SMAP.data(:,1)>=firstFrame,:);
    SMAP_txt = SMAP.textdata;
    
    s = scatter(SMAP_dat(:,2), SMAP_dat(:,3), 10, SMAP_dat(:,1), 'filled');
    axis equal
    title([num2str(movie) '/' num2str(file_count) '   "' files{movie} '"'], 'interpreter', 'none')
    
    proceed = 1;
    while proceed
        pause; % wait for a keypress
        currkey=get(gcf,'CurrentKey'); 
        if strcmp(currkey, 'return')
            proceed = 0;
        else
            continue
            
        end
    end
    delete(s);
end
close(gcf);