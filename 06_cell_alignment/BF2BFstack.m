[files, path] = uigetfile('*BF.tif','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

for BFnr = 1:file_count
    in_filename = [path files{BFnr}];
    
    load([in_filename(1:end-6) 'out_corr.mat'], 'shift_xy');
    
    im = imread(in_filename);
    im = imtranslate(im, shift_xy);
    imwrite(im, [in_filename(1:end-10) 'BF.tif'], 'writemode', 'append');
end