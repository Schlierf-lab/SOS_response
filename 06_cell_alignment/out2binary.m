% converts multiple outline files to a stack of cell-indexed binaries,
% filters for selector1 and 2 and deletes overlapping cells

[files, path] = uigetfile('*out_corr.mat','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

for BFnr = 1:file_count
    in_filename = [path files{BFnr}];
    
    raw_outl = load(in_filename);
    outlines = raw_outl.cellList.meshData{1, 1};
    div = raw_outl.cellList.selected1;
    excl = raw_outl.cellList.selected2;
    
    binary = zeros(512, 512);

    iter = 0;
    for i = 1:length(outlines)
        if ~isempty(outlines{1, i}.model) && div(i) && excl(i)
            cell = roipoly(binary, outlines{1, i}.model(:,1), outlines{1, i}.model(:,2));
            if binary(cell) == 0
                iter = iter + 1;
                binary(cell) = iter;
            end
        end
    end
   
    binary(binary>0) = 1./binary(binary>0);
    imwrite(binary, [in_filename(1:end-16) 'BN.tif'], 'writemode', 'append');
    
    binary = zeros(512, 512);

    iter = 0;
    for i = 1:length(outlines)
        if ~isempty(outlines{1, i}.model) && ~div(i) && excl(i)
            cell = roipoly(binary, outlines{1, i}.model(:,1), outlines{1, i}.model(:,2));
            if binary(cell) == 0
                iter = iter + 1;
                binary(cell) = iter;
            end
        end
    end
   
    binary(binary>0) = 1./binary(binary>0);
    imwrite(binary, [in_filename(1:end-16) 'BN_div.tif'], 'writemode', 'append');
    
end
