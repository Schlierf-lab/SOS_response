%% INPUT

dt = 0.02; % (s)
maxTrackLen = 50;
maxNumSegments = 3;
lenSegment = 4; % localizations
um_px = 0.1067;

%% load swift files

[files, path] = uigetfile('*filt.tracked.loc.txt','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end


%% convert to nicer format
trajectories = [];
tracklength = [];
ncells = 0;
all_params = [];
models = {};
for rep = 1:file_count
    in_filename = [path files{rep}];    
    tracks = importdata(in_filename, ',', 1);
    tracks = tracks.data;
    load([in_filename(1:end-20) 'out_corr.mat']);

    temptra = cell(max(tracks(:,18)),1);
    tempdisp = zeros(max(tracks(:,18)),1);
    templength = zeros(max(tracks(:,18)),1);
    
    for i = 1:max(tracks(:,18))
        if sum(tracks(:,18)==i) >= 2
            x_track = tracks(tracks(:,18)==i, 13);
            y_track = tracks(tracks(:,18)==i, 14);
            t_frame = tracks(tracks(:,18)==i, 1);
            t_time  = [0; cumsum(diff(t_frame))].*dt;
            cell_id = tracks(tracks(:,18)==i, 17) + ncells;
            temptra(i) = mat2cell([x_track, y_track, t_time, cell_id], length(x_track));
        end        
    end
    temptra = temptra(~cellfun('isempty',temptra(:,1)),:);
    templength = nonzeros(templength);
    trajectories = [trajectories; temptra];
    tracklength = [tracklength; templength];
    
    params = zeros(cellListN, 4); % col1: not dividing (1) or dividing (0), col2: included or excluded, col3: length, col4: area, col5: contour length, col6: movieNR
    for i = 1:cellListN
        model = cellList.meshData{1, 1}{1, i}.model;
        mesh = cellList.meshData{1, 1}{1, i}.mesh;        
        if mesh(1) ~= 0
            midx = mesh(:,1) + (mesh(:,3)-mesh(:,1)) ./ 2;
            midy = mesh(:,2) + (mesh(:,4)-mesh(:,2)) ./ 2;
            d = diff([midx, midy]);
            e = diff([model(:,2), model(:,2)]);
            celllength = sum(sqrt(sum(d.*d,2)))*um_px;
            contlength = sum(sqrt(sum(e.*e,2)))*um_px;
            area = polyarea(model(:,1).*um_px, model(:,2).*um_px);
            cellNr = i + ncells;
            params(i,:) = [cellNr celllength contlength area];
            models{cellNr} = model;
        end
    end
    all_params = [all_params; params];
    ncells = ncells + max(tracks(:,17));
end

%% extract D*

arrD=[];
mean_disp = [];
ntracks = 0;
for iter=1:length(trajectories)

    lenData=size(trajectories{iter},1);

    if (lenSegment<=lenData)&&(lenData<=maxTrackLen)
        ntracks = ntracks + 1;
        
        % MJD calculation
        displacement = zeros(lenData-1, 1);
        for step = 2:lenData
            displacement(step-1) = pdist([trajectories{iter}(step-1, 1) trajectories{iter}(step-1, 2); trajectories{iter}(step, 1) trajectories{iter}(step, 2)], 'euclidean');
        end
        displacement = displacement./(diff(trajectories{iter}(:,3)./dt));
        mean_disp(end+1, 1) = mean(displacement);
        mean_disp(end, 2) = length(displacement);
        
        % appD calculation
        numPossible=floor(lenData/lenSegment);
        numPossible=min([numPossible maxNumSegments]);

        for iterD=1:numPossible

            idx1=(iterD-1)*lenSegment+1;
            idx2=iterD*lenSegment;

            MSD=(trajectories{iter}(idx1:idx2,1)-trajectories{iter}(idx1,1)).^2+(trajectories{iter}(idx1:idx2,2)-trajectories{iter}(idx1,2)).^2;

            dy=MSD(2:end).*1e-6; % (µm^2)
            dx=trajectories{iter}(idx1+1:idx2,3)-trajectories{iter}(idx1,3); % (s)
            arrD(end+1,:)=[mean(dy./(4.*dx)) trajectories{iter}(1,4)];
        end
    end
end

%% single-cell

arrD_mean = zeros(max(arrD(:,2)), 2);
for s = 1:max(arrD(:,2))
    arrD_mean(s,:) = [mean(arrD(arrD(:,2) == s, 1)) length(arrD(arrD(:,2) == s, 1))];
end

%% plot 
minnroftracks = 10;
arrD_mean_a = arrD_mean(arrD_mean(:,2) >= minnroftracks,:);
cell_length = all_params(arrD_mean(:,2) >= minnroftracks,2);
models_a = models(arrD_mean(:,2) >= minnroftracks);

[~, order] = sort(arrD_mean_a(:,1));
% [~, order] = sort(cell_length);
arrD_mean_a = arrD_mean_a(order,:);
cell_length = cell_length(order);
models_a = models_a(order);

figure;
scatter(cell_length, arrD_mean_a(:,1));

figure;
histogram(arrD_mean_a(:,1))
xlabel('mean D* per cell')

rows = ceil(sqrt(length(models_a)));
columns = ceil(sqrt(length(models_a)));
spacing = 40;
% create grid for figure
ro = fliplr(1:rows);
co = 1:columns;
grid = combvec(co,ro)'.*spacing;

% assemble the figure
figure('renderer', 'painters'); hold on; axis equal;
map = flipud(brewermap(256, 'spectral'));
for k = 1:length(models_a)
    if ~isempty(models_a{k})
        xy = models_a{k};
        CoM = [mean([max(xy(:,1)), min(xy(:,1))]), mean([max(xy(:,2)), min(xy(:,2))])];
        xy_grid = xy - CoM + grid(k,:);
        color = map(1+round((arrD_mean_a(k,1)-min(arrD_mean_a(:,1)))/(max(arrD_mean_a(:,1))-min(arrD_mean_a(:,1)))*255),:);
        plot(xy_grid(:,1)*um_px, xy_grid(:,2)*um_px, 'color', color, 'LineWidth', 1.5)
    end
end
set(gca, 'color', 'none');
colormap(map);
bar = colorbar;
title(bar, 'D* [nm]')
