%% INPUT

dt = 0.02; % (s)
maxTrackLen = 50;
maxNumSegments = 3;
lenSegment = 4; % localizations

boolsave = 0;
savepath = 'A:\Leo\AcrB_TolC\TolC\20ms\';
savefile = 'Kan10_060.mat';

%% load swift files

[files, path] = uigetfile('*filt.csv','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end


datall = struct('tracks', []);
for movie = 1:file_count
    in_filename = [path files{movie}];    
    tracks = importdata([in_filename(1:end-4) '.tracked.loc.txt'], ',', 1);
    datall(movie).tracks = tracks.data;
end


%% convert to nicer format
    
xcol = 13;
ycol = 14;
ncol = 18;

trajectories = [];
tracklength = [];
ncells = 0;
for rep = 1:length(datall)

    temptra = cell(max(datall(rep).tracks(:,ncol)),1);
    tempdisp = zeros(max(datall(rep).tracks(:,ncol)),1);
    templength = zeros(max(datall(rep).tracks(:,ncol)),1);
    for i = 1:max(datall(rep).tracks(:,ncol))
        if sum(datall(rep).tracks(:,ncol)==i) >= 2
            x_track = datall(rep).tracks(datall(rep).tracks(:,ncol)==i, xcol);
            y_track = datall(rep).tracks(datall(rep).tracks(:,ncol)==i, ycol);
            t_frame = datall(rep).tracks(datall(rep).tracks(:,ncol)==i, 1);
            t_time  = [0; cumsum(diff(t_frame))].*dt;
            temptra(i) = mat2cell([x_track, y_track, t_time], length(x_track));
            templength(i) = length(x_track);
        end        
    end
    temptra = temptra(~cellfun('isempty',temptra(:,1)),:);
    templength = nonzeros(templength);
    trajectories = [trajectories; temptra];
    tracklength = [tracklength; templength];
    ncells = ncells + max(datall(rep).tracks(:,17));
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
            arrD(end+1)=mean(dy./(4.*dx));
        end
    end
end


%% plot 

% plot D*
binsize = 0.1;
binsize_MJD = 5;
 
eappD=(-5:binsize:1.5);
happD=histc(log10(arrD),eappD);
figure('Position' ,[356.2000  533.8000  196.0000  125.6000]); 
bar(eappD+0.5*binsize,happD./(sum(happD)),1,'linestyle', 'none');
%stairs(eappD,happD./(sum(happD)));
xlim([-3.5 1.0])
ylim([0 0.09])
xlabel('log_{10}(D*) (µm^2s^{-1})')
%ylabel('normalized amplitude');
title(['n_{tracks} = ' num2str(ntracks) ', n_{cells} = ' num2str(ncells)])

MJDw = mean_disp(:,1);
for i = 1:length(mean_disp)
    if mean_disp(i,2) > 1
        weight = repmat(mean_disp(i,1), mean_disp(i,2)-1, 1);
        MJDw = [MJDw; weight];
    end
end

figure('position', [356.2000  533.8000  196.0000  125.6000]);
hold on;
edges = (0:binsize_MJD:700);
hMJD = histc(MJDw, edges);

stairs(edges, hMJD./sum(hMJD))
xlim([0 200]);
ylim([0 0.14]);

xlabel('Mean Displacement (nm)')
ylabel('frequency')
title(['weighted, n = ' num2str(length(mean_disp))])

%plot track length
etracklength = (0:1:100);
htracklength = histc(tracklength, etracklength);
figure; 
bar(etracklength,htracklength);
xlabel('track length (localizations)')
ylabel('# of tracks')
xlim([0 100])
text(0.8, 0.8, ['median: ' num2str(median(tracklength),'%.1f') ' locs'], 'Units', 'normalized')
set(gca, 'yscale', 'log')


% plot localization precision
% tempall = vertcat(datall(:).tracks);
% figure;
% histogram(tempall(:,16), 300, 'BinLimits', [0 300])
% xlim([0 100])
% xlabel('Cramér Rao lower bound (nm)')
% text(0.8, 0.8, [num2str(mean(tempall(:,16)),'%.1f') ' nm'], 'Units', 'normalized')

%% save results
if boolsave
    save([savepath savefile], 'arrD', 'ntracks', 'ncells', 'trajectories', 'mean_disp');
end
