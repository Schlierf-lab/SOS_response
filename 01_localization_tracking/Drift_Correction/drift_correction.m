% Performs drift correction based on fiducial markers. Functions are taken 
% from PALMsiever (https://github.com/PALMsiever/palm-siever).
% --Leonard Schärfen, 4/27/2019

%%%INPUT%%%
noSplines = 60; 
firstframe = 501;
displayframe1 = 450;
saveIT = 1;
nm_px = 106.7;
dt = 0.02;
%%%INPUT%%%

addpath('bfmatlab');
[files, path] = uigetfile('*.csv', 'Select localization table...', 'MultiSelect','on');

% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

for movie = 1:file_count
      
    raw_output = importdata([path files{movie}]);
    Frame = raw_output.data(:,1);
    XPosition = raw_output.data(:,13);
    YPosition = raw_output.data(:,14);
    photons = raw_output.data(:,5);
    displayframe2 = max(Frame)-100;

    reader = bfGetReader([path files{movie}(1:end-3) 'nd2']);
    fluo1 = bfGetPlane(reader, displayframe1);
    fluo2 = bfGetPlane(reader, displayframe2);

    sel = Frame >= firstframe;

    f = figure('units','normalized','outerposition',[0 0 1 1]);
    s1 = subplot(2,2,2);
    imshow(fluo1, [])
    set(gca, 'YDir', 'normal');
    title(['frame ' num2str(displayframe1)])

    s2 = subplot(2,2,4);
    imshow(fluo2, [])
    set(gca, 'YDir', 'normal');
    title(['frame ' num2str(displayframe2)])
    colormap('jet');

    s3 = subplot(2,2,[1 3]);
    hold on
    axis equal
    axis off
    dscatter(XPosition(sel), YPosition(sel));
    title('select fiducial ROI...')

    roiList = {};
    proceed = true;
    % select fiducials: press 'return' to add another one, press 'space' to end
    while proceed
        roi = drawrectangle();
        roiList = [roiList, roi.Position];
        title('\rmpress \bfreturn \rmto add another fiducial, press \bfspace \rmto finish')
        pause; % wait for a keypress
        currkey=get(gcf,'CurrentKey'); 
        if strcmp(currkey, 'return')
            proceed=1;
            title('select fiducial ROI...')
        else
            proceed=0;
        end
    end
    close(f);

    
    driftTracks = getDriftTracks(roiList, XPosition, YPosition, Frame);
    [xc,yc,tDrift,xDrift,yDrift] = correctDrift(driftTracks,noSplines,Frame, XPosition,YPosition);

    % remove fiducials from localization table
    kickout = zeros(length(raw_output.data),1);
    for jj = 1:length(roiList)
        xmin = roiList{jj}(1);
        ymin = roiList{jj}(2);
        xmax = xmin + roiList{jj}(3);
        ymax = ymin + roiList{jj}(4);         
        roiout = inpolygon(raw_output.data(:,13), raw_output.data(:,14), [xmin xmin xmax xmax xmin], [ymin ymax ymax ymin ymin]);
        kickout = roiout | kickout;
    end
    
    % save results
    if saveIT
        movefile([path files{movie}], [path files{movie}(1:end-4) '_raw.csv'])
        disp('saving...');
        raw_output.data(:,13) = xc;
        raw_output.data(:,14) = yc;
        raw_output.data(:,2) = xc./nm_px;
        raw_output.data(:,3) = yc./nm_px;
        raw_output.data = raw_output.data(~kickout, :);
        fid = fopen([path files{movie}], 'w');
        fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', raw_output.textdata{:});
        for k=1:length(raw_output.data)
            fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',raw_output.data(k,:));
        end
        fclose(fid);
    end

    % plot
    CM = lines(numel(driftTracks));
    figure();
    legendCell={};
    hold on;
    for ii =1: numel(driftTracks)
       iiStr = num2str(ii);
       %X
       plot(driftTracks{ii}.t.*dt,driftTracks{ii}.x-driftTracks{ii}.x(1),'-', 'color', CM(ii,:));
       legendCell = [legendCell,['X ',iiStr]];
       %Y
       plot(driftTracks{ii}.t.*dt,driftTracks{ii}.y-driftTracks{ii}.y(1),'-', 'color', CM(ii,:));
       legendCell = [legendCell,['Y ',iiStr]];
    end

    plot(tDrift.*dt,xDrift,'r', 'Linewidth', 2);
    plot(tDrift.*dt,yDrift,'g', 'Linewidth', 2);
    legendCell = [legendCell,'X drift calculated','Y drift calculated'];

    legend(legendCell);
    xlabel('Time (s)');
    ylabel('Drift (nm)');
    title(files{movie}(1:end-4), 'interpreter', 'none');

    figure; hold on
    axis equal
    scatter(XPosition(sel), YPosition(sel), 10, 'filled')
    scatter(xc(sel & ~kickout), yc(sel & ~kickout), 10, 'filled')
    legend('uncorrected', 'corrected')

    % figure();
    % imshow(fluo1, []);
    % set(gca, 'YDir', 'normal');
    % hold on
    % for p = 1:numel(roiList)
    %     xmin = roiList{p}(1)/nm_px;
    %     ymin = roiList{p}(2)/nm_px;
    %     width = roiList{p}(3)/nm_px;
    %     height = roiList{p}(4)/nm_px;
    %     xmax = xmin + width;
    %     ymax = ymin + height;
    %     plot([xmin xmax], [ymin ymin], '-r', 'Linewidth', 1.5);
    %     plot([xmin xmin], [ymax ymin], '-r', 'Linewidth', 1.5);
    %     plot([xmin xmax], [ymax ymax], '-r', 'Linewidth', 1.5);
    %     plot([xmax xmax], [ymax ymin], '-r', 'Linewidth', 1.5);
    %     text(xmin, ymin-9, ['ROI ' num2str(p)], 'color', 'r', 'fontsize', 8);
    % end
end

disp('Done.');



function driftTracks = getDriftTracks(roiList, XPosition, YPosition, Frame)
disp('extracting traces...');
nTrack = numel(roiList);
driftTracks = cell(nTrack,1);
for ii = 1:nTrack
   xyLim = roiList{ii};
   driftTracks{ii} = getTrack(XPosition, YPosition, Frame, xyLim);
end
end

function driftTrack = getTrack(XPosition, YPosition, Frame, xyLim)

xmin = xyLim(1);
ymin = xyLim(2);
width = xyLim(3);
height = xyLim(4);
xmax = xmin + width;
ymax = ymin + height;

isTrack = XPosition>xmin & XPosition<xmax & YPosition>ymin & YPosition<ymax;
driftTrack.x = XPosition(isTrack);
driftTrack.y = YPosition(isTrack);
driftTrack.t = Frame(isTrack);

driftTrack = makeSingleValued(driftTrack);
end

function driftTrackSV = makeSingleValued(driftTrack)
%delete frames with multiple values

isBadRow =zeros(size(driftTrack.t));
tUn = unique(driftTrack.t);
for ii = 1:numel(tUn)
   isCurFrame = (driftTrack.t==tUn(ii));
   if sum(isCurFrame)>1
      isBadRow(isCurFrame) = 1;
   end
end

driftTrackSV.x = driftTrack.x(~isBadRow);
driftTrackSV.y = driftTrack.y(~isBadRow);
driftTrackSV.t = driftTrack.t(~isBadRow);
end