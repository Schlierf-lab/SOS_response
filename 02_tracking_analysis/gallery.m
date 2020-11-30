function gallery(swiftdata,rows,columns,spacing,minlength,maxlength)
    
    % throw all tracks and MJDs together
    alltracks = [];
    ind = 0;
    add = 0;
    
    
    for i = 1:length(swiftdata)
        for j = 1:length(swiftdata(i).tracks)
            ind = ind + 1;
            if ~isempty(swiftdata(i).segs(swiftdata(i).segs(:,1)==swiftdata(i).tracks(j,18),3)) && sum(swiftdata(i).segs(:,1)==swiftdata(i).tracks(j,18)) == 1 && swiftdata(i).segs(swiftdata(i).segs(:,1)==swiftdata(i).tracks(j,18),5) >= minlength && swiftdata(i).segs(swiftdata(i).segs(:,1)==swiftdata(i).tracks(j,18),4) <= 100
                alltracks(ind,:) = [swiftdata(i).tracks(j,18)+add, swiftdata(i).tracks(j,13), swiftdata(i).tracks(j,14), swiftdata(i).segs(swiftdata(i).segs(:,1)==swiftdata(i).tracks(j,18),3)];
            end 
        end
        add = max(alltracks(:,1));
    end
    
    alltracks = alltracks(any(alltracks,2),:);
    
    % randomly select tracks
    trackcount = rows*columns;
    seltracks = zeros(trackcount, 2);
    j = 1;
    loopcount = 0;
    while j <= trackcount
        sel = round(rand(1)*max(alltracks(:,1)));
        loopcount = loopcount + 1;
        if ismember(sel, alltracks(:,1)) && ~ismember(sel, seltracks(:,1)) && sum(alltracks(:,1)==sel) >= minlength && sum(alltracks(:,1)==sel) <= maxlength
            seltracks(j,:) = [sel, mean(alltracks(alltracks(:,1)==sel,4))];
            j = j + 1;
        end
        if loopcount > trackcount + 100000
            error('not enough tracks')
        end
    end
    
    % sort by MJD
    seltracks = sortrows(seltracks, 2);
    
    % create grid for figure
    ro = fliplr(1:rows);
    co = 1:columns;
    grid = combvec(co,ro)'.*spacing;
    
    % assemble the figure
    figure('renderer', 'painters'); hold on; axis equal;
    title(['Gallery of ' num2str(trackcount) ' tracks']);
    ylabel('[nm]');
    xlabel('[nm]');
    %map = brewermap(256, 'RdYlGn');
    map = flipud(brewermap(256, 'RdBu'));
    for k = 1:trackcount
        xy = alltracks(alltracks(:,1)==seltracks(k),2:3);
        CoM = [mean([max(xy(:,1)), min(xy(:,1))]), mean([max(xy(:,2)), min(xy(:,2))])];
        xy_grid = xy - CoM + grid(k,:);
        color = map(1+round((seltracks(k,2)-min(seltracks(:,2)))/(max(seltracks(:,2))-min(seltracks(:,2)))*255),:);
        plot(xy_grid(:,1), xy_grid(:,2), 'color', color, 'LineWidth', 1.5)
    end
    set(gca, 'color', 'none');
    colormap(map);
    bar = colorbar;
    title(bar, 'MJD [nm]')
    bar.Ticks = [0 1];
    bar.TickLabels = num2cell([round(min(seltracks(:,2))), round(max(seltracks(:,2)))]);
    scalebar('ScaleLengthRatio', 0.1, 'Colour', [0 0 0], 'Bold', 1, 'Unit', 'nm', 'Location', 'southeast');
end
