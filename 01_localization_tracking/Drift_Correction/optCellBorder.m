% Corrects brightfield - PALM movie offset by optimizing the number of localizations falling inside cells.
% --Andreas Hartmann

function [cellList, shift] = optCellBorder(localizations,cellList,cellListN,minShift,maxShift)
 
    % minShift, maxShift ... Pixel range for shifting the borders in x,y
    % direction
    % cellListN ... number of cells
    % cellList ... strcuture of cell borders

    % load localization data
    iniframe=localizations(:,1); % frame number of location i
    inix_Pix=localizations(:,2); % (Pixel) x adress of location i
    iniy_Pix=localizations(:,3); % (Pixel) y adress of location i

    arrOpt=[];

    hW=waitbar(0,'Identification of cell localizations...');

    arrShift=(minShift:0.25:maxShift);

    for iterX=1:length(arrShift)

        for iterY=1:length(arrShift)

            currN=0;
            for iterC=1:cellListN
                if ~isempty(cellList.meshData{1,1}{1,iterC}.model)
                    currN=currN+sum(inpolygon(inix_Pix,iniy_Pix,cellList.meshData{1,1}{1,iterC}.model(:,1)+arrShift(iterX),cellList.meshData{1,1}{1,iterC}.model(:,2)+arrShift(iterY)));
                end
            end

            arrOpt=[arrOpt;arrShift(iterX) arrShift(iterY) currN];
        end

        waitbar(iterX/length(arrShift));
    end

    close(hW);

    [~,indOpt]=max(arrOpt(:,3));

    opt_shiftX=arrOpt(indOpt,1);
    opt_shiftY=arrOpt(indOpt,2);
    
    shift = [opt_shiftX, opt_shiftY];
    
    for iterC=1:cellListN
        if ~isempty(cellList.meshData{1,1}{1,iterC}.model)
            cellList.meshData{1,1}{1,iterC}.model(:,1)=cellList.meshData{1,1}{1,iterC}.model(:,1)+opt_shiftX;
            cellList.meshData{1,1}{1,iterC}.model(:,2)=cellList.meshData{1,1}{1,iterC}.model(:,2)+opt_shiftY;
        end
    end
end