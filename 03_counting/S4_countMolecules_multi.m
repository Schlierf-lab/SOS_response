% Molecules counting using the extracted optimal blinking tolerance time
% for multiple files
% --Andreas Hartmann
% --Leonard Schärfen

clear all;

%% parameters
[FILE,PATH]=uigetfile('*.loc.txt','Localization File', 'multiselect', 'on');

DT_FRAME=20; % (ms) time per frame
max_displacement=55; % (nm)
NM_PERPIXEL=106.7;
FIRST_FRAME=500; % first active frame
boolBrdRegion=0; % select only molecules in the border region
BRD_REGION=-350; % (nm) thickness of the border region for molecule selection, see boolBrdRegion

addpath(genpath('scripts'));

%% load tauC Fit
if ~iscell(FILE)
    file_count = 1;
    FILE = {FILE};
else
    file_count = length(FILE);
end

%[FILE_FIT,PATH_FIT]=uigetfile('*SIM_Fit.mat','tauC-Fit');
load([PATH FILE{1}(1:end-24) 'SIM_Fit.mat'],'params');

optTauC=@(x) params(1)*(x-params(2))^params(3)+params(4);

%% load SWIFT data

for movie = 1:file_count

    locs=importdata([PATH FILE{movie}]);

    molID=locs.data(1:end-1,end-1);
    frames=locs.data(1:end-1,1);
    posX=locs.data(1:end-1,13); % (nm)
    posY=locs.data(1:end-1,14); % (nm)
    cellID=locs.data(1:end-1,17);

    molID=molID(frames>FIRST_FRAME);
    frames=frames(frames>FIRST_FRAME);
    posX=posX(frames>FIRST_FRAME);
    posY=posY(frames>FIRST_FRAME);
    cellID=cellID(frames>FIRST_FRAME);

    %% removal of close localizations per frame < max_displacement

    h=waitbar(0,'Removing localization twins...');

    redSWIFTDATA=[];

    for iterF=unique(frames)'

        sub_molID=molID(frames==iterF);
        sub_posX=posX(frames==iterF);
        sub_posY=posY(frames==iterF);
        sub_frames=frames(frames==iterF);
        sub_cellID=cellID(frames==iterF);

        if length(sub_molID)>1

            currSEL=[sub_frames(1) sub_posX(1) sub_posY(1) sub_molID(1) sub_cellID(1)];

            for iterID=2:length(sub_molID)

                distR=sqrt((sub_posX(iterID)-currSEL(:,2)).^2+(sub_posY(iterID)-currSEL(:,3)).^2);

                if sum(distR<=max_displacement)==0

                    currSEL(end+1,:)=[sub_frames(iterID) sub_posX(iterID) sub_posY(iterID) sub_molID(iterID) sub_cellID(iterID)];
                end
            end
        else

            currSEL=[sub_frames sub_posX sub_posY sub_molID sub_cellID];
        end

        redSWIFTDATA=[redSWIFTDATA;currSEL];

        waitbar(find(unique(frames)==iterF,1)/length(unique(frames)),h);
    end

    close(h);

    %% calculation of mean localization positions

    frames=redSWIFTDATA(:,1);
    posX=redSWIFTDATA(:,2); % (nm)
    posY=redSWIFTDATA(:,3); % (nm)
    molID=redSWIFTDATA(:,4);

    locMean=[];

    h=waitbar(0,'Calculating mean positions...');

    for iterL=sort(unique(molID))'

        meanXL=mean(posX(molID==iterL));
        meanYL=mean(posY(molID==iterL));
        framesL=frames(molID==iterL); 

        locMean(end+1,:)=[iterL meanXL meanYL min(framesL) max(framesL)];

        waitbar(find(sort(unique(molID))'==iterL)/length(sort(unique(molID))'),h);
    end
    close(h);

    %% analysis of reduced SWIFT data

    finalSWIFTDATA=redSWIFTDATA;

    h=waitbar(0,'Linking sections...');

    for iterM=sort(unique(molID))'

        if sum(finalSWIFTDATA(:,4)==iterM)>0

            meanXC=locMean(locMean(:,1)==iterM,2);
            meanYC=locMean(locMean(:,1)==iterM,3);
            frameMinC=locMean(locMean(:,1)==iterM,4);
            frameMaxC=locMean(locMean(:,1)==iterM,5);

            compMolID=locMean(locMean(:,1)~=iterM,1);
            distPIX=sqrt((locMean(locMean(:,1)~=iterM,2)-meanXC).^2+(locMean(locMean(:,1)~=iterM,3)-meanYC).^2); % (nm)

            % searching for molecules within max_displacement and temporal
            sameMol=compMolID(distPIX<=max_displacement);

            if ~isempty(sameMol)

                % merging molecule IDs
                for iterC=1:length(sameMol)

                    finalSWIFTDATA(finalSWIFTDATA(:,4)==sameMol(iterC),4)=iterM;
                end

                % recalculation of molecule information
                framesN=finalSWIFTDATA(finalSWIFTDATA(:,4)==iterM,1);
                meanXN=mean(finalSWIFTDATA(finalSWIFTDATA(:,4)==iterM,2));
                meanYN=mean(finalSWIFTDATA(finalSWIFTDATA(:,4)==iterM,3));

                locMean(locMean(:,1)==iterM,:)=[iterM meanXN meanYN min(framesN) max(framesN)];
            end
        end

        waitbar(find(sort(unique(molID))'==iterM)/length(sort(unique(molID))'),h);
    end

    close(h);

    %% Selection of molecules in the border region

    % parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellborders_unCorr=matfile([PATH FILE{movie}(1:end-20) 'out.mat']);
    cellborders_Corr=matfile([PATH FILE{movie}(1:end-20) 'out_corr.mat']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    finalFrames=finalSWIFTDATA(:,1);
    finalPosX=finalSWIFTDATA(:,2); % (nm)
    finalPosY=finalSWIFTDATA(:,3); % (nm)
    finalMolID=finalSWIFTDATA(:,4);
    finalCellID=finalSWIFTDATA(:,5);

    if boolBrdRegion==1

        finalSWIFTDATA(:,6)=0; % boolean for molecules inside the border region

        h=waitbar(0,'Searching for molecules in the border region...');

        % load borders
        for iterI=sort(unique(finalCellID))'

            loaddata_Corr=cellborders_Corr.cellList;
            border=loaddata_Corr.meshData{1}{iterI}.model.*NM_PERPIXEL; % (nm)

            % calculate shift vector
            diffVect=[border(2:end,1)-border(1:end-1,1) border(2:end,2)-border(1:end-1,2) zeros(length(border(:,1))-1,1)];
            diffVect=[diffVect;border(1,1)-border(end,1) border(1,2)-border(end,2) 0];
            orthVect=cross(diffVect,[zeros(length(border(:,1)),1) zeros(length(border(:,1)),1) ones(length(border(:,1)),1)]);
            shiftVect=BRD_REGION.*orthVect./sqrt(sum(orthVect.^2,2));

            % new cell border
            new_border=[(border(:,1)+shiftVect(:,1)) (border(:,2)+shiftVect(:,2))];

            idx_Brd=inpolygon(locMean(:,2),locMean(:,3),border(:,1),border(:,2));
            idx_newBrd=inpolygon(locMean(:,2),locMean(:,3),new_border(:,1),new_border(:,2));

            sel_molID=locMean(idx_Brd&~idx_newBrd,1);

            for iterID=sel_molID'

                finalSWIFTDATA(finalMolID==iterID,6)=1;
            end

            waitbar(find(sort(unique(finalCellID))'==iterI)/length(sort(unique(finalCellID))'),h);
        end

        close(h);
    else

        finalSWIFTDATA(:,6)=1;
    end

    % % Show localizations in cells
    % finalBoolBrd=finalSWIFTDATA(:,6);
    % 
    % % load borders
    % for iterI=sort(unique(finalCellID))'
    % 
    %     loaddata_Corr=cellborders_Corr.cellList;
    %     border=loaddata_Corr.meshData{1}{iterI}.model.*NM_PERPIXEL; % (nm)
    % 
    %     % calculate shift vector
    %     diffVect=[border(2:end,1)-border(1:end-1,1) border(2:end,2)-border(1:end-1,2) zeros(length(border(:,1))-1,1)];
    %     diffVect=[diffVect;border(1,1)-border(end,1) border(1,2)-border(end,2) 0];
    %     orthVect=cross(diffVect,[zeros(length(border(:,1)),1) zeros(length(border(:,1)),1) ones(length(border(:,1)),1)]);
    %     shiftVect=BRD_REGION.*orthVect./sqrt(sum(orthVect.^2,2));
    % 
    %     % new cell border
    %     new_border=[(border(:,1)+shiftVect(:,1)) (border(:,2)+shiftVect(:,2))];
    %     
    %     figure();
    %     plot(border(:,1),border(:,2),'-k');
    %     hold on;
    %     
    %     if boolBrdRegion==1
    %         
    %         plot(new_border(:,1),new_border(:,2),'--k');
    %     end
    %         
    %     sel_molID=finalMolID(finalCellID==iterI&finalBoolBrd==1);
    %         
    %     for iterID=sel_molID'
    %         
    %         pp1=plot(finalPosX(finalMolID==iterID),finalPosY(finalMolID==iterID),'-');
    % %         pp2=plot(locMean(locMean(:,1)==iterID,2),locMean(locMean(:,1)==iterID,3),'+');
    % %         set(pp2,'Color',get(pp1,'Color'));
    %     end
    %     
    %     axis equal;
    % end


    %% Calculation of number of molecules

    finalFrames=finalSWIFTDATA(:,1);
    finalPosX=finalSWIFTDATA(:,2); % (nm)
    finalPosY=finalSWIFTDATA(:,3); % (nm)
    finalMolID=finalSWIFTDATA(:,4);
    finalCellID=finalSWIFTDATA(:,5);
    finalBoolBrd=finalSWIFTDATA(:,6);

    countArr=zeros(length(unique(finalMolID)),4);
    idxCt=1;

    for iterK=unique(finalMolID)'

        subframes=sort(finalFrames(finalMolID==iterK));
        subframes=unique(subframes);

        gaps=setxor(subframes,(subframes(1):1:subframes(end)));

        if ~isempty(gaps)

            [bStart bLength]=burstLoc(gaps,1);

            tauC=0;

            for iterNB=1:100

                Nb=sum(bLength*DT_FRAME>tauC)+1;
                tauC=optTauC(Nb);
            end

            Nmol=Nb;
        else

            Nmol=1;
        end

        countArr(idxCt,:)=[locMean(locMean(:,1)==iterK,2) locMean(locMean(:,1)==iterK,3) Nmol sum(finalBoolBrd(finalMolID==iterK))==sum(finalMolID==iterK)];

        idxCt=idxCt+1;
    end

    %% Show result

    % parameters
    %%%%%%%%%%%
    minN=1;
    maxN=10;
    %%%%%%%%%%%

    % load borders
    loaddata_unCorr=cellborders_unCorr.cellList;
    loaddata_Corr=cellborders_Corr.cellList;

    % plot Result
    selX=countArr(minN<countArr(:,3)&countArr(:,3)<maxN&countArr(:,4)==1,1);
    selY=countArr(minN<countArr(:,3)&countArr(:,3)<maxN&countArr(:,4)==1,2);
    selC=countArr(minN<countArr(:,3)&countArr(:,3)<maxN&countArr(:,4)==1,3);

    cmap=colormap(jet);
    scatter(selX./NM_PERPIXEL,selY./NM_PERPIXEL,selC.*10,cmap(round(selC/maxN.*64),:),'filled');
    colorbar;
    caxis([minN maxN]);
    colormap(jet);
    hold on;

    for iterO=1:length(loaddata_unCorr.meshData{:})

        border=loaddata_Corr.meshData{1}{iterO}.model;

        plot(border(:,1),border(:,2),'-k','LineWidth',2);
    end

    set(gca,'Box','on','Color',[0.75 0.75 0.75]);
    xlabel('X');
    ylabel('Y');
    axis equal;

    title(FILE(1:end-20));

    save([PATH FILE{movie}(1:end-20) 'count.mat'],'countArr');
end