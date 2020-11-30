% Extraction of the on-time, off-time, blinking time and bleaching time distributions
% --Andreas Hartmann

clear all;

%% parameters
PATH='Photophysics_Dendra\';
FILE='Dd2_fermi_30_001_filt.tracked.loc.txt';
OUT='Dd2_Result_tc_INF.mat';
DT_FRAME=20; % (ms) time per frame
max_displacement=125; % (nm)
tc=Inf; % (ms) cut-off time

addpath(genpath('scripts'));

%% load SWIFT data
locs=importdata([PATH FILE]);

molID=locs.data(1:end-1,end-1);
frames=locs.data(1:end-1,1);
pixX=locs.data(1:end-1,2); % (Pixel)
pixY=locs.data(1:end-1,3); % (Pixel)
posX=locs.data(1:end-1,13); % (nm)
posY=locs.data(1:end-1,14); % (nm)

%% Calculation of apparent kinetics

app_ton=[];
app_toff=[];
app_Nblink=[];
app_tbleach=[];

h=waitbar(0,'Calculating apparent fluorescence kinetics...');

for iterK=unique(molID)'
           
    subframes=sort(frames(molID==iterK));
    subframes=unique(subframes);
    
    app_tbleach=[app_tbleach;subframes(end)-subframes(1)+1];
    
    [bStart bLength]=burstLoc(subframes,1);
    app_ton=[app_ton;bLength];
    
    gaps=setxor(subframes,(subframes(1):1:subframes(end)));
    
    if ~isempty(gaps)
        
        [bStart bLength]=burstLoc(gaps,1);
        app_toff=[app_toff;bLength];
        
        app_Nblink=[app_Nblink;length(bStart)];
    else
        
        app_Nblink=[app_Nblink;0];
        app_toff=[app_toff;0];
    end
    
    waitbar(find(unique(molID)==iterK)/length(unique(molID)),h);
end

close(h);

%% removal of close localizations per frame < max_displacement

h=waitbar(0,'Removing localization twins...');

redSWIFTDATA=[];

for iterF=unique(frames)'
    
    sub_molID=molID(frames==iterF);
    sub_posX=posX(frames==iterF);
    sub_posY=posY(frames==iterF);
    sub_frames=frames(frames==iterF);
        
    if length(sub_molID)>1
        
        currSEL=[sub_frames(1) sub_posX(1) sub_posY(1) sub_molID(1)];
        
        for iterID=2:length(sub_molID)
            
            distR=sqrt((sub_posX(iterID)-currSEL(:,2)).^2+(sub_posY(iterID)-currSEL(:,3)).^2);
            
            if sum(distR<=max_displacement)==0
                
                currSEL(end+1,:)=[sub_frames(iterID) sub_posX(iterID) sub_posY(iterID) sub_molID(iterID)];
            end
        end
    else
        
        currSEL=[sub_frames sub_posX sub_posY sub_molID];
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

h=waitbar(0,'Linking molecules...');

for iterM=sort(unique(molID))'
    
    if sum(finalSWIFTDATA(:,4)==iterM)>0
    
        meanXC=locMean(locMean(:,1)==iterM,2);
        meanYC=locMean(locMean(:,1)==iterM,3);
        frameMinC=locMean(locMean(:,1)==iterM,4);
        frameMaxC=locMean(locMean(:,1)==iterM,5);

        compMolID=locMean(locMean(:,1)~=iterM,1);
        distPIX=sqrt((locMean(locMean(:,1)~=iterM,2)-meanXC).^2+(locMean(locMean(:,1)~=iterM,3)-meanYC).^2); % (nm)
        distT=min([abs(locMean(locMean(:,1)~=iterM,4)-frameMaxC) abs(locMean(locMean(:,1)~=iterM,5)-frameMinC)],[],2).*DT_FRAME; % (ms)

        % searching for molecules within max_displacement and temporal
        % distance
        sameMol=compMolID(distPIX<=max_displacement&distT<=tc);

        if ~isempty(sameMol)

            % merging molecule IDs
            for iterC=1:length(sameMol)

                finalSWIFTDATA(finalSWIFTDATA(:,4)==sameMol(iterC),4)=iterM;
            end

            % recalculation of molecule information
            meanXN=mean(finalSWIFTDATA(finalSWIFTDATA(:,4)==iterM,2));
            meanYN=mean(finalSWIFTDATA(finalSWIFTDATA(:,4)==iterM,3));
            framesN=finalSWIFTDATA(finalSWIFTDATA(:,4)==iterM,1);

            locMean(locMean(:,1)==iterM,:)=[iterM meanXN meanYN min(framesN) max(framesN)];
        end
    end
    
    waitbar(find(sort(unique(molID))'==iterM)/length(sort(unique(molID))'),h);
end

close(h);

%% Show Result

figure();
hold on;

posX=finalSWIFTDATA(:,2); % (nm)
posY=finalSWIFTDATA(:,3); % (nm)
molID=finalSWIFTDATA(:,4);

for iterP=unique(molID)'
    
    plot(posX(molID==iterP),posY(molID==iterP),'.-');
end

axis equal

%% Calculation of corrected kinetics

finalFrames=finalSWIFTDATA(:,1);
finalMolID=finalSWIFTDATA(:,4);

ton=[];
toff=[];
Nblink=[];
tbleach=[];
total_ton=[];
total_toff=[];

h=waitbar(0,'Calculating corrected fluorescence kinetics...');

for iterK=unique(finalMolID)'
           
    subframes=sort(finalFrames(finalMolID==iterK));
    subframes=unique(subframes);
    
    tbleach=[tbleach;subframes(end)-subframes(1)+1];
    
    [bStart bLength]=burstLoc(subframes,1);
    ton=[ton;bLength];
    total_ton=[total_ton;sum(bLength)];
    
    gaps=setxor(subframes,(subframes(1):1:subframes(end)));
    
    if ~isempty(gaps)
        
        [bStart bLength]=burstLoc(gaps,1);
        toff=[toff;bLength];
        total_toff=[total_toff;sum(bLength)];
        
        Nblink=[Nblink;length(bStart)];
    else
        
        Nblink=[Nblink;0];
        toff=[toff;0];
        total_toff=[total_toff;0];
    end
    
    waitbar(find(unique(finalMolID)==iterK)/length(unique(finalMolID)),h);
end

close(h);

%% Show results

binsize=1;
bincenterB=(0:binsize:20);
bincenter=(1:binsize:500);
edgesB=[bincenterB bincenterB(end)+1]-0.5*binsize;
edges=[bincenter bincenter(end)+1]-0.5*binsize;

% number of blinking events
hNblink=histcounts(Nblink,edgesB);
app_hNblink=histcounts(app_Nblink,edgesB);
pNblink=hNblink./sum(hNblink);
app_pNblink=app_hNblink./sum(app_hNblink);

[cfunNb,~,~]=fit(bincenterB',pNblink','(1-a)*a^x','Lower',0);
ciNB=confint(cfunNb);

figure();
bar(bincenterB,pNblink,'hist');
hold on;
plot(bincenterB,app_pNblink,'sk');
plot((0:0.01:100),cfunNb((0:0.01:100)),'-r','LineWidth',1.5);
xlabel('N_{blink}');
ylabel('Probability');
xlim([-1 7]);
legend('corrected Data','apparent Data',['\eta^n(1-\eta) -> \eta=' num2str(cfunNb.a) '\pm' num2str((ciNB(2)-ciNB(1))/4)]);

% on-time histogram
hTon=histcounts(ton,edges);
app_hTon=histcounts(app_ton,edges);
pTon=hTon./sum(hTon);
app_pTon=app_hTon./sum(app_hTon);

[cfunOn,~,~]=fit(bincenter',pTon','a*exp(-b*x)','Weights',1./(hTon+1));
ciOn=confint(cfunOn);

figure();
plot(bincenter.*DT_FRAME./1000,app_pTon./DT_FRAME,'sb','Color',[0.5 0.5 0.5]);
hold on;
plot(bincenter.*DT_FRAME./1000,pTon./DT_FRAME,'ob');
plot((0.5:0.01:50).*DT_FRAME./1000,cfunOn((0.5:0.01:50))./DT_FRAME,'-r','LineWidth',1.5);
xlabel('T_{on} (s)');
ylabel('PDF (s^{-1})');
title('A*exp(-(k_b+k_d)*x)');
legend('apparent Data','corrected Data',['k_b+k_d=(' num2str(cfunOn.b/DT_FRAME*1000) '\pm' num2str((ciOn(2,2)-ciOn(1,2))/4/DT_FRAME*1000) ')s^{-1}']);
xlim([0 0.4]);

% off-time histogram
hToff=histcounts(toff,edges);
app_hToff=histcounts(app_toff,edges);
pToff=hToff./sum(hToff);
app_pToff=app_hToff./sum(app_hToff);

% [cfunOff,~,~]=fit(bincenter',pToff','a*exp(-b*x)','Weights',1./(hToff+1));
[cfunOff2,~,~]=fit(bincenter',pToff','a*exp(-b*x)+c*exp(-d*x)','Weights',1./(hToff+1),'Start',[0.9/0.7 0.9 0.7*0.14 0.14]);
% ciOff=confint(cfunOff);
ciOff2=confint(cfunOff2);

figure();
plot(bincenter.*DT_FRAME./1000,app_pToff./DT_FRAME,'sb','Color',[0.5 0.5 0.5]);
hold on;
plot(bincenter.*DT_FRAME./1000,pToff./DT_FRAME,'ob');
% plot((0.5:0.01:50).*DT_FRAME,cfunOff((0.5:0.01:50))./DT_FRAME,'--r','LineWidth',1.5);
plot((0.5:0.01:50).*DT_FRAME./1000,cfunOff2((0.5:0.01:50))./DT_FRAME,'-r','LineWidth',1.5);
xlabel('T_{off} (s)');
ylabel('PDF (s^{-1})');
% legend('Data',['A*exp(-k_r*x) -> k_r=(' num2str(cfunOff.b) '\pm' num2str((ciOff(2,2)-ciOff(1,2))/4) ')ms^{-1}']);
title('(k_{r1}*exp(-k_{r1}*x)+\alphak_{r2}*exp(-k_{r2}*x))/(1+\alpha)');
legend('apparent Data','corrected Data',['k_{r1}=(' num2str(cfunOff2.b/DT_FRAME*1000) '\pm' num2str((ciOff2(2,2)-ciOff2(1,2))/4/DT_FRAME*1000) ')s^{-1} | k_{r2}=(' num2str(cfunOff2.d/DT_FRAME*1000) '\pm' num2str((ciOff2(2,4)-ciOff2(1,4))/4/DT_FRAME*1000) ')s^{-1} | \alpha=' num2str((cfunOff2.b-cfunOff2.a)/cfunOff2.a)]);
xlim([0 0.75]);

% bleaching time histogram
hTbleach=histcounts(tbleach,edges);
app_hTbleach=histcounts(app_tbleach,edges);
pTbleach=hTbleach./sum(hTbleach);
app_pTbleach=app_hTbleach./sum(app_hTbleach);

figure();
plot(bincenter.*DT_FRAME./1000,app_pTbleach./DT_FRAME,'-sb','Color',[0.5 0.5 0.5]);
hold on;
plot(bincenter.*DT_FRAME./1000,pTbleach./DT_FRAME,'-ob');
xlabel('T_{bleach} (s)');
ylabel('PDF (s^{-1})');
xlim([0 0.5]);
legend('apparent Data','corrected Data');

%% save data

output.textdata={'app_ton','app_toff','app_Nblink','app_tbleach','corr_ton','corr_toff','corr_Nblink','corr_tbleach','total_ton','total_toff'};
output.data{1}=app_ton;
output.data{2}=app_toff;
output.data{3}=app_Nblink;
output.data{4}=app_tbleach;
output.data{5}=ton;
output.data{6}=toff;
output.data{7}=Nblink;
output.data{8}=tbleach;
output.data{9}=total_ton;
output.data{10}=total_toff;

save([PATH OUT],'output');