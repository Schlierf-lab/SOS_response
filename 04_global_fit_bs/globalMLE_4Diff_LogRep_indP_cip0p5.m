function globalMLE_4Diff_LogRep_indP_cip0p5()

    %% parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numSteps=3; % number of jumping steps for fit equation
    boolShow=1; % show fit result
    numbins=90; % number of histogram bins
    cond_x=[0 10 25 45 60 90 120 180]; % x-Axis of conditions
    
    % diffusion coefficients
    optD1=0.024714; % (µm^2/s)
    optD2=0.15202; % (µm^2/s)
    optD3=0.67744; % (µm^2/s)
    optD4=2.6645; % (µm^2/s)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    addpath('script');

    path='E:\04_LexA_project\00_live-cell_imaging\results\';
    %% load data   
    cond1_1=load([path 'results_cip20\cip00_000_1.mat'],'arrD');
    cond1_2=load([path 'results_cip20\cip00_000_2.mat'],'arrD');
    cond1_3=load([path 'results_cip20\cip00_000_3.mat'],'arrD');
    cond1_4=load([path 'results_cip20\cip00_000_4.mat'],'arrD');
    cond1_5=load([path 'results_cip20\cip00_000_5.mat'],'arrD');
    cond1_6=load([path 'results_cip20\cip00_000_6.mat'],'arrD');
    
    cond{1,1}=cond1_1.arrD;
    cond{1,2}=cond1_2.arrD;
    cond{1,3}=cond1_3.arrD;
    cond{1,4}=cond1_4.arrD;
    cond{1,5}=cond1_5.arrD;
    cond{1,6}=cond1_6.arrD;  
    numSub{1}=6;
    
    cond2_1=load([path 'results_cip05\cip05_010_1.mat'],'arrD');
    cond2_2=load([path 'results_cip05\cip05_010_2.mat'],'arrD');
    cond2_3=load([path 'results_cip05\cip05_010_3.mat'],'arrD');
    cond2_4=load([path 'results_cip05\cip05_010_4.mat'],'arrD');
    
    cond{2,1}=cond2_1.arrD;
    cond{2,2}=cond2_2.arrD;
    cond{2,3}=cond2_3.arrD;    
    cond{2,4}=cond2_4.arrD; 
    numSub{2}=4;

    cond3_1=load([path 'results_cip05\cip05_025_1.mat'],'arrD');
    cond3_2=load([path 'results_cip05\cip05_025_2.mat'],'arrD');
    cond3_3=load([path 'results_cip05\cip05_025_3.mat'],'arrD');    
        
    cond{3,1}=cond3_1.arrD;
    cond{3,2}=cond3_2.arrD;
    cond{3,3}=cond3_3.arrD;  
    numSub{3}=3;
    
    cond4_1=load([path 'results_cip05\cip05_045_1.mat'],'arrD');
    cond4_2=load([path 'results_cip05\cip05_045_2.mat'],'arrD');
    cond4_3=load([path 'results_cip05\cip05_045_3.mat'],'arrD');    
    
    cond{4,1}=cond4_1.arrD;
    cond{4,2}=cond4_2.arrD;
    cond{4,3}=cond4_3.arrD;
    numSub{4}=3;

    cond5_1=load([path 'results_cip05\cip05_060_1.mat'],'arrD');
    cond5_2=load([path 'results_cip05\cip05_060_2.mat'],'arrD');
    cond5_3=load([path 'results_cip05\cip05_060_3.mat'],'arrD');
    cond5_4=load([path 'results_cip05\cip05_060_4.mat'],'arrD');
    cond5_5=load([path 'results_cip05\cip05_060_5.mat'],'arrD');
    cond5_6=load([path 'results_cip05\cip05_060_6.mat'],'arrD');
    
    cond{5,1}=cond5_1.arrD;
    cond{5,2}=cond5_2.arrD;
    cond{5,3}=cond5_3.arrD; 
    cond{5,4}=cond5_4.arrD;
    cond{5,5}=cond5_5.arrD;
    cond{5,6}=cond5_6.arrD; 
    numSub{5}=6;
    
    cond6_1=load([path 'results_cip05\cip05_090_1.mat'],'arrD');
    cond6_2=load([path 'results_cip05\cip05_090_2.mat'],'arrD');
    cond6_3=load([path 'results_cip05\cip05_090_3.mat'],'arrD');    
    cond6_4=load([path 'results_cip05\cip05_090_4.mat'],'arrD');    
        
    cond{6,1}=cond6_1.arrD;
    cond{6,2}=cond6_2.arrD;
    cond{6,3}=cond6_3.arrD;    
    cond{6,4}=cond6_4.arrD; 
    numSub{6}=4;

    cond7_1=load([path 'results_cip05\cip05_120_1.mat'],'arrD');
    cond7_2=load([path 'results_cip05\cip05_120_2.mat'],'arrD');
    cond7_3=load([path 'results_cip05\cip05_120_3.mat'],'arrD');    
        
    cond{7,1}=cond7_1.arrD;
    cond{7,2}=cond7_2.arrD;
    cond{7,3}=cond7_3.arrD;
    numSub{7}=3;
    
    cond8_1=load([path 'results_cip05\cip05_180_1.mat'],'arrD');
    cond8_2=load([path 'results_cip05\cip05_180_2.mat'],'arrD');
    cond8_3=load([path 'results_cip05\cip05_180_3.mat'],'arrD');    
    
    cond{8,1}=cond8_1.arrD;
    cond{8,2}=cond8_2.arrD;
    cond{8,3}=cond8_3.arrD;  
    numSub{8}=3;

    %% bootstrapping fits
    
    % fit settings
    binstart=-3;
    binend=2;
    binsize=(binend-binstart)/numbins;
    binsizeF=(binend-binstart)/10000;
    edges=(binstart:binsize:binend);
    edgesF=(binstart:binsizeF:binend);
    
    % NEW NEW NEW
    %%%%%%%%%%%%%%%%%%%%%%%
    singleFun=@(y,a,b) a.*(((numSteps/b)^numSteps)*log(10)/factorial(numSteps-1)).*(10.^(numSteps.*y)).*exp(-(10.^y).*(numSteps/b));
    fitFun=@(y,a,b,c,d,e,f,g) a.*(((numSteps/d)^numSteps)*log(10)/factorial(numSteps-1)).*(10.^(numSteps.*y)).*exp(-(10.^y).*(numSteps/d))+b.*(((numSteps/e)^numSteps)*log(10)/factorial(numSteps-1)).*(10.^(numSteps.*y)).*exp(-(10.^y).*(numSteps/e))+c.*(((numSteps/f)^numSteps)*log(10)/factorial(numSteps-1)).*(10.^(numSteps.*y)).*exp(-(10.^y).*(numSteps/f))+(1-a-b-c).*(((numSteps/g)^numSteps)*log(10)/factorial(numSteps-1)).*(10.^(numSteps.*y)).*exp(-(10.^y).*(numSteps/g));
    pdf_2Diff=@(y,a,b,c,d,e,f,g) fitFun(y,a,b,c,d,e,f,g)./trapz(edgesF,fitFun(edgesF,a,b,c,d,e,f,g));
    %%%%%%%%%%%%%%%%%%%%%%%
    % NEW NEW NEW
    
%     singleFun=@(x,a,b) a.*((numSteps./b).^numSteps).*((10.^x).^(numSteps-1)).*exp(-numSteps.*(10.^x)./b)./factorial(numSteps-1);
%     fitFun=@(x,a,b,c,d,e,f,g) (a.*((numSteps./d).^numSteps).*((10.^x).^(numSteps-1)).*exp(-numSteps.*(10.^x)./d)./factorial(numSteps-1)+b.*((numSteps./e).^numSteps).*((10.^x).^(numSteps-1)).*exp(-numSteps.*(10.^x)./e)./factorial(numSteps-1)+c.*((numSteps./f).^numSteps).*((10.^x).^(numSteps-1)).*exp(-numSteps.*(10.^x)./f)./factorial(numSteps-1)+(1-a-b-c).*((numSteps./g).^numSteps).*((10.^x).^(numSteps-1)).*exp(-numSteps.*(10.^x)./g)./factorial(numSteps-1));
%     pdf_2Diff=@(x,a,b,c,d,e,f,g) fitFun(x,a,b,c,d,e,f,g)./trapz(edgesF,fitFun(edgesF,a,b,c,d,e,f,g));

    options=optimset('TolFun',45000);
    
    % boundary conditions   
    x0P=[0.003;0.064;0.20];
    lbP=[0;0;0];
    ubP=[1;1;1];

    % bottstrapping iterations for P fitting
    arrP1=-ones(size(cond,1),max(cell2mat(numSub)));
    arrP2=-ones(size(cond,1),max(cell2mat(numSub)));
    arrP3=-ones(size(cond,1),max(cell2mat(numSub)));
    allP1=-ones(size(cond,1),1);
    allP2=-ones(size(cond,1),1);
    allP3=-ones(size(cond,1),1);

    stdP1=zeros(size(cond,1),1);
    stdP2=zeros(size(cond,1),1);
    stdP3=zeros(size(cond,1),1);
    stdP4=zeros(size(cond,1),1);
    
    arr_hD=zeros(size(cond,1),length(edges));
   
    for iterC=1:size(cond,1)
        
        condAll=[];
        
        for iterS=1:numSub{iterC}
            
            curr_cond=cell2mat(cond(iterC,iterS));
            curr_cond=curr_cond(curr_cond~=0);
            
            condAll=[condAll;curr_cond(:)];
            
            curr_estMLE=fminsearchbnd(@MLE_P,x0P,lbP,ubP,options);
            
            arrP1(iterC,iterS)=curr_estMLE(1);
            arrP2(iterC,iterS)=curr_estMLE(2);
            arrP3(iterC,iterS)=curr_estMLE(3);
        end
        
        curr_cond=condAll(condAll~=0); 
        arr_hD(iterC,:)=histc(log10(condAll(condAll~=0)),edges)'; 
        
        estMLEAll=fminsearchbnd(@MLE_P,x0P,lbP,ubP,options);
        
        allP1(iterC)=estMLEAll(1);
        allP2(iterC)=estMLEAll(2);
        allP3(iterC)=estMLEAll(3);

        stdP1(iterC)=std(arrP1(iterC,1:numSub{iterC}));
        stdP2(iterC)=std(arrP2(iterC,1:numSub{iterC}));
        stdP3(iterC)=std(arrP3(iterC,1:numSub{iterC}));
        stdP4(iterC)=std(1-arrP1(iterC,1:numSub{iterC})-arrP2(iterC,1:numSub{iterC})-arrP3(iterC,1:numSub{iterC}));       
        
        disp([num2str(iterC) '/' num2str(size(cond,1))]);
    end
   
    if boolShow == 1 % not normalized

        for iterC=1:size(cond,1)
             
            sumFitAll=trapz(edgesF,fitFun(edgesF,allP1(iterC),allP2(iterC),allP3(iterC),optD1,optD2,optD3,optD4));
            sumhDAll=trapz(edges,arr_hD(iterC,:));
            facAll=sumhDAll/sumFitAll;
            
            Fit1=singleFun(edgesF,allP1(iterC),optD1);
            upperFit1=singleFun(edgesF,allP1(iterC)+stdP1(iterC),optD1);
            lowerFit1=singleFun(edgesF,allP1(iterC)-stdP1(iterC),optD1);
            Fit2=singleFun(edgesF,allP2(iterC),optD2);
            upperFit2=singleFun(edgesF,allP2(iterC)+stdP2(iterC),optD2);
            lowerFit2=singleFun(edgesF,allP2(iterC)-stdP2(iterC),optD2);
            Fit3=singleFun(edgesF,allP3(iterC),optD3);
            upperFit3=singleFun(edgesF,allP3(iterC)+stdP3(iterC),optD3);
            lowerFit3=singleFun(edgesF,allP3(iterC)-stdP3(iterC),optD3);
            Fit4=singleFun(edgesF,(1-allP1(iterC)-allP2(iterC)-allP3(iterC)),optD4);
            upperFit4=singleFun(edgesF,(1-allP1(iterC)-allP2(iterC)-allP3(iterC))+stdP4(iterC),optD4);
            lowerFit4=singleFun(edgesF,(1-allP1(iterC)-allP2(iterC)-allP3(iterC))-stdP4(iterC),optD4);
            
            figure();
            bar(edges,arr_hD(iterC,:),'BarWidth', 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'LineStyle', 'none');
            hold on;
            f1=fill([edgesF fliplr(edgesF)],[upperFit1 fliplr(lowerFit1)].*facAll,[0.6400 0.0800 0.1800],'EdgeColor','none');
            alpha(f1,0.5);
            f2=fill([edgesF fliplr(edgesF)],[upperFit2 fliplr(lowerFit2)].*facAll,[0.2700 0.7300 0.3600],'EdgeColor','none');
            alpha(f2,0.5);
            f3=fill([edgesF fliplr(edgesF)],[upperFit3 fliplr(lowerFit3)].*facAll,[0.3100 0.6200 0.8900],'EdgeColor','none');
            alpha(f3,0.5);
            f4=fill([edgesF fliplr(edgesF)],[upperFit4 fliplr(lowerFit4)].*facAll,[0.9569 0.5882 0.1922],'EdgeColor','none');
            alpha(f4,0.5);
            plot(edgesF,Fit1.*facAll,'Color',[0.6400 0.0800 0.1800],'LineWidth',1.5);
            plot(edgesF,Fit2.*facAll,'Color',[0.2700 0.7300 0.3600],'LineWidth',1.5);
            plot(edgesF,Fit3.*facAll,'Color',[0.3100 0.6200 0.8900],'LineWidth',1.5);
            plot(edgesF,Fit4.*facAll,'Color',[0.9569 0.5882 0.1922],'LineWidth',1.5); 
            plot(edgesF,(Fit1+Fit2+Fit3+Fit4).*facAll,'-k','LineWidth',2);
            xlim([binstart  binend]);
            LimY=get(gca,'YLim');
            ylim([0 LimY(2)]);
            xlabel('D* [µm^2s^-^1]');
            ylabel('number of molecules');
            title(['Cond' num2str(iterC) ]);
                        
            text(0.7, 0.8, {['P_1=(' num2str(allP1(iterC),'%.3f') '\pm' num2str(std(arrP1(iterC,:)),'%.3f') ')'];['P_2=(' num2str(allP2(iterC),'%.3f') '\pm' num2str(std(arrP2(iterC,:)),'%.3f') ')'];['P_3=(' num2str(allP3(iterC),'%.3f') '\pm' num2str(std(arrP3(iterC,:)),'%.3f') ')'];['P_4=(' num2str(mean(1-allP1(iterC)-allP2(iterC)-allP3(iterC)),'%.3f') '\pm' num2str(std(1-arrP1(iterC,:)-arrP2(iterC,:)-arrP3(iterC,:)),'%.3f') ')']}, 'Units', 'normalized');            
        end
    end  
    
    figure();
    errorbar(cond_x,allP1,stdP1,'-o','Color',[0.6400 0.0800 0.1800],'LineWidth',1.5);
    hold on;
    errorbar(cond_x,allP2,stdP2,'-o','Color',[0.2700 0.7300 0.3600],'LineWidth',1.5);
    errorbar(cond_x,allP3,stdP3,'-o','Color',[0.3100 0.6200 0.8900],'LineWidth',1.5);
    errorbar(cond_x,1-allP1-allP2-allP3,stdP4,'-o','Color',[0.9569 0.5882 0.1922],'LineWidth',1.5);
    xlabel('Time (s)');
    ylabel('Probability');
    legend('P1','P2','P3','P4');  
           
    function logP=MLE_P(params)
        
        P1=params(1);
        P2=params(2);
        P3=params(3);
     
        if P1+P2+P3>=1            
            logP=NaN;
        else
            logP=-sum(log(pdf_2Diff(log10(curr_cond),P1,P2,P3,optD1,optD2,optD3,optD4)));           
        end         
    end
end