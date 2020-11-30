function Fit_Tracking()    

    % - numerical infinitesimal calculation of populations
    % - logistic response function for new synthesis with global minimum 
    % rate ks01, individual response time ts, response speed w and 
    % response amplitude k02(i)>k01
    % - global initial amount of visible molecules v0
    % - without time point t=180min
    % - kxEff replaced by pf*kx
    % - fitted with Chi2
    %
    % Andreas Hartmann | 06.03.2020
    
    %% parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataT{1}=importdata('data\cip0p5_track.txt');
    dataT{2}=importdata('data\cip3_track.txt');
    dataT{3}=importdata('data\cip20_track.txt');
    dataC{1}=importdata('data\cip0p5.txt');
    dataC{2}=importdata('data\cip3.txt');
    dataC{3}=importdata('data\cip20_wo180.txt');
    fitData_cip0p5=load('data\ncip05_res.mat');
    fitData_cip3=load('data\ncip3_res.mat');
    fitData_cip20=load('data\ncip20_res.mat');

    ks01_opt=22.132;
    ks02_opt=[553.34363;830.42351;1527.0447];
    ts_opt=[102.1987;87.70705;56.06444];
    w_opt=[0.027723;0.039766;0.069146];
    kx02_opt=[3.5781;8.84254;11.1949];
    tx_opt=59.67; 
    kx01_opt=3.5781;

    aF=[0.1528;0.1417;0.1282];
    kF=[0.0322;0.0373;0.1001];
    pF0=0.0443; 
    
    detFrac=0.47933;
    dt=0.01; % (min)
    numBS=30; % number of data points for boot strapping    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    addpath('scripts');
    
    %% tracking data
    t=(0:dt:300); % (min)
    tdataT1=dataT{1}(:,1); % (min)
    tdataT2=dataT{2}(:,1); % (min)
    tdataT3=dataT{3}(:,1); % (min)  
    
    % CIP 0.5
    TdataT1=fitData_cip0p5.allP1;
    DdataT1=fitData_cip0p5.allP2;
    CdataT1=fitData_cip0p5.allP3;
    FdataT1=1-fitData_cip0p5.allP1-fitData_cip0p5.allP2-fitData_cip0p5.allP3;
    
    TdataT1Err=zeros(length(tdataT1),1);
    DdataT1Err=zeros(length(tdataT1),1);
    CdataT1Err=zeros(length(tdataT1),1);
    FdataT1Err=zeros(length(tdataT1),1);
    
    for iter=1:length(tdataT1)
        
        TdataT1Err(iter)=std(fitData_cip0p5.arrP1(iter,fitData_cip0p5.arrP1(iter,:)~=-1));
        DdataT1Err(iter)=std(fitData_cip0p5.arrP2(iter,fitData_cip0p5.arrP1(iter,:)~=-1));
        CdataT1Err(iter)=std(fitData_cip0p5.arrP3(iter,fitData_cip0p5.arrP1(iter,:)~=-1));
        FdataT1Err(iter)=std(fitData_cip0p5.arrP1(iter,fitData_cip0p5.arrP1(iter,:)~=-1)-fitData_cip0p5.arrP2(iter,fitData_cip0p5.arrP2(iter,:)~=-1)-fitData_cip0p5.arrP3(iter,fitData_cip0p5.arrP3(iter,:)~=-1));
    end    
    
    % CIP 3
    TdataT2=fitData_cip3.allP1;
    DdataT2=fitData_cip3.allP2;
    CdataT2=fitData_cip3.allP3;
    FdataT2=1-fitData_cip3.allP1-fitData_cip3.allP2-fitData_cip3.allP3;
    
    TdataT2Err=zeros(length(tdataT1),1);
    DdataT2Err=zeros(length(tdataT1),1);
    CdataT2Err=zeros(length(tdataT1),1);
    FdataT2Err=zeros(length(tdataT1),1);
    
    for iter=1:length(tdataT1)
        
        TdataT2Err(iter)=std(fitData_cip3.arrP1(iter,fitData_cip3.arrP1(iter,:)~=-1));
        DdataT2Err(iter)=std(fitData_cip3.arrP2(iter,fitData_cip3.arrP1(iter,:)~=-1));
        CdataT2Err(iter)=std(fitData_cip3.arrP3(iter,fitData_cip3.arrP1(iter,:)~=-1));
        FdataT2Err(iter)=std(fitData_cip3.arrP1(iter,fitData_cip3.arrP1(iter,:)~=-1)-fitData_cip3.arrP2(iter,fitData_cip3.arrP2(iter,:)~=-1)-fitData_cip3.arrP3(iter,fitData_cip3.arrP3(iter,:)~=-1));
    end
    
    % CIP 20
    TdataT3=fitData_cip20.allP1;
    DdataT3=fitData_cip20.allP2;
    CdataT3=fitData_cip20.allP3;
    FdataT3=1-fitData_cip20.allP1-fitData_cip20.allP2-fitData_cip20.allP3;
    
    TdataT3Err=zeros(length(tdataT1),1);
    DdataT3Err=zeros(length(tdataT1),1);
    CdataT3Err=zeros(length(tdataT1),1);
    FdataT3Err=zeros(length(tdataT1),1);
    
    for iter=1:length(tdataT1)
        
        TdataT3Err(iter)=std(fitData_cip20.arrP1(iter,fitData_cip20.arrP1(iter,:)~=-1));
        DdataT3Err(iter)=std(fitData_cip20.arrP2(iter,fitData_cip20.arrP1(iter,:)~=-1));
        CdataT3Err(iter)=std(fitData_cip20.arrP3(iter,fitData_cip20.arrP1(iter,:)~=-1));
        FdataT3Err(iter)=std(fitData_cip20.arrP1(iter,fitData_cip20.arrP1(iter,:)~=-1)-fitData_cip20.arrP2(iter,fitData_cip20.arrP2(iter,:)~=-1)-fitData_cip20.arrP3(iter,fitData_cip20.arrP3(iter,:)~=-1));
    end
    
    %% counting data
    tdataC1=dataC{1}(:,1); % (min)
    tdataC2=dataC{2}(:,1); % (min)
    tdataC3=dataC{3}(:,1); % (min) 
    
    VdataC1=dataC{1}(:,2)./detFrac;
    VdataC1Err=dataC{1}(:,3)./detFrac;
    VdataC2=dataC{2}(:,2)./detFrac;
    VdataC2Err=dataC{2}(:,3)./detFrac;
    VdataC3=dataC{3}(:,2)./detFrac;
    VdataC3Err=dataC{3}(:,3)./detFrac;
    
    [arrV1,arrV2,arrV3]=myV([ks01_opt;ks02_opt;ts_opt;w_opt;kx02_opt;tx_opt;kx01_opt]);
    
    interpV1=interp1(t,arrV1,tdataT1);
    interpV2=interp1(t,arrV2,tdataT1);
    interpV3=interp1(t,arrV3,tdataT1);
    
    %% synthesis and degradation rates
        
    arrTX=[0 tx_opt 210];
        
    arrKX_1=[kx01_opt kx02_opt(1) kx02_opt(1)];
    arrKX_2=[kx01_opt kx02_opt(2) kx02_opt(2)];
    arrKX_3=[kx01_opt kx02_opt(3) kx02_opt(3)];
        
    ks_opt{1}=(ks02_opt(1)-ks01_opt)./(1+exp(-w_opt(1)*(t-ts_opt(1))))+ks01_opt;
    ks_opt{2}=(ks02_opt(2)-ks01_opt)./(1+exp(-w_opt(2)*(t-ts_opt(2))))+ks01_opt;
    ks_opt{3}=(ks02_opt(3)-ks01_opt)./(1+exp(-w_opt(3)*(t-ts_opt(3))))+ks01_opt;
        
    kx_opt{1}=interp1(arrTX,arrKX_1,t);
    kx_opt{2}=interp1(arrTX,arrKX_2,t);
    kx_opt{3}=interp1(arrTX,arrKX_3,t);
    
    %% fitting
    % START VALUES

    ktc_=0.01;
    kdc_=0.01;
    kcd_=0.01;
    ktd_=0.01;
    kdt_=0.01;
    kf01_=0.01;
    kf02_=[0.01;0.01;0.01];
    tf_=50; 
    kct_=0.01;
        
    x0=[ktc_;kdc_;kcd_;ktd_;kdt_;kf01_;kf02_;tf_;kct_];
    
    % LOWER BOUNDS
    lb=1e-12.*ones(11,1);

    % UPPER BOUNDS
    ub=0.4; % ktc(1)
    ub=[ub;2]; % kdc(2)
    ub=[ub;2]; % kcd(3)
    ub=[ub;1]; % ktd(4)
    ub=[ub;1]; % kdt(5)
    ub=[ub;1]; % kf01(6)
    ub=[ub;10.*ones(3,1)]; % kf02(7:9)
    ub=[ub;150]; % tf(10)
    ub=[ub;1]; % kct(11)
    
    % CONSTRAINTS
    % inequalities
    Aineq=[];bineq=[];
%     Aineq=[Aineq;1,zeros(1,10)];bineq=[bineq;0]; % ktc=0 
    Aineq=[Aineq;0,0,0,1,zeros(1,7)];bineq=[bineq;0]; % ktd=0    

    % equalities
    A=[];b=[];   

    options=optimoptions(@fmincon,'Algorithm','sqp','Display','iter');    
    problem=createOptimProblem('fmincon','objective',@myFun,'x0',x0,'Aeq',A,'beq',b,'Aineq',Aineq,'bineq',bineq,'lb',lb,'ub',ub,'options',options);
    problem.options.MaxIterations=1000;
    problem.options.MaxFunctionEvaluations=45000;
    
    % multiple fitting for boot strapping
    arrParams=zeros(numBS,11);
    arrChi2=zeros(numBS,1);
    
    arrRPC1=zeros(numBS,length(tdataT1));
    arrRPC2=zeros(numBS,length(tdataT2));
    arrRPC3=zeros(numBS,length(tdataT3));    

    for iterBS=1:numBS
        
        currTdataT1_=zeros(length(tdataT1),1);
        currTdataT2_=zeros(length(tdataT2),1);
        currTdataT3_=zeros(length(tdataT3),1);
        
        currDdataT1_=zeros(length(tdataT1),1);
        currDdataT2_=zeros(length(tdataT2),1);
        currDdataT3_=zeros(length(tdataT3),1);        
        
        currCdataT1_=zeros(length(tdataT1),1);
        currCdataT2_=zeros(length(tdataT2),1);
        currCdataT3_=zeros(length(tdataT3),1);        
        
        currFdataT1_=zeros(length(tdataT1),1);
        currFdataT2_=zeros(length(tdataT2),1);
        currFdataT3_=zeros(length(tdataT3),1);
        
        for iterTS=1:length(tdataT1)
            
            randRPT=randi(sum(fitData_cip0p5.arrP1(iterTS,:)~=-1),1);
            currTdataT1_(iterTS)=fitData_cip0p5.arrP1(iterTS,randRPT);
            currDdataT1_(iterTS)=fitData_cip0p5.arrP2(iterTS,randRPT);
            currCdataT1_(iterTS)=fitData_cip0p5.arrP3(iterTS,randRPT);
            currFdataT1_(iterTS)=1-fitData_cip0p5.arrP1(iterTS,randRPT)-fitData_cip0p5.arrP2(iterTS,randRPT)-fitData_cip0p5.arrP3(iterTS,randRPT);
            arrRPC1(iterBS,iterTS)=randRPT;
        end
        
        for iterTS=1:length(tdataT2)
            
            randRPT=randi(sum(fitData_cip3.arrP1(iterTS,:)~=-1),1);
            currTdataT2_(iterTS)=fitData_cip3.arrP1(iterTS,randRPT);
            currDdataT2_(iterTS)=fitData_cip3.arrP2(iterTS,randRPT);
            currCdataT2_(iterTS)=fitData_cip3.arrP3(iterTS,randRPT);
            currFdataT2_(iterTS)=1-fitData_cip3.arrP1(iterTS,randRPT)-fitData_cip3.arrP2(iterTS,randRPT)-fitData_cip3.arrP3(iterTS,randRPT); 
            arrRPC2(iterBS,iterTS)=randRPT;
        end
        
        for iterTS=1:length(tdataT3)
            
            randRPT=randi(sum(fitData_cip20.arrP1(iterTS,:)~=-1),1);
            currTdataT3_(iterTS)=fitData_cip20.arrP1(iterTS,randRPT);
            currDdataT3_(iterTS)=fitData_cip20.arrP2(iterTS,randRPT);
            currCdataT3_(iterTS)=fitData_cip20.arrP3(iterTS,randRPT);
            currFdataT3_(iterTS)=1-fitData_cip20.arrP1(iterTS,randRPT)-fitData_cip20.arrP2(iterTS,randRPT)-fitData_cip20.arrP3(iterTS,randRPT);
            arrRPC3(iterBS,iterTS)=randRPT;
        end        

        currTdataT1=currTdataT1_.*interpV1;
        currDdataT1=currDdataT1_.*interpV1;
        currCdataT1=currCdataT1_.*interpV1;
        currFdataT1=currFdataT1_.*interpV1;

        currTdataT2=currTdataT2_.*interpV2;
        currDdataT2=currDdataT2_.*interpV2;
        currCdataT2=currCdataT2_.*interpV2;
        currFdataT2=currFdataT2_.*interpV2;        
        
        currTdataT3=currTdataT3_.*interpV3;
        currDdataT3=currDdataT3_.*interpV3;
        currCdataT3=currCdataT3_.*interpV3;
        currFdataT3=currFdataT3_.*interpV3;        
        
        estParams=fmincon(problem);
        
        arrChi2(iterBS)=myFun(estParams);
        
        arrParams(iterBS,:)=estParams(:)';
        
        disp([num2str(iterBS) '\' num2str(numBS)]);
    end    
    
    % fit of the combining data
    
    currTdataT1=TdataT1.*interpV1;
    currDdataT1=DdataT1.*interpV1;
    currCdataT1=CdataT1.*interpV1;
    currFdataT1=FdataT1.*interpV1;
    
    currTdataT2=TdataT2.*interpV2;
    currDdataT2=DdataT2.*interpV2;
    currCdataT2=CdataT2.*interpV2;
    currFdataT2=FdataT2.*interpV2;
    
    currTdataT3=TdataT3.*interpV3;
    currDdataT3=DdataT3.*interpV3;
    currCdataT3=CdataT3.*interpV3;
    currFdataT3=FdataT3.*interpV3;
    
    estParams=fmincon(problem);   
     
    [optT1,optT2,optT3,optC1,optC2,optC3,optD1,optD2,optD3,optF1,optF2,optF3]=myPlot(estParams);

    % optimal parameters
    opt_ktc=estParams(1);
    opt_kdc=estParams(2);
    opt_kcd=estParams(3);
    opt_ktd=estParams(4);
    opt_kdt=estParams(5);
    opt_kf01=estParams(6);
    opt_kf02=estParams(7:9);
    opt_tf=estParams(10);
    opt_kct=estParams(11);

    % bootstrapping results
    conf_ktc=lrSD(arrParams(:,1)');
    conf_kdc=lrSD(arrParams(:,2)');
    conf_kcd=lrSD(arrParams(:,3)');
    conf_ktd=lrSD(arrParams(:,4)');
    conf_kdt=lrSD(arrParams(:,5)'); 
    conf_kf01=lrSD(arrParams(:,6)');    
    conf_kf02=lrSD(arrParams(:,7:9)');    
    conf_tf=lrSD(arrParams(:,10)');  
    conf_kct=lrSD(arrParams(:,11)');  
    
    
	arrTFO=[0 conf_tf(2) 300];
        
    arrKFO_1=[conf_kf01(2) conf_kf02(1,2) conf_kf02(1,2)];
    arrKFO_2=[conf_kf01(2) conf_kf02(2,2) conf_kf02(2,2)];
    arrKFO_3=[conf_kf01(2) conf_kf02(3,2) conf_kf02(3,2)];
        
    kfO{1}=interp1(arrTFO,arrKFO_1,t);
    kfO{2}=interp1(arrTFO,arrKFO_2,t);
    kfO{3}=interp1(arrTFO,arrKFO_3,t);
    
    % show average result 
    [meanT1,meanT2,meanT3,meanC1,meanC2,meanC3,meanD1,meanD2,meanD3,meanF1,meanF2,meanF3]=myPlot(mean(arrParams));    
    
    fig1=figure();
    subplot(2,2,1);
    plot([0.5 3 20],arrParams(:,6).*ones(1,3),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_kf01(:,2).*ones(1,3),conf_kf01(:,1).*ones(1,3),conf_kf01(:,3).*ones(1,3),'-or');
    ylabel('\itk\rm_{f01} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    xlabel('[cip]');  
    
    subplot(2,2,2);
    plot([0.5 3 20],arrParams(:,7:9),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_kf02(:,2),conf_kf02(:,1),conf_kf02(:,3),'-or');
    ylabel('\itk\rm_{f02} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    xlabel('[cip]');  
    
    subplot(2,2,3);
    plot([0.5 3 20],arrParams(:,10).*ones(1,3),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_tf(:,2).*ones(1,3),conf_tf(:,1).*ones(1,3),conf_tf(:,3).*ones(1,3),'-or');
    ylabel('\itt\rm_{f} (min)');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    xlabel('[cip]');
    
    subplot(2,2,4);
    plot(t,kfO{1});
    hold on;  
    plot(t,kfO{2});
    plot(t,kfO{3});
    ylabel('\itk\rm_{f} (min^{-1})');
    xlim([0 300]); 
    xlabel('t (min)');  

    fig2=figure();
    subplot(2,3,1);
    plot([0.5 3 20],arrParams(:,1).*ones(1,3),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_ktc(:,2).*ones(1,3),conf_ktc(:,1).*ones(1,3),conf_ktc(:,3).*ones(1,3),'-or');
    ylabel('\itk\rm_{tc} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    xlabel('[cip]');     
   
    subplot(2,3,2);
    plot([0.5 3 20],arrParams(:,2).*ones(1,3),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_kdc(:,2).*ones(1,3),conf_kdc(:,1).*ones(1,3),conf_kdc(:,3).*ones(1,3),'-or');
    ylabel('\itk\rm_{dc} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    xlabel('[cip]');

    subplot(2,3,3);
    plot([0.5 3 20],arrParams(:,3).*ones(1,3),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_kcd(:,2).*ones(1,3),conf_kcd(:,1).*ones(1,3),conf_kcd(:,3).*ones(1,3),'-or');
    ylabel('\itk\rm_{cd} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    xlabel('[cip]'); 
    
    subplot(2,3,4);
    plot([0.5 3 20],arrParams(:,4).*ones(1,3),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_ktd(:,2).*ones(1,3),conf_ktd(:,1).*ones(1,3),conf_ktd(:,3).*ones(1,3),'-or');
    ylabel('\itk\rm_{td} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    xlabel('[cip]');     
    
    subplot(2,3,5);
    plot([0.5 3 20],arrParams(:,5).*ones(1,3),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_kdt(:,2).*ones(1,3),conf_kdt(:,1).*ones(1,3),conf_kdt(:,3).*ones(1,3),'-or');
    ylabel('\itk\rm_{dt} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    xlabel('[cip]');     
    
    subplot(2,3,6);
    plot([0.5 3 20],arrParams(:,11).*ones(1,3),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_kct(:,2).*ones(1,3),conf_kct(:,1).*ones(1,3),conf_kct(:,3).*ones(1,3),'-or');
    ylabel('\itk\rm_{ct} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    xlabel('[cip]');   

    fig4=figure();
    subplot(2,2,1);
    errorbar(tdataT1,TdataT1.*interpV1,TdataT1Err.*interpV1,'ok');
    hold on;
    plot(t,optT1,'-b');
    plot(t,meanT1,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Target bound');
    
    subplot(2,2,2);
    errorbar(tdataT1,DdataT1.*interpV1,DdataT1Err.*interpV1,'ok');
    hold on;
    plot(t,optD1,'-b');
    plot(t,meanD1,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('DNA bound');
    
    subplot(2,2,3);
    errorbar(tdataT1,CdataT1.*interpV1,CdataT1Err.*interpV1,'ok');
    hold on;
    plot(t,optC1,'-b');
    plot(t,meanC1,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Cytoplasmic');
    
    subplot(2,2,4);
    errorbar(tdataT1,FdataT1.*interpV1,FdataT1Err.*interpV1,'ok');
    hold on;
    plot(t,optF1,'-b');
    plot(t,meanF1,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Cleaved Fragments');

    % figureX(fig4);
    
    fig5=figure();
    subplot(2,2,1);
    errorbar(tdataT2,TdataT2.*interpV2,TdataT2Err.*interpV2,'ok');
    hold on;
    plot(t,optT2,'-b');
    plot(t,meanT2,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Target bound');
    
    subplot(2,2,2);
    errorbar(tdataT2,DdataT2.*interpV2,DdataT2Err.*interpV2,'ok');
    hold on;
    plot(t,optD2,'-b');
    plot(t,meanD2,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('DNA bound');
    
    subplot(2,2,3);
    errorbar(tdataT2,CdataT2.*interpV2,CdataT2Err.*interpV2,'ok');
    hold on;
    plot(t,optC2,'-b');
    plot(t,meanC2,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Cytoplasmic');
    
    subplot(2,2,4);
    errorbar(tdataT2,FdataT2.*interpV2,FdataT2Err.*interpV2,'ok');
    hold on;
    plot(t,optF2,'-b');
    plot(t,meanF2,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Cleaved Fragments');    

    % figureX(fig5);
    
    fig6=figure();
    subplot(2,2,1);
    errorbar(tdataT3,TdataT3.*interpV3,TdataT3Err.*interpV3,'ok');
    hold on;
    plot(t,optT3,'-b');
    plot(t,meanT3,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Target bound');
    
    subplot(2,2,2);
    errorbar(tdataT3,DdataT3.*interpV3,DdataT3Err.*interpV3,'ok');
    hold on;
    plot(t,optD3,'-b');
    plot(t,meanD3,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('DNA bound');
    
    subplot(2,2,3);
    errorbar(tdataT3,CdataT3.*interpV3,CdataT3Err.*interpV3,'ok');
    hold on;
    plot(t,optC3,'-b');
    plot(t,meanC3,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Cytoplasmic');
    
    subplot(2,2,4);
    errorbar(tdataT3,FdataT3.*interpV3,FdataT3Err.*interpV3,'ok');
    hold on;
    plot(t,optF3,'-b');
    plot(t,meanF3,'--b');    
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Cleaved Fragments');  
    
    str=['ktc=' num2str(opt_ktc')];
    disp(str);
    str=['kdc=' num2str(opt_kdc')];
    disp(str);
    str=['kcd=' num2str(opt_kcd')];
    disp(str);
    str=['ktd=' num2str(opt_ktd')];
    disp(str);
    str=['kdt=' num2str(opt_kdt')];
    disp(str);
    str=['kct=' num2str(opt_kct')];
    disp(str);    
    str=['kf01=' num2str(opt_kf01')];
    disp(str);    
    str=['kf02=' num2str(opt_kf02')];
    disp(str);    
    str=['tf=' num2str(opt_tf')];
    disp(str);    
    
    
    % saving results
    save('Result_combCurve.mat','estParams');
    save('Result_Bootstrapping_Curves.mat','arrParams');
    save('Result_Bootstrapping_Chi2.mat','arrChi2');
    save('Result_selected_Replicate_Cip0p5.mat','arrRPC1');
    save('Result_selected_Replicate_Cip3.mat','arrRPC2');
    save('Result_selected_Replicate_Cip20.mat','arrRPC3');  

   function [arrV1,arrV2,arrV3]=myV(params)
        
        % counting parameters
        ks01=params(1);
        ks02=params(2:4);
        ts=params(5:7);
        w=params(8:10);
        kx02=params(11:13);
        tx=params(14);
        kx01=params(15);

        v0=VdataC1(1);
        
        % synthesis function
        arrTX=[0 tx 210];
        
        arrKX_1=[kx01 kx02(1) kx02(1)];
        arrKX_2=[kx01 kx02(2) kx02(2)];
        arrKX_3=[kx01 kx02(3) kx02(3)];        
        
        % synthesis function
        ks{1}=(ks02(1)-ks01)./(1+exp(-w(1)*(t-ts(1))))+ks01;
        ks{2}=(ks02(2)-ks01)./(1+exp(-w(2)*(t-ts(2))))+ks01;
        ks{3}=(ks02(3)-ks01)./(1+exp(-w(3)*(t-ts(3))))+ks01;
        
        kx{1}=interp1(arrTX,arrKX_1,t);
        kx{2}=interp1(arrTX,arrKX_2,t);
        kx{3}=interp1(arrTX,arrKX_3,t);       

        pF{1}=aF(1).*(1-exp(-kF(1).*t))+pF0;
        pF{2}=aF(2).*(1-exp(-kF(2).*t))+pF0;
        pF{3}=aF(3).*(1-exp(-kF(3).*t))+pF0;
        
        arrV1=zeros(length(t),1);
        arrV2=zeros(length(t),1);
        arrV3=zeros(length(t),1);        
        
        for iterT=1:length(t)
            
            if iterT==1
                
                arrV1(iterT)=v0;
                arrV2(iterT)=v0;
                arrV3(iterT)=v0;
            else  

                
                dV1=(ks{1}(iterT-1)-kx{1}(iterT-1)*pF{1}(iterT-1)*arrV1(iterT-1))*dt;
                dV2=(ks{2}(iterT-1)-kx{2}(iterT-1)*pF{2}(iterT-1)*arrV2(iterT-1))*dt;
                dV3=(ks{3}(iterT-1)-kx{3}(iterT-1)*pF{3}(iterT-1)*arrV3(iterT-1))*dt;               
                
                arrV1(iterT)=arrV1(iterT-1)+dV1;
                arrV2(iterT)=arrV2(iterT-1)+dV2;
                arrV3(iterT)=arrV3(iterT-1)+dV3;                          
            end
        end                       
   end    

    function chi2=myFun(params)
        
        % response function of SOS triggered RecA-LexA interaction
        v0=VdataC1(1);        
        c0=[CdataT1(1);CdataT2(1);CdataT3(1)].*v0;
        f0=[FdataT1(1);FdataT2(1);FdataT3(1)].*v0;
        t0=[TdataT1(1);TdataT2(1);TdataT3(1)].*v0;
        d0=[DdataT1(1);DdataT2(1);DdataT3(1)].*v0;         

        ktc=params(1);
        kdc=params(2);
        kcd=params(3);
        ktd=params(4);
        kdt=params(5);
        kf01=params(6);
        kf02=params(7:9);
        tf=params(10);
        kct=params(11);
        
        % synthesis function
        arrTF=[0 tf 300];
        
        arrKF_1=[kf01 kf02(1) kf02(1)];
        arrKF_2=[kf01 kf02(2) kf02(2)];
        arrKF_3=[kf01 kf02(3) kf02(3)];        
        
        kf{1}=interp1(arrTF,arrKF_1,t);
        kf{2}=interp1(arrTF,arrKF_2,t);
        kf{3}=interp1(arrTF,arrKF_3,t);          
       
        arrT1=zeros(length(t),1);
        arrT2=zeros(length(t),1);
        arrT3=zeros(length(t),1);
        arrD1=zeros(length(t),1);
        arrD2=zeros(length(t),1);
        arrD3=zeros(length(t),1);
        arrC1=zeros(length(t),1);
        arrC2=zeros(length(t),1);
        arrC3=zeros(length(t),1);
        arrF1=zeros(length(t),1);
        arrF2=zeros(length(t),1);
        arrF3=zeros(length(t),1);
        
        for iterT=1:length(t)
            
            if iterT==1
                
                arrC1(iterT)=c0(1);
                arrC2(iterT)=c0(2);
                arrC3(iterT)=c0(3);
                arrD1(iterT)=d0(1);
                arrD2(iterT)=d0(2);
                arrD3(iterT)=d0(3);
                arrT1(iterT)=t0(1);
                arrT2(iterT)=t0(2);
                arrT3(iterT)=t0(3);
                arrF1(iterT)=f0(1);
                arrF2(iterT)=f0(2);
                arrF3(iterT)=f0(3);
            else
               
                dT1=(kct*arrC1(iterT-1)+kdt*arrD1(iterT-1)-ktd*arrT1(iterT-1)-ktc*arrT1(iterT-1))*dt;
                dT2=(kct*arrC2(iterT-1)+kdt*arrD2(iterT-1)-ktd*arrT2(iterT-1)-ktc*arrT2(iterT-1))*dt;
                dT3=(kct*arrC2(iterT-1)+kdt*arrD3(iterT-1)-ktd*arrT3(iterT-1)-ktc*arrT3(iterT-1))*dt;
                
                dD1=(kcd*arrC1(iterT-1)+ktd*arrT1(iterT-1)-kdc*arrD1(iterT-1)-kdt*arrD1(iterT-1))*dt;
                dD2=(kcd*arrC2(iterT-1)+ktd*arrT2(iterT-1)-kdc*arrD2(iterT-1)-kdt*arrD2(iterT-1))*dt;
                dD3=(kcd*arrC3(iterT-1)+ktd*arrT3(iterT-1)-kdc*arrD3(iterT-1)-kdt*arrD3(iterT-1))*dt;
                
                dC1=(ks_opt{1}(iterT-1)+kdc*arrD1(iterT-1)+ktc*arrT1(iterT-1)-kct*arrC1(iterT-1)-kcd*arrC1(iterT-1)-kf{1}(iterT-1)*arrC1(iterT-1))*dt;
                dC2=(ks_opt{2}(iterT-1)+kdc*arrD2(iterT-1)+ktc*arrT2(iterT-1)-kct*arrC2(iterT-1)-kcd*arrC2(iterT-1)-kf{2}(iterT-1)*arrC2(iterT-1))*dt;
                dC3=(ks_opt{3}(iterT-1)+kdc*arrD3(iterT-1)+ktc*arrT3(iterT-1)-kct*arrC3(iterT-1)-kcd*arrC3(iterT-1)-kf{3}(iterT-1)*arrC3(iterT-1))*dt;
                                
                dF1=(kf{1}(iterT-1)*arrC1(iterT-1)-kx_opt{1}(iterT-1)*arrF1(iterT-1))*dt;
                dF2=(kf{2}(iterT-1)*arrC2(iterT-1)-kx_opt{2}(iterT-1)*arrF2(iterT-1))*dt;
                dF3=(kf{3}(iterT-1)*arrC3(iterT-1)-kx_opt{3}(iterT-1)*arrF3(iterT-1))*dt;

                arrT1(iterT)=arrT1(iterT-1)+dT1;
                arrT2(iterT)=arrT2(iterT-1)+dT2;
                arrT3(iterT)=arrT3(iterT-1)+dT3;
                               
                arrD1(iterT)=arrD1(iterT-1)+dD1;
                arrD2(iterT)=arrD2(iterT-1)+dD2;
                arrD3(iterT)=arrD3(iterT-1)+dD3;
                
                arrC1(iterT)=arrC1(iterT-1)+dC1;
                arrC2(iterT)=arrC2(iterT-1)+dC2;
                arrC3(iterT)=arrC3(iterT-1)+dC3;
                
                arrF1(iterT)=arrF1(iterT-1)+dF1;
                arrF2(iterT)=arrF2(iterT-1)+dF2;
                arrF3(iterT)=arrF3(iterT-1)+dF3;
                               
            end
        end

        interpT1=interp1(t,arrT1,tdataT1);
        interpT2=interp1(t,arrT2,tdataT2);
        interpT3=interp1(t,arrT3,tdataT3);
                
        interpD1=interp1(t,arrD1,tdataT1);
        interpD2=interp1(t,arrD2,tdataT2);
        interpD3=interp1(t,arrD3,tdataT3);
        
        interpC1=interp1(t,arrC1,tdataT1);
        interpC2=interp1(t,arrC2,tdataT2);
        interpC3=interp1(t,arrC3,tdataT3);
        
        interpF1=interp1(t,arrF1,tdataT1);        
        interpF2=interp1(t,arrF2,tdataT2);        
        interpF3=interp1(t,arrF3,tdataT3);        
              
        chi2T_T1=sum((interpT1(currTdataT1~=0)-currTdataT1(currTdataT1~=0)).^2./interpT1(currTdataT1~=0));
        chi2T_T2=sum((interpT2(currTdataT2~=0)-currTdataT2(currTdataT2~=0)).^2./interpT2(currTdataT2~=0));
        chi2T_T3=sum((interpT3(currTdataT3~=0)-currTdataT3(currTdataT3~=0)).^2./interpT3(currTdataT3~=0));
        chi2T_D1=sum((interpD1(currDdataT1~=0)-currDdataT1(currDdataT1~=0)).^2./interpD1(currDdataT1~=0));
        chi2T_D2=sum((interpD2(currDdataT2~=0)-currDdataT2(currDdataT2~=0)).^2./interpD2(currDdataT2~=0));
        chi2T_D3=sum((interpD3(currDdataT3~=0)-currDdataT3(currDdataT3~=0)).^2./interpD3(currDdataT3~=0));
        chi2T_C1=sum((interpC1(currCdataT1~=0)-currCdataT1(currCdataT1~=0)).^2./interpC1(currCdataT1~=0));
        chi2T_C2=sum((interpC2(currCdataT2~=0)-currCdataT2(currCdataT2~=0)).^2./interpC2(currCdataT2~=0));
        chi2T_C3=sum((interpC3(currCdataT3~=0)-currCdataT3(currCdataT3~=0)).^2./interpC3(currCdataT3~=0));
        chi2T_F1=sum((interpF1(currFdataT1~=0)-currFdataT1(currFdataT1~=0)).^2./interpF1(currFdataT1~=0));
        chi2T_F2=sum((interpF2(currFdataT2~=0)-currFdataT2(currFdataT2~=0)).^2./interpF2(currFdataT2~=0));
        chi2T_F3=sum((interpF3(currFdataT3~=0)-currFdataT3(currFdataT3~=0)).^2./interpF3(currFdataT3~=0));

        chi2=chi2T_T1+chi2T_T2+chi2T_T3+chi2T_D1+chi2T_D2+chi2T_D3+chi2T_C1+chi2T_C2+chi2T_C3+chi2T_F1+chi2T_F2+chi2T_F3;
        
        if (chi2<0)||isnan(chi2)
            
            disp('function created NaN value!');
        end
    end

    function [arrT1,arrT2,arrT3,arrC1,arrC2,arrC3,arrD1,arrD2,arrD3,arrF1,arrF2,arrF3]=myPlot(params)
        
        % response function of SOS triggered RecA-LexA interaction
        v0=VdataC1(1);        
        c0=[CdataT1(1);CdataT2(1);CdataT3(1)].*v0;
        f0=[FdataT1(1);FdataT2(1);FdataT3(1)].*v0;
        t0=[TdataT1(1);TdataT2(1);TdataT3(1)].*v0;
        d0=[DdataT1(1);DdataT2(1);DdataT3(1)].*v0;         

        ktc=params(1);
        kdc=params(2);
        kcd=params(3);
        ktd=params(4);
        kdt=params(5);
        kf01=params(6);
        kf02=params(7:9);
        tf=params(10);
        kct=params(11);
        
        % synthesis function
        arrTF=[0 tf 300];
        
        arrKF_1=[kf01 kf02(1) kf02(1)];
        arrKF_2=[kf01 kf02(2) kf02(2)];
        arrKF_3=[kf01 kf02(3) kf02(3)];        
        
        kf{1}=interp1(arrTF,arrKF_1,t);
        kf{2}=interp1(arrTF,arrKF_2,t);
        kf{3}=interp1(arrTF,arrKF_3,t);          
       
        arrT1=zeros(length(t),1);
        arrT2=zeros(length(t),1);
        arrT3=zeros(length(t),1);
        arrD1=zeros(length(t),1);
        arrD2=zeros(length(t),1);
        arrD3=zeros(length(t),1);
        arrC1=zeros(length(t),1);
        arrC2=zeros(length(t),1);
        arrC3=zeros(length(t),1);
        arrF1=zeros(length(t),1);
        arrF2=zeros(length(t),1);
        arrF3=zeros(length(t),1);
        
        for iterT=1:length(t)
            
            if iterT==1
                
                arrC1(iterT)=c0(1);
                arrC2(iterT)=c0(2);
                arrC3(iterT)=c0(3);
                arrD1(iterT)=d0(1);
                arrD2(iterT)=d0(2);
                arrD3(iterT)=d0(3);
                arrT1(iterT)=t0(1);
                arrT2(iterT)=t0(2);
                arrT3(iterT)=t0(3);
                arrF1(iterT)=f0(1);
                arrF2(iterT)=f0(2);
                arrF3(iterT)=f0(3);
            else
               
                dT1=(kct*arrC1(iterT-1)+kdt*arrD1(iterT-1)-ktd*arrT1(iterT-1)-ktc*arrT1(iterT-1))*dt;
                dT2=(kct*arrC2(iterT-1)+kdt*arrD2(iterT-1)-ktd*arrT2(iterT-1)-ktc*arrT2(iterT-1))*dt;
                dT3=(kct*arrC2(iterT-1)+kdt*arrD3(iterT-1)-ktd*arrT3(iterT-1)-ktc*arrT3(iterT-1))*dt;
                
                dD1=(kcd*arrC1(iterT-1)+ktd*arrT1(iterT-1)-kdc*arrD1(iterT-1)-kdt*arrD1(iterT-1))*dt;
                dD2=(kcd*arrC2(iterT-1)+ktd*arrT2(iterT-1)-kdc*arrD2(iterT-1)-kdt*arrD2(iterT-1))*dt;
                dD3=(kcd*arrC3(iterT-1)+ktd*arrT3(iterT-1)-kdc*arrD3(iterT-1)-kdt*arrD3(iterT-1))*dt;
                
                dC1=(ks_opt{1}(iterT-1)+kdc*arrD1(iterT-1)+ktc*arrT1(iterT-1)-kct*arrC1(iterT-1)-kcd*arrC1(iterT-1)-kf{1}(iterT-1)*arrC1(iterT-1))*dt;
                dC2=(ks_opt{2}(iterT-1)+kdc*arrD2(iterT-1)+ktc*arrT2(iterT-1)-kct*arrC2(iterT-1)-kcd*arrC2(iterT-1)-kf{2}(iterT-1)*arrC2(iterT-1))*dt;
                dC3=(ks_opt{3}(iterT-1)+kdc*arrD3(iterT-1)+ktc*arrT3(iterT-1)-kct*arrC3(iterT-1)-kcd*arrC3(iterT-1)-kf{3}(iterT-1)*arrC3(iterT-1))*dt;
                                
                dF1=(kf{1}(iterT-1)*arrC1(iterT-1)-kx_opt{1}(iterT-1)*arrF1(iterT-1))*dt;
                dF2=(kf{2}(iterT-1)*arrC2(iterT-1)-kx_opt{2}(iterT-1)*arrF2(iterT-1))*dt;
                dF3=(kf{3}(iterT-1)*arrC3(iterT-1)-kx_opt{3}(iterT-1)*arrF3(iterT-1))*dt;

                arrT1(iterT)=arrT1(iterT-1)+dT1;
                arrT2(iterT)=arrT2(iterT-1)+dT2;
                arrT3(iterT)=arrT3(iterT-1)+dT3;
                               
                arrD1(iterT)=arrD1(iterT-1)+dD1;
                arrD2(iterT)=arrD2(iterT-1)+dD2;
                arrD3(iterT)=arrD3(iterT-1)+dD3;
                
                arrC1(iterT)=arrC1(iterT-1)+dC1;
                arrC2(iterT)=arrC2(iterT-1)+dC2;
                arrC3(iterT)=arrC3(iterT-1)+dC3;
                
                arrF1(iterT)=arrF1(iterT-1)+dF1;
                arrF2(iterT)=arrF2(iterT-1)+dF2;
                arrF3(iterT)=arrF3(iterT-1)+dF3;                             
            end
        end
    end 
end