function RespTimeDelay_vs_CIP()    

    % - numerical infinitesimal calculation of populations
    % - logistic response function for new synthesis with global minimum 
    % rate ks01, individual response time ts, response speed w and 
    % response amplitude k02(i)>k01
    % - without time point t=180min
    % - kxEff replaced by pf*kx
    %
    % Andreas Hartmann | 03.06.2020
    
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
    numBS_C=30;
    numBS=30;
    
    pathCount='data\Result_Bootstrapping_Curves.mat';
    pathCountOPT='data\Result_combCurve.mat';
    
    FitCount=load(pathCount);
    FitCountOPT=load(pathCountOPT);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    addpath('scripts');
    
    %% tracking data
    t=(0:dt:1000); % (min)
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
    VdataC1=dataC{1}(:,2)./detFrac;
    
    [arrV1,arrV2,arrV3]=myV([ks01_opt;ks02_opt;ts_opt;w_opt;kx02_opt;tx_opt;kx01_opt]);
    
    interpV1=interp1(t,arrV1,tdataT1);
    interpV2=interp1(t,arrV2,tdataT1);
    interpV3=interp1(t,arrV3,tdataT1);
    
    %% synthesis and degradation rates
        
    arrTX=[0 tx_opt t(end)];
        
    arrKX_1=[kx01_opt kx02_opt(1) kx02_opt(1)];
    arrKX_2=[kx01_opt kx02_opt(2) kx02_opt(2)];
    arrKX_3=[kx01_opt kx02_opt(3) kx02_opt(3)];
        
    ks_opt{1}=(ks02_opt(1)-ks01_opt)./(1+exp(-w_opt(1)*(t-ts_opt(1))))+ks01_opt;
    ks_opt{2}=(ks02_opt(2)-ks01_opt)./(1+exp(-w_opt(2)*(t-ts_opt(2))))+ks01_opt;
    ks_opt{3}=(ks02_opt(3)-ks01_opt)./(1+exp(-w_opt(3)*(t-ts_opt(3))))+ks01_opt;
        
    kx_opt{1}=interp1(arrTX,arrKX_1,t);
    kx_opt{2}=interp1(arrTX,arrKX_2,t);
    kx_opt{3}=interp1(arrTX,arrKX_3,t);
     

    
    %% fitting of T1  
    
        arrV1EC=zeros(length(tdataT1),100);
    arrV2EC=zeros(length(tdataT2),100);
    arrV3EC=zeros(length(tdataT3),100);
    
    % error evaluation
    for iterER=1:100
        
        [rndV1,rndV2,rndV3]=myV(FitCount.arrParams(randi(numBS_C,1),:)');        
        
        arrV1EC(:,iterER)=interp1(t,rndV1,tdataT1);
        arrV2EC(:,iterER)=interp1(t,rndV2,tdataT1);
        arrV3EC(:,iterER)=interp1(t,rndV3,tdataT1);                    
    end
    
    stdV1=std(arrV1EC,[],2);
    stdV2=std(arrV2EC,[],2);
    stdV3=std(arrV3EC,[],2); 
    
    errT1=sqrt((interpV1.*TdataT1Err).^2+(TdataT1.*stdV1).^2);
    errT2=sqrt((interpV2.*TdataT2Err).^2+(TdataT2.*stdV2).^2);
    errT3=sqrt((interpV3.*TdataT3Err).^2+(TdataT3.*stdV3).^2);
       
    
    arrTC_T1=zeros(numBS,1);
    arrTC_T2=zeros(numBS,1);
    arrTC_T3=zeros(numBS,1);
    
    for iterBS=1:numBS
        
        currTdataT1_=zeros(length(tdataT1),1);
        currTdataT2_=zeros(length(tdataT2),1);
        currTdataT3_=zeros(length(tdataT3),1);
        
        for iterTS=1:length(tdataT1)
            
            randRPT=randi(sum(fitData_cip0p5.arrP1(iterTS,:)~=-1),1);
            currTdataT1_(iterTS)=fitData_cip0p5.arrP1(iterTS,randRPT);
        end
        
        for iterTS=1:length(tdataT2)
            
            randRPT=randi(sum(fitData_cip3.arrP1(iterTS,:)~=-1),1);
            currTdataT2_(iterTS)=fitData_cip3.arrP1(iterTS,randRPT);
        end
        
        for iterTS=1:length(tdataT3)
            
            randRPT=randi(sum(fitData_cip20.arrP1(iterTS,:)~=-1),1);
            currTdataT3_(iterTS)=fitData_cip20.arrP1(iterTS,randRPT);
        end        

        [rndV1,rndV2,rndV3]=myV(FitCount.arrParams(randi(numBS_C,1),:)');        
              
        interpRNDV1=interp1(t,rndV1,tdataT1);
        interpRNDV2=interp1(t,rndV2,tdataT1);
        interpRNDV3=interp1(t,rndV3,tdataT1);
        
        currTdataT1=currTdataT1_.*interpRNDV1;
        currTdataT2=currTdataT2_.*interpRNDV2;
        currTdataT3=currTdataT3_.*interpRNDV3;
     
        [cfunT1,~,~]=fit(tdataT1(1:4),currTdataT1(1:4),'a*exp(-b*x)+c','Start',[currTdataT1(1) 1/tdataT1(2) currTdataT1(4)],'Weights',1./errT1(1:4),'Lower',[0 0 0]);
        [cfunT2,~,~]=fit(tdataT2(1:4),currTdataT2(1:4),'a*exp(-b*x)+c','Start',[currTdataT2(1) 1/tdataT2(2) currTdataT2(4)],'Weights',1./errT2(1:4),'Lower',[0 0 0]);
        [cfunT3,~,~]=fit(tdataT3(1:4),currTdataT3(1:4),'a*exp(-b*x)+c','Start',[currTdataT3(1) 1/tdataT3(2) currTdataT3(4)],'Weights',1./errT3(1:4),'Lower',[0 0 0]);        
        
        arrTC_T1(iterBS)=1/cfunT1.b;
        arrTC_T2(iterBS)=1/cfunT2.b;
        arrTC_T3(iterBS)=1/cfunT3.b;
        
        disp([num2str(iterBS) '\' num2str(numBS)]);
    end    
    
    currTdataT1=TdataT1.*interpV1;    
    currTdataT2=TdataT2.*interpV2;
    currTdataT3=TdataT3.*interpV3;
    
    [cfunT1,~,~]=fit(tdataT1(1:4),currTdataT1(1:4),'a*exp(-b*x)+c','Start',[currTdataT1(1) 1/tdataT1(2) currTdataT1(4)],'Weights',1./errT1(1:4),'Lower',[0 0 0]);
    [cfunT2,~,~]=fit(tdataT2(1:4),currTdataT2(1:4),'a*exp(-b*x)+c','Start',[currTdataT2(1) 1/tdataT2(2) currTdataT2(4)],'Weights',1./errT2(1:4),'Lower',[0 0 0]);
    [cfunT3,~,~]=fit(tdataT3(1:4),currTdataT3(1:4),'a*exp(-b*x)+c','Start',[currTdataT3(1) 1/tdataT3(2) currTdataT3(4)],'Weights',1./errT3(1:4),'Lower',[0 0 0]);

    conf_TC_T1=lrSD(arrTC_T1');
    conf_TC_T2=lrSD(arrTC_T2');
    conf_TC_T3=lrSD(arrTC_T3');
    
    figure();
    subplot(2,3,1);
    errorbar(tdataT1,currTdataT1,sqrt((interpV1.*TdataT1Err).^2+(TdataT1.*stdV1).^2),'ok');
    hold on;
    plot(cfunT1);
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Target bound [CIP]=0.5');    
    
    subplot(2,3,2);
    errorbar(tdataT2,currTdataT2,sqrt((interpV2.*TdataT2Err).^2+(TdataT2.*stdV2).^2),'ok');
    hold on;
    plot(cfunT2);
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Target bound [CIP]=3');    
    
    subplot(2,3,3);
    errorbar(tdataT3,currTdataT3,sqrt((interpV3.*TdataT3Err).^2+(TdataT3.*stdV3).^2),'ok');
    hold on;
    plot(cfunT3);
    xlabel('time (min)');
    xlim([0 200]);
    ylabel('molecules per µm^{2}');
    title('Target bound [CIP]=20');    
    
    tresp=FitCountOPT.estParams(5:7);
    conf_ts=lrSD(FitCount.arrParams(:,5:7)');
    
    diffTS=tresp(:)-[1/cfunT1.b 1/cfunT2.b 1/cfunT3.b]';
    err_diffTS=sqrt([std(FitCount.arrParams(:,5)) std(FitCount.arrParams(:,6)) std(FitCount.arrParams(:,7))].^2+[std(arrTC_T1) std(arrTC_T2) std(arrTC_T3)].^2);
    
    subplot(2,3,4);
    errorbar([1/cfunT1.b 1/cfunT2.b 1/cfunT3.b],tresp,[std(FitCount.arrParams(:,5)) std(FitCount.arrParams(:,6)) std(FitCount.arrParams(:,7))],[std(FitCount.arrParams(:,5)) std(FitCount.arrParams(:,6)) std(FitCount.arrParams(:,7))],[std(arrTC_T1) std(arrTC_T2) std(arrTC_T3)],[std(arrTC_T1) std(arrTC_T2) std(arrTC_T3)],'ok')
    xlabel('characteristic decay time of T (min)');
    ylabel('response time of \itv\rm_s(t) (min)');
    xlim([0 20]);
    ylim([0 150]);

    subplot(2,3,5);
    errorbar([0.5 3 20],diffTS,err_diffTS,'ok')
    set(gca,'XScale','log');
    xlabel('[CIP]');
    ylabel('response time delay (min)');
    xlim([0.1 100]);
    ylim([0 150]);
    
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
        arrTX=[0 tx t(end)];
        
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
        arrTF=[0 tf t(end)];
        
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