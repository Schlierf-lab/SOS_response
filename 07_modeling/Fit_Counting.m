% Extraction of the on-time, off-time, blinking time and bleaching time distributions
% --Andreas Hartmann

function Fit_Counting()    

    % - numerical infinitesimal calculation of populations
    % - logistic response function for new synthesis with global minimum 
    % rate ks01, individual response time ts, response speed w and 
    % response amplitude k02(i)>k01
    % - global rate kx01<kx02(i)
    % - global initial amount of visible molecules v0
    % - without time point t=180min
    % - kxEff replaced by pf*kx
    % - fitted with Chi2
    %
    % Andreas Hartmann | 26.05.2020
    
    %% parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataC{1}=importdata('data\cip0p5.txt');
    dataC{2}=importdata('data\cip3.txt');
    dataC{3}=importdata('data\cip20_wo180.txt');
    
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
    
    %% fitting
    % START VALUES
    % counting
    ks01_=21.62;
    ks02_=[553.4;830.4;1527];
    ts_=[103.4;88.98;56.42];
    w_=[0.02775;0.04005;0.06981];
    kx02_=[3.52;8.737;11.14];
    tx_=61.95; 
    kx01_=3.436;    

    x0=[ks01_;ks02_;ts_;w_;kx02_;tx_;kx01_];
    
    % LOWER BOUNDS
    lb=1e-12.*ones(15,1);
    
    % UPPER BOUNDS
    ub=200; % ks01 (1)
    ub=[ub;750;1200;1750]; % ks02 (2:4) CHNAGE THIS TO 2000 FROM 4000
    ub=[ub;200.*ones(3,1)]; % ts (5:7)
    ub=[ub;0.05;0.06;0.1]; % w (8:10)
    ub=[ub;2000.*ones(3,1)]; % kx02 (11:13)
    ub=[ub;75]; % tx (14)   
    ub=[ub;10]; % kx01 (15)   
    
    % CONSTRAINTS
    % inequalities
    Aineq=[];bineq=[];
    Aineq=[Aineq;1,-1,zeros(1,13)];bineq=[bineq;0]; % ks01<=ks02(1),
    Aineq=[Aineq;0,1,-1,zeros(1,12)];bineq=[bineq;0]; % ks02(1)<=ks02(2)    
    Aineq=[Aineq;0,0,1,-1,zeros(1,11)];bineq=[bineq;0]; % ks02(2)<=ks02(3)
    Aineq=[Aineq;0,0,0,0,-1,1,zeros(1,9)];bineq=[bineq;0]; % ts(2)<=ts(1)
    Aineq=[Aineq;0,0,0,0,0,-1,1,zeros(1,8)];bineq=[bineq;0]; % ts(3)<=ts(2)
    Aineq=[Aineq;zeros(1,7),1,-1,zeros(1,6)];bineq=[bineq;0]; % w(1)<=w(2)    
    Aineq=[Aineq;zeros(1,8),1,-1,zeros(1,5)];bineq=[bineq;0]; % w(2)<=w(3)  

    Aineq=[Aineq;zeros(1,10),-1,zeros(1,3),1];bineq=[bineq;0]; % kx01<=kx02(1)  
    Aineq=[Aineq;zeros(1,10),1,-1,zeros(1,3)];bineq=[bineq;0]; % kx02(1)<=kx02(2)    
    Aineq=[Aineq;zeros(1,11),1,-1,0,0];bineq=[bineq;0]; % kx02(2)<=kx02(3)

    % equalities
    A=[];b=[];  

    options=optimoptions(@fmincon,'Algorithm','sqp','Display','iter');    
    problem=createOptimProblem('fmincon','objective',@myfun,'x0',x0,'Aeq',A,'beq',b,'Aineq',Aineq,'bineq',bineq,'lb',lb,'ub',ub,'options',options);
    problem.options.MaxIterations=1000;
    problem.options.MaxFunctionEvaluations=45000;
    
    % multiple fitting for boot strapping
    arrParams=zeros(numBS,15);
    arrChi2=zeros(numBS,1);  

    for iterBS=1:numBS       

        currVdataC1=randn(length(VdataC1),1).*VdataC1Err+VdataC1;
        currVdataC2=randn(length(VdataC2),1).*VdataC2Err+VdataC2;
        currVdataC3=randn(length(VdataC3),1).*VdataC3Err+VdataC3;
        
        estParams=fmincon(problem);
        
        arrChi2(iterBS)=myfun(estParams);
        
        arrParams(iterBS,:)=estParams(:)';
        
        disp([num2str(iterBS) '\' num2str(numBS)]);
    end   
    
    currVdataC1=VdataC1;
    currVdataC2=VdataC2;
    currVdataC3=VdataC3;
    
    estParams=fmincon(problem);   
     
    [optV1,optV2,optV3]=myPlot(estParams);
    red_chi2=myRedC(estParams);
    
    % counting result of combining fit
    opt_ks01=estParams(1);
    opt_ks02=estParams(2:4);
    opt_ts=estParams(5:7);
    opt_w=estParams(8:10);
    opt_kx02=estParams(11:13);
    opt_tx=estParams(14);
    opt_kx01=estParams(15);

    % counting SD results of bootstrapping
    conf_ks01=lrSD(arrParams(:,1)');
    conf_ks02=lrSD(arrParams(:,2:4)');
    conf_ts=lrSD(arrParams(:,5:7)');
    conf_w=lrSD(arrParams(:,8:10)');
    conf_kx02=lrSD(arrParams(:,11:13)');
    conf_tx=lrSD(arrParams(:,14)');
    conf_kx01=lrSD(arrParams(:,15)');

    % show average result      
    [meanV1,meanV2,meanV3]=myPlot(mean(arrParams));
    
    conf_arrTX=[0 conf_tx(:,2) 210];
        
    conf_arrKX_1=[conf_kx01(1,2) conf_kx02(1,2) conf_kx02(1,2)];
    conf_arrKX_2=[conf_kx01(1,2) conf_kx02(2,2) conf_kx02(2,2)];
    conf_arrKX_3=[conf_kx01(1,2) conf_kx02(3,2) conf_kx02(3,2)];

    opt_arrTX=[0 opt_tx 210];
        
    opt_arrKX_1=[opt_kx01 opt_kx02(1) opt_kx02(1)];
    opt_arrKX_2=[opt_kx01 opt_kx02(2) opt_kx02(2)];
    opt_arrKX_3=[opt_kx01 opt_kx02(3) opt_kx02(3)];    

    fig1=figure();
    subplot(2,2,1);
    plot([0.5 3 20],arrParams(:,2:4),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_ks02(:,2),conf_ks02(:,1),conf_ks02(:,3),'-or');
    plot([0.5 3 20],[arrParams(:,1) arrParams(:,1) arrParams(:,1)],'.k');
    errorbar([0.5 3 20],[conf_ks01(:,2);conf_ks01(:,2);conf_ks01(:,2)],[conf_ks01(:,1);conf_ks01(:,1);conf_ks01(:,1)],[conf_ks01(:,3);conf_ks01(:,3);conf_ks01(:,3)],'-or');
    ylabel('\itv\rm_{s01} & \itv\rm_{s02} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    xlabel('[cip]');
    title('Response curve \itv\rm_s(t)');
    
    subplot(2,2,2);
    plot([0.5 3 20],arrParams(:,5:7),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_ts(:,2),conf_ts(:,1),conf_ts(:,3),'-or');
    ylabel('\itt\rm_{s} (min)');
    set(gca,'XScale','log');
    xlim([0.1 100]);    
    xlabel('[cip]'); 
     
    subplot(2,2,3);
    plot([0.5 3 20],arrParams(:,8:10),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_w(:,2),conf_w(:,1),conf_w(:,3),'-or');
    ylabel('\itw');
    set(gca,'XScale','log');
    xlim([0.1 100]);    
    xlabel('[cip]');     
    
    subplot(2,2,4);
    plot(t,((conf_ks02(1,2)-conf_ks01(2))./(1+exp(-conf_w(1,2)*(t-conf_ts(1,2))))+conf_ks01(2)));
    hold on;
    plot(t,((conf_ks02(2,2)-conf_ks01(2))./(1+exp(-conf_w(2,2)*(t-conf_ts(2,2))))+conf_ks01(2)));
    plot(t,((conf_ks02(3,2)-conf_ks01(2))./(1+exp(-conf_w(3,2)*(t-conf_ts(3,2))))+conf_ks01(2)));
    xlim([0 200]);
    ylabel('\itv\rm_s(t) (min^{-1})');
    xlabel('time (min)');     

    figure();
    subplot(2,2,1);
    plot([0.5 3 20],arrParams(:,15).*ones(1,3),'.k');
    hold on;
    errorbar([0.5 3 20],conf_kx01(1,2).*ones(1,3),conf_kx01(1,1).*ones(1,3),conf_kx01(1,3).*ones(1,3),'-ob');
    ylabel('\itk\rm_{x01} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    ylim([0 15]); 
    xlabel('[cip]'); 
    title('Response curve \itk\rm_x(t)');
    
    subplot(2,2,2);
    plot([0.5 3 20],arrParams(:,11:13),'.k');
    hold on;   
    errorbar([0.5 3 20],conf_kx02(:,2),conf_kx02(:,1),conf_kx02(:,3),'-ob');
    ylabel('\itk\rm_{x02} (min^{-1})');
    set(gca,'XScale','log');
    xlim([0.1 100]); 
    ylim([0 15]); 
    xlabel('[cip]'); 
    
    subplot(2,2,3);
    plot([0.5 3 20],arrParams(:,14).*ones(1,3),'.k');
    hold on;
    errorbar([0.5 3 20],conf_tx(:,2).*ones(1,3),conf_tx(:,1).*ones(1,3),conf_tx(:,3).*ones(1,3),'-ob');
    ylabel('\itt\rm_{x} (min)');
    set(gca,'XScale','log');
    xlim([0.1 100]);    
    ylim([0 150]); 
    xlabel('[cip]');
    
    subplot(2,2,4);
    plot(conf_arrTX,conf_arrKX_1,'-');
    hold on;
    plot(conf_arrTX,conf_arrKX_2,'-');
    plot(conf_arrTX,conf_arrKX_3,'-');
    xlim([0 200]);
    ylim([0 15]); 
    ylabel('\itk\rm_x(t) (min^{-1})');
    xlabel('time (min)'); 
    
    % figureX(fig1);

    fig2=figure();
    annotation(gcf,'textbox',...
    [0.0589285714285714 0.832333334539618 0.176785709549274 0.0642857130794299],...
    'String',{'\bfk rate results of bootstrapping'},'LineStyle','none','HorizontalAlignment','center');
    
    subplot(2,2,2);
    errorbar(tdataC1,VdataC1,VdataC1Err,'ok');
    hold on;
    plot(t,optV1,'-b');
    plot(t,meanV1,'--b');    
    xlabel('time (min)');
    xlim([0 250]);
    ylabel('LexA concentration (molecules/µm^2)');
    title('cip 0.5');

    subplot(2,2,3);
    errorbar(tdataC2,VdataC2,VdataC2Err,'ok');
    hold on;
    plot(t,optV2,'-b');
    plot(t,meanV2,'--b');    
    xlabel('time (min)');
    xlim([0 250]);
    ylabel('LexA concentration (molecules/µm^2)');
    title('cip 3');
    
    subplot(2,2,4);
    errorbar(tdataC3,VdataC3,VdataC3Err,'ok');
    hold on;
    plot(t,optV3,'-b');
    plot(t,meanV3,'--b');    
    xlabel('time (min)');
    xlim([0 250]);
    ylabel('LexA concentration (molecules/µm^2)');
    title('cip 20');    
    
    str=['ks01=' num2str(opt_ks01)];
    disp(str);
    str=['ks02=' num2str(opt_ks02')];
    disp(str);
    str=['ts=' num2str(opt_ts')];
    disp(str);
    str=['w=' num2str(opt_w')];
    disp(str);
    str=['kx01=' num2str(opt_kx01')];
    disp(str);    
    str=['kx02=' num2str(opt_kx02')];
    disp(str);
    str=['tx=' num2str(opt_tx')];
    disp(str); 
    str=['red.Chi2=' num2str(red_chi2)];
    disp(str); 
    
    % saving results
    save('Result_combCurve.mat','estParams');
    save('Result_Bootstrapping_Curves.mat','arrParams');
    save('Result_Bootstrapping_Chi2.mat','arrChi2'); 
    
    function chi2=myfun(params)
        
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
        
        interpV1=interp1(t,arrV1,tdataC1);
        interpV2=interp1(t,arrV2,tdataC2);
        interpV3=interp1(t,arrV3,tdataC3);
              
        chi2C1=sum((interpV1(currVdataC1~=0)-currVdataC1(currVdataC1~=0)).^2./interpV1(currVdataC1~=0))./sum(currVdataC1(currVdataC1~=0));
        chi2C2=sum((interpV2(currVdataC2~=0)-currVdataC2(currVdataC2~=0)).^2./interpV2(currVdataC2~=0))./sum(currVdataC2(currVdataC2~=0));
        chi2C3=sum((interpV3(currVdataC3~=0)-currVdataC3(currVdataC3~=0)).^2./interpV3(currVdataC3~=0))./sum(currVdataC3(currVdataC3~=0));
        chi2=chi2C1+chi2C2+chi2C3;
        
        if (chi2<0)||isnan(chi2)
            
            disp('function created NaN value!');
        end
    end

    function [arrV1,arrV2,arrV3]=myPlot(params)
        
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

    function red_chi2=myRedC(params)

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

        interpV1=interp1(t,arrV1,tdataC1);
        interpV2=interp1(t,arrV2,tdataC2);
        interpV3=interp1(t,arrV3,tdataC3);

        chi2C1=sum(((interpV1(currVdataC1~=0)-currVdataC1(currVdataC1~=0))./VdataC1Err(currVdataC1~=0)).^2);
        chi2C2=sum(((interpV2(currVdataC2~=0)-currVdataC2(currVdataC2~=0))./VdataC1Err(currVdataC2~=0)).^2);
        chi2C3=sum(((interpV3(currVdataC3~=0)-currVdataC3(currVdataC3~=0))./VdataC1Err(currVdataC3~=0)).^2);


        red_chi2=(chi2C1+chi2C2+chi2C3)/(sum(currVdataC1~=0)+sum(currVdataC2~=0)+sum(currVdataC3~=0)-length(params));
    end
end
