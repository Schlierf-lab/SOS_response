function result=S2_simOptTc_measAct_varyMol_fast()

    addpath(genpath('scripts'));
    
    % parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DT_FRAME=20; % (ms) per frame
    NUM_TRACE=5000; % number of traces
    MAX_MOLECULES=100; % maximal number of molecules per trace
    STEP_MOLECULES=1;
    FIRST_FRAME=501; % first active frame
    
    % Toff
    kr1=12.8194; % (s^-1)
    kr2=1.7197; % (s^-1)
    prob_r1=0.6151;
    prob_r2=0.1895;
    
    % Ton
    k1=40.1297; % (s^-1)
    k2=6.9079; % (s^-1)
    prob1=0.9034;
   
    % blinking
    n=0.3048;
    tauC_blink=0.15; % (s)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    meanK=1/(prob1/k1+(1-prob1)/k2);
    p_tc=prob_r1/(prob_r1+prob_r2)*exp(-kr1*tauC_blink)+prob_r2/(prob_r1+prob_r2)*exp(-kr2*tauC_blink);
    corr_n=n/(p_tc+n*(1-p_tc));
    meanKd=meanK*corr_n; % (s^-1) 
    meanKb=meanK*(1-corr_n); % (s^-1) 

    % A -> state 1 | D -> state 2 | B -> state 3    
    result=[];
    
    % load activation distribution from location file
    [FILE,PATH]=uigetfile('*loc.txt','MultiSelect','on');
    
    if ~iscell(FILE)
        
       fileArr{1}=FILE; 
    else
        fileArr=FILE;        
    end
    
    numLocs=zeros(300000,1);
    maxFrame=300000;
    
    h=waitbar(0,'Loading activation distribution...');
    
    for iterFL=1:length(fileArr)
        
       locs=importdata([PATH fileArr{iterFL}]);
       frames=locs.data(1:end-1,1);        
       
        for iterF=1:max(frames)

            numLocs(iterF)=numLocs(iterF)+sum(frames==iterF);
        end     
        
        if max(frames)<maxFrame
            
            maxFrame=max(frames);
        end
        
        waitbar(iterFL/length(fileArr));
    end
    
    close(h);
    
    numLocs=numLocs(FIRST_FRAME:maxFrame);
    
    [cfunAD,gof,output]=fit((1:1:length(numLocs))',numLocs,'poly9');
        
    figure();
    subplot(1,2,1);
    plot(numLocs,'-k');
    hold on;
    plot((1:1:length(numLocs)),cfunAD((1:1:length(numLocs))),'LineWidth',2);
    xlabel('Frame');
    xlim([0 length(numLocs)]);
    ylabel('Number of active molecules');
    pause(0.1);
    
    framesFine=(1:1:length(numLocs));   
    
    % iterate number of molecules 
    for iterNM=2:STEP_MOLECULES:MAX_MOLECULES
    
        % Create random traces
        % drawing the activation frame       
        actF=randDISTR(framesFine,cfunAD(framesFine),NUM_TRACE*iterNM)+rand(NUM_TRACE*iterNM,1);

        % draw dwell times of state 1 (A) and 2 (D)
        dwell=zeros(100,NUM_TRACE*iterNM,2);
        
        % draw Ton from double exponential
        state1=binornd(1,prob1,100,NUM_TRACE*iterNM);
        mu1=exprnd(1/k1,100,NUM_TRACE*iterNM);
        mu2=exprnd(1/k2,100,NUM_TRACE*iterNM);
        mu=mu2;
        mu(state1==1)=mu1(state1==1);
        dwell(:,:,1)=mu;
%         dwell(:,:,1)=exprnd(1/(kd+kb),100,NUM_TRACE*iterNM);
        
        % draw Toff from double exponential
        state1=binornd(1,prob_r1/(prob_r1+prob_r2),100,NUM_TRACE*iterNM);
        mu1=exprnd(1/kr1,100,NUM_TRACE*iterNM);
        mu2=exprnd(1/kr2,100,NUM_TRACE*iterNM);       
        mu=mu2;
        mu(state1==1)=mu1(state1==1); 
        dwell(:,:,2)=mu;
%         dwell(:,:,2)=exprnd(1/kr,100,NUM_TRACE*iterNM); 

        clear state1 mu1 mu2 mu;
        
        % draw state transitions
        transS1=3-binornd(1,meanKd/(meanKd+meanKb),NUM_TRACE*iterNM*100,1);

        % save results in temporal arrays
        rndTon=1./zeros(NUM_TRACE*iterNM,100); % (frames)
        rndToff=1./zeros(NUM_TRACE*iterNM,100); % (frames)
       
        idxTS=1;
        
        for iterM=1:NUM_TRACE*iterNM

            idxOn=1;
            idxOff=1;
            state=1; % starting in the active state

            % assign dwell times to molecule
            while state~=3

                if state==1

                    rndTon(iterM,idxOn)=dwell(idxOn,iterM,1)*1000/DT_FRAME;
                    idxOn=idxOn+1;
                    
                    state=transS1(idxTS); % take transition from array
                    idxTS=idxTS+1;

                elseif state==2

                    rndToff(iterM,idxOff)=dwell(idxOff,iterM,2)*1000/DT_FRAME;
                    idxOff=idxOff+1;
                    
                    state=1; % only possible transition D -> A
                end
            end
        end
        
        % Get frames -> half exposed frames are active molecules
        molF=cell(NUM_TRACE*iterNM,1); 
        
        parfor iterM=1:NUM_TRACE*iterNM
            
            varON=rndTon(iterM,:);
            varOFF=rndToff(iterM,:);
                        
            currF=dwell2frames(actF(iterM),varON(varON~=Inf),varOFF(varOFF~=Inf)); 
                
            % draw again if frames are empty
            while isempty(currF)
                
                rndTonX=1./zeros(100,1); % (frames)
                rndToffX=1./zeros(100,1); % (frames)
                
                idxOn=1;
                idxOff=1;
                state=1; % starting in the active state
                
                % assign dwell times to molecule
                while state~=3
                    
                    if state==1

                        if binornd(1,prob1)==1
                            
                            rndTonX(idxOn)=exprnd(1/k1)*1000/DT_FRAME;
                        else
                            rndTonX(idxOn)=exprnd(1/k2)*1000/DT_FRAME;
                        end
%                         rndTonX(idxOn)=exprnd(1/(kd+kb))*1000/DT_FRAME;
                        idxOn=idxOn+1;                        
                        
                        state=3-binornd(1,meanKd/(meanKd+meanKb));
                        
                    elseif state==2
                        
                        if binornd(1,prob_r1)==1
                            
                            rndToffX(idxOff)=exprnd(1/kr1)*1000/DT_FRAME;
                        else
                            rndToffX(idxOff)=exprnd(1/kr2)*1000/DT_FRAME;
                        end

                        idxOff=idxOff+1;
                        
                        state=1; % only possible transition D -> A
                    end
                end
                
                varONX=rndTonX;
                varOFFX=rndToffX;
                
                currF=dwell2frames(actF(iterM),varONX(varONX~=Inf),varOFFX(varOFFX~=Inf));
            end
            
            molF{iterM}=currF;
        end
        
        % merge molecules
        assNum=reshape((1:1:NUM_TRACE*iterNM),[],iterNM);
        
        molF_NEW=cell(NUM_TRACE,1);
        
        for iterS=1:NUM_TRACE
            
            arrFS=[];
            
            for iterSL=1:iterNM
                
                arrFS=[arrFS;molF{assNum(iterS,iterSL)}];
            end
            
            molF_NEW{iterS}=sort(arrFS);
        end    

        % find optimal tC
        options=optimset('Display','iter','TolFun',47000);
        estimate=fminsearchbnd(@tauC2Nb,200,DT_FRAME,2000,options);
        
        % save results
        result=[result;iterNM estimate];
        
        save([PATH FILE{1}(1:end-24) 'SIM_TauC.mat'],'result');
        
        disp(['Num. Mol. ' num2str(iterNM) '| tauC=' num2str(estimate) 'ms']);
        
        subplot(1,2,2);
        plot(result(:,1),result(:,2),'-ok');
        xlabel('Number of molecules per trace');
        ylabel('\tau_{c,opt} (ms)');
        xlim([0 MAX_MOLECULES+1]);
        ylim([0 2000]);
        
        pause(0.1);
    end 
    
    function diffNb=tauC2Nb(params)
        
        tauC=params(1);
        
        arrNb=zeros(NUM_TRACE,1);
        idxNb=1;
        
        for iterN=1:NUM_TRACE
            
            subframes=unique(molF_NEW{iterN});
            
            gaps=setxor(subframes,(subframes(1):1:subframes(end)));
            
            if ~isempty(gaps)
            
                [gStart gLength]=burstLoc(gaps,1);
            else
                gLength=0;
            end
            
            arrNb(idxNb)=sum(gLength.*DT_FRAME>tauC)+1;
            
            idxNb=idxNb+1;
        end
        
        diffNb=abs(mean(arrNb)-iterNM);
    end
end


