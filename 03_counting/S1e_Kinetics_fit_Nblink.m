clear all;

%% parameters
DT_FRAME=20; % (ms) time per frame

%% Load data
[FILE,PATH]=uigetfile('*.mat','MultiSelect','on');

if ~iscell(FILE)
    
    fileArr{1}=FILE;
else
    fileArr=FILE;
end

Nblink=[];

for iterFL=1:length(fileArr)

    load([PATH fileArr{iterFL}]);

    Nblink=[Nblink;output.data{7}];
end

%% Show result
binsize=1;
bincenter=(0:binsize:10);
edges=[bincenter bincenter(end)+1]-0.5*binsize;

% on-time histogram
hNblink=histcounts(Nblink,edges);
pNblink=hNblink./sum(hNblink);

% [cfunNB,~,~]=fit(bincenter',pNblink','(1-a)*(a^x)','Start',0.3,'Lower',0,'Weights',1./hNblink);
[cfunNB,~,~]=fit(bincenter',pNblink','(1-a)*(a^x)','Lower',0);
ciNB=confint(cfunNB);

figure();
s1=bar(bincenter,pNblink,'hist');
set(s1,'facecolor',[0.1 0.5 1]);
hold on;
plot((0:0.1:10),cfunNB((0:0.1:10)),'-r','LineWidth',2);
xlim([-1 11]);
xlabel('N_{blink}');
ylabel('Probability');
title('\eta^n(1-\eta)');
legend('Data',['\eta=(' num2str(cfunNB.a) '\pm' num2str((ciNB(2,1)-ciNB(1,1))/4) ')']);

disp(['n=(' num2str(cfunNB.a) '\pm' num2str((ciNB(2,1)-ciNB(1,1))/4) ')']);