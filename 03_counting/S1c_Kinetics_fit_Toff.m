% Analysis of the off-time distribution
% --Andreas Hartmann

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

toff=[];

for iterFL=1:length(fileArr)

    load([PATH fileArr{iterFL}]);

    toff=[toff;output.data{6}];
end

%% Show result
binsize=1;
bincenter=(1:binsize:500);
edges=[bincenter bincenter(end)+1]-0.5*binsize;

% off-time histogram
hToff=histcounts(toff,edges);
pToff=hToff./sum(hToff);

[cfunOff,~,~]=fit(bincenter(hToff~=0)',log(hToff(hToff~=0)'),'log(a*exp(-b*x)+c*exp(-d*x)+f*exp(-g*x))','Start',[0.7 0.01 0.2 0.1 0.1 1],'Lower',[0 0 0 0 0 0]);
ciOff=confint(cfunOff);

figure();
s1=semilogy(bincenter.*DT_FRAME./1000,hToff,'ok');
set(s1, 'markerfacecolor', [0 0 0],'MarkerSize',4);
hold on;
plot((0.5:0.01:500).*DT_FRAME./1000,exp(cfunOff((0.5:0.01:500))),'-r','LineWidth',2);
xlabel('T_{off} (s)');
ylabel('Number of events');
title('(a_1*exp(-k_{r1}*x)+a_2*exp(-k_{r2}*x)+a_3*exp(-k_{r3}*x)');
legend('Data',['k_{r1}=(' num2str(cfunOff.b/DT_FRAME*1000) '\pm' num2str((ciOff(2,2)-ciOff(1,2))/4/DT_FRAME*1000) ')s^{-1} | k_{r2}=(' num2str(cfunOff.d/DT_FRAME*1000) '\pm' num2str((ciOff(2,4)-ciOff(1,4))/4/DT_FRAME*1000) ')s^{-1} | k_{r3}=(' num2str(cfunOff.g/DT_FRAME*1000) '\pm' num2str((ciOff(2,6)-ciOff(1,6))/4/DT_FRAME*1000) ')s^{-1}']);

a_=cfunOff.a/cfunOff.b;
c_=cfunOff.c/cfunOff.d;
f_=cfunOff.f/cfunOff.g;

P1=a_/(a_+c_+f_);
P2=c_/(a_+c_+f_);
P3=f_/(a_+c_+f_);

% error propagation
err_a=(ciOff(2,1)-ciOff(1,1))/4;
err_b=(ciOff(2,2)-ciOff(1,2))/4;
err_c=(ciOff(2,3)-ciOff(1,3))/4;
err_d=(ciOff(2,4)-ciOff(1,4))/4;
err_f=(ciOff(2,5)-ciOff(1,5))/4;
err_g=(ciOff(2,6)-ciOff(1,6))/4;

randA=(err_a.*randn(1000,1)+cfunOff.a);
randB=(err_b.*randn(1000,1)+cfunOff.b)./DT_FRAME.*1000;
randC=(err_c.*randn(1000,1)+cfunOff.c);
randD=(err_d.*randn(1000,1)+cfunOff.d)./DT_FRAME.*1000;
randF=(err_f.*randn(1000,1)+cfunOff.f);
randG=(err_g.*randn(1000,1)+cfunOff.g)./DT_FRAME.*1000;

randP1=(randA./randB)./((randA./randB)+(randC./randD)+(randF./randG));
errP1=std(randP1);
randP2=(randC./randD)./((randA./randB)+(randC./randD)+(randF./randG));
errP2=std(randP2);
randP3=(randF./randG)./((randA./randB)+(randC./randD)+(randF./randG));
errP3=std(randP3);

disp(['kr1=(' num2str(cfunOff.b/DT_FRAME*1000) '+-' num2str((ciOff(2,2)-ciOff(1,2))/4/DT_FRAME*1000) ')s^-1']);
disp(['kr2=(' num2str(cfunOff.d/DT_FRAME*1000) '+-' num2str((ciOff(2,4)-ciOff(1,4))/4/DT_FRAME*1000) ')s^-1']);
disp(['kr3=(' num2str(cfunOff.g/DT_FRAME*1000) '+-' num2str((ciOff(2,6)-ciOff(1,6))/4/DT_FRAME*1000) ')s^-1']);
disp(['P1=(' num2str(P1) '+-' num2str(errP1) ')']);
disp(['P2=(' num2str(P2) '+-' num2str(errP2) ')']);
disp(['P3=(' num2str(P3) '+-' num2str(errP3) ')']);