clear all;

%% parameters
DT_FRAME=20; % (ms) time per frame

%% Load data
[FILE,PATH]=uigetfile('*.mat','MultiSelect','on');

addpath('scripts');

if ~iscell(FILE)
    
    fileArr{1}=FILE;
else
    fileArr=FILE;
end

ton=[];
total_ton=[];
total_toff=[];

for iterFL=1:length(fileArr)

    load([PATH fileArr{iterFL}]);

    ton=[ton;output.data{5}];
    total_ton=[total_ton;output.data{9}];
    total_toff=[total_toff;output.data{10}];
end

%% Show result
binsize=1;
bincenter=(1:binsize:37);
edges=[bincenter bincenter(end)+1]-0.5*binsize;

% on-time histogram
hTon=histcounts(ton,edges);
pTon=hTon./sum(hTon);

[cfunOn2,~,~]=fit(bincenter(hTon~=0)',log(hTon(hTon~=0)'),'log(a*exp(-b*x)+c*exp(-d*x))','Start',[0.7 0.01 0.3 1],'Lower',[0 0 0 0]);
ciOn2=confint(cfunOn2);
meanTau=mle(ton.*DT_FRAME./1000,'distribution','exp');
errmeanTau=std(ton.*DT_FRAME./1000)/sqrt(length(ton));
funcMT=exp(-bincenter.*DT_FRAME./1000./meanTau);
funcMT=funcMT./sum(funcMT).*sum(hTon);

% Show result
figure();
s1=semilogy(bincenter.*DT_FRAME./1000,hTon,'ok');
set(s1, 'markerfacecolor', [0 0 0],'MarkerSize',4);
hold on;
plot((0.5:0.01:37).*DT_FRAME./1000,exp(cfunOn2((0.5:0.01:37))),'-r','LineWidth',2);
xlim([0 0.7]);
xlabel('T_{on} (s)');
ylabel('Number of events');

a_=cfunOn2.a/cfunOn2.b;
c_=cfunOn2.c/cfunOn2.d;

P1=a_/(a_+c_);
P2=c_/(a_+c_);

% error propagation
K1=cfunOn2.b/DT_FRAME*1000; % (s^-1)
errK1=(ciOn2(2,2)-ciOn2(1,2))/4./DT_FRAME.*1000; % (s^-1)
err_a=(ciOn2(2,1)-ciOn2(1,1))/4;
K2=cfunOn2.d/DT_FRAME*1000; % (s^-1)
errK2=(ciOn2(2,4)-ciOn2(1,4))/4./DT_FRAME.*1000; % (s^-1)
err_c=(ciOn2(2,3)-ciOn2(1,3))/4;

randA=(err_a.*randn(10000,1)+cfunOn2.a);
randK1=errK1.*randn(10000,1)+K1;
randC=(err_c.*randn(10000,1)+cfunOn2.c);
randK2=errK2.*randn(10000,1)+K2;

randP1=(randA./randK1)./((randA./randK1)+(randC./randK2));
errP1=std(randP1);
randP2=(randC./randK2)./((randA./randK1)+(randC./randK2));
errP2=std(randP2);

meanK=1/(P1/K1+P2/K2); % (s^-1)

arr_meanK=zeros(100,1);

for iter=1:length(arr_meanK)

    curr_ton=randDISTR(bincenter,hTon,10000);
    
    arr_meanK(iter)=1/mean(curr_ton.*DT_FRAME./1000);
end

err_meanK=std(arr_meanK);

disp(['k_1=(' num2str(K1) '+-' num2str(errK1) ')s^-1 | k_2=(' num2str(K2) '+-' num2str(errK2) ')s^-1']);
disp(['P1=(' num2str(P1) '+-' num2str(errP1)]);
disp(['P2=(' num2str(P2) '+-' num2str(errP2)]);
disp(['<k>=(' num2str(meanK) '+-' num2str(err_meanK) ')s^-1']);

