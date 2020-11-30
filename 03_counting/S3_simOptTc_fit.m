% Fit of the optimal blinking tolerance time vs. underlying number of molecules
% --Andreas Hartmann

clear all;

% load file
[FILE,PATH]=uigetfile('*SIM_TauC.mat');
load([PATH FILE]);

% fit distribution
x0=[2000 0 -1 0];
% lb=[0 0 0 0];
 
% [cfun,gof,output]=fit(result(:,1),result(:,2),'a*exp(-b*x)+c*exp(-d*x)+f*exp(-g*x)','Start',x0,'Lower',lb,'Weights',1./result(:,2));
% [cfun,gof,output]=fit(result(:,1),result(:,2),'a*(x-b)^c+d','Start',x0,'Lower',lb,'Weights',1./result(:,2));
% [cfun,gof,output]=fit(result(:,1),result(:,2),'a*(x-b)^c+d','Start',x0,'Lower',lb);
[cfun,gof,output]=fit(result(:,1),result(:,2),'a*(x-b)^c+d','Start',x0);

figure();
s1=subplot(2,1,1);
plot(result(:,1),result(:,2),'ok');
hold on;
p1=plot(cfun);
set(p1,'LineWidth',2);
xlim([0 100]);
ylim([0 2000]);
xlabel('');
ylabel('\it t\rm_{c,opt} (ms)');
set(s1,'XTickLabels',{''});
legend('data','a*(x-b)^c+d');

subplot(2,1,2);
plot(result(:,1),(result(:,2)-cfun(result(:,1)))./sqrt(result(:,2)),'-k');
xlim([0 100]);
xlabel('Number of molecules per trace');
ylabel('w. residual');

params=[cfun.a cfun.b cfun.c cfun.d];
save([PATH FILE(1:end-12) 'SIM_Fit.mat']);

disp(cfun);
