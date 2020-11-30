%% LB
x = [0.01 0.1 0.5 0.75 1.5 3 4 20];

OD_wt = [6.86	6.22	6.11	4.64	2.36	0.728	0.085	0.099];
OD_tg = [6.23	5.01	5.48	4.69	1.83	0.26	0.09	0.084];
OD_mt = [6.53	4.67	0.87	0.459	0.221	0.118	0.066	0.055];

OD_wt = OD_wt./OD_wt(1);
OD_tg = OD_tg./OD_tg(1);
OD_mt = OD_mt./OD_mt(1);

figure('position', [693.8000  598.6000  278.4000  156.8000]); hold on
plot(x, OD_wt, '-x', 'linewidth', 2)
plot(x, OD_tg, '-x', 'linewidth', 2)
plot(x, OD_mt, '-x', 'linewidth', 2)
legend('LexA', 'LexA-PAmCherry', 'LexAS119A-PAmCherry')
set(gca, 'xscale', 'log')
title('LB')
ylim([0 1])
xlabel('cipro concentration')
ylabel('OD relative to untreated')

%% M9

OD_wt = [1.110	1.047	0.499	0.269	0.092	0.010	0.006	0.000];
OD_tg = [1.183	1.094	0.651	0.308	0.045	0.022	0.020	0.004];
OD_mt = [0.799	0.694	0.122	0.026	0.024	0.022	0.020	0.002];

OD_wt = OD_wt./OD_wt(1);
OD_tg = OD_tg./OD_tg(1);
OD_mt = OD_mt./OD_mt(1);

figure('position', [693.8000  598.6000  278.4000  156.8000]); hold on
plot(x, OD_wt, '-x', 'linewidth', 2)
plot(x, OD_tg, '-x', 'linewidth', 2)
plot(x, OD_mt, '-x', 'linewidth', 2)
%legend('LexA', 'LexA-PAmCherry', 'LexAS119A-PAmCherry')
set(gca, 'xscale', 'log')
title('M9')
ylim([0 1])
xlabel('cipro concentration')
ylabel('OD relative to untreated')