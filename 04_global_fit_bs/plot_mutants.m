
x = [0 10 20 30 40 50 60 70 80 90 100 110 120];
max_ost = 1.5;


file = 'nmutants_res.mat';
% cip     
color = colormap('lines');
load(file);

figure(1); hold on
plot(x(1:size(arrP1, 1)), allP1, 's', 'color', 'r', 'markerfacecolor', 'r',  'linewidth', 2)

figure(2); hold on
plot(x(1:size(arrP1, 1)), allP2, 's', 'color', 'r', 'markerfacecolor', 'r',  'linewidth', 2)

figure(3); hold on
plot(x(1:size(arrP1, 1)), allP3, 's', 'color', 'r', 'markerfacecolor', 'r',  'linewidth', 2)

figure(4); hold on
plot(x(1:size(arrP1, 1)), 1-allP1-allP2-allP3, 's', 'color', 'r', 'markerfacecolor', 'r',  'linewidth', 2)

figure(5); hold on
plot(x(1:size(arrP1, 1)), allP1, 's', 'color', [0 0 0], 'markerfacecolor', [0 0 0],  'linewidth', 2)
plot(x(1:size(arrP1, 1)), allP2, 's', 'color', [1 0 0], 'markerfacecolor', [1 0 0],  'linewidth', 2)
plot(x(1:size(arrP1, 1)), allP3, 's', 'color', [0 1 0], 'markerfacecolor', [0 1 0],  'linewidth', 2)
plot(x(1:size(arrP1, 1)), 1-allP1-allP2-allP3, 's', 'color', [0 0 1], 'markerfacecolor', [0 0 1],  'linewidth', 2)

for i = 1:size(arrP1, 1)
    for j = 1:length(arrP1(i,arrP1(i,:)~=-1))
        figure(1)
        plot(x(i)+(rand*2*max_ost-max_ost), arrP1(i, j), 'o', 'color', [0 0 0], 'markersize', 3, 'HandleVisibility', 'off')
        figure(2)
        plot(x(i)+(rand*2*max_ost-max_ost), arrP2(i, j), 'o', 'color', [0 0 0], 'markersize', 3, 'HandleVisibility', 'off')
        figure(3)
        plot(x(i)+(rand*2*max_ost-max_ost), arrP3(i, j), 'o', 'color', [0 0 0], 'markersize', 3, 'HandleVisibility', 'off')
        figure(4)
        plot(x(i)+(rand*2*max_ost-max_ost), 1-arrP1(i, j)-arrP2(i, j)-arrP3(i, j), 'o', 'color', [0 0 0], 'markersize', 3, 'HandleVisibility', 'off')
        figure(5)
        plot(x(i)+(rand*2*max_ost-max_ost), arrP1(i, j), 'o', 'color', [0 0 0], 'markersize', 3, 'HandleVisibility', 'off')
        plot(x(i)+(rand*2*max_ost-max_ost), arrP2(i, j), 'o', 'color', [1 0 0], 'markersize', 3, 'HandleVisibility', 'off')
        plot(x(i)+(rand*2*max_ost-max_ost), arrP3(i, j), 'o', 'color', [0 1 0], 'markersize', 3, 'HandleVisibility', 'off')
        plot(x(i)+(rand*2*max_ost-max_ost), 1-arrP1(i, j)-arrP2(i, j)-arrP3(i, j), 'o', 'color', [0 0 1], 'markersize', 3, 'HandleVisibility', 'off')
    end
end




figure(1)
set(gcf, 'position', [359.4000  509.8000  468.0000  252.0000])
xlim([-5 75])
xlabel('incubation time (min)')
ylabel('probability')
title('target-bound')
xticklabels({'wt', 'S119A', 'E45K', '3x', 'wtst', 'S119Ast', 'E45Kst', '3xst'})
xticks(x(1:size(arrP1, 1)))

figure(2)
set(gcf, 'position', [359.4000  509.8000  468.0000  252.0000])
xlim([-5 75])
title('DNA-bound')
xlabel('incubation time (min)')
ylabel('probability')
xticklabels({'wt', 'S119A', 'E45K', '3x', 'wtst', 'S119Ast', 'E45Kst', '3xst'})
xticks(x(1:size(arrP1, 1)))

figure(3)
set(gcf, 'position', [359.4000  509.8000  468.0000  252.0000])
xlim([-5 75])
title('cytoplasmic')
xlabel('incubation time (min)')
ylabel('probability')
xticklabels({'wt', 'S119A', 'E45K', '3x', 'wtst', 'S119Ast', 'E45Kst', '3xst'})
xticks(x(1:size(arrP1, 1)))

figure(4)
set(gcf, 'position', [359.4000  509.8000  468.0000  252.0000])
xlim([-5 75])
title('cleaved')
xlabel('incubation time (min)')
ylabel('probability')
xticklabels({'wt', 'S119A', 'E45K', '3x', 'wtst', 'S119Ast', 'E45Kst', '3xst'})
xticks(x(1:size(arrP1, 1)))

figure(5)
set(gcf, 'position', [359.4000  509.8000  468.0000  252.0000])
xlim([-5 75])
title('all states')
xlabel('incubation time (min)')
ylabel('probability')
xticklabels({'wt', 'S119A', 'E45K', '3x', 'wtst', 'S119Ast', 'E45Kst', '3xst'})
xticks(x(1:size(arrP1, 1)))
% set(gca, 'yscale', 'log')

