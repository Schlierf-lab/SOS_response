%% load data

[files, path] = uigetfile('.mat');
file_location = [path files];
load(file_location);

%% scatter plot of Dvals vs weight fraction for the most probable model

states = mode(out.L(round(length(out.L)/2):end));
states = 4;

model_sel = out.L == states;
model_sel(1:round(length(out.L)/2)) = 0;
appDs = out.Dvals(model_sel);
probs = out.Pi(model_sel);

meanappD = mean([appDs{:}], 2);
[meanappDs, order] = sort(meanappD); %sort D and save the order

appD_final = zeros(length(appDs), 2*states);
for kk = 1:length(appDs)
    appD_final(kk,1:states) = appDs{kk}';
    appD_final(kk,states+1:end) = nonzeros(probs{kk})';
end

figure; hold on;
stat_data = zeros(states, 4);
for pp = order'
   scatter(appD_final(:,pp), appD_final(:,pp+states), 20, 'filled', 'MarkerFaceAlpha', 0.45)
   stat_data(pp,1) = mean(appD_final(:,pp));
   stat_data(pp,2) = std(appD_final(:,pp));
   stat_data(pp,3) = mean(appD_final(:,pp+states));
   stat_data(pp,4) = std(appD_final(:,pp+states));
end

title(files(1:end-4), 'Interpreter', 'none')
leg = legend(num2str(stat_data(order,1), '%.3f'), 'location', 'northwest');
legend('boxoff')
title(leg, 'mean \itD (\mum^2 s^{-1})')
set(gca,'xscale','log')
xlim([0.001 2.5])
ylim([0 0.7])
xlabel('Diffusion Coeffecient (\mum^2 s^{-1})')
ylabel('Weight Fraction')

%% transition matrix
figure;
allTM = cat(3, out.TransMat{model_sel}); %put all TMs into 3D array
meanTM = mean(allTM, 3); %calculate mean for all TMs in 3rd dimension
meanTM = meanTM(:,1:end-1); % remove weird zeros
meanTMs = meanTM(order, order); %order according to the sorted D array
imagesc(meanTMs)
title(['Transition Matrix, ' files(1:end-4)], 'Interpreter', 'none')
xlabel('transition to \itD (\mum^2 s^{-1})')
xticks(1:states)
xticklabels(num2str(meanappDs, '%.3f'));
ylabel('transition from \itD (\mum^2 s^{-1})')
yticks(1:states)
yticklabels(num2str(meanappDs, '%.3f'));
axis square
map = importdata('TM_map.txt');
colormap(map)
c = colorbar;
c.Label.String = 'transition probability';
caxis([0 1])
set(gca,'YDir','normal')

%% plot the # of states per iteration
figure
plot(out.L,'b','linewidth',2)
ylim([0 12])
title(['Mobility Sates, ' files(1:end-4)], 'Interpreter', 'none')
ylabel('Number of Sates')
xlabel('Iteration')


%% scatter plot of Dvals vs iteration
figure;
c = colormap('lines');
c = c(1:15, :);
plot_arr = [];
for ii=1:length(out.L)-1
    Dsort=sort(out.Dvals{ii});
    plot_arr = [plot_arr; [ii*ones(length(Dsort), 1), Dsort, c(1:length(Dsort), :)]];
end

scatter(plot_arr(:,1), plot_arr(:,2), 7, plot_arr(:,3:end), 'filled')
set(gca,'yscale','log')
axis tight
title(files(1:end-4), 'Interpreter', 'none')
ylabel('Diffusion Coeffecient, µm^2/s')
xlabel('Iterations')
ylim([0.001 2.5])

clear out;