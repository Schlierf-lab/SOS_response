% plot raw data from D* and MSD analysis

[files, path] = uigetfile('*.mat','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

arrD_tot = [];
MJD_tot = [];
ncells_tot = 0;
ntracks_tot = 0;
for rep = 1:file_count
    load([path files{rep}]);
    arrD_tot = [arrD_tot arrD];
    MJD_tot = [MJD_tot; mean_disp];
    ncells_tot = ncells_tot + ncells;
    ntracks_tot = ntracks_tot + ntracks;
end

%% plot 

% plot D*
binsize = 0.1;
binsize_MJD = 5;
 
eappD=(-5:binsize:1.5);
happD=histc(log10(arrD_tot),eappD);
figure('Position' ,[356.2000  533.8000  196.0000  125.6000]); 
%bar(eappD+0.5*binsize,happD./(sum(happD)),1,'linestyle', 'none');
stairs(eappD,happD./(sum(happD)));
set(gca,'Color','none','FontName','Arial', 'box', 'off');
xlim([-3.0 1.5])
ylim([0 0.1])
%xlabel('log_{10}(D*) (µm^2s^{-1})')
%ylabel('normalized amplitude');
disp(['n_{tracks} = ' num2str(ntracks_tot) ', n_{cells} = ' num2str(ncells_tot)])

MJDw = zeros(length(MJD_tot)*10, 1);
MJDw(1:length(MJD_tot)) = MJD_tot(:,1);
currpos = length(MJD_tot) + 1;
for i = 1:length(MJD_tot)
    if MJD_tot(i,2) > 1
        weight = repmat(MJD_tot(i,1), MJD_tot(i,2)-1, 1);
        MJDw(currpos:currpos+MJD_tot(i,2)-2) = weight;
        currpos = currpos + MJD_tot(i,2)-2;
    end
end
MJDw = nonzeros(MJDw);

figure('Position' ,[554.6000  533.8000  196.0000  125.6000]);
hold on;
edges = (0:binsize_MJD:700);
hMJD = histc(MJDw, edges);

stairs(edges, hMJD./sum(hMJD), 'color', 'r')
set(gca,'Color','none','FontName','Arial', 'box', 'off');
%plot(edges + 0.5*binsize_MJD, hMJD./sum(hMJD))
xlim([0 350]);
ylim([0 0.06]);
% xlabel('Mean Displacement (nm)')
% ylabel('frequency')
% title(['weighted, n = ' num2str(length(mean_disp))])