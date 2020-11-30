
[files, path] = uigetfile('*.txt','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

all = zeros(118,1);
for condit = 1:file_count
    rdist = importdata([path files{condit}]);
    
    all(:,end+1) = rdist(2,:)'./sum(rdist(2,:)');
    
end

all(:,1) = rdist(1,:)';
all([1, end],:) = [];

figure; hold on
for i = 2:size(all, 2)
    plot(all(:,1), all(:,i), 'linewidth', 3)
end
xlim([0 1])
legend('show');

x = [0 10 25 45 60 90 120 180];
y = all(59,2:end);