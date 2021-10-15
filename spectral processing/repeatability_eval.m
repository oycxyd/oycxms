%% repeatability evaluation with pairwise correlation between data from same cell line
% load csv
[filename,filepath] = uigetfile;
file = readcell([filepath,filename]);

sample_data = file(2:end,6:end);
sample_data = cell2mat(sample_data);
sample_labels = file(2:end,1);

list = unique(sample_labels);
pearsons = [];

% calculate pairwise correlations
for i = 1:length(list)
    select = strcmp(sample_labels, list{i});
    specs_select = sample_data(select,:);
    dist = 1-pdist(specs_select,'correlation');
    pearsons(i,1) = mean(dist);
    pearsons(i,2) = std(dist);
end
pearsons = cat(2,list,num2cell(pearsons));
writecell(pearsons,[filepath,'repeatabilities.csv'])
disp('results saved.')
