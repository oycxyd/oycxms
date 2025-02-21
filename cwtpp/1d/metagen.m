csv1 = readcell('iKnife Metadata pending iA.xlsx');
csv2 = readcell('metadata.csv');

%% match ID and assign metadata
for i = 2:length(csv2(2:end,2))+1
% for i = 2:3
    % ID_i = split(csv2{i,2},'_');
    ID_i = csv2{i,2};
    match = find(strcmpi(csv1(:,3),ID_i));
    try
        csv2(i,1) = csv1(match,1);
        csv2(i,6) = csv1(match,2);
        csv2(i,7) = csv1(match,3);
        csv2(i,8) = csv1(match,4);
        csv2(i,9) = csv1(match,5);
        csv2(i,10) = csv1(match,6);
        csv2(i,11) = csv1(match,7);
        csv2(i,12) = csv1(match,8);
        csv2(i,13) = csv1(match,9);
    catch
        warning('metadata not found.')
        continue
    end
end

%% clean up missing values
for j = 1:size(csv2,2)
mask = cellfun(@ismissing, csv2(:,j), 'UniformOutput', false);
    for i = 1:length(mask)
        if cell2mat(mask(i)) == 1
            csv2{i,j} = 'TBD';
        end
    end
end

writecell(csv2,['metadata.csv'])
