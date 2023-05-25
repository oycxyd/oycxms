function [output] = h5toCSV()
try
    dname = uigetdir();
    cd (dname);
    filenames=dir('*.raw');

    if exist('metadata.csv','file') > 0
            metadata = readcell('metadata.csv');
            use_metadata = 1;
            disp('metadata found in folder.')
    end

    data_matrix= [];
    for i = 1:length(filenames)
        filename = filenames(i).name;
        [dum,~,mz]=h5toMat([filename,'\datacube_aligned.h5']);
        data_matrix = cat(1,data_matrix,mean(dum));
        clear dum
    end
    data_matrix = cat(1,mz,data_matrix);

    output = cat(2,metadata(:,1:3),num2cell(data_matrix));
    writecell(output,'aligned_for_postprocessing.csv')
catch
    warning('no preprocessed files found?');
end
% 
% cells = unique(metadata(:,3));
% GD = metadata(:,3);
% GD{1} = 'GD2 level';
% 
% i =5;
% ind = find(contains(metadata(:,3),cells(i)));
% disp(cells(i))
% GD(ind) = cellstr('negative');
% 
% output = cat(2,GD,output);