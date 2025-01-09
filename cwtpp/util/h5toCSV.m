function [output] = h5toCSV()
try
    dname = uigetdir();
    cd (dname);
    filenames=dir('*.raw');

    if exist('metadata.csv','file') > 0
            metadata = readcell('metadata.csv');
            disp('metadata found in folder.')
    end

    data_matrix= [];
    labels = {};labels{1} = 'Class';labels = labels';
    % labels2 = {};labels2{1} = 'Tumour type';labels2 = labels2';
    files = {};files{1} = 'File'; files = files';
    start_scan = {};start_scan{1} = 'Start_scan'; start_scan = start_scan';
    end_scan = {};end_scan{1} = 'End_scan'; end_scan = end_scan';
    
    for i = 1:length(filenames)
        disp(i)
        filename = filenames(i).name;
        id = find(strcmp(metadata(:,2), filename));
        [dum,~,mz]=h5toMat([filename,'\datacube_aligned.h5']);
        scan_num = 1;
        for j = 1:length(id)
            num_scans = cell2mat(metadata(id(j),4))-cell2mat(metadata(id(j),3));
            start_scan = cat(1,start_scan,metadata(id(j),3));
            end_scan = cat(1,end_scan,metadata(id(j),4));
            labels = cat(1,labels,char(metadata(id(j),1)));
            % labels2 = cat(1,labels2,char(metadata(id(j),6)));
            files = cat(1,files,char(metadata(id(j),2)));
            if num_scans == 0
                data_matrix = cat(1,data_matrix,(dum(scan_num,:)));
            else
                data_matrix = cat(1,data_matrix,mean(dum(scan_num:scan_num+num_scans,:)));
            end
            scan_num = scan_num+num_scans;
        end
        clear dum
    end
    data_matrix = cat(1,mz,data_matrix);
    
    output = cat(2,end_scan,num2cell(data_matrix));
    output = cat(2,start_scan,output);
    output = cat(2,files,output);
    output = cat(2,labels,output);
    % output = cat(2,labels2,output);
    
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