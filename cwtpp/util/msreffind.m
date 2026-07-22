function [mzs_output, ms_references, ppm_dis_ini] = msreffind(mode, options)
    arguments
        mode (1,:) char
        options.path (1,:) = pwd
        options.threshold (1,:) {mustBeNumeric,mustBeReal} = 100
        options.wsdatasets = []
        options.wsmzs = []
    end
    
    datasets = {};mzs = {};
    if isempty(nargin)
        dname = uigetdir();
        cd (dname);
        filenames=dir('*.raw');
        datasets = {};mzs = {};
        for i = 1:length(filenames)
            filename = filenames(i).name;
            [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
%             [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
%             datasets = {};
        end
    else
        if strcmp(mode,'raw')
            cd (options.path);
            filenames=dir('*.raw');
            for i = 1:length(filenames)
                filename = filenames(i).name;
                [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
    %             [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
    %             datasets = {};
            end
        elseif strcmp(mode,'recal')
            disp('finding reference with recalibrated data.')
            filenames=[dir('*.raw');dir('*_recal.h5')];
            % filenames=[dir('*.raw');dir('*_filt.h5')];
            for i = 1:length(filenames)
                filename = filenames(i).name;
                [~,~,ext] = fileparts(filename);
                if strcmp(ext,'.raw')
                    [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'\datacube_recal.h5']);
                else
                    [datasets{i},dimes{i},mzs{i}]=h5toMat([filename]);
                end
    %             [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
    %             datasets = {};
            end
        elseif strcmp(mode,'workspace')
            datasets = options.wsdatasets;
            mzs = options.wsmzs;
        else
            filenames=dir('*.mz5');
            % disp('finding reference in mz5 files.')
            for i = 1:length(filenames)
                filename = filenames(i).name;filename = filename(1:end-4);
                % [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
                [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'.h5']);
            end
        end
    end
  
    for i = 1:length(datasets)
       lens(i) = size(datasets{i},2);
    end
    [~,I] = sort(lens,'descend');
    datasets = datasets(I);
    mzs = mzs(I);
    mzs_recal = mzs;

%     mz_recal = mzs{1};
    mz_recal = mzs_recal{1};
%     mzs_recal = mzs;
    ms_references = [];
    % if size(datasets{1},1) ~= 1
    %     [~,locs] = findpeaks(mean(datasets{1}),'MinPeakProminence',mean(mean(datasets{1}))/3);
    % else
        locs = 1:length(mz_recal);
    % end
    ms_references(:,1) = mz_recal(locs);

%     mzs_recal = {};
    for n = 2:length(mzs_recal)
        % disp(n)
        mz_raw = mzs_recal{n};
%         mz_raw = mzs{n};
        % if size(datasets{n},1) ~= 1
        %     [~,locs] = findpeaks(mean(datasets{n}),'MinPeakProminence',mean(mean(datasets{n}))/3);
        % else
            % [~,locs] = findpeaks((datasets{n}),'MinPeakProminence',mean(mean(datasets{n}))/3);
            locs = 1:length(mz_raw);
        % end
        mz_raw_p = mz_raw(locs);
        for m = 1:size(ms_references,1)
                % disp(m)
                [diff, ind] = min( abs(mz_raw_p-ms_references(m,1)) );
                ppm = diff/ms_references(m,1)*10^6;
                if ppm <= options.threshold
                    [~, ind] = min( abs(mz_raw-mz_raw_p(ind)) );
                    ms_references(m,n) = mz_raw(ind);
                else
                    ms_references(m,:)=0;
                end
        end     
    end
    ms_references(ms_references(:,1)==0,:)=[];
    [~,ia] = unique(ms_references(:,2));ms_references = ms_references(ia,:);
    ms_references = cat(2,ms_references,mean(ms_references,2));
    
    ini = zeros(length(mzs_recal),1);
    for i = 1:length(mzs_recal)
        % disp(i)
        dum = mzs_recal{i};
        [~,dumInd] = ismember(ms_references(:,i), dum);
        dum(dumInd) = ms_references(:,end);
        mzs_recal{i} = dum;
        dif = abs(mzs{i}-mzs_recal{i});
        xs = sqrt(mzs_recal{i});xs = xs(~(dif==0))';dif = dif(~(dif==0))';
        ini(i)  = (xs' * dif) / (xs' * xs);
        clear dum
    end
    mzs_output(I) = mzs_recal;
    ppm_dis_ini(I) = ini;
    % ms_references = ms_references(:,1);
end