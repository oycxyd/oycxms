function [data_aligned, mz_recal, I] = mia(path,mode,options)
    arguments
        path (1,:) = pwd
        mode (1,:) char = 'raw'
        options.mzlow (1,:) {mustBeNumeric,mustBeReal} = 50
        options.mzhigh (1,:) {mustBeNumeric,mustBeReal} = 1200
        options.matchppm (1,:) {mustBeNumeric,mustBeReal} = 100
        options.freq (1,:) {mustBeNumeric,mustBeReal} = 0.4
        options.save = 1
        options.autoparams = 1
        options.wsdatasets = []
        options.wsmzs = []
        options.verbose {mustBeNumericOrLogical} = true;
    end

%% Data import
    if isempty(path)
        dname = uigetdir();
        cd (dname);
        filenames=[dir('*.raw')];
        datasets = {};mzs = {};
        if isempty(filenames)
            filenames=[dir('*.mz5')];
            for i = 1:length(filenames)
                filename = filenames(i).name;filename = filename(1:end-4);
                % [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
                [~,dimes{i},mzs{i}]=h5toMat([filename,'.h5']);
    %             datasets = {};
            end
        else
            for i = 1:length(filenames)
                filename = filenames(i).name;
                % [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
                [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
    %             datasets = {};
            end
        end
    else
        dname = path;
        try
            cd (dname);
        catch
            mode = 'workspace';
        end
        filenames=dir('*.raw');
        datasets = {};mzs = {};
        
        if strcmp(mode, 'msp')
            if options.verbose
                disp('using resampled data.')
            end
            for i = 1:length(filenames)
                filename = filenames(i).name;
                [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube_msp.h5']);
%             [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
%             datasets = {};
            end
        elseif strcmp(mode, 'recal')
            if options.verbose
                disp('using recalibrated data.')
            end
            options.matchppm = 30;
            filenames=[dir('*.raw');dir('*_recal.h5')];
            for i = 1:length(filenames)
                filename = filenames(i).name;
                [~,~,ext] = fileparts(filename);
                if strcmp(ext,'.raw')
                    [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube_recal.h5']);
                    
                else
                    [~,dimes{i},mzs{i}]=h5toMat([filename]);
                end
%             [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
%             datasets = {};
            end
        elseif strcmp(mode, 'workspace')
            if options.verbose
                disp('using defined workspace data.')
            end
            datasets = options.wsdatasets;
            mzs = options.wsmzs;
            options.save = 0;
        elseif strcmp(mode, 'raw')
            if options.verbose
                disp('detected raw files.')
            end
            for i = 1:length(filenames)
                filename = filenames(i).name;
                [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
%             datasets = {};
            end
        else
            if options.verbose
                disp('detected mz5 files.')
            end
            filenames=[dir('*.mz5')];
            for i = 1:length(filenames)
                filename = filenames(i).name;filename = filename(1:end-4);
                % [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
                [~,dimes{i},mzs{i}]=h5toMat([filename,'.h5']);
    %             datasets = {};
            end
        end
    end

%% virtual lock mass correction
    if options.autoparams
        disp('autoparams on, searching optimal ppm threshold by MI...')
        [ppmthresh] = migsearch_batch( ...
        path=options.path, ...
        mode=options.mode, ...
        thresholds=[5 10 20 30 40 50 60 70 80 100 150 300], ...
        mzlow=50, ...
        mzhigh=1200, ...
        minVLMfreq=0.9, ...
        miBins=5, ...
        verbose=true);
        options.matchppm = ppmthresh;
    end
    
    for i = 1:length(mzs)
        lens(i) = length(mzs{i});
    end    
    mean_length=median(lens);
    % mean_length=max(lens);
    [~,ind1]=min(abs(lens-mean_length));
    I = 1:length(mzs);
    if strcmp(mode, 'recal')
        [mzs_recal, VLM, VLM_ind, ppm_dis_ini] = msreffind('recal');
    else
        if strcmp(mode, 'workspace')
            [mzs_recal, VLM, VLM_ind, ppm_dis_ini] = msreffind('workspace', ...
                wsmzs=mzs);
        elseif strcmp(mode, 'raw')
            [mzs_recal, VLM, VLM_ind, ppm_dis_ini] = msreffind('raw', ...
                threshold = options.matchppm);
        else
            [mzs_recal, VLM, VLM_ind, ppm_dis_ini] = msreffind('mz5', ...
                threshold = options.matchppm);
        end
    end
    
    dum2 = mzs_recal;
    dum2{1} = mzs_recal{ind1};
    dum2{ind1} = mzs_recal{1};
    mzs_recal = dum2;
    clear dum2
    I(1) = ind1;
    I(ind1) = 1;
    if options.verbose
        disp('using the following reference masses: ')
        disp([num2str(VLM(:,end))])
    end

%% MI-based alignment
    if options.verbose
        disp(['using a threshold of ',num2str(options.matchppm), ' ppm'])
    end
    mz_recal = cwt2cmz(mzs_recal,'minFreq',ceil(length(mzs)*options.freq), ...
        'mzRange',[options.mzlow options.mzhigh], ...
        'ppm',options.matchppm);
    % [~,dum_ind] = ismember(VLM(:,end),mz_recal);dum_ind = dum_ind(dum_ind~=0);
    % idx_not1 = setdiff(1:numel(mz_recal), dum_ind);
    if strcmp(mode,'workspace')
        [a,b] = maxent(mzs{1},datasets{1});
        [data1] = a*randn(size(datasets{1},1),length(mz_recal))+b;
        % [~,indd1,indd2] = intersect(round(mzs{1},4),round(mz_recal,4));
        % idx_not2 = setdiff(1:numel(mzs_recal{1}), VLM_ind(:,1));
        if strcmp(mode,'recal')
            % [~,indd1,indd2] = intersect(round(mzs{1},1),round(mz_recal,1));
            [indd1,indd2] = fixedmatch(mzs{1}, mz_recal,10^-(log10(median(diff(mz_recal))/min(diff(mz_recal)))/2));
        else
            [indd1,indd2] = adaptmatch(mzs_recal{1}, mz_recal,ppm_dis_ini);
        end
        dum = datasets{1};data1(:,indd2) = dum(:,indd1);
        % data1(:,dum_ind) = dum(:,VLM_ind(1,:));
        % [mz_recal_f,data1] = merge_peaks(mz_recal,data1,0.001);
        % [data1]=mzmatchfill(mzs{1},datasets{1},mz_recal,options.matchppm);
        data_aligned{1} = data1;
    else 
        if strcmp(mode,'mz5')
            filename = filenames(I(1)).name(1:end-4);
            [data1]=h5toMat([filename,'.h5']);
            out_dir = pwd;
        elseif strcmp(mode,'recal')
            filename = filenames(I(1)).name;
            if contains(filename,'.h5')
                [data1]=h5toMat([filename]);
                out_dir = pwd;
            else
                [data1]=h5toMat([filenames(I(1)).name,'\datacube_recal.h5']);
                out_dir = filenames(I(1)).name;
            end
        elseif strcmp(mode,'raw')
            [data1]=h5toMat([filenames(I(1)).name,'\datacube.h5']);
            out_dir = filenames(I(1)).name;
        end
        [a,b] = maxent(mzs{1},data1);
        [data11] = a*randn(size(data1,1),length(mz_recal))+b;
        % [~,indd1,indd2] = intersect(round(mzs{1},1),round(mz_recal,1));
        % idx_not2 = setdiff(1:numel(mzs_recal{1}), VLM_ind(:,1));
        if strcmp(mode,'recal')
            % [~,indd1,indd2] = intersect(round(mzs{1},1),round(mz_recal,1));
            [indd1,indd2] = fixedmatch(mzs{1}, mz_recal,10^-(log10(median(diff(mz_recal))/min(diff(mz_recal)))/2));
        else
            [indd1,indd2] = adaptmatch(mzs_recal{1}, mz_recal,ppm_dis_ini);
        end
        data11(:,indd2) = data1(:,indd1);
        % data11(:,VLM_ind(:,1)) = data1(:,dum_ind);
        % [mz_recal_f,data11] = merge_peaks(mz_recal,data11,0.001);
        data1 = data11;
        % [data1]=mzmatchfill(mzs{1},data1,mz_recal,options.matchppm);
    end
    
    % save
    if options.save == 1
        if strcmp(mode,'mz5')
            if isfile([filename,'_aligned.h5']) > 0
                delete([filename,'_aligned.h5'])
            end
            h5create([filename,'_aligned.h5'],'/datacube',size(data1))
            h5create([filename,'_aligned.h5'],'/mz',size(mz_recal));
            h5create([filename,'_aligned.h5'],'/dims',size(dimes{I(1)}));
            h5write([filename,'_aligned.h5'],'/datacube',data1)
            h5write([filename,'_aligned.h5'],'/mz',mz_recal)
            h5write([filename,'_aligned.h5'],'/dims',dimes{I(1)})
        elseif strcmp(mode,'recal')
            if contains(filename,'h5')
                filename = filename(1:end-3);
                if isfile([filename,'_aligned.h5']) > 0
                    delete([filename,'_aligned.h5'])
                end
                h5create([filename,'_aligned.h5'],'/datacube',size(data1))
                h5create([filename,'_aligned.h5'],'/mz',size(mz_recal));
                h5create([filename,'_aligned.h5'],'/dims',size(dimes{I(1)}));
                h5write([filename,'_aligned.h5'],'/datacube',data1)
                h5write([filename,'_aligned.h5'],'/mz',mz_recal)
                h5write([filename,'_aligned.h5'],'/dims',dimes{I(1)})
            else
                if isfile([out_dir,'/datacube_aligned.h5']) > 0
                    delete([out_dir,'/datacube_aligned.h5'])
                end
                h5create([out_dir,'/datacube_aligned.h5'],'/datacube',size(data1))
                h5create([out_dir,'/datacube_aligned.h5'],'/mz',size(mz_recal));
                h5create([out_dir,'/datacube_aligned.h5'],'/dims',size(dimes{I(1)}));
                h5write([out_dir,'/datacube_aligned.h5'],'/datacube',data1)
                h5write([out_dir,'/datacube_aligned.h5'],'/mz',mz_recal)
                h5write([out_dir,'/datacube_aligned.h5'],'/dims',dimes{I(1)})
            end
        elseif strcmp(mode,'raw')
            if isfile([out_dir,'/datacube_aligned.h5']) > 0
                delete([out_dir,'/datacube_aligned.h5'])
            end
            h5create([out_dir,'/datacube_aligned.h5'],'/datacube',size(data1))
            h5create([out_dir,'/datacube_aligned.h5'],'/mz',size(mz_recal));
            h5create([out_dir,'/datacube_aligned.h5'],'/dims',size(dimes{I(1)}));
            h5write([out_dir,'/datacube_aligned.h5'],'/datacube',data1)
            h5write([out_dir,'/datacube_aligned.h5'],'/mz',mz_recal)
            h5write([out_dir,'/datacube_aligned.h5'],'/dims',dimes{I(1)})
        end
    end
    
    for n = 1:length(mzs)-1
        if options.verbose
            disp(['aligining file ',num2str(n),'/',num2str(length(mzs)-1)])
        end
        if strcmp(mode,'workspace')
            data2 = datasets{n+1};
        else
            if strcmp(mode,'mz5')
                filename = filenames(I(n+1)).name;filename = filename(1:end-4);
                data2 = h5toMat([filename,'.h5']);
            elseif strcmp(mode,'recal')
                filename = filenames(I(n+1)).name;
                if contains(filename,'.h5')
                    [data2]=h5toMat([filename]);
                else
                    data2 = h5toMat([filename,'\datacube_recal.h5']);
                    out_dir = filename;
                end
            elseif strcmp(mode,'raw')
                filename = filenames(I(n+1)).name;
                data2 = h5toMat([filename,'\datacube.h5']);
                out_dir = filename;
            end
        end

        mz2 = mzs_recal{n+1};
        [a,b] = maxent(mz2,data2);
        data_mia = a*randn(size(data2,1),length(mz_recal))+b;
        % [~,indd1,indd2] = intersect(round(mz2,4),round(mz_recal,4));
        % idx_not2 = setdiff(1:numel(mz2), VLM_ind(:,n+1));
        if strcmp(mode,'recal')
            % [~,indd1,indd2] = intersect(round(mz2,1),round(mz_recal,1));
            [indd1,indd2] = fixedmatch(mz2, mz_recal, 10^-(log10(median(diff(mz_recal))/min(diff(mz_recal)))/2));
        else
            [indd1,indd2] = adaptmatch(mz2, mz_recal,ppm_dis_ini);
        end
        data_mia(:,indd2) = data2(:,indd1);
        % data_mia(:,dum_ind) = data2(:,VLM_ind(:,n+1));
        % [~,data_mia] = merge_peaks(mz_recal,data_mia,0.001);
        
        % save
        if options.save == 1        
%                 disp(n)
            disp('saving...');
            if strcmp(mode,'mz5')
                if isfile([filename,'_aligned.h5']) > 0
                    disp('Found existing aligned file. Saving over.')
                    delete([filename,'_aligned.h5'])
                end
                h5create([filename,'_aligned.h5'],'/datacube',size(data_mia))
                h5create([filename,'_aligned.h5'],'/mz',size(mz_recal));
                h5create([filename,'_aligned.h5'],'/dims',size(dimes{I(n+1)}));
                h5write([filename,'_aligned.h5'],'/datacube',data_mia)
                h5write([filename,'_aligned.h5'],'/mz',mz_recal)
                h5write([filename,'_aligned.h5'],'/dims',dimes{I(n+1)})
            elseif strcmp(mode,'recal')
                if contains(filename,'.h5')
                    filename = filename(1:end-3);
                    if isfile([filename,'_aligned.h5']) > 0
                        disp('Found existing aligned file. Saving over.')
                        delete([filename,'_aligned.h5'])
                    end
                    h5create([filename,'_aligned.h5'],'/datacube',size(data_mia))
                    h5create([filename,'_aligned.h5'],'/mz',size(mz_recal));
                    h5create([filename,'_aligned.h5'],'/dims',size(dimes{I(n+1)}));
                    h5write([filename,'_aligned.h5'],'/datacube',data_mia)
                    h5write([filename,'_aligned.h5'],'/mz',mz_recal)
                    h5write([filename,'_aligned.h5'],'/dims',dimes{I(n+1)})
                else
                    if isfile([out_dir,'/datacube_aligned.h5']) > 0
                        disp('Found existing aligned file. Saving over.')
                        delete([out_dir,'/datacube_aligned.h5'])
                    end
                    h5create([out_dir,'/datacube_aligned.h5'],'/datacube',size(data_mia))
                    h5create([out_dir,'/datacube_aligned.h5'],'/mz',size(mz_recal));
                    h5create([out_dir,'/datacube_aligned.h5'],'/dims',size(dimes{I(n+1)}));
                    h5write([out_dir,'/datacube_aligned.h5'],'/datacube',data_mia)
                    h5write([out_dir,'/datacube_aligned.h5'],'/mz',mz_recal)
                    h5write([out_dir,'/datacube_aligned.h5'],'/dims',dimes{I(n+1)})
                end
            elseif strcmp(mode,'raw')
                if isfile([out_dir,'/datacube_aligned.h5']) > 0
                    disp('Found existing aligned file. Saving over.')
                    delete([out_dir,'/datacube_aligned.h5'])
                end
                h5create([out_dir,'/datacube_aligned.h5'],'/datacube',size(data_mia))
                h5create([out_dir,'/datacube_aligned.h5'],'/mz',size(mz_recal));
                h5create([out_dir,'/datacube_aligned.h5'],'/dims',size(dimes{I(n+1)}));
                h5write([out_dir,'/datacube_aligned.h5'],'/datacube',data_mia)
                h5write([out_dir,'/datacube_aligned.h5'],'/mz',mz_recal)
                h5write([out_dir,'/datacube_aligned.h5'],'/dims',dimes{I(n+1)})
            end
        else
            data_aligned{n+1} = data_mia;
            % disp('workspace check')
        end
        clear data_mia
    end
function [a,b,ent] = maxent(mz,data)
    ent = [];
    for i = 1:length(mz)
        ent(i) = entropy(data(:,i));
    end
    [~,ii] =  max(ent);
    a = sqrt(std(data(:,ii)))/100;b = mean(data(:,ii))/100;

function [ind1,ind2] = adaptmatch(mz, mz_recal,ini)
    for n = 1:length(mz_recal)
        [dif(n),ind1(n)] = min(abs(mz-mz_recal(n)));
    end
    ppm_dis = ini*sqrt(mz_recal);
    margin = 1.5;
    matched = (dif<= ppm_dis*margin);
    ind1 = ind1(matched);[ind1,iu] = unique(ind1);
    ind2 = (1:n);ind2 = ind2(matched);ind2 = ind2(iu);

function [ind1,ind2] = fixedmatch(mz, mz_recal, thresh)
    for n = 1:length(mz_recal)
        [dif(n),ind1(n)] = min(abs(mz-mz_recal(n)));
    end
    matched = (dif<=thresh);
    ind1 = ind1(matched);[ind1,iu] = unique(ind1);
    ind2 = (1:n);ind2 = ind2(matched);ind2 = ind2(iu);