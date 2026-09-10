function [ppmthresh,metric, feature_scale, MIscore] = migsearch(mzs, datasets, options)
    arguments
        mzs
        datasets
        options.filenames = []
        options.mzlow (1,:) {mustBeNumeric,mustBeReal} = 50
        options.mzhigh (1,:) {mustBeNumeric,mustBeReal} = 1200
        options.freq (1,:) {mustBeNumeric,mustBeReal} = 0.4
    end
    range = [5, 10, 20, 30, 40, 50, 60, 70, 80, 100, 150, 300];
    mean_specs = {};
    filenames = options.filenames;
    if isempty(datasets)
        for n = 1:length(filenames)
            try
                mode = 'mz5';
                filename = filenames(n).name;filename = filename(1:end-4);
                [dataset_n]=h5toMat([filename,'.h5']);
            catch
                mode = 'raw';
                filename = filenames(n).name;
                [dataset_n]=h5toMat([filename,'\datacube.h5']);
            end
            mean_specs{n} = mean(dataset_n);
        end
        clear dataset_n
    else
        mode = 'workspace';
        for n = 1:length(datasets)
            mean_specs{n} = mean(datasets{n});
        end
    end

    %% Find references as virtual lock masses  
    for i = 1:length(mzs)
        lens(i) = length(mzs{i});
    end    
    mean_length=median(lens);
    % mean_length=max(lens);
    [~,ind1]=min(abs(lens-mean_length));
    I = 1:length(mzs);

    % tic
    metric = zeros(length(mean_specs)-1,length(range));
    nfeatures = zeros(length(mean_specs)-1,length(range));
    for n = 1:length(range)
        % disp(range(n))
        if strcmp(mode, 'recal')
            [mzs_recal, ~, VLM_ind, ppm_dis_ini] = msreffind('recal');
        else
            if strcmp(mode, 'workspace')
                [mzs_recal, ~, VLM_ind, ppm_dis_ini] = msreffind('workspace', ...
                    wsdatasets = datasets,wsmzs=mzs, threshold=range(n));
            elseif strcmp(mode, 'raw')
                [mzs_recal, ~, VLM_ind, ppm_dis_ini] = msreffind('raw', ...
                    threshold = range(n));
            else
                [mzs_recal, ~, VLM_ind, ppm_dis_ini] = msreffind('workspace', ...
                     wsdatasets = mean_specs, wsmzs=mzs, threshold = range(n));
            end
        end
        if isempty(VLM_ind)
            warning('migsearch:NoVLM', ...
                'No virtual lock masses at %.3f ppm; skipping.', range(n));
        
            metric(:,n) = NaN;
            nfeatures(:,n) = 0;
            continue
        end
        dum1 = mean_specs;
        dum1{1} = mean_specs{ind1};
        dum1{ind1} = mean_specs{1};
        dum2 = mzs_recal;
        dum2{1} = mzs_recal{ind1};
        dum2{ind1} = mzs_recal{1};
        mean_specs = dum1;
        mzs_recal = dum2;
        clear dum1
        clear dum2
        I(1) = ind1;
        I(ind1) = 1;

        % mz_recal = VLM(:,end);
        mz_recal = cwt2cmz(mzs_recal,'minFreq',ceil(length(mzs)*options.freq), ...
    'mzRange',[options.mzlow options.mzhigh], ...
    'ppm',range(n));
        % [~,dum_ind] = ismember(VLM(:,end),mz_recal);
        % idx_not1 = setdiff(1:numel(mz_recal), dum_ind);% just match the ones not in VLM
        % idx_not2 = setdiff(1:numel(mzs_recal{1}), VLM_ind(:,1));% just match the ones not in VLM
        %%
        [indd1,indd2] = adaptmatch(mzs_recal{1}, mz_recal,ppm_dis_ini);
        ref_spec = ones(1,length(mz_recal));
        ref_spec(:,indd2) = mean_specs{1}(:,indd1);
        ref_spec(:,VLM_ind(1,:)) = mean_specs{1}(:,VLM_ind(1,:));
        for m = 2:length(mean_specs)
            other_spec = ones(1,length(mz_recal));
            % idx_not2 = setdiff(1:numel(mzs_recal{m}), VLM_ind(:,m));% just match the ones not in VLM
            [indd1,indd2] = adaptmatch(mzs_recal{m}, mz_recal,ppm_dis_ini);
            other_spec(:,indd2) = mean_specs{m}(:,indd1);
            other_spec(:,VLM_ind(m,:)) = mean_specs{m}(:,VLM_ind(m,:));
            MImeasure = mi_cont_cont(ref_spec, other_spec, 5)/mi_cont_cont(ref_spec, ref_spec, 5);
            metric(m-1,n) = MImeasure;
            nfeatures(m-1,n) = length(mz_recal);
            % [~,Mind] = max(metric);
            % ppmthresh()
        end
            % nfeatures(n) = length(mz_recal);
    end
    % toc
    % weight by number of features
    feature_scale = abs(nfeatures-mean(lens));
    if size(metric) > 1
        MIscore = mean(metric.*(1-(feature_scale/max(feature_scale(:)))));
    else
        MIscore = (metric.*(1-(feature_scale/max(feature_scale(:)))));
    end
    
    valid = isfinite(MIscore);
    
    if ~any(valid)
        error('migsearch:NoValidThreshold', ...
            'No tested ppm threshold generated usable virtual lock masses.');
    end
    
    [~,local_ind] = max(MIscore(valid));
    valid_ind = find(valid);
    Mind = valid_ind(local_ind);
    
    ppmthresh = range(Mind);
    % [~,Mind] = min(abs(round(MIscore,1)-round(max(MIscore),1)));
    % ppmthresh = range(Mind);
    
    % figure,plot(range, feature_scale/max(feature_scale))
    % figure,plot(range, 1-(feature_scale/max(feature_scale(:))).^2)
    % figure,plot(range,metric(2:end,:).*(1-(feature_scale/max(feature_scale(:))).^2))
    % figure,plot(range,nfeatures)

    %%
    function [ind1,ind2] = adaptmatch(mz, mz_recal,ini)
        for n = 1:length(mz_recal)
            [dif(n),ind1(n)] = min(abs(mz-mz_recal(n)));
        end
        ppm_dis = ini*sqrt(mz_recal);
        margin = 1.5;
        matched = (dif<= ppm_dis*margin);
        ind1 = ind1(matched);[ind1,iu] = unique(ind1);
        ind2 = (1:n);ind2 = ind2(matched);ind2 = ind2(iu);

