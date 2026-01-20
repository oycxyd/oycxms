%% using dtw to align (input is a cell of datasets to be aligned with their corresponding m/z vectors)
% matlab dtw function distance metric has a bug, not using atm while
% investigation ongoing
function [data_aligned, mz_recal, I] = dtwa(path,mode,options)
    arguments
        path (1,:) = pwd
        mode (1,:) char = 'raw'
        options.matchppm (1,:) {mustBeNumeric,mustBeReal} = 100
        options.save = 1
    end
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
            varargin{3} = 'workspace';
        end
        filenames=dir('*.raw');
        datasets = {};mzs = {};
        
        if strcmp(mode, 'msp')
            disp('using resampled data.')
            for i = 1:length(filenames)
                filename = filenames(i).name;
                [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube_msp.h5']);
%             [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
%             datasets = {};
            end
        elseif strcmp(mode, 'recal')
            disp('using recalibrated data.')
            options.matchppm = 50;
            filenames=[dir('*.raw');dir('*_recal.h5')];
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
        elseif strcmp(mode, 'workspace')
            disp('using defined workspace data.')
            datasets = varargin{1};
            mzs = varargin{2};
            options.save = 0;
        elseif strcmp(mode, 'raw')
            for i = 1:length(filenames)
            filename = filenames(i).name;
            [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
%             datasets = {};
            end
        else
            disp('detected mz5 files.')
            filenames=[dir('*.mz5')];
            for i = 1:length(filenames)
                filename = filenames(i).name;filename = filename(1:end-4);
                % [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
                [~,dimes{i},mzs{i}]=h5toMat([filename,'.h5']);
    %             datasets = {};
            end
        end
    end
    
    % if length(varargin)<3
    %     method = 'mean_spec';
    % else
    %     method = varargin{3};
    % end
    answer = questdlg('which method to use?', ...
    'Question','meanspec','refmass','cancel');
    switch answer
        case 'meanspec'
            method = 'meanspec';
            disp([answer ' OK.'])
        case 'refmass'
            method = 'mz';
            disp([answer ' OK.'])
        case 'cancel'
            disp([answer ' Aborted.'])
    end
    
    for i = 1:length(mzs)
        lens(i) = length(mzs{i});
    end
    
    if strcmp(method,'mz')
        disp('aligning using mz')
        mean_length=median(lens);
        % mean_length=max(lens);
        [~,ind1]=min(abs(lens-mean_length));
        I = 1:length(mzs);
        % dum1 = datasets;
        if strcmp(mode, 'recal')
            [mzs_recal, ms_references] = msreffind('recal');
        else
            if strcmp(mode, 'workspace')
                [mzs_recal, ms_references] = msreffind('workspace', ...
                    wsdatasets = datasets,wsmzs=mzs);
            elseif strcmp(mode, 'raw')
                [mzs_recal, ms_references] = msreffind('raw', ...
                    threshold = options.matchppm);
            else
                [mzs_recal, ms_references] = msreffind('mz5', ...
                    threshold = options.matchppm);
            end
        end
        % [mzs_recal, ms_references] = msreffind();
        mzs = mzs_recal;
        dum2 = mzs;
        % dum1{1} = datasets{ind1};
        dum2{1} = mzs{ind1};
        % dum1{ind1} = datasets{1};
        dum2{ind1} = mzs{1};
        % datasets = dum1;
        mzs = dum2;
        I(1) = ind1;
        I(ind1) = 1;
        % clear dum1
        clear dum2
        disp('using the following reference masses: ')
        disp([num2str(ms_references(:,1))])
    else
        [~,I] = sort(lens,'descend');
%         datasets = datasets(I);
        mzs = mzs(I);
    end

%% inter-data recalibration step (under development)
%     threshold = 200;
% %     mzs = {}; mzs{1}=mz1;mzs{2}=mz2;mzs{3}=mz3;mzs{4}=mz4;mzs{5}=mz51;
% %     ms_references = zeros(length(mz_recal52),length(mzs)+1);
%     ms_references = [];
%     ms_references(:,1) = mz_recal;
%     mzs_recal = {};
%     for n = 1:length(mzs)
% %     for n = 1:50
%         disp(n)
%         mz_raw = mzs{n};
%         [~,locs] = findpeaks(mean(datasets{n}),'MinPeakProminence',mean(mean(datasets{n})));
%         mz_raw_p = mz_raw(locs);
%         [mz_new_p] = MSrecal(mz_raw_p,references, 1000);
% %         for p = 1:3+1
% %             if p == 1
% %                 new_mz = coefs(end);
% %             else
% %                 new_mz = new_mz + coefs(end-p+1)*mz_raw.^(p-1);
% %             end
% %         end
%         mz_raw(locs) = mz_new_p;
%         [mz_new] = MSrecal(mz_raw,references, 1000);
%         mzs_recal{n} = mz_new;
% %         for m = 1:length(ms_references)
% %                 [diff, ind] = min( abs(mz_new-ms_references(m,1)) );
% %     %             [diff, ind] = min( abs(mz_raw-ms_references(m,1)) );
% %                 ppm = diff/ms_references(m,1)*10^6;
% %                 if ppm <= threshold
% %                     ms_references(m,n+1) = ind;
% %                 else
% %                     ms_references(m,:)=0;
% %                 end
% %         end     
%     end
%     ms_references(ms_references(:,1)==0,:)=[];
%     
%     for i = 1:length(mzs_recal)
% %     for i = 1:2
%         dum = mzs_recal{i};
%         dum(ms_references(:,i+1)) = ms_references(:,1);
%         mzs_recal{i} = dum;
%         clear dum
%     end
%% DTW-based alignment
    % mz_recal = mzs{1};
    % mz_low = min(unique(round(cell2mat(mzs_recal),4)));
    mz_low = 50;
    % mz_high = max(unique(round(cell2mat(mzs_recal),4)));
    mz_high = 1200;
    disp(['using a threshold of ',num2str(options.matchppm), ' ppm'])
    mz_recal = cwt2cmz(mzs,'minFreq',ceil(length(mzs)*0.4), ...
        'mzRange',[mz_low mz_high], ...
        'ppm',options.matchppm);
    % mz_recal = unique(round(cell2mat(mzs_recal),4));
    % test = diff(mz_recal);test = cat(2,1,test);indc = (test>0.01);
    % mz_recal = mz_recal(indc);
    if strcmp(mode,'workspace')
        % [data1] = datasets{1};
        if strcmp(method,'mz')
            [a,b] = maxent(mzs{1},datasets{1});
            [data1] = a*randn(size(datasets{1},1),length(mz_recal))+b;
            % [~,indd1,indd2] = intersect(round(mzs{1},4),round(mz_recal,4));
            [indd1,indd2] = adaptmatch(mzs{1}, mz_recal);
            if strcmp(mode,'recal')
                % [~,indd1,indd2] = intersect(round(mzs{1},1),round(mz_recal,1));
                [indd1,indd2] = fixedmatch(mzs{1}, mz_recal,10^-(log10(median(diff(mz_recal))/min(diff(mz_recal)))/2));
            else
                [indd1,indd2] = adaptmatch(mzs{1}, mz_recal);
            end
            dum = datasets{1};data1(:,indd2) = dum(:,indd1);
            % [mz_recal_f,data1] = merge_peaks(mz_recal,data1,0.001);
            % [data1]=mzmatchfill(mzs{1},datasets{1},mz_recal,options.matchppm);
        end
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
        else
            [data1]=h5toMat([filenames(I(1)).name,'\datacube.h5']);
            out_dir = filenames(I(1)).name;
        end
        if strcmp(method,'mz')
            [a,b] = maxent(mzs{1},data1);
            [data11] = a*randn(size(data1,1),length(mz_recal))+b;
            % [~,indd1,indd2] = intersect(round(mzs{1},1),round(mz_recal,1));
            [indd1,indd2] = adaptmatch(mzs{1}, mz_recal);
            if strcmp(mode,'recal')
                % [~,indd1,indd2] = intersect(round(mzs{1},1),round(mz_recal,1));
                [indd1,indd2] = fixedmatch(mzs{1}, mz_recal,10^-(log10(median(diff(mz_recal))/min(diff(mz_recal)))/2));
            else
                [indd1,indd2] = adaptmatch(mzs{1}, mz_recal);
            end
            data11(:,indd2) = data1(:,indd1);
            % [mz_recal_f,data11] = merge_peaks(mz_recal,data11,0.001);
            data1 = data11;
            % [data1]=mzmatchfill(mzs{1},data1,mz_recal,options.matchppm);
        end
    end
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
    end
    
%     mz_recal(ms_references(:,n+1))=ms_references(:,1);

    % data_length = zeros(length(mzs),1);
    % mz2_length = zeros(length(mzs),1);
    for n = 1:length(mzs)-1
        disp(['aligining file ',num2str(n),'/',num2str(length(mzs)-1)])
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
            else
                filename = filenames(I(n+1)).name;
                data2 = h5toMat([filename,'\datacube.h5']);
                out_dir = filename;
            end
        end
        % data_length(n) = size(data2,2);
        mz2 = mzs{n+1};
        % mz2_length(n) = length(mz2);
        if strcmp(method,'mean_spec')
            disp('aligning using mean spectra')
            if size(data1,1) > 1
                [~,i1,i2] = dtw(mean(data1), mean(data2),'symmkl');
            else
                [~,i1,i2] = dtw((data1), (data2),'symmkl');
            end
        else
            disp('aligning using mz')
            % [~,i1,i2] = dtw(mz_recal, mz2,'symmkl');
        end
        data_dtw = [];
        if strcmp(method,'mz')
            [a,b] = maxent(mz2,data2);
            data_dtw = a*randn(size(data2,1),length(mz_recal))+b;
            % [~,indd1,indd2] = intersect(round(mz2,4),round(mz_recal,4));
            [indd1,indd2] = adaptmatch(mz2, mz_recal);
            if strcmp(mode,'recal')
                % [~,indd1,indd2] = intersect(round(mz2,1),round(mz_recal,1));
                [indd1,indd2] = fixedmatch(mz2, mz_recal, 10^-(log10(median(diff(mz_recal))/min(diff(mz_recal)))/2));
            else
                [indd1,indd2] = adaptmatch(mz2, mz_recal);
            end
            data_dtw(:,indd2) = data2(:,indd1);
            % [~,data_dtw] = merge_peaks(mz_recal,data_dtw,0.001);
        else
            counter1 = parfor_wait(size(data2,1), 'Waitbar', true);
            parfor j = 1:size(data2,1)
                counter1.Send;
                spec_dtw = async(data2(j,:),i2, i1);
                data_dtw(j,:) = spec_dtw;
            end
            counter1.Destroy
        end
        
        if options.save == 1        
%                 disp(n)
            disp('saving...');
            if strcmp(mode,'mz5')
                if isfile([filename,'_aligned.h5']) > 0
                    disp('Found existing aligned file. Saving over.')
                    delete([filename,'_aligned.h5'])
                end
                h5create([filename,'_aligned.h5'],'/datacube',size(data_dtw))
                h5create([filename,'_aligned.h5'],'/mz',size(mz_recal));
                h5create([filename,'_aligned.h5'],'/dims',size(dimes{I(n+1)}));
                h5write([filename,'_aligned.h5'],'/datacube',data_dtw)
                h5write([filename,'_aligned.h5'],'/mz',mz_recal)
                h5write([filename,'_aligned.h5'],'/dims',dimes{I(n+1)})
            elseif strcmp(mode,'recal')
                if contains(filename,'.h5')
                    filename = filename(1:end-3);
                    if isfile([filename,'_aligned.h5']) > 0
                        disp('Found existing aligned file. Saving over.')
                        delete([filename,'_aligned.h5'])
                    end
                    h5create([filename,'_aligned.h5'],'/datacube',size(data_dtw))
                    h5create([filename,'_aligned.h5'],'/mz',size(mz_recal));
                    h5create([filename,'_aligned.h5'],'/dims',size(dimes{I(n+1)}));
                    h5write([filename,'_aligned.h5'],'/datacube',data_dtw)
                    h5write([filename,'_aligned.h5'],'/mz',mz_recal)
                    h5write([filename,'_aligned.h5'],'/dims',dimes{I(n+1)})
                else
                    if isfile([out_dir,'/datacube_aligned.h5']) > 0
                        disp('Found existing aligned file. Saving over.')
                        delete([out_dir,'/datacube_aligned.h5'])
                    end
                    h5create([out_dir,'/datacube_aligned.h5'],'/datacube',size(data_dtw))
                    h5create([out_dir,'/datacube_aligned.h5'],'/mz',size(mz_recal));
                    h5create([out_dir,'/datacube_aligned.h5'],'/dims',size(dimes{I(n+1)}));
                    h5write([out_dir,'/datacube_aligned.h5'],'/datacube',data_dtw)
                    h5write([out_dir,'/datacube_aligned.h5'],'/mz',mz_recal)
                    h5write([out_dir,'/datacube_aligned.h5'],'/dims',dimes{I(n+1)})
                end
            else
                if isfile([out_dir,'/datacube_aligned.h5']) > 0
                    delete([out_dir,'/datacube_aligned.h5'])
                end
                h5create([out_dir,'/datacube_aligned.h5'],'/datacube',size(data_dtw))
                h5create([out_dir,'/datacube_aligned.h5'],'/mz',size(mz_recal));
                h5create([out_dir,'/datacube_aligned.h5'],'/dims',size(dimes{I(n+1)}));
                h5write([out_dir,'/datacube_aligned.h5'],'/datacube',data_dtw)
                h5write([out_dir,'/datacube_aligned.h5'],'/mz',mz_recal)
                h5write([out_dir,'/datacube_aligned.h5'],'/dims',dimes{I(n+1)})
            end
        else
            data_aligned{n+1} = data_dtw;
            disp('workspace check')
        end
        clear data_dtw
    end
    
%% visualise results
%     answer = questdlg('Visualise results of alignment (using PCA)?', ...
%     'Question');
%         switch answer
%             case 'Yes'
%                 disp([answer ' OK.'])
%                 specs_aligned = [];
%                 labels = {};
%                 for m = 1:length(data_aligned)
%                     data_m = data_aligned{m};
%                     if length(data_aligned)> 5
%                         disp('a lot of data! taking 10% pixels per file only.')
%                         percentage = 10;
%                         sample_size = ceil(size(data_m,1)*percentage/100);
%                         idx = randi(sample_size,[1 sample_size]);
%                         data_m = data_m(idx,:);
%                         specs_aligned = cat(1,specs_aligned,data_m);
%                     else
%                         specs_aligned = cat(1,specs_aligned,data_aligned{m});
%                     end
%                     dum = {};
%                     for i = 1:size(data_m,1)
%                         dum{i}=['data',num2str(m)];
%                     end
%                     dum = dum';
%                     labels = cat(1,labels,dum);
%                     clear dum
%                 end
% 
%                 [output] = component_analysis(log_trans(TIC_norm(specs_aligned)),3,[0 0]);
%                 scores = output{2};
%                 figure,gscatter(scores(1,:),scores(2,:),labels)   
%             case 'No'
%                 disp([answer ' OK.'])
%             case 'Cancel'
%                 disp([answer ' Aborted.'])  
%         end
%     data_aligned = data_aligned(I);
    % [~, ~,~,coefs]=peakfit([mz2 mz_recal],0,0,1,28,n,1,0,0,0,1);
    % for p = 1:n+1
    %     if p == 1
    %         new_mz = coefs(end);
    %     else
    %         new_mz = new_mz + coefs(end-p+1)*mz_recal.^(p-1);
    %     end
    % end
    %%
    % [indd1,indd2] = adaptmatch(mzs{1}, mz_recal);
    % [a,b] = maxent(mzs{1},data1);

    %%
function [a,b,ent] = maxent(mz,data)
    ent = [];
    for i = 1:length(mz)
        ent(i) = entropy(data(:,i));
    end
    [~,ii] =  max(ent);
    a = sqrt(std(data(:,ii)))/100;b = mean(data(:,ii))/100;

function [ind1,ind2] = adaptmatch(mz, mz_recal)
    for n = 1:length(mz_recal)
        [dif(n),ind1(n)] = min(abs(mz-mz_recal(n)));
    end
    ini = mean(dif)/10;rmse = sqrt(abs(sum(dif-ini*sqrt(mz_recal))))/length(dif);
    sign = 1;
    while true
        step = ini/10;
        ini = ini+step*sign;
        rmse_new = sqrt(abs(sum(dif-ini*sqrt(mz_recal))))/length(dif);
        if rmse_new > rmse
            sign = sign*-1;
        end
        if rmse_new-rmse < rmse/100
            break
        end
        rmse = rmse_new;
    end
    ppm_dis = ini*sqrt(mz_recal);
    margin = 1.1;
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
% function [data_aligned]=mzmatchfill(mz_in,data_in,mz_recal,thresh)
%     [a,b] = maxent(mz_in,data_in);
%     [data_aligned] = a*randn(size(data_in,1),length(mz_recal))+b;    
%     for n = 1:length(mz_recal)
%         % [inds]=(abs(mz_in-mz_recal(n)))/mz_recal(n)*1e6 <= thresh;
%         [mdiff,mind]=min(abs(mz_in-mz_recal(n)));
%         if mdiff/mz_recal(n)*1e6 <= thresh
%         % if isempty(find(inds,1)) == 0
%             % data_aligned(:,n) = sum(data_in(2,inds));
%             data_aligned(:,n) = (data_in(2,mind));
%         else
%             disp(['no match at ',num2str(mz_recal(n))])
%         end
%     end
