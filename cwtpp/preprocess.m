%% new preprocessing workflow
%% MAIN FUNCTION
function [specs_interp, global_mz, cwtpeaks] = preprocess(varargin)

tic
%% load in raw data
    if isempty(varargin)
        dname = uigetdir();
        %my dir, it will just go to root if this folder doesn't exist
        cd (dname);
        filenames=[dir('*.raw');dir('*.imzml')];
    else
        filenames = varargin{1};
        if ischar(filenames)
            filenames = cellstr(filenames);
        end
    end
%% load in metadata if available
    if exist('metadata.csv','file') > 0
        metadata = readcell('metadata.csv');
        use_metadata = 1;
        disp('metadata found in folder.')
    else
        use_metadata = 0;
    end

    for n = 1:length(filenames)
        disp(['preprocessing file ',num2str(n),' of ',num2str(length(filenames))])
%% TIC image & ROI detection
% TIC_image = TIC(datacube1);
% TIC_image = reshape(TIC_image,dims_n1);
% TIC_image = TIC_image';
% figure,imagesc(TIC_image);colormap('magma');axis image
% 
% % edge detection-based
% % [mask, threshold] = edge(TIC_image,'sobel');
% 
% % clustering-based
% [clusters,kcaimage] = cluster_analysis(3,Data_log,dims_n1);
% mask = (clusters == 1);
% mask = reshape(mask,dims_n1);
% mask = mask';
% 
% figure,imagesc(TIC_image.*mask);colormap('magma');axis image

%% peak detection w CWTPP
    try
        filename = filenames(n).name;
        if use_metadata == 1
            disp('using metadata!')
            meta_i = find(contains(metadata(:,1),filename));
            data_select = {};
            for p = 1:length(meta_i)
                scans = cell2mat([metadata(meta_i(p),4),metadata(meta_i(p),5)]);             
                for q = scans(1):scans(2)
                    [mz,spectrum] = readraw2spec(filename,q);
                    scan = cat(1,mz,spectrum);
                    data_select = cat(1,data_select,scan);
                end
            end
        end
    catch
        if iscell(filenames)
            filenames = filenames{1};
        end
        filename = filenames;
    end
    try
        if use_metadata == 1
            [cwtpeaks,dims,mode] = cwtpp(data_select,100,use_metadata);
        else
            [cwtpeaks,dims,mode] = cwtpp(filename);
        end
        %% define a global axis (vector in HS data thats ~ mean/median)
        spectral_lens = [];
        specs = cwtpeaks;
%         specs = raw_specs;
        for p =1:length(specs)
            spectral_lens(p) = length(cell2mat(specs(p)));
        end
        mean_length=max(spectral_lens);
%         mean_length=median(spectral_lens);
        [~,ind1]=min(abs(spectral_lens-mean_length));

        dum=(cell2mat(specs(ind1)));
        global_mz=dum(1,:);
        clear dum

%% Interpolate all spectra to global axis
        specs_interp=interpSpec(specs,global_mz,length(specs) );
        specs_interp(isnan(specs_interp)|isinf(specs_interp))=0;
        % parallelised

        % check mean spec & TIC image
        if use_metadata ~= 1
        TIC_image = TICimg(specs_interp,dims,0);
        figure      
            subplot(1,2,1);
            imagesc(TIC_image);axis image;colormap('magma');
            title('TIC image', 'Interpreter', 'none')
        end
        
        subplot(1,2,2);
        stem(global_mz,mean(specs_interp),'Marker','none');
        title(filename, 'Interpreter', 'none')
    
%% save datacube to h5
        if strcmp(mode,'.raw')
            if isfile([filename,'/datacube.h5']) == 0
                h5create([filename,'/datacube.h5'],'/datacube',size(specs_interp))
                h5create([filename,'/datacube.h5'],'/mz',size(global_mz));
                h5create([filename,'/datacube.h5'],'/dims',size(dims));
                h5write([filename,'/datacube.h5'],'/datacube',specs_interp)
                h5write([filename,'/datacube.h5'],'/mz',global_mz)
                h5write([filename,'/datacube.h5'],'/dims',dims)
            else
                delete([filename,'/datacube.h5'])
                h5create([filename,'/datacube.h5'],'/datacube',size(specs_interp))
                h5create([filename,'/datacube.h5'],'/mz',size(global_mz));
                h5create([filename,'/datacube.h5'],'/dims',size(dims));
                h5write([filename,'/datacube.h5'],'/datacube',specs_interp)
                h5write([filename,'/datacube.h5'],'/mz',global_mz)
                h5write([filename,'/datacube.h5'],'/dims',dims)
            end
        else
            if isfile(['datacube.h5']) == 0
                h5create(['datacube.h5'],'/datacube',size(specs_interp))
                h5create(['datacube.h5'],'/mz',size(global_mz));
                h5create(['datacube.h5'],'/dims',size(dims));
                h5write(['datacube.h5'],'/datacube',specs_interp)
                h5write(['datacube.h5'],'/mz',global_mz)
                h5write(['datacube.h5'],'/dims',dims)
            else
                delete(['datacube.h5'])
                h5create(['datacube.h5'],'/datacube',size(specs_interp))
                h5create(['datacube.h5'],'/mz',size(global_mz));
                h5create(['datacube.h5'],'/dims',size(dims));
                h5write(['datacube.h5'],'/datacube',specs_interp)
                h5write(['datacube.h5'],'/mz',global_mz)
                h5write(['datacube.h5'],'/dims',dims)
            end
        end
    catch
        warning(['something went wrong with the file ',filename,' probably due to missing pixels. Skipping to next.'])
    end
    end
    disp('All Done!')
    
    if length(filenames) > 1
        answer = questdlg('Align all datasets to common m/z axis?', ...
    'Question');
        switch answer
            case 'Yes'
                disp([answer ' OK.'])
                datasets = {};dimes = {};mzs = {};
                for i = 1:length(filenames)
                    filename = filenames(i).name;
                    [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
                end
                [data_aligned,mz_new, I] = dtwa(datasets,mzs);
                for i = 1:length(filenames)
                    disp(i)
                    filename = filenames(I(i)).name;
                    if isfile([filename,'/datacube_aligned.h5']) == 0
                        h5create([filename,'/datacube_aligned.h5'],'/datacube',size(data_aligned{i}))
                        h5create([filename,'/datacube_aligned.h5'],'/mz',size(mz_new));
                        h5create([filename,'/datacube_aligned.h5'],'/dims',size(dimes{I(i)}));
                        h5write([filename,'/datacube_aligned.h5'],'/datacube',data_aligned{i})
                        h5write([filename,'/datacube_aligned.h5'],'/mz',mz_new)
                        h5write([filename,'/datacube_aligned.h5'],'/dims',dimes{I(i)})
                    else
                        delete([filename,'/datacube_aligned.h5'])
                        h5create([filename,'/datacube_aligned.h5'],'/datacube',size(data_aligned{i}))
                        h5create([filename,'/datacube_aligned.h5'],'/mz',size(mz_new));
                        h5create([filename,'/datacube_aligned.h5'],'/dims',size(dimes{I(i)}));
                        h5write([filename,'/datacube_aligned.h5'],'/datacube',data_aligned{i})
                        h5write([filename,'/datacube_aligned.h5'],'/mz',mz_new)
                        h5write([filename,'/datacube_aligned.h5'],'/dims',dimes{I(i)})
                    end
                end
            case 'No'
                disp([answer ' OK.'])
            case 'Cancel'
                disp([answer ' Aborted.'])  
        end
    end
toc
end