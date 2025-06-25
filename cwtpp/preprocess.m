%% new preprocessing workflow
%% MAIN FUNCTION
function [aligned_peaks, cmz, cwtpeaks] = preprocess(filenames,options)

arguments
    filenames string = 'default'
    options.TQswitch = 0
    options.mode string = '.raw'
end
% tic
%% load in raw data
    if strcmp(filenames,'default')
        dname = uigetdir();
        %my dir, it will just go to root if this folder doesn't exist
        cd (dname);
        filenames=[dir('*.raw');dir('*.imzml');dir('*.mz5')];
    else
        if strcmp(options.TQswitch,'TQ')
            dname = uigetdir();
            cd (dname);
            filenames=[dir('*.raw');dir('*.imzml');dir('*.mz5')];
        else
            filenames = cellstr(filenames);
        end
    end
%% load in metadata if available
    if exist('metadata.csv','file') == 2
        metadata = readcell('metadata.csv');
        use_metadata = 1;
        disp('metadata found in folder.')
    else
        use_metadata = 0;
    end

    logs = {};
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
            meta_i = find(contains(metadata(:,2),filename));
            data_select = {};
            for p = 1:length(meta_i)
                scans = cell2mat([metadata(meta_i(p),3),metadata(meta_i(p),4)]);             
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
        if strcmp(options.TQswitch,'TQ')
            disp('TQ data!')
            if contains(filename,'.raw')
                [cwtpeaks,dims] = raw2mat(filename);
                dims = [dims(1)+1,dims(2)-1];
                cwtpeaks = cwtpeaks(1:dims(1)*dims(2));
            end
            if contains(filename,'.imzML')
                disp('imzml')
                [cwtpeaks,dims] = load_imzml(filename);
            end
            if contains(filename,'.mz5')
                [cwtpeaks,dims] = mz5toMat(filename);
            end
            options.mode = '.raw';
        else
            if use_metadata == 1
                [cwtpeaks,dims,options.mode] = cwtpp(data_select,Imin = 100,use_metadata=1,dir=filename);
            else
                [cwtpeaks,dims,options.mode] = cwtpp(filename,Imin = 100);
            end
        end
%% define a global axis (vector in HS data thats ~ mean/median)
        % specs = cwtpeaks;
        % [~,ind1] = max(cellfun(@length, specs));
        % global_mz=specs{ind1}(1,:);

%% Interpolate all spectra to global axis
        % specs_interp=interpSpec(specs,global_mz,length(specs) );
        % specs_interp(isnan(specs_interp)|isinf(specs_interp))=0;
        % parallelised

%% NEW peak matching scheme w/o interpolation
        if strcmp(options.TQswitch,'TQ')
            cmz = cell2mat(cellfun(@(x) x(1,:), cwtpeaks(1,1), 'UniformOutput', false));
            aligned_peaks=interpSpec(cwtpeaks,cmz,length(cwtpeaks));
        else
            if use_metadata == 1
                [cmz,aligned_peaks] = matchSpec(cwtpeaks,freq=0.5);% ask for presence in >50% scans
            else
                [cmz,aligned_peaks] = matchSpec(cwtpeaks);
            end
        end
        
        
%% check mean spec & TIC image
        if use_metadata ~= 1
            try
                TIC_image = TICimg(aligned_peaks,dims,0);
                figure      
                subplot(1,2,1);
                imagesc(TIC_image);axis image;colormap('magma');
                title('TIC image', 'Interpreter', 'none')
                subplot(1,2,2);
                stem(cmz,mean(aligned_peaks),'Marker','none');
                title(filename, 'Interpreter', 'none')
            catch
                % continue
            end
        end
    
%% save datacube to h5
        if strcmp(options.mode,'.raw')
            if contains(filename,'.raw')
                if isfile([filename,'/datacube.h5']) == 0
                    h5create([filename,'/datacube.h5'],'/datacube',size(aligned_peaks))
                    h5create([filename,'/datacube.h5'],'/mz',size(cmz));
                    h5create([filename,'/datacube.h5'],'/dims',size(dims));
                    h5write([filename,'/datacube.h5'],'/datacube',aligned_peaks)
                    h5write([filename,'/datacube.h5'],'/mz',cmz)
                    h5write([filename,'/datacube.h5'],'/dims',dims)
                else
                    delete([filename,'/datacube.h5'])
                    h5create([filename,'/datacube.h5'],'/datacube',size(aligned_peaks))
                    h5create([filename,'/datacube.h5'],'/mz',size(cmz));
                    h5create([filename,'/datacube.h5'],'/dims',size(dims));
                    h5write([filename,'/datacube.h5'],'/datacube',aligned_peaks)
                    h5write([filename,'/datacube.h5'],'/mz',cmz)
                    h5write([filename,'/datacube.h5'],'/dims',dims)
                end
            else
                [~,file] = fileparts(filename);
                if isfile([file,'.h5']) == 0
                    h5create([file,'.h5'],'/datacube',size(aligned_peaks))
                    h5create([file,'.h5'],'/mz',size(cmz));
                    h5create([file,'.h5'],'/dims',size(dims));
                    h5write([file,'.h5'],'/datacube',aligned_peaks)
                    h5write([file,'.h5'],'/mz',cmz)
                    h5write([file,'.h5'],'/dims',dims)
                else
                    delete([file,'.h5'])
                    h5create([file,'.h5'],'/datacube',size(aligned_peaks))
                    h5create([file,'.h5'],'/mz',size(cmz));
                    h5create([file,'.h5'],'/dims',size(dims));
                    h5write([file,'.h5'],'/datacube',aligned_peaks)
                    h5write([file,'.h5'],'/mz',cmz)
                    h5write([file,'.h5'],'/dims',dims)
                end
            end
        else
            [~,file] = fileparts(filename);
            if isfile([file,'.h5']) == 0
                h5create([file,'.h5'],'/datacube',size(aligned_peaks))
                h5create([file,'.h5'],'/mz',size(cmz));
                h5create([file,'.h5'],'/dims',size(dims));
                h5write([file,'.h5'],'/datacube',aligned_peaks)
                h5write([file,'.h5'],'/mz',cmz)
                h5write([file,'.h5'],'/dims',dims)
            else
                delete([file,'.h5'])
                h5create([file,'.h5'],'/datacube',size(aligned_peaks))
                h5create([file,'.h5'],'/mz',size(cmz));
                h5create([file,'.h5'],'/dims',size(dims));
                h5write([file,'.h5'],'/datacube',aligned_peaks)
                h5write([file,'.h5'],'/mz',cmz)
                h5write([file,'.h5'],'/dims',dims)
            end
        end
    catch e
%         warning(['something went wrong with the file ',filename,' probably due to missing pixels. Skipping to next.'])
        fprintf(2,'There was an error! The message was:\n%s',e.message);
        logs{n} = filename;
    end
    end
    if ~isempty(logs)
        writecell(logs','error_logs.csv')
    end
    disp('All Done!')
    
    if n > 1
        answer = questdlg('Align all datasets to common m/z axis?', ...
    'Question');
        switch answer
            case 'Yes'
                disp([answer ' OK.'])
                dtwa(dname);
            case 'No'
                disp([answer ' OK.'])
            case 'Cancel'
                disp([answer ' Aborted.'])  
        end
    end
% toc
end