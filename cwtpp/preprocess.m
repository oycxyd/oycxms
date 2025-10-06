%% new preprocessing workflow
%% MAIN FUNCTION
function [aligned_peaks, cmz, cwtpeaks] = preprocess(filenames,options)

arguments
    filenames string = 'default'
    options.TQswitch = 0
    options.mode string = '.raw'
    options.T0 {mustBeNumeric} = 1000;% default parameters for simulated annealing, can be optimised
    options.ite {mustBeNumeric} = 100;
    options.Nparticles {mustBeNumeric} = 1/2;
    options.Imin {mustBeNumeric} = 100;
    options.G {mustBeNumeric} = 3;
    options.L {mustBeNumeric} = 2;
    options.pmthresh = [];
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
    try
        filename = filenames(n).name;
        [~,~,ext] = fileparts(filename);
        if use_metadata == 1
            disp('using metadata!')
            meta_i = find(contains(metadata(:,2),filename));
            data_select = {};
            for p = 1:length(meta_i)
                scans = cell2mat([metadata(meta_i(p),3),metadata(meta_i(p),4)]);             
                for q = scans(1):scans(2)
                    if strcmp(ext, '.raw')
                        [mz,spectrum] = readraw2spec(filename,q);
                    else
                        options.mode = '.mz5';
                        Sindex = h5read(filename,['/SpectrumIndex']);
                        if q == 1
                            spectrum = h5read(filename,['/SpectrumIntensity'],1,double(Sindex(q)))';
                            mz = h5read(filename,['/SpectrumMZ'],1,double(Sindex(q)))';
                        else
                            try
                                spectrum = h5read(filename,['/SpectrumIntensity'],double(Sindex(q-1))+1,double(Sindex(q)-Sindex(q-1)))';
                                mz = h5read(filename,['/SpectrumMZ'],double(Sindex(q-1))+1,double(Sindex(q)-Sindex(q-1)))';
                            catch
                                spectrum = h5read(filename,['/SpectrumIntensity'],double(Sindex(q-1))+1,1)';
                                mz = h5read(filename,['/SpectrumMZ'],double(Sindex(q-1))+1,1)';
                            end
                        end
                        for i=2:length(mz)
                            mz(i) = mz(i-1)+mz(i);
                        end
                    end
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
%% peak detection w CWTPP
    try
        if strcmp(options.TQswitch,'TQ')
            disp('TQ data!')
            if contains(filename,'.raw')
                [cwtpeaks,dims] = raw2mat(filename);
                dims = [dims(1)+1,dims(2)-1];
                cwtpeaks = cwtpeaks(1:dims(1)*dims(2));
                options.mode = '.raw';
            elseif contains(filename,'.imzML')
                disp('imzml')
                [cwtpeaks,dims] = load_imzml(filename);
                options.mode = '.imzml';
            else
                [cwtpeaks,dims] = mz5toMat(filename);
                options.mode = '.mz5';
            end
            
        else
            if use_metadata == 1
                    [cwtpeaks,dims,options.mode] = cwtpp(data_select, T0 = options.T0, ite = options.ite, ...
                    Imin = options.Imin,use_metadata=1,dir=filename);
            else
                [cwtpeaks,dims,options.mode] = cwtpp(filename,T0 = options.T0, ite = options.ite, ...
                    Imin = options.Imin);
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
                [cmz,aligned_peaks] = matchSpec(cwtpeaks,options.pmthresh,freq=0.5);% ask for presence in >50% scans
            else
                [cmz,aligned_peaks] = matchSpec(cwtpeaks,options.pmthresh);
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
        if isempty(dims)
            dims = [0; 0];
        end
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
                if strcmp(options.mode,'.mz5')
                    dtwa(dname,'mz5');
                else
                    dtwa(dname);
                end
            case 'No'
                disp([answer ' OK.'])
            case 'Cancel'
                disp([answer ' Aborted.'])  
        end
    end
% toc
end