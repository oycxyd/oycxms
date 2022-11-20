%% new preprocessing workflow
function [specs_interp, global_mz, cwtpeaks] = preprocess(varargin)
clear vars
tic
%% load in raw data
    if isempty(varargin)
        dname = uigetdir('D:\BOX\Box Sync\');
        %my dir, it will just go to root if this folder doesn't exist
        cd (dname);
        filenames=dir('*.raw');

%         filenames=specs(1).name;
%         [specs_raw,~,~,xy2D] = desiReadRaw(filenames,0);
%         dims = fliplr(size(xy2D));
    else
        filenames = varargin{1};
        if isstring(filenames)
            filenames = cellstr(filenames);
        end
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
    catch
        filename = filenames;
    end
    try
        [cwtpeaks,dims] = cwtpp(filename);
        save([filename,'\cwtpeaks'],'cwtpeaks')
        %% define a global axis (vector in HS data thats ~ mean/median)
        spectral_lens = [];
        raw_specs = cwtpeaks;
        for p =1:length(raw_specs)
            spectral_lens(p) = length(cell2mat(raw_specs(p)));
        end
        mean_length=max(spectral_lens);
        [~,ind1]=min(abs(spectral_lens-mean_length));

        dum=(cell2mat(raw_specs(ind1)));
        global_mz=dum(1,:);
        clear dum

%% Interpolate all spectra to global axis
        specs_interp=interpSpec(raw_specs,global_mz,length(raw_specs) );
        specs_interp(isnan(specs_interp)|isinf(specs_interp))=0;
        % parallelised

        % check mean spec
        figure, stem(global_mz,mean(specs_interp),'Marker','none');
        title(filename)
    
%% save datacube to h5
        if isfile([filename,'/datacube.h5']) == 0
            h5create([filename,'/datacube.h5'],'/datacube',size(specs_interp))
            h5create([filename,'/datacube.h5'],'/mz',size(global_mz));
            h5create([filename,'/datacube.h5'],'/dims',size(dims));
            h5write([filename,'/datacube.h5'],'/datacube',specs_interp)
            h5write([filename,'/datacube.h5'],'/mz',global_mz)
            h5write([filename,'/datacube.h5'],'/dims',dims)
        else
            h5write([filename,'/datacube.h5'],'/datacube',specs_interp)
            h5write([filename,'/datacube.h5'],'/mz',global_mz)
            h5write([filename,'/datacube.h5'],'/dims',dims)
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
                datasets = {};mzs = {};
                for i = 1:length(filenames)
                    filename = filenames(i).name;
                    [datasets{i},dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
                end
                [data_aligned,mz_new, I] = dtwa(datasets,mzs);
                for i = 1:length(filenames)
                    disp(i)
                    filename = filenames(I(i)).name;
                    h5create([filename,'/datacube_aligned_new.h5'],'/datacube',size(data_aligned{i}))
                    h5create([filename,'/datacube_aligned_new.h5'],'/mz',size(mz_new));
                    h5create([filename,'/datacube_aligned_new.h5'],'/dims',size(dimes{I(i)}));
                    h5write([filename,'/datacube_aligned_new.h5'],'/datacube',data_aligned{i})
                    h5write([filename,'/datacube_aligned_new.h5'],'/mz',mz_new)
                    h5write([filename,'/datacube_aligned_new.h5'],'/dims',dimes{I(i)})
                end
            case 'No'
                disp([answer ' OK.'])
            case 'Cancel'
                disp([answer ' Aborted.'])  
        end
    end
toc
end