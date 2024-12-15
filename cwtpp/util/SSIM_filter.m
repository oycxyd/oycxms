%% Use SSIM to identify tissue-specific peaks
%% use ion image of a reference mass (reference) to filter masses 
% that do not give images that resemble tissues on the slide. 
%% need the m/z axis of your data (mz), your data cube (data) and the dimension of your image (dims).

function [SS,filter,simmin,simmax]=SSIM_filter(reference,mz, data, dims)

    SS = [];
    
%     mz = mz_LD1;
%     data = specs_LD1_neg;
%     dims = dims_LD1;
%     reference = 174.0408;
    if length(reference) ==1
        [~,ref_index] = min(abs(mz-reference));
        ref_image =reshape(data(:,ref_index),[dims(1) dims(2)]);
        ref_image=rot90(ref_image,-1);
        ref_image=fliplr(ref_image);
        ref_image = (ref_image - min(min(ref_image)))/max(max(ref_image));
    else
        ref_image = reference;
        ref_image = (ref_image - min(min(ref_image)))/max(max(ref_image));
    end
    
%     figure;imagesc(ref_image);colorMap = jet(256);colormap(colorMap); colorbar;
%     ax=gca;
%     ax.PlotBoxAspectRatio=[dims(1) dims(2)  1];
    
    counter2 = waitbar(0,'Calculating structural similarities...');
    for i = 1:length(mz)
%         disp(i)
        waitbar(i/length(mz), counter2);
        image = ion_image(mz(i), mz, data, dims,0);
        if max(image(:))~= 0
            image_norm = (image - min(min(image)))/max(max(image));
        else
            image_norm = image;
        end
        SS(i,1) = ssim (image_norm, ref_image);
        SS(i,2) = multissim(image_norm, ref_image);
    end
    close(counter2);
    figure,boxplot(SS,'Notch','on','Labels',{'SSIM','MULTISSIM'})
    
    answer = questdlg('which metric to use?', ...
        'Question','SSIM','Multi-SSIM','Cancel');
    switch answer
        case 'SSIM'
            metric = SS(:,1);
        case 'Multi-SSIM'
            metric = SS(:,2);
        case 'Cancel'
            disp([answer ' Aborted.'])
    end
    % determine a threshold
    for i = 1:10
        simmin = min(metric);
        simmax = max(metric);
        interval = (simmax-simmin)/10;
        similarity = simmin + interval*(i-1);
        [~,index] = min(abs(metric-similarity));
        image = ion_image(mz(index), mz, data, dims,0);
        if i == 1
            fig = figure;
            fig.WindowState = 'maximized';
            subplot(2,5,i)
            imagesc(image);colorMap = jet(256);colormap(colorMap); colorbar;
            axis image
            title(['SSIM: ', num2str(metric(index))])
        else
            subplot(2,5,i)
            imagesc(image);colorMap = jet(256);colormap(colorMap); colorbar;
            axis image
            title(['SSIM: ', num2str(metric(index))])
        end
    end
    threshold = inputdlg({['Threshold ',num2str(i), ' :']}, 'choose a threshold', [1 50]); 
%     threshold = 0.7;
    filter = (metric > str2double(threshold{:}));
    close all
end
