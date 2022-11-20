function [output,save_masks, check_roi]=ROIselect(mz,Data)
%% ROI selection & return spectra of pixels in ROI

% *need to define mz axis and the raw hyperspectral datacube (or 'flattened'
% version, both Nx x Ny x Ns or NxNy x Ns dims are fine)



stop=0;
k=1;
while stop == 0
    roiname = sprintf('%s_%d','roi',k);
    
    
    
    eval([roiname, '= drawassisted;']);
%     eval([roiname, '= drawrectangle;']);

    answer = questdlg('select another ROI?', ...
        'Question');
    switch answer
        case 'Yes'
            disp([answer ' OK, draw another one.'])
            k=k+1;
            continue
        case 'No'
            disp([answer ' Finished.'])
            stop = 1;
        case 'Cancel'
            disp([answer ' Aborted.'])
            break
    end
    eval(['check_roi=',roiname,';']);
end

roi_dims = check_roi.Position;
disp(['the dims of the ROI is [',num2str(round(roi_dims(3))),',',num2str(round(roi_dims(4))),']']);

output=[];
label={};
ID={};
label{1}='Label';
ID{1}='ID';


% datanames1 = inputdlg({'raw hyperspectral matrix variable name:'}, 'Specify your variables', [1 35]);
% datanames2 = inputdlg({'m/z variable name:'}, 'Specify your variables', [1 35]);
% 
% eval(['global ',char(datanames1)])
% eval(['global ',char(datanames2)])

label_names={};
for i=1:k
label_names{i} = inputdlg({['Label name ',num2str(i), ' :']}, 'Label your ROI', [1 50]); 
end

save_masks = {};
for i=1:k
    eval(['mask=createMask(roi_',num2str(i),');']);
    mask=fliplr(mask);mask=rot90(mask,1);
    save_masks{i} = mask;
    ROIspectra=Data(mask,:);
    number_pixels=size(ROIspectra);
    number_pixels=number_pixels(1);
    if i==1
        n=2;m=number_pixels+1;
        for j=n:m
            label{j}=char(label_names{i});
            ID{j}=num2str(j-1);
        end
    else
        n=m+1;m=m+number_pixels;
        for j=n:m
            label{j}=char(label_names{i});
            ID{j}=num2str(j-1);
        end
    end
    output=cat(1,output, ROIspectra);
end


label=label';
ID=ID';

output=cat(1,mz,output);
output=cat(2,label,num2cell(output));
output=cat(2,ID,output);

answer = questdlg('save to a CSV file?', ...
        'Question');
    switch answer
        case 'Yes'
            filename=inputdlg({'Save as:'}, 'Name your file', [1 50]);
            writecell(output,[char(filename),'.csv'])
            disp('OK, saved.')
%             close all
        case 'No'
            disp([answer ' Finished.'])
%             close all
        case 'Cancel'
            disp([answer ' Aborted.'])
            close all
    end


%check ROI accuracy
% input=spectra;
% tic
% [vcatest,Up,ymean,sing_values] = mvsa(input',p);
% toc
% vc1=[];
% vcimage1=[];
% for i=1:p
%     vc1(i,:)=Up(i,:);vcimage1(:,:,i)=reshape(vc1(i,:),[xsize1, ysize1]);
% %     vcimage1(:,2:2:end,i)=flipud(vcimage1(:,2:2:end,i));
%     vcimage_dum=rot90(vcimage1(:,:,i),-1);
%     vcimage_dum=fliplr(vcimage_dum);
%     figure;imagesc(vcimage_dum);colorMap = jet(256);colormap(colorMap); colorbar
% %     ax=gca;
% %     ax.PlotBoxAspectRatio=[xsize ysize 1];
% end
end
