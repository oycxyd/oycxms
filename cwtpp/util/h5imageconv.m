function [peak_I] = h5imageconv(varargin)

% spatial filter for noisy images (using TIC image generally)
% TIC_image = TICimg(data,dims,0);

% [SS,filter] = SSIM_filter(ion1,mz,data,dims);

if isempty(varargin)
    [file,location] = uigetfile('*.h5');
    cd (location);
else
    file = varargin{1};
end

[data,dims,mz] = h5toMat(file);
counter1 = waitbar(0,'converting to images...');
peak_I = [];
if exist('output_images','dir')==0
    mkdir('output_images')
end
for i = 1:length(mz)
% for i = 1:top
    waitbar(i/size(data,2), counter1);
% %         waitbar(i/length(specs_LR), counter1);
% %         filename = ['liver1peak',num2str(i)];
    Img = ion_image(mz(i),mz,data,dims,0,'magma');
% %     Img = Img.*mask;
% %         image = specs_HR(i,:,:);
% %         image = reshape(image,[82 124]);
    peak_I(i) = max(Img(:));
%     I = image-min(image(:));
    Img = Img/max(Img(:));
    rgb = cat(3,Img,Img,Img);
    I_name = num2str(mz(i));
    %     I=im2uint16(image);
%         imwrite(I,['C:\Users\OYCX\Desktop\temp\training data for SSR\liver\test_tiff_new\HR\',filename,'.tif'])
    imwrite(rgb,['output_images\',I_name,'.tif'])
end
% save(['output_images\','Ipeak'],'peak_I')
save(['Ipeak'],'peak_I')
% clear data;clear dims;clear mz
close(counter1);

%% read restored images back in -> h5
% peak_I = cell2mat(readcell('brain4_maxs.csv'));
% 
% % files=dir('*.png');
% files=dir('*.tif');
% filenames = {files.name};
% filenames = natsortfiles(filenames);
% 
% counter2 = parfor_wait(length(filenames), 'Waitbar', true);
% specs_res = [];
% parfor i = 1:length(filenames)
%     counter2.Send;
% %     disp(i)
%     filename = filenames{i};
% %     [~,mz_name] = (fileparts(filename));
% %     mz_name = str2double(mz_name);
% %     mz(i) = mz_name;
% %     [~, mz_ind] = min( abs(mz-mz_name) );
% 
%     image = im2double(imread(filename));
%     image = sum(image,3)/3;
% %     image = imresize(image,size(Img));
% %     image(2:2:end,:)=fliplr(image(2:2:end,:));
% 
%     image = rot90(image,-1);
%     image = fliplr(image);
% 
% 
%     specs_res(:,i) = (image(:)-min(image(:)))*peak_I(i);
% %     specs_res(:,i) = image(:); 
% end
% counter2.Destroy
% % new_dims = (dims)*4;
% 
% % scores = output{2};
% % rgb = cat(3,scores(1,:),scores(2,:),scores(3,:));
% % rgb = reshape (rgb, [fliplr(dims1)*4,3]);
% % figure,imagesc((((rgb))));axis image
% 
% h5create('datacube_HyReS.h5','/datacube',size(specs_res))
% h5create('datacube_HyReS.h5','/mz',size(mz))
% h5create('datacube_HyReS.h5','/dims',size(new_dims))
% h5write('datacube_HyReS.h5','/datacube',specs_res)
% h5write('datacube_HyReS.h5','/mz',mz)
% h5write('datacube_HyReS.h5','/dims',new_dims)
