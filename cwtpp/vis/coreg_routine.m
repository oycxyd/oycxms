function [image_out, tforms] = coreg_routine(varargin)
image_in = varargin{1};
gt = varargin{2};
if nargin > 2
    params = varargin{3};
    param1 = params(1);param2 = params(2);
else
    param1 = 15; %radius 
    param2 = 500; %# iterations
end

% normalise first
image_in = (image_in-min(image_in(:)))/max(image_in(:));
gt = (gt-min(gt(:)))/max(gt(:));
f=figure;
subplot(1,2,1),
imshowpair(gt,image_in,'scaling','independent');
title('input')

opti_switch = 1;
tforms = {};
n = 1;
%% find unsupervised rigid transform
while opti_switch == 1
    [optimizer,metric] = imregconfig("multimodal");
    
    optimizer.InitialRadius = optimizer.InitialRadius/param1;
    optimizer.MaximumIterations = param2;
    % tformSimilarity = imregtform(image_in(:,:,3),gt(:,:,2),'rigid',optimizer,metric);
    
    % unsupervised
    tic
    tform0 = imregtform(image_in, gt,"affine",optimizer,metric);
    image_out= imwarp(image_in,tform0,"OutputView",imref2d(size(gt)));
    toc
    subplot(1,2,2),
    imshowpair(gt,image_out,'scaling','independent');
    f.WindowState = 'maximized';
    title('aligned')
    tforms{n} = tform0;
    [~,~,moc] = manders(image_in,image_out);
    
     answer = questdlg(['overlap score= ',num2str(moc),'. happy? If No, then optimisation routine will go again.'], ...
'Question');
    switch answer
        case 'Yes'
            opti_switch = 0;
            break
        case 'No'
            image_in = image_out;
            n = n+1;
            continue
        case 'Cancel'
            disp([answer ' Aborted.'])
            break
    end
end
% tform = imregtform(gt(:,:,2),image_in(:,:,1), "affine",optimizer,metric,...
% 'InitialTransformation',tformSimilarity);

%% manual rigid transform with control points %% TBD
% [mp,fp] = cpselect(image_in,gt,'WAIT',true);
% mp_adj = cpcorr(mp,fp,image_in,gt);
% tform2 = fitgeotrans(mp_adj,fp,"projective");

%% apply transform to warp original image(s) TBD
% aligned1 = imwarp(gt(:,:,1),tform,"OutputView",imref2d(size(test)));
% aligned2 = test;
% aligned3 = imwarp(gt(:,:,3),tform,"OutputView",imref2d(size(test)));
% 
% original = test;
% aligned1 = imwarp(original(:,:,1),tform,"OutputView",imref2d(size(image_out)));
% aligned2 = imwarp(original(:,:,2),tform,"OutputView",imref2d(size(image_out)));
% aligned3 = imwarp(original(:,:,3),tform,"OutputView",imref2d(size(image_out)));
% % 
% aligned = cat(3,aligned1,aligned2,aligned3);
% 
% % aligned = imwarp(image_in,tform,"OutputView",imref2d(size(image_out)));
% aligned2 = imwarp(image_in,tform2,"OutputView",imref2d(size(image_out)));
% % figure,imshowpair(aligned,image_in,'montage');
% % figure,imshowpair(gt,aligned,'montage');
% figure,imshowpair(aligned2,gt,'scaling','independent');
% 
% data_pos_new = [];
% for i = 1:length(mz_pos)
% % for i = 1:2
%     disp(i)
%     image_in = ion_image(mz_pos(i),mz_pos,data_pos,dims_pos,0,'magma');
%     aligned = imwarp(image_in,tform0,"OutputView",imref2d(size(gt)));
%     aligned = rot90(aligned,1);
%     aligned = fliplr(aligned);
%     data_pos_new(:,i) = flipud(reshape(aligned,[size(data_pos,1) 1]));
% end
% 
% data_neg_new = [];
% for i = 1:length(mz_neg)
%     image_in = ion_image(mz_neg(i),mz_neg,data_neg,dims_neg,0,'magma');
%     aligned = imresize(image_in,size(gt));
%     aligned = rot90(aligned,1);
%     aligned = fliplr(aligned);
%     data_neg_new(:,i) = flipud(reshape(aligned,[size(data_pos,1) 1]));
% end

end