function output_images = RGBimage(varargin)
    if nargin < 3
        images = varargin{1};
        dims = varargin{2};
%         rgb = cat(3,(images(:,:,1)),(images(:,:,2)),(images(:,:,3)));
        rgb_ind = perms(1:3);
        for i = 1:length(rgb_ind)
            ind = rgb_ind(i,:);
            rgb = cat(3,imadjust(images(:,:,ind(1))),imadjust(images(:,:,ind(2))),imadjust(images(:,:,ind(3))));
    %         rgb = reshape (rgb, [dims(1) dims(2)]);
            rgb = (rgb-min(rgb(:)))/max(rgb(:));
            rgb = im2uint8(rgb);
            figure,imagesc(rgb);axis image
            output_images{i} = rgb;
        end
    end
end