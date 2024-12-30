function [kcadata,Imagekca,C]=cluster_analysis(input, dims, k, plot)
% start=[temp_sub.data(10,:);temp_sub.data(400,:);];


tic
if nargin < 3
    [kcadata,C,~,k]=kmeans_opt(input,10);
    plot = 1;
    disp(['the optimal number of clusters is ',num2str(k)])
else
    if nargin < 4
        plot = 1;
    end
    [kcadata,C] = kmeans(input,k,'MaxIter',100,'Distance','sqeuclidean','Replicates',4);
end
% [kcadata,C] = kmeans(input,k,'MaxIter',100,'Distance','cityblock','Replicates',4);

%             MSI
xsize = dims(1);
ysize = dims(2);

if size(input,2) > 3
    Imagekca=reshape(kcadata,[xsize ysize]);
    Image_dum = Imagekca;
    Image_dum=rot90(Image_dum,-1);
    Image_dum=fliplr(Image_dum);
    Imagekca = Image_dum;
else
    Imagekca=reshape(kcadata,[ysize xsize]);
end

if plot == 1
    figure;imagesc(0:xsize,0:ysize,Imagekca);colorMap = jet(256);colormap(colorMap); colorbar;
    title(num2str(['KCA image with k = ', num2str(k)]));
    axis image
end
% ax.PlotBoxAspectRatio=[xsize ysize  1];
% Imagekca(:,2:2:end)=flipud(Imagekca(:,2:2:end));
% Imagekca=rot90(Imagekca);
% figure;imagesc(Imagekca);colorMap = jet(256);colormap(colorMap); colorbar;
toc