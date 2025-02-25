function [network,pairs,results] = coloc2(mz,data,dims,mask,save)

x = 1:length(mz);
[X,Y] = meshgrid(x);
Xu = triu(X,1);
Yu = triu(Y,1);
pairs = [nonzeros(Xu(:)),nonzeros(Yu(:))];

if mask == 1
    [file,path] = uigetfile('*.h5');
    h5mask=h5read(fullfile(path,file),['/mask']);
    % h5mask = h5mask-1;
    % h5mask = logical(h5mask);
    figure,imagesc(h5mask);axis image
    data_r = data(logical(h5mask),:);
end
% new_mask = strcmp (h5mask,'TRUE');

tic
results = zeros(length(pairs),5);
counter1 = parfor_wait(length(pairs), 'Waitbar', true); 
parfor i = 1:ceil(length(pairs))
    counter1.Send;
% parfor i = 1:10    
%     disp(i)
%     waitbar(i/length(pairs), counter1);
%     waitbar(i/100, counter1);
    if mask
        image1 = data_r(:,pairs(i,1));
        image2 = data_r(:,pairs(i,2));
    else
        image1 = data(:,pairs(i,1));
        image2 = data(:,pairs(i,2));
    end
    image1 = image1/sum(image1(:));  
    image2 = image2/sum(image2(:));
    [rho,pval] = corr(image1,image2, 'Type','Spearman');
    [M1, M2, MOC]=manders(image1, image2);
    results(i,:) = [rho, pval, M1, M2, MOC];
end
counter1.Destroy
toc

figure,boxplot([results(:,1) results(:,5)],'Notch','on','Labels',{'spearman','MOC'})
figure,boxplot([results(:,2)],'Notch','on','Labels',{'pvals'})
[~, ~, ~, adj_p]=fdr_bh(results(:,2));%p-value correction
results(:,2) = adj_p;

condition = (abs(results(:,1))>0.9) & (results(:,2)<0.05) & (results(:,5)>0.9);
network = pairs(condition,:);
results = results(condition,:);
% network = pairs;

% save coloc images & network plot
if save
    labels = {};
    for i = 1:length(network)
        if i==1
            if not(isfolder('network'))
                disp('saving.')
                mkdir('network')
                cd ('network');
            else
                cd ('network');
            end
        end
    % for i = 1:length(network)
        pair = network(i,:);
        red = ion_image(mz(pair(1)),mz,data,dims,0);
        green = ion_image(mz(pair(2)),mz,data,dims,0);
        if mask
            red = red.*h5mask;
            green = green.*h5mask;
        end
        % rgb = cat(3,red/sum(red(:)),green/sum(green(:)),red*0);
        rgb = imfuse(imadjust(red/sum(red(:))),imadjust(green/sum(green(:))),'falsecolor','Scaling','independent','ColorChannels',[1 2 0]);
        % hsv = rgb2hsv(rgb);
       
        gfilter = imgaussfilt(rgb,0.6);
        % for j = 1:3
        % %     mfilter(:,:,i) = medfilt2(gfilter(:,:,i),[2 2]);
        %     mfilter(:,:,j) = medfilt2(rgb(:,:,j),[3 3]);
        % end
        % figure,imagesc(rgb);axis image
        figure('visible','off'),imagesc(gfilter);axis image
        pair_name = [num2str(mz(pair(1))),'-',num2str(mz(pair(2))),', corr=',num2str(results(i,1))];
        labels{i} = pair_name;
        title(pair_name)
        saveas(gcf,[pair_name,'.png'])
        close all
    end
end

end

%visualisation of networks_under dev
% s={};
% t={};
% for i = 1:length(network)
%     s{i} = num2str(mz(network(i,1)));
%     t{i} = num2str(mz(network(i,2)));
% end
% G = digraph(s,t);
% G = digraph(s,t,[],G.Nodes,'omitselfloops');
% figure,plot(G,'Layout','force','NodeLabel',G.Nodes.Name,'NodeColor',...
%     'k','EdgeAlpha',0.5,'UseGravity',true)
% figure,plot(G,'Layout','circle','NodeLabel',G.Nodes.Name,...
%     'NodeColor','k','EdgeAlpha',1.0,'EdgeColor',[0.4940 0.1840 0.5560])
% title('colocalisation networks')
% fig=gcf;
% set(gcf, 'PaperPosition', [0 0 20 20])    % can be bigger than screen 
% set(gcf, 'PaperSize', [20 20]) 
% print(gcf, 'networks.png','-dpng');
