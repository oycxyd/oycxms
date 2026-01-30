function [mask] = KCAmask (data,dims)
    [~,Imagekca]= cluster_analysis((data),dims,2,0);
    Imagekca = medfilt2(Imagekca);
    se90 = strel('line',2,90);
    se0 = strel('line',2,0);
    BW2dil = imdilate(Imagekca,[se90 se0]);
    % figure,imagesc(BW2dil);axis image
    BW2dil(BW2dil==2) = 0;
    BW2dil(BW2dil==1) = 1;
    BW2 = imfill(BW2dil);
    % BW2(BW2==1) = 0;
    % BW2(BW2==2) = 1;
    mask = BW2;
    figure,imagesc(BW2);axis image
    
    answer = questdlg('keep/save this mask?', ...
'Question');
    switch answer
        case 'Yes'
            disp([answer ' OK.'])
                if isfile(['mask.h5']) == 0
                    h5create(['mask.h5'],'/mask',size(mask))
                    h5write(['mask.h5'],'/mask',mask)
                else
                    delete(['mask.h5'])
                    h5create(['mask.h5'],'/mask',size(mask))
                    h5write(['mask.h5'],'/mask',mask)
                end
        case 'No'
            disp([answer ' OK.'])
        case 'Cancel'
            disp([answer ' Aborted.'])
    end
    close all
end
    

