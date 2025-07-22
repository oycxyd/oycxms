function batch_recal(varargin)
    dname = uigetdir();
    %my dir, it will just go to root if this folder doesn't exist
    cd (dname);
    filenames=[dir('*.raw');dir('*.mz5')];


    mzs_recal = {};
    ppms = [];
    if isempty(varargin)
        threshold = 300;%ppm
        references = 255.2330;
    else
        threshold = varargin{1};
        references = varargin{2};
    end
    

    for i = 1:length(filenames)
        disp(['recalibrating ', num2str(i),'/',num2str(length(filenames))])
        filename = filenames(i).name;
        [~,file,ext] = fileparts(filename);
        if strcmp(ext,'.mz5')
            [~,dimes{i},mzs{i}]=h5toMat([file,'.h5']);
            [mzs_recal{i},ppm] = MSrecal(mzs{i},references, threshold);
            copyfile([file,'.h5'],[file,'_recal.h5'])
            h5write([file,'_recal.h5'],'/mz',mzs_recal{i})
        else
            [~,dimes{i},mzs{i}]=h5toMat([filename,'\datacube.h5']);
            [mzs_recal{i},ppm] = MSrecal(mzs{i},references, threshold);
            copyfile([filename,'/datacube.h5'],[filename,'/datacube_recal.h5'])
            h5write([filename,'/datacube_recal.h5'],'/mz',mzs_recal{i})
        end
        ppms(i) = mean(ppm);
    end
    disp (['the maximum average ppm after recalibration is ',num2str(max(ppms))])
end
