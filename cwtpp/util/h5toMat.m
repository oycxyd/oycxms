function [varargout]=h5toMat(varargin)
if ischar(varargin{1}) == 1
    filename = varargin{1};
else
    dname = uigetdir('D:\BOX\Box Sync\');%office, it will just go to root if this folder doesn't exist
    cd (dname);
    specs=[dir('*.h5');dir('*.hdf5')];
    if length(specs) > 1
        k = varargin{1};
    else
        k = 1;
    end
    filename=specs(k).name;
end
%% load h5 files in a directory %%
% works for any h5 of any size, define varibales accordingly
% e.g. [mz,Data,size]=h5toMat; 

try
    info=h5info(filename);
    datasets=info.Datasets;
    if isempty(datasets)==1
        %just for the old GUI, which saves data into Groups
         warning('Data not found in Datasets, will try Groups');
         datasets=info.Groups(2).Datasets;
         number_variables=length(datasets);
         for i=1:number_variables
            dataname=[datasets(i).Name];
            temp=h5read(filename,['/Brillouin data/',dataname]);
            varargout{i} = double(temp);
         end
    else
        number_variables=length(datasets);
        for i=1:number_variables
            dataname=[datasets(i).Name];
            temp=h5read(filename,['/',dataname]);
            varargout{i} = double(temp);
        end
%     eval([dataname,'=double(','data',num2str(i),');'])
    end
catch
   warning('Error! No/multiple files found.');
end
end