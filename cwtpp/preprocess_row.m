function [aligned_peaks, cmz, cwtpeaks] = preprocess_row(filenames,options)

arguments
    filenames string = 'default'
    options.TQswitch = 0
    options.mode string = '.mz5'
    options.T0 {mustBeNumeric} = 1000;
    options.ite {mustBeNumeric} = 100;
    options.Nparticles {mustBeNumeric} = 1/2;
    options.Imin {mustBeNumeric} = 100;
    options.G {mustBeNumeric} = 3;
    options.L {mustBeNumeric} = 2;
    options.pmthresh = [];
end

%% load mz5 file list
if strcmp(filenames,'default')
    dname = uigetdir();
else
    dname = filenames;
end

cd(dname);

files = dir(fullfile(dname, '*.mz5'));

if isempty(files)
    error('No .mz5 files found in the folder: %s', dname);
end

filelist = fullfile(dname, {files.name});


%% peak detection using mz5toMat_row
disp('Loading mz5 image composed of multiple files...')
[cwtpeaks, dims] = mz5toMat_row(filelist);
options.mode = '.mz5';

%% peak matching
[cmz, aligned_peaks] = matchSpec(cwtpeaks, options.pmthresh);

%% TIC image
try
    TIC_image = TICimg(aligned_peaks, dims, 0);
    figure
    subplot(1,2,1);
    imagesc(TIC_image); axis image; colormap('magma');
    title('TIC image', 'Interpreter', 'none')
    subplot(1,2,2);
    stem(cmz, mean(aligned_peaks), 'Marker', 'none');
    title('Mean Spectrum', 'Interpreter', 'none')
catch
end

%% save datacube
[~, file] = fileparts(dname);
h5name = [file '.h5'];

if isfile(h5name)
    delete(h5name)
end

h5create(h5name,'/datacube',size(aligned_peaks))
h5create(h5name,'/mz',size(cmz));
h5create(h5name,'/dims',size(dims));
h5write(h5name,'/datacube',aligned_peaks)
h5write(h5name,'/mz',cmz)
h5write(h5name,'/dims',dims)

disp('All Done!')

end
