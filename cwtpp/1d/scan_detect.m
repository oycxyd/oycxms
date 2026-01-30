dname = uigetdir();
%my dir, it will just go to root if this folder doesn't exist
cd (dname);
filenames=[dir('*.raw');dir('*.imzml');dir('*.mz5')];

labels = {};labels{1} = 'Class';labels = labels';
samples = {};samples{1} = 'Sample';samples = samples';
files = {};files{1} = 'File'; files = files';
start_scan = {};start_scan{1} = 'Start_scan'; start_scan = start_scan';
end_scan = {};end_scan{1} = 'End_scan'; end_scan = end_scan';
errors = {};errors = errors';

counter1 = parfor_wait(length(filenames), 'Waitbar', true);
disp('detecting scans from raw file(s)')

for i = 1:length(filenames)
    counter1.Send;
    disp(['file ',num2str(i),'/',num2str(length(filenames))])
    filename = filenames(i).name;
    [~,~,ext] = fileparts(filename);
    try
        if strcmp(ext, '.raw')
            disp('detected raw files')
            [raw_specs] = raw2mat(filename);
            chromo = TIC(raw_specs);
        else
            disp('detected mz5 files')
            [raw_specs,~,~,chromo] = mz5toMat(filename);
        end
        [pks,locs] = findpeaks(chromo,'MinPeakProminence',mean(chromo)/2,'MinPeakDistance',3);
        % figure,plot(chromo)
        % hold on
        % stem(locs,pks)
        peak_width = 3;
        if ~isempty(locs)
            for n = 1:length(locs)
                files = cat(1,files,filename);
                if locs(n)-peak_width<=0
                    start_scan = cat(1,start_scan,1);
                else
                    start_scan = cat(1,start_scan,locs(n)-peak_width);
                end
                if locs(n)+peak_width>length(chromo)
                    end_scan = cat(1,end_scan,length(chromo));
                else
                    end_scan = cat(1,end_scan,locs(n)+peak_width);
                end
                labels = cat(1,labels,'TBD');
                samples = cat(1,samples,i);
            end
        else
            disp(['no scans detected in ',filename,' using all scans!'])
            errors = cat(1,errors,filename);
            files = cat(1,files,filename);
            start_scan = cat(1,start_scan,1);
            end_scan = cat(1,end_scan,length(chromo));
            labels = cat(1,labels,'TBD');
            samples = cat(1,samples,i);
            continue
        end
    catch
        warning(['cannot open ',filename,'!']);
        errors = cat(1,errors,filename);
        continue
    end
end
counter1.Destroy

output = cat(2,start_scan,end_scan);
output = cat(2,files,output);
output = cat(2,labels,output);
output = cat(2,output,samples);
writecell(output,['metadata.csv'])