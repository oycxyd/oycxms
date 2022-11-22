function [mz, spectrum,xy] = readraw2spec(filename,scan)
% tic
[p1,p2] = watersPackages(filename);

scan_size = calllib('MassLynxRaw','getScanSize',p2,1,scan);

% Read spectrum from .raw
mzp  = libpointer('singlePtr',single(zeros(1,scan_size)));    
intp = libpointer('singlePtr',single(zeros(1,scan_size)));    
calllib('MassLynxRaw','readSpectrum',p2,1,scan,mzp,intp);
xp = libpointer('singlePtr',0);   
yp = libpointer('singlePtr',0);    
calllib('MassLynxRaw','getXYCoordinates',p1,1,scan,xp,yp); 

mz = double(mzp.Value); 
spectrum = double(intp.Value);
xy = [xp.Value,yp.Value];
% toc

% need to free up allocated space to avoid memory leak!!!
calllib('MassLynxRaw','delCMassLynxRawReader',p1);
calllib('MassLynxRaw','delCMassLynxRawScanReader',p2);

end