function [pks] = cwt2cmz(cwt,varargin)
%cwt2cmz Convert a series of disparate peaks into an m/z vector
%   This vector serves as a reference vector against which the individual
%   spectra can be aligned. Default values are provided, but should be
%   adjusted for your dataset, especially MINFREQ.
%
%   INPUTS
%   cwt - structure of peak lists output from Yuchen's function
%   ppm - width for smoothing
%   minFreq - scalar value, minimum number of times a peak must appear
%   mzRange - [1 x 2] of low and high m/z values
%   mzRes - fractional value for m/z resolution
%   plot - true/false value for plotting figure on complet
%
%   OUTPUT
%   pks - a vector of peaks that satisfy the requirements
%
%   James McKenzie, 2025.


% Required inputs and default values
p = inputParser;
p.addParameter('ppm',20,@isnumeric); % wider accomodates larger shifts
p.addParameter('minFreq',5,@isnumeric); % a count of 5 over ppm window
p.addParameter('mzRange',[100 1000],@isnumeric);
p.addParameter('mzRes',0.001,@isnumeric);
p.addParameter('plot',false,@islogical);
p.parse(varargin{:});
% p.Results

% Standard error checking
assert(numel(p.Results.mzRange) == 2,'Specify low and high m/z values');
assert(p.Results.mzRange(1) < p.Results.mzRange(2),'Specify low m/z value first');
% assert(isfield(cwt,'list_of_peaks'),'Incorrect `cwt` input format');

% Create mz vector
mzVec = p.Results.mzRange(1):p.Results.mzRes:p.Results.mzRange(2);

% Determine frequency of peak occurence
[fq] = peakFreq(cwt,mzVec,p.Results.mzRes);

% Variable width smoothing
[sm2,sm1] = varSmooth(mzVec,fq,p.Results.mzRes,p.Results.ppm);

% Quick function for local maxima
lm = @(x) x > [x(2:end) NaN] & x > [NaN x(1:end-1)];

% These are the peaks with a frequency above the cut off
pp = sm2 >= p.Results.minFreq & lm(sm2);
% disp(['Number of peaks = ' int2str(sum(pp))]);

% Output the final vector of picked peaks
pks = mzVec(pp);

% Plot a figure?
if p.Results.plot
    [~] = plotFigure(mzVec,fq,sm1,sm2,pp,p.Results.minFreq);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fq] = peakFreq(cwt,mzVec,mzRes)
% Determine the frequency of occurence for each peak. m/z values are
% rounded using the supplied resolution, typically 0.001.

% Convert m/z values to indices
f1 = @(x) round((x-mzVec(1)) / mzRes) + 1;

% Combine all peaks
h = horzcat(cwt{:});
h = sort(h(1,:))';

% Trim out m/z values outside the range
h = h(h >= mzVec(1) & h <= mzVec(end));

% Convert m/z values to indices
hidx = f1(h);

% Sum frequently occuring peaks, and transpose final result
fq = accumarray(hidx,ones(size(hidx)),[f1(mzVec(end)) 1])';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sm2,sm1] = varSmooth(mzVec,fq,mzRes,ppm)
% Not entirely sure if this is necessary because there appears to be little
% variability in the peak vectors from the cwt output. However, the
% principal of variable width smoothing will work to coalesce peaks which
% deviate. Select the window size as appropriate.


% Convert to integer data point widths
w = round((mzVec * ppm / 1e6) / mzRes);

% Ensure an odd number
isEven = mod(w,2) == 0;
w(isEven) = w(isEven) + 1;

% Unique values
unqW = unique(w);

sm1 = zeros(size(mzVec));
sm2 = zeros(size(mzVec));

% Loop through all window sizes
for n = 1:numel(unqW)

    tic;

    % All values to use window size unqW(n)
    fx = w == unqW(n);

    % Sum over this window size
    sm1(fx) = movsum(fq(fx),unqW(n));
    
    % Create a Gauss window of same size...
    gw = gausswin(unqW(n),3);
    gw = gw / sum(gw); % normalise sum to 1
    tmp = filter(gw,1,sm1(fx));

    % Shift signal backwards
    hg = ceil(unqW(n) / 2);
    sm2(fx) = [tmp(hg:end) zeros(1,hg-1)];

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = plotFigure(mzVec,fq,sm1,sm2,pp,minFreq)

fig = figure; 
hold on;

stem(mzVec,fq,'DisplayName','Peak Frequency');

plot(mzVec,sm1,'DisplayName','Sum Smoothed'); 
plot(mzVec,sm2,'DisplayName','Gauss Smoothed');


scatter(mzVec(pp),sm2(pp),120,'r','o',...
    'DisplayName','Picked Peaks');

plot(mzVec([1 end]),[minFreq minFreq],...
    'Color','k','LineStyle','--',...
    'DisplayName','Minimum Frequency')

legend
box on;
grid on;
axis tight;
xlabel('m/z');
ylabel('Frequency');
set(gca,'FontSize',14);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
