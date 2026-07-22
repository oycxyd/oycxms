function ridgetable = CrazyClimber(ww, T, Nparticles, varargin)
% function rr = CrazyClimber(ww, T, Nparticles, varargin)
%
%   Performs a simulated-annealing ridge location using the algorithm
%   described in Carmona et. al (1999).  Returns a matrix of the size of
%   the field passed to it that contains the mean residence times of all
%   particles at each scale/translation pair.
%
%   Required:
%       ww - a RxRxR+ field
%       T - a cooling schedule
%       Nparticles - the number of particles to use
%
%   Optional:
%       TimingData - if set to 1, uses 'tic' and 'toc' to print timing
%           data to the MATLAB command window.  Defaults to 0.


p = inputParser;
p.addRequired('ww', @isnumeric);
p.addRequired('T', @isnumeric);
p.addRequired('Nparticles', @isnumeric);
p.addOptional('TimingData', 0, @(x)(x==1)||(x==1));
p.parse(ww,T,Nparticles,varargin{:});

opts = p.Results;

ww = opts.ww;
T = opts.T;
Nparticles = opts.Nparticles;


W = abs(ww);    
[r c] = size(W);
MC = zeros(r,c);
P.sc = floor(rand(1, Nparticles)*r)+1;
P.tr = floor(rand(1, Nparticles)*c)+1;

Niter = length(T);
if opts.TimingData
    tic
end
for lc = 1:Niter

    dx = (rand(1, Nparticles)>0.5)*2-1;
    dy = (rand(1, Nparticles)>0.5)*2-1;

    dy(P.sc == 1) = 1;
    dy(P.sc == r) = -1;
    
    proposedA = P.sc+dy;

    P.tr = mod(P.tr + dx, c) + 1;

    
    Wstart = W( ( (P.tr-1).*r) + (P.sc));
    Wend = W( ( (P.tr-1).*r) + (proposedA));

    DM = Wend-Wstart;

    switchers = (rand(1, Nparticles)<exp(DM/T(lc)))|(DM>0);
    stayers = 1-switchers;

    P.sc = sum([switchers;stayers].*[proposedA;P.sc]);
    
    MC( ((P.tr-1).*r) + (P.sc) ) =MC( ((P.tr-1).*r) + (P.sc)) + W( ((P.tr-1).*r) + (P.sc));
end
if opts.TimingData
    toc
end
MC = MC./Niter;
ridgetable = 1-MC;
