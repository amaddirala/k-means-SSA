function [allmode] = eemd(Y,Nstd,NE)

%   FUNCTION allmode = eemd(Y,Nstd,NE)
%
% INPUT:
%       Y:      Inputted data;
%       Nstd:   ratio of the standard deviation of the added noise and that of Y;
%       NE:     Ensemble number for the EEMD : Number of trials
% OUTPUT:
%       A matrix of N*(m+1) matrix, where N is the length of the input
%       data Y, and m=fix(log2(N))-1. Column 1 is the original data, columns 2, 3, ...
%       m are the IMFs from high to low frequency, and comlumn (m+1) is the
%       residual (over all trend).
%
% NOTE:
%       It should be noted that when Nstd is set to zero and NE is set to 1, the
%       program degenerates to a EMD program.
%
% References can be found in the "Reference" section.
%
% The original code was prepared by Zhaohua Wu. For questions, please read the "Q&A" section of his website or
% contact
%   zwu@fsu.edu
%
% K. Sweeney 2013



% Variables
if size(Y,1) > size(Y,2)
    Y = Y';
end

xsize   = length(Y);
dd      = 1:xsize;
Ystd    = std(Y);
Y       = Y/Ystd;

% Dont know why this is... but determines the number of modes for later
TNM     = fix(log2(xsize))-1;           % Total number of modes
TNM2    = TNM+2;                        % Used to store the residual position
mode    = zeros(xsize, TNM2);           % Stores all the determined IMFs

% initialises the allmode matrix : where all the ensembles are stored
allmode = zeros(xsize, TNM2);

fprintf('NE : 00');
for iii = 1 : NE    % Runs through the number of Ensembles
    if iii < 10
        eval(['fprintf(''\b' num2str(iii) ''');']);
    else
        eval(['fprintf(''\b\b' num2str(iii) ''');']);
    end
    % For EEMD - Adds random white noise to the signal
    temp       = randn(1,xsize) * Nstd;    % Create the white noise to add to each ensemble.
    X1         = Y + temp;                 % Original signal plus the white noise
    mode(:,1)  = Y;                        % Sets the original signal to be first Column of output
    
    % Sets Variables for use later
    x_orig      = X1;
    xend        = x_orig;
    
    % Initialises nmode to 1 at start
    nmode       = 1;
    monotonic   = 0;
    
    %Runs through and calculates the IMF
    while nmode <= TNM && monotonic ~=1
        % Sets the start vector as the end of the last determined vector
        xstart = xend;
        %resets Iterations
        iter   = 1;
        
        while iter <= 10
            % The locations and values of the maxima minima are returned
            % flag is used to show error
            [spmax, spmin] = extrema(xstart);
            
            % Calculates the upper and lower envelope
            upper   = spline(spmax(:,1),spmax(:,2),dd);
            lower   = spline(spmin(:,1),spmin(:,2),dd);
            % Determines the new mean
            mean_ul = (upper + lower)/2;
            % The next "start" signal is the original signal minus the mean
            xstart  = xstart - mean_ul;
            iter    = iter +1;
        end
        % Sets the final signal as the current last signal minus the
        % detemined previous signal
        xend    = xend - xstart;
        
        % Increments nmode so program will finish when nmode > TNM
        nmode   = nmode+1;
        
        % Adds the latest Signal to the mode matrix
        mode(:,nmode) = xstart;
        
        % check for monotonicity
        [spmax, spmin] = extrema(xend);
        if size(spmax,1) == size(spmin,1)
            if spmax == spmin
                mode(:,nmode +1) = xend;
                monotonic        = 1;
            end
        end
    end
    
    % If this is the last run through the code. Set the residual to the
    % last Column
    if nmode > TNM && monotonic ~=1
        mode(:,TNM2) = xend;
    end
    
    % Updates the output by adding the current calculates IMF to previous
    % Ensemble IMF calculations
    allmode = allmode + mode;
    
end

blank_cols              = intersect(find(max(mode,[],1) == 0 ),find(min(mode,[],1) == 0 ));
allmode(:,blank_cols)   = [];


if NE == 100
    fprintf('\b\b\b\b\b\b\b\b');
else
    fprintf('\b\b\b\b\b\b\b');
end

% Output
allmode = allmode/NE;       % IMF's : average of all ensembles
allmode = allmode*Ystd;


