function [SPSD] = fun_Calculate_SPSD(Data,Fs,win,noverlap,band)

%% Inputs
% Data: data array (n channels x n timepoints)
% Fs: Sampling Frequency (Hz)
% win: Length of window (samples)
% noverlap: window overlap (samples)
% band: [lowfreq highfreq]

%% Outputs
% SPSD: Power Spectral Density within each window (n channels x n timepoints)

%% Calculate sigma power spectrum
% "Power was calculated in the frequency range of 7–15 Hz using windows 500-ms long and
% moving steps of 100 ms"


numchan = size(Data,1);
numwins = floor((size(Data,2)-win)/(win-noverlap))+1; % Power spectrum will have that many samples (timepoints)
SPSD = zeros(numwins,numchan-2);

notNaN = ~all(isnan(Data),2);
for i=1:numwins
    start_smp = (i-1)*(win-noverlap)+1;
    if i==numwins, end_smp = size(Data,2); else end_smp = start_smp + win; end;
    SPSD(i,notNaN) = bandpower(Data(notNaN,start_smp:end_smp)',Fs,band);
end
SPSD(:,~notNaN) = NaN;
SPSD = SPSD';

