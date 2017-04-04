function [spindle_locs] = fun_ndehghani_detect_spindles_v1_Nima(Data,Fs,type)

%% Input
% spsd_data: [timepoints]x[channels]
% type : string  e.g., 'EEG','MEG'

numchan = size(Data,1);

% %% Normalize raw MEG/EEG channels ???
% switch type
%     case 'MEG'
%         Data = Data./repmat(max(abs(Data),[],2),1,size(Data,2)); % normalize each channel independently
%     case 'EEG'
%         Data = Data/max(max(abs(Data))); % normalize all channels together
% end

% Mark bad channels
bad_ch = zeros(numchan,1);
for ch=1:numchan
    bad_ch(ch) = any(isnan(Data(ch,:)));
end
bad_ch = logical(bad_ch);
chans = ~bad_ch;

%% Calculate sigma power spectrum
% "Power was calculated in the frequency range of 7–15 Hz using windows 500-ms long and
% moving steps of 100 ms"

win = 0.5*Fs;   % 500ms window
noverlap = 0.4*Fs;  % number of samples to overlap between windows (400ms overlap)
band = [7 15];
numwins = floor((size(Data,2)-win)/(win-noverlap))+1; % Power spectrum will have that many samples (timepoints)

% % % Power = zeros(numwins,numchan);
% % % for i=1:numwins
% % %     start_smp = (i-1)*(win-noverlap)+1;
% % %     if i==numwins, end_smp = size(Data,2); else end_smp = start_smp + win; end;
% % %     Power(i,chans) = bandpower(Data(chans,start_smp:end_smp)',Fs,[7 15]);
% % % end
% % % Power = Power';

[Power] = fun_Calculate_SPSD(Data,Fs,win,noverlap,band);

Power(bad_ch,:) = NaN;

size(Power)

% Normalize SPSD MEG/EEG channels
switch type
    case 'EEG'
        Power = Power/max(max(abs(Power))); % normalize all channels together
    otherwise
        Power = Power./repmat(max(abs(Power),[],2),1,size(Power,2)); % normalize each channel independently
end

[spindle_locs] = SPSD_tally_GD(Power',type);