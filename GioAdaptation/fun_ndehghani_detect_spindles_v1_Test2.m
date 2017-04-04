function [spindle_locs] = fun_ndehghani_detect_spindles_v1_Test2(Data,Fs,type)
% Data: EEG or MEG data from a 30sec epoch, (channels) x (timepoints)
% Fs: data sampling frequency
% type: specify 'EEG' for EEG data or 'MEG' for MEG data (differs only in
% channel normalization method)

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

%% Normalize SPSD MEG/EEG channels
% Implemented as in Figure 2

switch type
    case 'MEG'
        Power = Power./repmat(max(abs(Power),[],2),1,size(Power,2)); % normalize each channel independently
    case 'EEG'
        Power = Power/max(max(abs(Power))); % normalize all channels together
end

%% Smooth with 5 passes of 3-point moving average filter
% "each sensor’s SPSD data were smoothed by five sequential passes
% of a three point smoothing filter."

mov_ave = 1/3*[1 1 1];
num_passes = 5;

for ch=1:numchan  
        for i=1:num_passes
            Power(ch,:) = conv(Power(ch,:),mov_ave,'same'); % keep central part so as to avoid lag
        end
end

%% Detect local maxima and std across channels
% "All resulting local maxima were found, and the standard
% deviation of these maxima was calculated across all sensors"

sigma_max = diff(Power,[],2);

sigma_max_locs = []; % keep points where some channel (at least one) has maximum
sigma_max_locs_ch = cell(numchan,1);
Power_max = [];
flag = 0;
for i=2:size(sigma_max,2)
    for ch=1:numchan
        if and(sigma_max(ch,i)<0,sigma_max(ch,i-1)>0) % zero-crossings (points of local maxima)
            sigma_max_locs_ch{ch} = [sigma_max_locs_ch{ch} i];
            Power_max = [Power_max Power(ch,i)]; % stack maxima for entire epoch across all channels
            flag = 1;
        end
    end
    if flag, sigma_max_locs = [sigma_max_locs i]; end
end
sigma_max_locs = unique(sigma_max_locs);

sigma_max_std = std(Power_max); % calculate standard deviation across all maxima


%% Calculate channel count and retain channels >2*std
% "The number of channels with maxima greater than two standard deviations
%  above zero occurring within a 100-ms window was counted." 

max_channels = cell(numwins,1);
chan_count = zeros(numwins,1); % channel count for every frame (100ms-point in the power spectrum)

for t=1:numwins % for every 100ms-point
    
    if ~ismember(t,sigma_max_locs), continue; end % if not in the location of a maximum leave zero channel count and go to next 100ms-point
    
    for ch=1:numchan % look at every channel: does it have a maximum in this 100ms-point?
        ch_max = ismember(t,sigma_max_locs_ch{ch});
        
        if ~ch_max, continue; end % if no maximum go to next channel

        if Power(ch,t) > 2*sigma_max_std % else check if max > 2*std of all maxima
            max_channels{t} = [max_channels{t} ch]; % keep which channels have maxima at each 100ms-point
            chan_count(t) = chan_count(t) + 1;
        end
    end
end


%% Square, smooth (10 passes of 3-point moving average) and normalize
% "This count was squared, then convolved ten times with a three point smoothing
% filter, then normalized"

mov_ave = 1/3*[1 1 1];
num_passes = 10;

chan_count_sq = chan_count.^2;

for i=1:num_passes
%     chan_count_sq = filtfilt(mov_ave,1,chan_count_sq);
    chan_count_sq = conv(chan_count_sq,mov_ave,'same');
end

chan_count_sq = chan_count_sq/numchan;


% "the largest 80% were retained as spindles" ???

%% ========== THIS? ==========
% % Find maxima
% chan_count_sq_max = zeros(length(chan_count_sq)-1,1);
% chan_count_sq_max(diff(chan_count_sq)<0 & diff([0; chan_count_sq(1:end-1)])>0) = 1;
% % for i=1:length(chan_count_sq_max)-1
% %     if and(chan_count_sq_max(i)==1,chan_count_sq_max(i+1)==1), chan_count_sq_max(i+1)=0; end
% % end
% 
% % Extract spindles: keep 80% largest maxima
% keep = round(0.8*sum(chan_count_sq_max));
% 
% spindle_timeline = zeros(size(chan_count_sq));
% spindle_timeline(chan_count_sq_max==1) = chan_count_sq(chan_count_sq_max==1);
% 
% [~,ind] = sort(spindle_timeline,'descend');
% 
% spindle_timeline(ind(keep+1:end)) = 0;
% spindle_locs = ind(1:keep);

%% ========= OR THIS? ===========
% Extract spindles: keep 80% largest values
keep = round(0.8*length(chan_count_sq));

spindle_timeline = chan_count_sq;
[~,ind] = sort(spindle_timeline,'descend');

spindle_timeline(ind(keep+1:end)) = 0;
spindle_locs = ind(1:keep);

%% Plot
figure; plot(linspace(0,30,length(chan_count_sq)),chan_count_sq)
hold on; plot(linspace(0,30,length(chan_count_sq)),spindle_timeline,'r')
