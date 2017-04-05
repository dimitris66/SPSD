% clear all; close all; clc

addpath('/cluster/manoach/Dimitris/Matlab_Scripts/Preprocessing_Pipeline');
addpath(genpath('/cluster/manoach/K24_MEG_EEG_Study/Code/Spindle_Detection_Algorithms/depedencies'));
addpath('/cluster/manoach/K24_MEG_EEG_Study/Code/');

addpath(genpath('~/gitcode/SPSD'));

% cd '/cluster/manoach/K24_MEG_EEG_Study/Code/Nimas_algo'

root_dir = '/cluster/manoach/K24_MEG_EEG_Study/K24_data_Preprocessed/';
subj_list = {'KC_01','KC_02','KC_03','KC_04','KC_05','KC_06','KC_07','KC_08','KS_01','KS_08','KS_09'};
visits_list={'MST','Baseline'};
ss=3;vv=2;

Fs = 200; % !! DOWNSAMPLED !!
tthresh = '400ms';
%--------------------------------------------------

tthreshD = str2double(tthresh(1:3));

source_dir = fullfile(root_dir,subj_list{ss},visits_list{vv});

if strcmp(visits_list{vv},'MST')
    re = '(N\d|0[13])_ArtRej_raw.fif';
elseif strcmp(visits_list{vv},'Baseline')
    re = '(N\d|01)_ArtRej_raw.fif';
end
source_files = dir(fullfile(source_dir, '/K*'));
source_files2 = struct2cell(source_files);
source_files2 = source_files2(1,:);
source_files = source_files(~cellfun(@isempty,regexp(source_files2,re,'once')));
clear source_files2

for ff = 1%:length(source_files)
    
    clear Data data
    
    [raw] = fiff_setup_read_raw(fullfile(source_dir,source_files(ff).name));
    Data = fiff_read_raw_segment(raw,raw.first_samp,raw.last_samp);
    [Data,~,~,~] = fun_activate_projections(raw,Data);
    
    [~,fname,~] = fileparts(source_files(ff).name);
    if strcmp(visits_list{vv},'MST'), expr = 'K[CS]{1}_\d\d_(N\d|0[13])'; else expr = 'K[CS]{1}_\d\d_(N\d|01)'; end
    [a,aa] = regexpi(fname,expr);
    fname = fname(a:aa);
    
    % find MEG-EEG channels of interest
    load(fullfile('/cluster/manoach/K24_MEG_EEG_Study/Code', 'K24_channel_names.mat'));
    [~,channels] = ismember(channel_names,raw.info.ch_names);
    
    [chans_eeg,~] = fun_Pick_EEG_MEG_Channels(raw,'EEG');
    [chans_meg,~] = fun_Pick_EEG_MEG_Channels(raw,'MEG');
    [chans_mag,~] = fun_Pick_EEG_MEG_Channels(raw,'MAG');
    [chans_grad,~] = fun_Pick_EEG_MEG_Channels(raw,'GRAD');
    
    % exclude bad channels
    if ~isempty(raw.info.bads)
        [~,b] = ismember(raw.info.bads,raw.info.ch_names);
        Data(b,:) = NaN;
    end
    
    display(source_files(ff).name)
    display(strcat('size of data = ', num2str(size(Data))))
    
    numEpochs = ceil(length(Data)/(Fs*30));
    data_segms = cell(numEpochs,1);
    for i=1:numEpochs-1
        data_segms{i} = Data(:,(i-1)*Fs*30+1:i*Fs*30);
    end
    data_segms{end} = Data(:,i*Fs*30+1:end);
    data_segms{end} = fun_padEpoch_signal(data_segms{end},Fs,30);
    
    numEpochs = ceil(length(Data)/(Fs*30));
    data_segms = cell(numEpochs,1);
    for i=1:numEpochs-1
        data_segms{i} = Data(:,(i-1)*Fs*30+1:i*Fs*30);
    end
    data_segms{end} = Data(:,i*Fs*30+1:end);
    data_segms{end} = fun_padEpoch_signal(data_segms{end},Fs,30);
    
    display('Running Sequential Power Spectral Density algorithm for spindle detection.');
    
    spindle_locs = cell(numEpochs,3);
    
    % Gio version
    for ee=1:numEpochs
            [spindle_locs{ee,1}] = fun_ndehghani_detect_spindles_v1_Test2(data_segms{ee}(chans_eeg,:),Fs,'EEG');
            [spindle_locs{ee,2}] = fun_ndehghani_detect_spindles_v1_Test2(data_segms{ee}(chans_mag,:),Fs,'MEG');
            [spindle_locs{ee,3}] = fun_ndehghani_detect_spindles_v1_Test2(data_segms{ee}(chans_grad,:),Fs,'MEG');
    end
    
    % Nima version
    for ee=1:numEpochs
            [spindle_locs{ee,1}] = fun_ndehghani_detect_spindles_v1_Nima(data_segms{ee}(chans_eeg,:),Fs,'EEG');
            [spindle_locs{ee,2}] = fun_ndehghani_detect_spindles_v1_Nima(data_segms{ee}(chans_mag,:),Fs,'MEG');
            [spindle_locs{ee,3}] = fun_ndehghani_detect_spindles_v1_Nima(data_segms{ee}(chans_grad,:),Fs,'MEG');
    end
    
    dest_dir = '/cluster/manoach/K24_MEG_EEG_Study/Code/Spindle_Detection_Algorithms/SPSD/scratch';
    save(fullfile(dest_dir,[fname '_Nimas_N2']),'spindle_locs','-v7.3');
    display('succesfully finished spindle detection')
    
end



