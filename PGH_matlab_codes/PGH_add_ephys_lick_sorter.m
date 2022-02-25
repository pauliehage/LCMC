%% function add_ephys_lick_sorter
function [LICKS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS] = PGH_add_ephys_lick_sorter(path_to_rec, current_unit, params, funcs)
% load EPHYS sorted DATA
if ~strcmp(path_to_rec(end), filesep);path_to_rec = [path_to_rec filesep];end
path_to_unit = [path_to_rec, 'analyzed_data', filesep,...
    'units', filesep,current_unit,filesep];
file_name_ = dir([path_to_unit,current_unit,'_sorted_*.mat']);
file_name = file_name_(1).name;
fprintf(['Loading ', file_name, ' ... ']);
load([path_to_unit file_name],'ss_index','cs_index','sample_rate',...
    't_start','type','waveform');
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = path_to_unit;
fprintf(' --> Completed. \n')

% load EPHYS EVENT DATA
path_to_tongue = [path_to_rec,'analyzed_data',filesep,'behavior_data', filesep,...
    'tongue',filesep];
file_name_ = dir([path_to_tongue '*_EVE1_aligned.mat']);
file_name = file_name_(1).name;
EPHYS.CH_EVE = load([path_to_tongue file_name]);
EPHYS.CH_EVE.EPHYS_time_1K  = EPHYS.CH_EVE.EPHYS_time_1K(:);
EPHYS.CH_EVE.VID_time_1K = EPHYS.CH_EVE.VID_time_1K(:);

% load BEHAVE DATA
file_name_ = dir([path_to_tongue '*_ANALYZED.mat']);
file_name = file_name_(1).name;
BEHAVE = load([path_to_tongue file_name]);
BEHAVE.EXPERIMENT_PARAMS.EPHYS_file_name = EPHYS.CH_sorted_file_name;
BEHAVE.EXPERIMENT_PARAMS.EPHYS_file_path = EPHYS.CH_sorted_file_path;

% build EPHYS.CH_sorted from DATA_PSORT
ch_time = (t_start + (0:length(ss_index)-1)/sample_rate)';
SS_index = ss_index;
CS_index = cs_index;
SS_time = ch_time(SS_index);
CS_time = ch_time(CS_index);

EPHYS.CH_sorted.SS_data.SS_ind = SS_index;
EPHYS.CH_sorted.CS_data.CS_ind = CS_index;
EPHYS.CH_sorted.SS_data.SS_time = SS_time;
EPHYS.CH_sorted.CS_data.CS_time = CS_time;
EPHYS.CH_sorted.waveform = waveform;
EPHYS.CH_sorted.duration = ch_time(end) - ch_time(1);

% build SSxSS AUTO PROBABILITY
clearvars -except EPHYS BEHAVE type params funcs
SS_time   = EPHYS.CH_sorted.SS_data.SS_time;
CS_time   = EPHYS.CH_sorted.CS_data.CS_time;
Corr_data = ESN_correlogram(SS_time, CS_time);
EPHYS.CH_sorted.Corr_data = Corr_data;

% SS & CS train_aligned
clearvars -except EPHYS BEHAVE type params funcs
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K;
length_time_ = length(EPHYS_time_1K);
CS_time = EPHYS.CH_sorted.CS_data.CS_time;
if isempty(CS_time)
    CS_time = EPHYS_time_1K(1);
end
CS_time(end+1) = max([EPHYS_time_1K(end), CS_time(end)])+1;
SS_time = EPHYS.CH_sorted.SS_data.SS_time;
if isempty(SS_time)
    SS_time = EPHYS_time_1K(1);
end
SS_time(end+1) = max([EPHYS_time_1K(end), SS_time(end)])+1;
EPHYS_CS_train_1K = false(size(EPHYS_time_1K));
EPHYS_SS_train_1K = false(size(EPHYS_time_1K));
counter_CS = find(CS_time >= EPHYS_time_1K(1), 1, 'first');
counter_SS = find(SS_time >= EPHYS_time_1K(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = EPHYS_time_1K(counter_time_point);
    if time_ponit_>=CS_time(counter_CS)
        EPHYS_CS_train_1K(counter_time_point) = true;
        counter_CS = counter_CS + 1;
    end
    if time_ponit_>=SS_time(counter_SS)
        EPHYS_SS_train_1K(counter_time_point) = true;
        counter_SS = counter_SS + 1;
    end
end
EPHYS.CH_EVE.EPHYS_CS_train_1K = EPHYS_CS_train_1K;
EPHYS.CH_EVE.EPHYS_SS_train_1K = EPHYS_SS_train_1K;


% extract alignment
clearvars -except EPHYS BEHAVE type params funcs

data_type_tongue_list    = {'tongue_dm', 'tongue_vm', 'tongue_ang'};
data_type_BEHAVE_list = {'tongue_dm_stream', 'tongue_vm_stream', 'tongue_ang_stream'};
data_type_neuro_list  = {'neuro_SS', 'neuro_CS'};
data_type_EPHYS_list  = {'EPHYS_SS_train_1K', 'EPHYS_CS_train_1K'};
event_type_list       = {'onset', 'vmax', 'dmax', 'vmin', 'offset'};

length_trace = params.lick.length_trace;
inds_span    = params.lick.inds_span;


if isempty(data_type_tongue_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
fprintf(['Analyzing ', EPHYS.CH_sorted_file_name, ' ... ']);
REFRENCE_TIME = EPHYS.CH_EVE.VID_time_1K;
length_time_  = length(REFRENCE_TIME);
num_licks      = length(BEHAVE.LICKS_ALL_DATA.tag);

% Init variables
% use nan for tongue data and logical false for neuro data
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_data_type_tongue = 1 : length(data_type_tongue_list)
        data_type_tongue_name = data_type_tongue_list{counter_data_type_tongue};
        BEHAVE.LICKS_ALL_DATA.([data_type_tongue_name '_' event_type_name]) = nan(length_trace, num_licks);
    end
    for counter_data_type_neuro = 1 : length(data_type_neuro_list)
        data_type_neuro_name = data_type_neuro_list{counter_data_type_neuro};
        BEHAVE.LICKS_ALL_DATA.([data_type_neuro_name '_' event_type_name]) = false(length_trace, num_licks);
    end
end

% Build stream dataset.
% concatenate the data using the cell2mat and then interpolate the potential missing data
BEHAVE.stream = struct;
BEHAVE.stream.('time_1K') = REFRENCE_TIME;
time_1K = BEHAVE.LICKS_ALL_DATA.time_1K_stream;
for counter_variable = 1 : 1 : length(data_type_BEHAVE_list)
    variable_name = data_type_BEHAVE_list{counter_variable};
    variable = BEHAVE.LICKS_ALL_DATA.(variable_name)(:)';
    BEHAVE.stream.(variable_name) = interp1(time_1K, variable, REFRENCE_TIME, 'nearest', 'extrap');
end

% Loop over licks and add the kinematic traces to the LICKS_ALL_DATA
for counter_lick = 1 : num_licks
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        event_time = BEHAVE.LICKS_ALL_DATA.(['time' '_' event_type_name])(1,counter_lick);
        if isnan(event_time)
            % if the event_time is nan, then skip the event.
            continue;
        end
        event_ind = find(REFRENCE_TIME >= event_time, 1, 'first');
        if isempty(event_ind)
            % if the event_ind is empty, then skip the event.
            BEHAVE.LICKS_ALL_DATA.validity(1, counter_lick) = false;
            continue;
        end
        event_inds = repmat( event_ind, 1, length(inds_span)) + repmat(inds_span(:)', length(event_ind), 1);
        event_inds( event_inds < 1 ) = 1;
        event_inds( event_inds > length_time_ ) = length_time_;
        for counter_data_type_tongue = 1 : length(data_type_tongue_list)
            data_type_tongue_name = data_type_tongue_list{counter_data_type_tongue};
            data_type_BEHAVE_name = data_type_BEHAVE_list{counter_data_type_tongue};
            BEHAVE.LICKS_ALL_DATA.([data_type_tongue_name '_' event_type_name])(:,counter_lick) = ...
                reshape(BEHAVE.stream.(data_type_BEHAVE_name)(event_inds), length_trace, 1);
        end

    end
end

% Loop over licks and add neuro_SS & neuro_CS to the LICKS_ALL_DATA
REFRENCE_TIME = EPHYS.CH_EVE.align_LED.VID_LED_xcorr_time_1K(:); % the initial ind will be drawn from BEHAVE_EB_xcorr_time_1K
length_time_  = length(EPHYS.CH_EVE.EPHYS_time_1K); % the converted ind will be applied to EPHYS_time_1K
for counter_lick = 1 : num_licks
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        event_time = BEHAVE.LICKS_ALL_DATA.(['time' '_' event_type_name])(1,counter_lick);
        if isnan(event_time)
            % if the event_time is nan, then skip the event.
            continue;
        end
        event_ind = find(REFRENCE_TIME >= event_time, 1, 'first');
        if isempty(event_ind)
            % if the event_ind is empty, then skip the event.
            BEHAVE.LICKS_ALL_DATA.validity(1, counter_lick) = false;
            continue;
        end
        % set the event_ind_converted based on align_states
        event_ind_converted = EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_1K(event_ind);
        if isempty(event_ind_converted)
            % if the event_ind_converted is empty, then skip the event.
            BEHAVE.LICKS_ALL_DATA.validity(1, counter_lick) = false;
            continue;
        end
        event_inds_converted = repmat( event_ind_converted, 1, length(inds_span)) + repmat(inds_span(:)', length(event_ind_converted), 1);
        event_inds_converted( event_inds_converted < 1 ) = 1;
        event_inds_converted( event_inds_converted > length_time_ ) = length_time_;
        for counter_data_type_neuro_list = 1 : length(data_type_neuro_list)
            data_type_neuro_name = data_type_neuro_list{counter_data_type_neuro_list};
            data_type_EPHYS_name = data_type_EPHYS_list{counter_data_type_neuro_list};
            BEHAVE.LICKS_ALL_DATA.([data_type_neuro_name '_' event_type_name])(:,counter_lick) = ...
                logical(reshape(EPHYS.CH_EVE.(data_type_EPHYS_name)(event_inds_converted), length_trace, 1));
        end
    end
end

% Add neuro_CS_count_onset_offset BEHAVE.LICKS_ALL_DATA
% count the number of the the CS from lick onset to offset
BEHAVE.LICKS_ALL_DATA.neuro_CS_count_onset_offset   = zeros(size(BEHAVE.LICKS_ALL_DATA.tag));
BEHAVE.LICKS_ALL_DATA.neuro_CS_count_onset_offset_period   = zeros(size(BEHAVE.LICKS_ALL_DATA.tag));

for counter_lick = 1 : num_licks
    ind_start     = (length_trace/2);
    time_onset    = BEHAVE.LICKS_ALL_DATA.time_onset(   1, counter_lick);
    time_offset   = BEHAVE.LICKS_ALL_DATA.time_offset(   1, counter_lick);
    duration_lick = round((time_offset - time_onset)*1000.0);

    if (duration_lick <= 250)
        ind_end = ind_start + duration_lick;
        neuro_CS_count_onset_offset_period = duration_lick;
    elseif (duration_lick > 250)
        ind_end = ind_start + 250;
        neuro_CS_count_onset_offset_period = 250;
%     elseif (duration_lick <= 0)
%         ind_end = ind_start;
%         neuro_CS_count_onset_offset_period = nan;
    else
        % this condition covers the nan values
%         ind_end = ind_start;
%         neuro_CS_count_onset_offset_period = nan; 
        disp('Error');
    end

    BEHAVE.LICKS_ALL_DATA.neuro_CS_count_onset_offset_period(1,counter_lick) = neuro_CS_count_onset_offset_period;
    BEHAVE.LICKS_ALL_DATA.neuro_CS_count_onset_offset(  1, counter_lick) = nansum(BEHAVE.LICKS_ALL_DATA.neuro_CS_onset(  ind_start:ind_end,   counter_lick));
end

fprintf(' --> Completed. \n')

% Add Neural_Properties
clearvars -except EPHYS BEHAVE type params funcs
EPHYS.Neural_Properties = struct();
EPHYS.Neural_Properties.type = type;
EPHYS.Neural_Properties.waveform = EPHYS.CH_sorted.waveform;

if length(EPHYS.CH_sorted.SS_data.SS_time)>1
    EPHYS.Neural_Properties.SS_num = length(EPHYS.CH_sorted.SS_data.SS_time);
    EPHYS.Neural_Properties.SS_duration = EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.SS_firing_rate = EPHYS.Neural_Properties.SS_num / EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.SS_time = EPHYS.CH_sorted.SS_data.SS_time;
else
    EPHYS.Neural_Properties.SS_num = 0;
    EPHYS.Neural_Properties.SS_duration = EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.SS_firing_rate = 0;
    EPHYS.Neural_Properties.SS_time = [];
end

if length(EPHYS.CH_sorted.CS_data.CS_time)>1
    EPHYS.Neural_Properties.CS_num = length(EPHYS.CH_sorted.CS_data.CS_time);
    EPHYS.Neural_Properties.CS_firing_rate = EPHYS.Neural_Properties.CS_num / EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.CS_time = EPHYS.CH_sorted.CS_data.CS_time;
else
    EPHYS.Neural_Properties.CS_num = 0;
    EPHYS.Neural_Properties.CS_firing_rate = 0;
    EPHYS.Neural_Properties.CS_time = [];
end
EPHYS.Neural_Properties.Corr_data_CS_inds_span     = nanmean(EPHYS.CH_sorted.Corr_data.CS_inds_span);
EPHYS.Neural_Properties.Corr_data_CS_bin_size_time = nanmean(EPHYS.CH_sorted.Corr_data.CS_bin_size_time);
EPHYS.Neural_Properties.Corr_data_SS_inds_span     = nanmean(EPHYS.CH_sorted.Corr_data.SS_inds_span);
EPHYS.Neural_Properties.Corr_data_SS_bin_size_time = nanmean(EPHYS.CH_sorted.Corr_data.SS_bin_size_time);
EPHYS.Neural_Properties.Corr_data_SS_SSxSS_AUTO    = nanmean(EPHYS.CH_sorted.Corr_data.SS_SSxSS_AUTO);
EPHYS.Neural_Properties.Corr_data_CS_CSxSS_AUTO    = nanmean(EPHYS.CH_sorted.Corr_data.CS_CSxSS_AUTO);

% remove fields
% BEHAVE.LICKS_ALL_DATA = rmfield(BEHAVE.LICKS_ALL_DATA, {'time_1K_stream', 'tongue_dm_stream', 'tongue_vm_stream', ...
%     'tongue_ang_stream', 'tongue_dm', 'tongue_vm', 'tongue_ang', 'tongue_tip_px', 'tongue_tip_py', 'tongue_r_px', ...
%     'tongue_r_py', 'tongue_l_px', 'tongue_l_py', 'tongue_mid_px', 'tongue_mid_py', 'nose_r_px', 'nose_r_py', ...
%      'nose_l_px', 'nose_l_py', 'rew_r_px', 'rew_r_py', 'rew_l_px', 'rew_l_py', 'rtube_r_px', 'rtube_r_py', ...
%      'rtube_l_px', 'rtube_l_py', 'ltube_r_px', 'ltube_r_py', 'ltube_l_px', 'ltube_l_py', 'rew_capacity_r', 'rew_capacity_l' });

BEHAVE.LICKS_ALL_DATA = rmfield(BEHAVE.LICKS_ALL_DATA, {'time_1K_stream', 'tongue_dm_stream', 'tongue_vm_stream', ...
    'tongue_ang_stream', 'tongue_dm', 'tongue_vm', 'tongue_ang', 'rew_capacity_r', 'rew_capacity_l' });


% outputs
LICKS_ALL_DATA = BEHAVE.LICKS_ALL_DATA;
Neural_Properties = EPHYS.Neural_Properties;
EXPERIMENT_PARAMS = BEHAVE.EXPERIMENT_PARAMS;
end

