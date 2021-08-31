function ESN_event_alignment_v4_rand_corr_post201911_COMBINED
%% clear
clear;

%% load EPHYS EVENT DATA
file_path = pwd;
[file_name,file_path] = uigetfile([file_path filesep 'all_channels.events'], 'Select all_channels.events file');
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
[ch_data, ch_time, ch_info] = load_open_ephys_data([file_path file_name]);
EPHYS.CH_EVE.ch_data = ch_data;
EPHYS.CH_EVE.ch_time = ch_time;
EPHYS.CH_EVE.ch_info = ch_info;
EPHYS.file_name_CH_EVE = file_name;
EPHYS.file_path_CH_EVE = file_path;
fprintf(' --> Completed. \n')

EPHYS.debug_figures = false;

%% load ADC1 and ADC2
[file_name,file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep '*.continuous'], 'Select .continuous file');
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
[ch_data, ch_time, ch_info] = load_open_ephys_data([file_path file_name]);
EPHYS.ADC1.ch_data = ch_data;
EPHYS.ADC1.ch_time = ch_time;
EPHYS.ADC1.ch_info = ch_info;
EPHYS.time_30K   = double(EPHYS.ADC1.ch_time);
fprintf(' --> Completed. \n')

%% load BEHAVE DATA
[file_name,file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep '*.mat'], 'Select _corrective_saccades file');
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
data = load([file_path file_name],'data');
BEHAVE = data.data;
fprintf(' --> Completed. \n')

%% BEHAVE Data
clearvars -except EPHYS BEHAVE
fprintf(['Filter BEHAVE Eye Data', ' ... ']);

min_length = min([ ...
    length(BEHAVE.eyelink_time),...
    length(BEHAVE.t),...
    length(BEHAVE.state_value),...
    ]);
% timeseries
inds_invalid   = false(min_length, 1);
time_eyelink   = double(BEHAVE.eyelink_time(1:min_length)');          inds_invalid = isnan(time_eyelink) | inds_invalid;
time_tgt       = double(BEHAVE.t(1:min_length)');                     inds_invalid = isnan(time_tgt)     | inds_invalid;
state          = double(BEHAVE.state_value(1:min_length)');           inds_invalid = isnan(state)        | inds_invalid;
% correct for the bias between time_eyelink and time_tgt
time_eyelink   = time_eyelink .* (time_tgt(end)-time_tgt(1)) ./ (time_eyelink(end)-time_eyelink(1));
time_eyelink   = time_eyelink - time_eyelink(1) + time_tgt(1);
time_1K        = (time_eyelink(1) : 0.001 : time_eyelink(end))';
time_100        = (time_eyelink(1) : 0.01 : time_eyelink(end))';

% make non unique points of eye traces invalid
inds_invalid = ([false; (diff(time_eyelink)==0)]) | inds_invalid;
% remove invalid values
time_eyelink(inds_invalid) = [];
state(       inds_invalid) = [];

% reconstruct data
state_1K    = interp1(time_eyelink, state,    time_1K, 'nearest','extrap');

% reconstruct data 100
state_100    = interp1(time_eyelink, state,    time_100, 'nearest','extrap');

BEHAVE.time_100  = time_100;
BEHAVE.time_1K  = time_1K;
BEHAVE.state_1K    = state_1K;
BEHAVE.state_100    = state_100;

fprintf(' --> Completed. \n')

%% extract event data in ephys
clearvars -except EPHYS BEHAVE
fprintf(['Building EPHYS Event Data', ' ... ']);
% ch_data
% 2 : STR_TARGET_PURSUIT 0
% 3 : STR_TARGET_FIXATION 1
% - : DETECT_SACCADE_START 2
% 4 : DETECT_SACCADE_END 3
% - : END_TARGET_FIXATION 4
% 5 : 1Hz ossilation 5
% 6 : photodiode: STR_TARGET_FIXATION+DETECT_SACCADE_END
% eventId
% 1 : rising
% 0 : falling
EPHYS.CH_EVE.data = [EPHYS.CH_EVE.ch_time(:) EPHYS.CH_EVE.ch_data(:) EPHYS.CH_EVE.ch_info.eventId(:)];
fprintf(' --> Completed. \n')
%{
enum state {
0   INIT = 0,
1	STR_TARGET_PURSUIT,
2	STR_TARGET_PRESENT,
3	STR_TARGET_FIXATION,
4	CUE_TARGET_PRESENT,
5	DETECT_SACCADE_START,
6	SACCADE,
7	DETECT_SACCADE_END,
8	DELIVER_REWARD,
9	END_TARGET_FIXATION,
10	INCORRECT_SACCADE,
11	ITI };
%}

%% 1K
    %% Build EPHYS Alignment events
clearvars -except EPHYS BEHAVE
fprintf(['Building EPHYS Alignment events', ' ... ']);

time_reference      = ( ESN_Round(EPHYS.time_30K(1),0.001) : 0.001 : ESN_Round(EPHYS.time_30K(end),0.001) )';
length_time         = length(time_reference);
time_STR_TARGET_PURSUIT_rise = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 2) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_STR_TARGET_PURSUIT_fall = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 2) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
time_STR_TARGET_FIXATION_rise = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 3) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_STR_TARGET_FIXATION_fall = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 3) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
time_DETECT_SACCADE_END_rise = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_DETECT_SACCADE_END_fall = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
time_photodiode_rise      = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 6) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_photodiode_fall      = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 6) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
variable_list = {...
    '_STR_TARGET_PURSUIT_rise' ,'_STR_TARGET_PURSUIT_fall', ...
    '_STR_TARGET_FIXATION_rise','_STR_TARGET_FIXATION_fall', ...
    '_DETECT_SACCADE_END_rise' ,'_DETECT_SACCADE_END_fall', ...
    '_photodiode_rise','_photodiode_fall'};

if isempty(time_photodiode_rise)
    time_photodiode_rise = [time_STR_TARGET_FIXATION_rise; time_DETECT_SACCADE_END_rise];
    time_photodiode_rise = sort(time_photodiode_rise);
end

if isempty(time_photodiode_fall)
    time_photodiode_fall = [time_STR_TARGET_FIXATION_fall; time_DETECT_SACCADE_END_fall];
    time_photodiode_fall = sort(time_photodiode_fall);
end

EPHYS.Alignment.time_1K = time_reference;
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'EPHYS.Alignment.time' variable_name ' = ' 'time' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'time_temp_' ' = ' 'time' variable_name ';']);
    time_temp_(end+1) = max([time_reference(end), time_temp_(end)])+1;
    event_temp_       = false(length_time, 1);
    counter_temp_     = find(time_temp_ >= time_reference(1), 1, 'first');
    eval([ 'time'    variable_name ' = ' 'time_temp_'    ';']);
    eval([ 'event'   variable_name ' = ' 'event_temp_'   ';']);
    eval([ 'counter' variable_name ' = ' 'counter_temp_' ';']);
end
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    %{
    % Do not use 'eval', it is very slow, instead use the actual variables
    for counter_variable = 1 : 1 : length(variable_list)
        variable_name = variable_list{counter_variable};
        eval([ ...
            'if ( time_ponit_ >= time_' variable_name '(  counter_' variable_name ') );', ...
                'event_' variable_name '(    counter_time_point) = true;', ...
                'counter_' variable_name '   = counter_' variable_name '   + 1;', ...
            'end;', ...
        ]);
    end
    % Here is the template for actual variables
    if time_ponit_ >= time_VARIABLE(  counter_VARIABLE)
        event_VARIABLE(    counter_time_point) = true;
        counter_VARIABLE   = counter_VARIABLE   + 1;
    end
    %}
    % STR_TARGET_PURSUIT_rise
    if time_ponit_ >= time_STR_TARGET_PURSUIT_rise(  counter_STR_TARGET_PURSUIT_rise)
        event_STR_TARGET_PURSUIT_rise(    counter_time_point) = true;
        counter_STR_TARGET_PURSUIT_rise   = counter_STR_TARGET_PURSUIT_rise   + 1;
    end
    % STR_TARGET_PURSUIT_fall
    if time_ponit_ >= time_STR_TARGET_PURSUIT_fall(  counter_STR_TARGET_PURSUIT_fall)
        event_STR_TARGET_PURSUIT_fall(    counter_time_point) = true;
        counter_STR_TARGET_PURSUIT_fall   = counter_STR_TARGET_PURSUIT_fall   + 1;
    end
    % STR_TARGET_FIXATION_rise
    if time_ponit_ >= time_STR_TARGET_FIXATION_rise(  counter_STR_TARGET_FIXATION_rise)
        event_STR_TARGET_FIXATION_rise(    counter_time_point) = true;
        counter_STR_TARGET_FIXATION_rise   = counter_STR_TARGET_FIXATION_rise   + 1;
    end
    %STR_TARGET_FIXATION_fall
    if time_ponit_ >= time_STR_TARGET_FIXATION_fall(  counter_STR_TARGET_FIXATION_fall)
        event_STR_TARGET_FIXATION_fall(    counter_time_point) = true;
        counter_STR_TARGET_FIXATION_fall   = counter_STR_TARGET_FIXATION_fall   + 1;
    end
    % DETECT_SACCADE_END_rise
    if time_ponit_ >= time_DETECT_SACCADE_END_rise(  counter_DETECT_SACCADE_END_rise)
        event_DETECT_SACCADE_END_rise(    counter_time_point) = true;
        counter_DETECT_SACCADE_END_rise   = counter_DETECT_SACCADE_END_rise   + 1;
    end
    % DETECT_SACCADE_END_fall
    if time_ponit_ >= time_DETECT_SACCADE_END_fall(  counter_DETECT_SACCADE_END_fall)
        event_DETECT_SACCADE_END_fall(    counter_time_point) = true;
        counter_DETECT_SACCADE_END_fall   = counter_DETECT_SACCADE_END_fall   + 1;
    end
    % photodiode_rise
    if time_ponit_ >= time_photodiode_rise(     counter_photodiode_rise)
        event_photodiode_rise(       counter_time_point) = true;
        counter_photodiode_rise      = counter_photodiode_rise      + 1;
    end
    % photodiode_fall
    if time_ponit_ >= time_photodiode_fall(     counter_photodiode_fall)
        event_photodiode_fall(       counter_time_point) = true;
        counter_photodiode_fall      = counter_photodiode_fall      + 1;
    end
end

event_STR_TARGET_PURSUIT  = false(length_time, 1);
event_STR_TARGET_FIXATION = false(length_time, 1);
event_DETECT_SACCADE_END  = false(length_time, 1);
event_photodiode          = false(length_time, 1);
flag_STR_TARGET_PURSUIT   = false;
flag_STR_TARGET_FIXATION  = false;
flag_DETECT_SACCADE_END   = false;
flag_photodiode           = false;
for counter_time_point = 1 : length_time
    flag_STR_TARGET_PURSUIT  = flag_STR_TARGET_PURSUIT  ||   event_STR_TARGET_PURSUIT_rise( counter_time_point);
    flag_STR_TARGET_PURSUIT  = flag_STR_TARGET_PURSUIT  && (~event_STR_TARGET_PURSUIT_fall( counter_time_point));
    flag_STR_TARGET_FIXATION = flag_STR_TARGET_FIXATION ||   event_STR_TARGET_FIXATION_rise(counter_time_point);
    flag_STR_TARGET_FIXATION = flag_STR_TARGET_FIXATION && (~event_STR_TARGET_FIXATION_fall(counter_time_point));
    flag_DETECT_SACCADE_END  = flag_DETECT_SACCADE_END  ||   event_DETECT_SACCADE_END_rise( counter_time_point);
    flag_DETECT_SACCADE_END  = flag_DETECT_SACCADE_END  && (~event_DETECT_SACCADE_END_fall( counter_time_point));
    flag_photodiode          = flag_photodiode          ||   event_photodiode_rise(         counter_time_point);
    flag_photodiode          = flag_photodiode          && (~event_photodiode_fall(         counter_time_point));
    event_STR_TARGET_PURSUIT( counter_time_point) = flag_STR_TARGET_PURSUIT;
    event_STR_TARGET_FIXATION(counter_time_point) = flag_STR_TARGET_FIXATION;
    event_DETECT_SACCADE_END( counter_time_point) = flag_DETECT_SACCADE_END;
    event_photodiode(         counter_time_point) = flag_photodiode;
end

variable_list = {...
    '_STR_TARGET_PURSUIT_rise' ,'_STR_TARGET_PURSUIT_fall', '_STR_TARGET_PURSUIT', ...
    '_STR_TARGET_FIXATION_rise','_STR_TARGET_FIXATION_fall', '_STR_TARGET_FIXATION', ...
    '_DETECT_SACCADE_END_rise' ,'_DETECT_SACCADE_END_fall', '_DETECT_SACCADE_END', ...
    '_photodiode_rise','_photodiode_fall', '_photodiode'};
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'EPHYS.Alignment.event' variable_name ' = ' 'event' variable_name ';']);
end

event_state_combined = ...
    double(event_STR_TARGET_PURSUIT)  .* 1 + ...
    double(event_STR_TARGET_FIXATION) .* 3 + ...
    double(event_DETECT_SACCADE_END)  .* 7 ;

event_photodiode_combined = ...
    double(event_photodiode)  .* 1;

EPHYS.Alignment.event_state_combined      = event_state_combined;
EPHYS.Alignment.event_photodiode_combined = event_photodiode_combined;

state_description = [
'STR_TARGET_PURSUIT: 1 , ', ...
'STR_TARGET_FIXATION: 3 , ', ...
'DETECT_SACCADE_END: 7 , '];

photodiode_description = [
'STR_TARGET_FIXATION: 1 , ', ...
'DETECT_SACCADE_END: 1 , ', ...
];
EPHYS.Alignment.state_description      = state_description;
EPHYS.Alignment.photodiode_description = photodiode_description;

fprintf(' --> Completed. \n');
    %% Build BEHAVE Alignment events
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE Alignment events', ' ... ']);
time_reference      = BEHAVE.time_1K;
BEHAVE.Alignment.time_1K = time_reference;
BEHAVE_state = BEHAVE.state_1K(:);
event_STR_TARGET_PURSUIT  = (BEHAVE_state == 1);
event_STR_TARGET_FIXATION = (BEHAVE_state == 3);
event_DETECT_SACCADE_END =  (BEHAVE_state == 7);
event_photodiode = event_STR_TARGET_FIXATION | event_DETECT_SACCADE_END;

event_state_combined = ...
    double(event_STR_TARGET_PURSUIT)  .* 1 + ...
    double(event_STR_TARGET_FIXATION) .* 3 + ...
    double(event_DETECT_SACCADE_END)  .* 7 ;

event_photodiode_combined = ...
    double(event_photodiode)  .* 1;

BEHAVE.Alignment.event_state_combined      = event_state_combined;
BEHAVE.Alignment.event_photodiode_combined = event_photodiode_combined;

state_description = [
'STR_TARGET_PURSUIT: 1 , ', ...
'STR_TARGET_FIXATION: 3 , ', ...
'DETECT_SACCADE_END: 7 , '];

photodiode_description = [
'STR_TARGET_FIXATION: 1 , ', ...
'DETECT_SACCADE_END: 1 , ', ...
];
BEHAVE.Alignment.state_description      = state_description;
BEHAVE.Alignment.photodiode_description = photodiode_description;

fprintf(' --> Completed. \n');
    %% ALIGN EPHYS and BEHAVE state_combined through xcorr and dtw
clearvars -except EPHYS BEHAVE
fprintf(['Aligning EPHYS and BEHAVE state_combined', ' ... ']);
EPHYS_time_1K              = EPHYS.Alignment.time_1K;
EPHYS_time_30K             = EPHYS.time_30K;
BEHAVE_time_1K             = BEHAVE.Alignment.time_1K;
EPHYS_state_combined       = EPHYS.Alignment.event_state_combined;
BEHAVE_state_combined      = BEHAVE.Alignment.event_state_combined;

% state_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(EPHYS_state_combined, BEHAVE_state_combined); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);

if  sample_diff > 0
    EPHYS_EB_xcorr_state_combined_1K  = EPHYS_state_combined(  abs(sample_diff):end);
    EPHYS_EB_xcorr_time_1K            = EPHYS_time_1K(            abs(sample_diff):end);
    BEHAVE_EB_xcorr_state_combined_1K = BEHAVE_state_combined;
    BEHAVE_EB_xcorr_time_1K           = BEHAVE_time_1K;
elseif sample_diff < 0
    EPHYS_EB_xcorr_state_combined_1K  = EPHYS_state_combined;
    EPHYS_EB_xcorr_time_1K            = EPHYS_time_1K;
    BEHAVE_EB_xcorr_state_combined_1K = BEHAVE_state_combined( abs(sample_diff):end);
    BEHAVE_EB_xcorr_time_1K           = BEHAVE_time_1K(           abs(sample_diff):end);
end
% state_combined: make the vectors the same size
if length(BEHAVE_EB_xcorr_state_combined_1K) ~= length(EPHYS_EB_xcorr_state_combined_1K)
    min_length = min([ length(BEHAVE_EB_xcorr_state_combined_1K),  length(EPHYS_EB_xcorr_state_combined_1K)]);
    EPHYS_EB_xcorr_state_combined_1K  = EPHYS_EB_xcorr_state_combined_1K(  1:min_length);
    EPHYS_EB_xcorr_time_1K            = EPHYS_EB_xcorr_time_1K(            1:min_length);
    BEHAVE_EB_xcorr_state_combined_1K = BEHAVE_EB_xcorr_state_combined_1K( 1:min_length);
    BEHAVE_EB_xcorr_time_1K           = BEHAVE_EB_xcorr_time_1K(           1:min_length);
end

% %% state_combined: find the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
sampling_freq = 1000.0;
cutoff_freq = 25.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
EPHYS_EB_xcorr_state_combined_1K_filt  = filtfilt(b_butter,a_butter,EPHYS_EB_xcorr_state_combined_1K);
BEHAVE_EB_xcorr_state_combined_1K_filt = filtfilt(b_butter,a_butter,BEHAVE_EB_xcorr_state_combined_1K);
% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(EPHYS_EB_xcorr_time_1K ) / 500);
ind_edges = round(linspace(1, length(EPHYS_EB_xcorr_time_1K ), ind_edge_width));
ind_edges(1) = 0;
% init and loop over chunks
EPHYS_EB_inds_DTW  = cell((length(ind_edges)-1), 1);
BEHAVE_EB_inds_DTW = cell((length(ind_edges)-1), 1);
for counter_chunk = 1 : 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    EPHYS_EB_state_combined_chunk  = EPHYS_EB_xcorr_state_combined_1K_filt(inds_chunk);
    BEHAVE_EB_state_combined_chunk = BEHAVE_EB_xcorr_state_combined_1K_filt(inds_chunk);
    [~,ix,iy] = dtw(EPHYS_EB_state_combined_chunk,BEHAVE_EB_state_combined_chunk, 15, 'absolute');  % allow upto 15ms warp
    EPHYS_EB_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    BEHAVE_EB_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end
EPHYS_EB_inds_DTW  = cell2mat(EPHYS_EB_inds_DTW);
BEHAVE_EB_inds_DTW = cell2mat(BEHAVE_EB_inds_DTW);
EPHYS_EB_inds      = ( 1 : 1 : length(EPHYS_EB_xcorr_time_1K ) )';
BEHAVE_EB_inds     = ( 1 : 1 : length(BEHAVE_EB_xcorr_time_1K) )';
% dtw works by replicating the inds to match the two signals, here we
% reverse the replicated inds to generate two matched signals but with the
% size of original signals.
EB_ind_convert_from_EPHYS_to_BEHAVE = nan(size(EPHYS_EB_inds));
EB_ind_convert_from_BEHAVE_to_EPHYS = nan(size(BEHAVE_EB_inds));
for counter_ind = 1 : 1 : length(EPHYS_EB_inds_DTW)
    ind_EPHYS_EB_DTW  = EPHYS_EB_inds_DTW(counter_ind);
    ind_BEHAVE_EB_DTW = BEHAVE_EB_inds_DTW(counter_ind);
    EB_ind_convert_from_EPHYS_to_BEHAVE(ind_BEHAVE_EB_DTW) = ind_EPHYS_EB_DTW;
    EB_ind_convert_from_BEHAVE_to_EPHYS(ind_EPHYS_EB_DTW)  = ind_BEHAVE_EB_DTW;
end

time_reference      = EPHYS_time_30K(:);
length_time         = length(time_reference);
time_EPHYS_EB_xcorr_1K = EPHYS_EB_xcorr_time_1K;
time_EPHYS_EB_xcorr_1K(end+1) = max([time_reference(end), time_EPHYS_EB_xcorr_1K(end)])+1;
event_EPHYS_EB_xcorr_30K       = nan(length(EPHYS_EB_xcorr_time_1K), 1);
counter_EPHYS_EB_xcorr     = find(time_EPHYS_EB_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_EB_xcorr_1K(  counter_EPHYS_EB_xcorr)
        event_EPHYS_EB_xcorr_30K(    counter_EPHYS_EB_xcorr) = counter_time_point;
        counter_EPHYS_EB_xcorr   = counter_EPHYS_EB_xcorr   + 1;
    end
end

time_reference      = EPHYS_time_1K(:);
length_time         = length(time_reference);
time_EPHYS_EB_xcorr_1K = EPHYS_EB_xcorr_time_1K;
time_EPHYS_EB_xcorr_1K(end+1) = max([time_reference(end), time_EPHYS_EB_xcorr_1K(end)])+1;
event_EPHYS_EB_xcorr_1K       = nan(length(EPHYS_EB_xcorr_time_1K), 1);
counter_EPHYS_EB_xcorr     = find(time_EPHYS_EB_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_EB_xcorr_1K(  counter_EPHYS_EB_xcorr)
        event_EPHYS_EB_xcorr_1K(    counter_EPHYS_EB_xcorr) = counter_time_point;
        counter_EPHYS_EB_xcorr   = counter_EPHYS_EB_xcorr   + 1;
    end
end

time_reference      = BEHAVE_time_1K(:);
length_time         = length(time_reference);
time_BEHAVE_EB_xcorr_1K = BEHAVE_EB_xcorr_time_1K;
time_BEHAVE_EB_xcorr_1K(end+1) = max([time_reference(end), time_BEHAVE_EB_xcorr_1K(end)])+1;
event_BEHAVE_EB_xcorr_1K       = nan(length(BEHAVE_EB_xcorr_time_1K), 1);
counter_BEHAVE_EB_xcorr     = find(time_BEHAVE_EB_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_BEHAVE_EB_xcorr_1K(  counter_BEHAVE_EB_xcorr)
        event_BEHAVE_EB_xcorr_1K(    counter_BEHAVE_EB_xcorr) = counter_time_point;
        counter_BEHAVE_EB_xcorr   = counter_BEHAVE_EB_xcorr   + 1;
    end
end

EPHYS_EB_xcorr_ind_30K   = event_EPHYS_EB_xcorr_30K;
EPHYS_EB_xcorr_ind_1K    = event_EPHYS_EB_xcorr_1K;
BEHAVE_EB_xcorr_ind_1K   = event_BEHAVE_EB_xcorr_1K;
EPHYS_EB_aligned_ind_30K = EPHYS_EB_xcorr_ind_30K(EB_ind_convert_from_EPHYS_to_BEHAVE);
EPHYS_EB_aligned_ind_1K  = EPHYS_EB_xcorr_ind_1K( EB_ind_convert_from_EPHYS_to_BEHAVE);
BEHAVE_EB_aligned_ind_1K = BEHAVE_EB_xcorr_ind_1K(EB_ind_convert_from_BEHAVE_to_EPHYS);

EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_30K          = EPHYS_EB_aligned_ind_30K;
EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K           = EPHYS_EB_aligned_ind_1K;
EPHYS.CH_EVE.align_states.BEHAVE_EB_aligned_ind_1K          = BEHAVE_EB_aligned_ind_1K;
EPHYS.CH_EVE.align_states.EB_ind_convert_from_BEHAVE_to_EPHYS = EB_ind_convert_from_BEHAVE_to_EPHYS;
EPHYS.CH_EVE.align_states.EB_ind_convert_from_EPHYS_to_BEHAVE = EB_ind_convert_from_EPHYS_to_BEHAVE;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_time_1K            = EPHYS_EB_xcorr_time_1K;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_time_1K           = BEHAVE_EB_xcorr_time_1K;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_ind_30K            = EPHYS_EB_xcorr_ind_30K;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_ind_1K             = EPHYS_EB_xcorr_ind_1K;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K            = BEHAVE_EB_xcorr_ind_1K;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_state_combined_1K  = EPHYS_EB_xcorr_state_combined_1K;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_state_combined_1K = BEHAVE_EB_xcorr_state_combined_1K;
% EPHYS.CH_EVE.align_states.state_description            = BEHAVE.Alignment.state_description;
fprintf(' --> Completed. \n');

if EPHYS.debug_figures
clf(figure(1));subplot(2,1,1);plot(EPHYS_state_combined);subplot(2,1,2);plot(BEHAVE_state_combined);
clf(figure(2));subplot(2,1,1);plot(EPHYS_EB_xcorr_state_combined_1K);subplot(2,1,2);plot(BEHAVE_EB_xcorr_state_combined_1K);
clf(figure(3));subplot(2,1,1);plot(EPHYS_EB_xcorr_state_combined_1K);subplot(2,1,2);plot(BEHAVE_EB_xcorr_state_combined_1K(EB_ind_convert_from_BEHAVE_to_EPHYS));
end
disp(['length EPHYS: ' num2str( EPHYS.Alignment.time_1K(end)-EPHYS.Alignment.time_1K(1) )])
disp(['length BEHAVE: ' num2str( BEHAVE.Alignment.time_1K(end)-BEHAVE.Alignment.time_1K(1) )])
disp(['xcorr diff: ' num2str( sample_diff )])
    %% ALIGN EPHYS and BEHAVE photodiode_combined through xcorr and dtw
clearvars -except EPHYS BEHAVE
fprintf(['Aligning EPHYS and BEHAVE photodiode_combined', ' ... ']);
EPHYS_time_1K              = EPHYS.Alignment.time_1K;
EPHYS_time_30K             = EPHYS.time_30K;
BEHAVE_time_1K             = BEHAVE.Alignment.time_1K;
EPHYS_photodiode_combined  = EPHYS.Alignment.event_photodiode_combined;
BEHAVE_photodiode_combined = BEHAVE.Alignment.event_photodiode_combined;

% photodiode_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(EPHYS_photodiode_combined+1, BEHAVE_photodiode_combined+1); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);
if sample_diff > 0
    EPHYS_PD_xcorr_photodiode_combined_1K  = EPHYS_photodiode_combined(  abs(sample_diff):end);
    EPHYS_PD_xcorr_time_1K                 = EPHYS_time_1K(                 abs(sample_diff):end);
    BEHAVE_PD_xcorr_photodiode_combined_1K = BEHAVE_photodiode_combined;
    BEHAVE_PD_xcorr_time_1K                = BEHAVE_time_1K;
elseif sample_diff < 0
    EPHYS_PD_xcorr_photodiode_combined_1K  = EPHYS_photodiode_combined;
    EPHYS_PD_xcorr_time_1K                 = EPHYS_time_1K;
    BEHAVE_PD_xcorr_photodiode_combined_1K = BEHAVE_photodiode_combined( abs(sample_diff):end);
    BEHAVE_PD_xcorr_time_1K                = BEHAVE_time_1K(                abs(sample_diff):end);
end
% photodiode_combined: make the vectors the same size
if length(BEHAVE_PD_xcorr_photodiode_combined_1K) ~= length(EPHYS_PD_xcorr_photodiode_combined_1K)
    min_length = min([ length(BEHAVE_PD_xcorr_photodiode_combined_1K),  length(EPHYS_PD_xcorr_photodiode_combined_1K)]);
    EPHYS_PD_xcorr_photodiode_combined_1K  = EPHYS_PD_xcorr_photodiode_combined_1K(  1:min_length);
    EPHYS_PD_xcorr_time_1K                 = EPHYS_PD_xcorr_time_1K(                 1:min_length);
    BEHAVE_PD_xcorr_photodiode_combined_1K = BEHAVE_PD_xcorr_photodiode_combined_1K( 1:min_length);
    BEHAVE_PD_xcorr_time_1K                = BEHAVE_PD_xcorr_time_1K(                1:min_length);
end

% photodiode_combined: find the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
sampling_freq = 1000.0;
cutoff_freq = 25.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
EPHYS_PD_xcorr_photodiode_combined_1K_filt  = filtfilt(b_butter,a_butter,EPHYS_PD_xcorr_photodiode_combined_1K);
BEHAVE_PD_xcorr_photodiode_combined_1K_filt = filtfilt(b_butter,a_butter,BEHAVE_PD_xcorr_photodiode_combined_1K);
% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(EPHYS_PD_xcorr_time_1K ) / 500);
ind_edges = round(linspace(1, length(EPHYS_PD_xcorr_time_1K ), ind_edge_width));
ind_edges(1) = 0;
% init and loop over chunks
EPHYS_PD_inds_DTW  = cell((length(ind_edges)-1), 1);
BEHAVE_PD_inds_DTW = cell((length(ind_edges)-1), 1);
for counter_chunk = 1 : 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    EPHYS_PD_photodiode_combined_chunk  = EPHYS_PD_xcorr_photodiode_combined_1K_filt(inds_chunk);
    BEHAVE_PD_photodiode_combined_chunk = BEHAVE_PD_xcorr_photodiode_combined_1K_filt(inds_chunk);
    [~,ix,iy] = dtw(EPHYS_PD_photodiode_combined_chunk,BEHAVE_PD_photodiode_combined_chunk, 50, 'absolute'); % allow upto 50ms warp
    EPHYS_PD_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    BEHAVE_PD_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end
EPHYS_PD_inds_DTW  = cell2mat(EPHYS_PD_inds_DTW);
BEHAVE_PD_inds_DTW = cell2mat(BEHAVE_PD_inds_DTW);
EPHYS_PD_inds      = ( 1 : 1 : length(EPHYS_PD_xcorr_time_1K ) )';
BEHAVE_PD_inds     = ( 1 : 1 : length(BEHAVE_PD_xcorr_time_1K) )';
% dtw works by replicating the inds to match the two signals, here we
% reverse the replicated inds to generate two matched signals but with the
% size of original signals.
PD_ind_convert_from_EPHYS_to_BEHAVE = nan(size(EPHYS_PD_inds));
PD_ind_convert_from_BEHAVE_to_EPHYS = nan(size(BEHAVE_PD_inds));
for counter_ind = 1 : 1 : length(EPHYS_PD_inds_DTW)
    ind_EPHYS_PD_DTW  = EPHYS_PD_inds_DTW(counter_ind);
    ind_BEHAVE_PD_DTW = BEHAVE_PD_inds_DTW(counter_ind);
    PD_ind_convert_from_EPHYS_to_BEHAVE(ind_BEHAVE_PD_DTW) = ind_EPHYS_PD_DTW;
    PD_ind_convert_from_BEHAVE_to_EPHYS(ind_EPHYS_PD_DTW)  = ind_BEHAVE_PD_DTW;
end

time_reference      = EPHYS_time_30K(:);
length_time         = length(time_reference);
time_EPHYS_PD_xcorr_1K = EPHYS_PD_xcorr_time_1K;
time_EPHYS_PD_xcorr_1K(end+1) = max([time_reference(end), time_EPHYS_PD_xcorr_1K(end)])+1;
event_EPHYS_PD_xcorr_30K       = nan(length(EPHYS_PD_xcorr_time_1K), 1);
counter_EPHYS_PD_xcorr     = find(time_EPHYS_PD_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_PD_xcorr_1K(  counter_EPHYS_PD_xcorr)
        event_EPHYS_PD_xcorr_30K(    counter_EPHYS_PD_xcorr) = counter_time_point;
        counter_EPHYS_PD_xcorr   = counter_EPHYS_PD_xcorr   + 1;
    end
end

time_reference      = EPHYS_time_1K(:);
length_time         = length(time_reference);
time_EPHYS_PD_xcorr_1K = EPHYS_PD_xcorr_time_1K;
time_EPHYS_PD_xcorr_1K(end+1) = max([time_reference(end), time_EPHYS_PD_xcorr_1K(end)])+1;
event_EPHYS_PD_xcorr_1K       = nan(length(EPHYS_PD_xcorr_time_1K), 1);
counter_EPHYS_PD_xcorr     = find(time_EPHYS_PD_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_PD_xcorr_1K(  counter_EPHYS_PD_xcorr)
        event_EPHYS_PD_xcorr_1K(    counter_EPHYS_PD_xcorr) = counter_time_point;
        counter_EPHYS_PD_xcorr   = counter_EPHYS_PD_xcorr   + 1;
    end
end

time_reference      = BEHAVE_time_1K(:);
length_time         = length(time_reference);
time_BEHAVE_PD_xcorr_1K = BEHAVE_PD_xcorr_time_1K;
time_BEHAVE_PD_xcorr_1K(end+1) = max([time_reference(end), time_BEHAVE_PD_xcorr_1K(end)])+1;
event_BEHAVE_PD_xcorr_1K       = nan(length(BEHAVE_PD_xcorr_time_1K), 1);
counter_BEHAVE_PD_xcorr     = find(time_BEHAVE_PD_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_BEHAVE_PD_xcorr_1K(  counter_BEHAVE_PD_xcorr)
        event_BEHAVE_PD_xcorr_1K(    counter_BEHAVE_PD_xcorr) = counter_time_point;
        counter_BEHAVE_PD_xcorr   = counter_BEHAVE_PD_xcorr   + 1;
    end
end

EPHYS_PD_xcorr_ind_30K   = event_EPHYS_PD_xcorr_30K;
EPHYS_PD_xcorr_ind_1K    = event_EPHYS_PD_xcorr_1K;
BEHAVE_PD_xcorr_ind_1K   = event_BEHAVE_PD_xcorr_1K;
EPHYS_PD_aligned_ind_30K = EPHYS_PD_xcorr_ind_30K(PD_ind_convert_from_EPHYS_to_BEHAVE);
EPHYS_PD_aligned_ind_1K  = EPHYS_PD_xcorr_ind_1K( PD_ind_convert_from_EPHYS_to_BEHAVE);
BEHAVE_PD_aligned_ind_1K = BEHAVE_PD_xcorr_ind_1K(PD_ind_convert_from_BEHAVE_to_EPHYS);

EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_30K               = EPHYS_PD_aligned_ind_30K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_1K                = EPHYS_PD_aligned_ind_1K;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_aligned_ind_1K               = BEHAVE_PD_aligned_ind_1K;
EPHYS.CH_EVE.align_photodiode.PD_ind_convert_from_BEHAVE_to_EPHYS = PD_ind_convert_from_BEHAVE_to_EPHYS;
EPHYS.CH_EVE.align_photodiode.PD_ind_convert_from_EPHYS_to_BEHAVE = PD_ind_convert_from_EPHYS_to_BEHAVE;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_time_1K                 = EPHYS_PD_xcorr_time_1K;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_time_1K                = BEHAVE_PD_xcorr_time_1K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_ind_30K                 = EPHYS_PD_xcorr_ind_30K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_ind_1K                  = EPHYS_PD_xcorr_ind_1K;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_ind_1K                 = BEHAVE_PD_xcorr_ind_1K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_photodiode_combined_1K  = EPHYS_PD_xcorr_photodiode_combined_1K;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_photodiode_combined_1K = BEHAVE_PD_xcorr_photodiode_combined_1K;
EPHYS.CH_EVE.align_photodiode.photodiode_description              = BEHAVE.Alignment.photodiode_description;
fprintf(' --> Completed. \n');
    %% Save EPHYS EVENT DATA
clearvars -except EPHYS BEHAVE
EPHYS_time_1K    = EPHYS.Alignment.time_1K;
EPHYS_time_30K   = EPHYS.time_30K;
BEHAVE_time_1K   = BEHAVE.Alignment.time_1K;
align_photodiode = EPHYS.CH_EVE.align_photodiode;
align_states     = EPHYS.CH_EVE.align_states;

parts = strsplit(EPHYS.file_path_CH_EVE, ["/", "\"]);
name = erase(parts(length(parts) - 2), '-');
file_name = char(extractBetween(name, 3, 15));
file_path = EPHYS.file_path_CH_EVE;
file_name = [file_name '_EVE1_aligned.mat'];
fprintf([file_name ': Saving EPHYS Event Data ...'])
save([file_path '../analyzed_data/' file_name], '-v7.3');
fprintf(' --> Completed. \n')
%% 100
    %% Build EPHYS Alignment events
clearvars -except EPHYS BEHAVE
fprintf(['Building EPHYS Alignment events', ' ... ']);

time_reference      = ( ESN_Round(EPHYS.time_30K(1),0.001) : 0.01 : ESN_Round(EPHYS.time_30K(end),0.001) )';
length_time         = length(time_reference);
time_STR_TARGET_PURSUIT_rise = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 2) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_STR_TARGET_PURSUIT_fall = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 2) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
time_STR_TARGET_FIXATION_rise = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 3) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_STR_TARGET_FIXATION_fall = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 3) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
time_DETECT_SACCADE_END_rise = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_DETECT_SACCADE_END_fall = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
time_photodiode_rise      = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 6) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1) , 0.001);
time_photodiode_fall      = ESN_Round( EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 6) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1) , 0.001);
variable_list = {...
    '_STR_TARGET_PURSUIT_rise' ,'_STR_TARGET_PURSUIT_fall', ...
    '_STR_TARGET_FIXATION_rise','_STR_TARGET_FIXATION_fall', ...
    '_DETECT_SACCADE_END_rise' ,'_DETECT_SACCADE_END_fall', ...
    '_photodiode_rise','_photodiode_fall'};

if isempty(time_photodiode_rise)
    time_photodiode_rise = [time_STR_TARGET_FIXATION_rise; time_DETECT_SACCADE_END_rise];
    time_photodiode_rise = sort(time_photodiode_rise);
end

if isempty(time_photodiode_fall)
    time_photodiode_fall = [time_STR_TARGET_FIXATION_fall; time_DETECT_SACCADE_END_fall];
    time_photodiode_fall = sort(time_photodiode_fall);
end

EPHYS.Alignment.time_100 = time_reference;
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'EPHYS.Alignment.time' variable_name ' = ' 'time' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'time_temp_' ' = ' 'time' variable_name ';']);
    time_temp_(end+1) = max([time_reference(end), time_temp_(end)])+1;
    event_temp_       = false(length_time, 1);
    counter_temp_     = find(time_temp_ >= time_reference(1), 1, 'first');
    eval([ 'time'    variable_name ' = ' 'time_temp_'    ';']);
    eval([ 'event'   variable_name ' = ' 'event_temp_'   ';']);
    eval([ 'counter' variable_name ' = ' 'counter_temp_' ';']);
end
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    %{
    % Do not use 'eval', it is very slow, instead use the actual variables
    for counter_variable = 1 : 1 : length(variable_list)
        variable_name = variable_list{counter_variable};
        eval([ ...
            'if ( time_ponit_ >= time_' variable_name '(  counter_' variable_name ') );', ...
                'event_' variable_name '(    counter_time_point) = true;', ...
                'counter_' variable_name '   = counter_' variable_name '   + 1;', ...
            'end;', ...
        ]);
    end
    % Here is the template for actual variables
    if time_ponit_ >= time_VARIABLE(  counter_VARIABLE)
        event_VARIABLE(    counter_time_point) = true;
        counter_VARIABLE   = counter_VARIABLE   + 1;
    end
    %}
    % STR_TARGET_PURSUIT_rise
    if time_ponit_ >= time_STR_TARGET_PURSUIT_rise(  counter_STR_TARGET_PURSUIT_rise)
        event_STR_TARGET_PURSUIT_rise(    counter_time_point) = true;
        counter_STR_TARGET_PURSUIT_rise   = counter_STR_TARGET_PURSUIT_rise   + 1;
    end
    % STR_TARGET_PURSUIT_fall
    if time_ponit_ >= time_STR_TARGET_PURSUIT_fall(  counter_STR_TARGET_PURSUIT_fall)
        event_STR_TARGET_PURSUIT_fall(    counter_time_point) = true;
        counter_STR_TARGET_PURSUIT_fall   = counter_STR_TARGET_PURSUIT_fall   + 1;
    end
    % STR_TARGET_FIXATION_rise
    if time_ponit_ >= time_STR_TARGET_FIXATION_rise(  counter_STR_TARGET_FIXATION_rise)
        event_STR_TARGET_FIXATION_rise(    counter_time_point) = true;
        counter_STR_TARGET_FIXATION_rise   = counter_STR_TARGET_FIXATION_rise   + 1;
    end
    %STR_TARGET_FIXATION_fall
    if time_ponit_ >= time_STR_TARGET_FIXATION_fall(  counter_STR_TARGET_FIXATION_fall)
        event_STR_TARGET_FIXATION_fall(    counter_time_point) = true;
        counter_STR_TARGET_FIXATION_fall   = counter_STR_TARGET_FIXATION_fall   + 1;
    end
    % DETECT_SACCADE_END_rise
    if time_ponit_ >= time_DETECT_SACCADE_END_rise(  counter_DETECT_SACCADE_END_rise)
        event_DETECT_SACCADE_END_rise(    counter_time_point) = true;
        counter_DETECT_SACCADE_END_rise   = counter_DETECT_SACCADE_END_rise   + 1;
    end
    % DETECT_SACCADE_END_fall
    if time_ponit_ >= time_DETECT_SACCADE_END_fall(  counter_DETECT_SACCADE_END_fall)
        event_DETECT_SACCADE_END_fall(    counter_time_point) = true;
        counter_DETECT_SACCADE_END_fall   = counter_DETECT_SACCADE_END_fall   + 1;
    end
    % photodiode_rise
    if time_ponit_ >= time_photodiode_rise(     counter_photodiode_rise)
        event_photodiode_rise(       counter_time_point) = true;
        counter_photodiode_rise      = counter_photodiode_rise      + 1;
    end
    % photodiode_fall
    if time_ponit_ >= time_photodiode_fall(     counter_photodiode_fall)
        event_photodiode_fall(       counter_time_point) = true;
        counter_photodiode_fall      = counter_photodiode_fall      + 1;
    end
end

event_STR_TARGET_PURSUIT  = false(length_time, 1);
event_STR_TARGET_FIXATION = false(length_time, 1);
event_DETECT_SACCADE_END  = false(length_time, 1);
event_photodiode          = false(length_time, 1);
flag_STR_TARGET_PURSUIT   = false;
flag_STR_TARGET_FIXATION  = false;
flag_DETECT_SACCADE_END   = false;
flag_photodiode           = false;
for counter_time_point = 1 : length_time
    flag_STR_TARGET_PURSUIT  = flag_STR_TARGET_PURSUIT  ||   event_STR_TARGET_PURSUIT_rise( counter_time_point);
    flag_STR_TARGET_PURSUIT  = flag_STR_TARGET_PURSUIT  && (~event_STR_TARGET_PURSUIT_fall( counter_time_point));
    flag_STR_TARGET_FIXATION = flag_STR_TARGET_FIXATION ||   event_STR_TARGET_FIXATION_rise(counter_time_point);
    flag_STR_TARGET_FIXATION = flag_STR_TARGET_FIXATION && (~event_STR_TARGET_FIXATION_fall(counter_time_point));
    flag_DETECT_SACCADE_END  = flag_DETECT_SACCADE_END  ||   event_DETECT_SACCADE_END_rise( counter_time_point);
    flag_DETECT_SACCADE_END  = flag_DETECT_SACCADE_END  && (~event_DETECT_SACCADE_END_fall( counter_time_point));
    flag_photodiode          = flag_photodiode          ||   event_photodiode_rise(         counter_time_point);
    flag_photodiode          = flag_photodiode          && (~event_photodiode_fall(         counter_time_point));
    event_STR_TARGET_PURSUIT( counter_time_point) = flag_STR_TARGET_PURSUIT;
    event_STR_TARGET_FIXATION(counter_time_point) = flag_STR_TARGET_FIXATION;
    event_DETECT_SACCADE_END( counter_time_point) = flag_DETECT_SACCADE_END;
    event_photodiode(         counter_time_point) = flag_photodiode;
end

variable_list = {...
    '_STR_TARGET_PURSUIT_rise' ,'_STR_TARGET_PURSUIT_fall', '_STR_TARGET_PURSUIT', ...
    '_STR_TARGET_FIXATION_rise','_STR_TARGET_FIXATION_fall', '_STR_TARGET_FIXATION', ...
    '_DETECT_SACCADE_END_rise' ,'_DETECT_SACCADE_END_fall', '_DETECT_SACCADE_END', ...
    '_photodiode_rise','_photodiode_fall', '_photodiode'};
for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'EPHYS.Alignment.event' variable_name ' = ' 'event' variable_name ';']);
end

event_state_combined = ...
    double(event_STR_TARGET_PURSUIT)  .* 1 + ...
    double(event_STR_TARGET_FIXATION) .* 3 + ...
    double(event_DETECT_SACCADE_END)  .* 7 ;

event_photodiode_combined = ...
    double(event_photodiode)  .* 1;

EPHYS.Alignment.event_state_combined      = event_state_combined;
EPHYS.Alignment.event_photodiode_combined = event_photodiode_combined;

state_description = [
'STR_TARGET_PURSUIT: 1 , ', ...
'STR_TARGET_FIXATION: 3 , ', ...
'DETECT_SACCADE_END: 7 , '];

photodiode_description = [
'STR_TARGET_FIXATION: 1 , ', ...
'DETECT_SACCADE_END: 1 , ', ...
];
EPHYS.Alignment.state_description      = state_description;
EPHYS.Alignment.photodiode_description = photodiode_description;

fprintf(' --> Completed. \n');
    %% Build BEHAVE Alignment events
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE Alignment events', ' ... ']);
time_reference      = BEHAVE.time_100;
BEHAVE.Alignment.time_100 = time_reference;
BEHAVE.Alignment.time_1K = BEHAVE.time_1K;
BEHAVE_state = BEHAVE.state_100(:);
event_STR_TARGET_PURSUIT  = (BEHAVE_state == 1);
event_STR_TARGET_FIXATION = (BEHAVE_state == 3);
event_DETECT_SACCADE_END =  (BEHAVE_state == 7);
event_photodiode = event_STR_TARGET_FIXATION | event_DETECT_SACCADE_END;

event_state_combined = ...
    double(event_STR_TARGET_PURSUIT)  .* 1 + ...
    double(event_STR_TARGET_FIXATION) .* 3 + ...
    double(event_DETECT_SACCADE_END)  .* 7 ;

event_photodiode_combined = ...
    double(event_photodiode)  .* 1;

BEHAVE.Alignment.event_state_combined      = event_state_combined;
BEHAVE.Alignment.event_photodiode_combined = event_photodiode_combined;

state_description = [
'STR_TARGET_PURSUIT: 1 , ', ...
'STR_TARGET_FIXATION: 3 , ', ...
'DETECT_SACCADE_END: 7 , '];

photodiode_description = [
'STR_TARGET_FIXATION: 1 , ', ...
'DETECT_SACCADE_END: 1 , ', ...
];
BEHAVE.Alignment.state_description      = state_description;
BEHAVE.Alignment.photodiode_description = photodiode_description;

fprintf(' --> Completed. \n');
    %% ALIGN EPHYS and BEHAVE state_combined through xcorr and dtw
clearvars -except EPHYS BEHAVE
fprintf(['Aligning EPHYS and BEHAVE state_combined', ' ... ']);
EPHYS_time_100              = EPHYS.Alignment.time_100;
EPHYS_time_30K             = EPHYS.time_30K;
BEHAVE_time_100             = BEHAVE.Alignment.time_100;
EPHYS_state_combined       = EPHYS.Alignment.event_state_combined;
BEHAVE_state_combined      = BEHAVE.Alignment.event_state_combined;

% state_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(EPHYS_state_combined, BEHAVE_state_combined); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);

if  sample_diff > 0
    EPHYS_EB_xcorr_state_combined_100  = EPHYS_state_combined(  abs(sample_diff):end);
    EPHYS_EB_xcorr_time_100            = EPHYS_time_100(            abs(sample_diff):end);
    BEHAVE_EB_xcorr_state_combined_100 = BEHAVE_state_combined;
    BEHAVE_EB_xcorr_time_100           = BEHAVE_time_100;
elseif sample_diff < 0
    EPHYS_EB_xcorr_state_combined_100  = EPHYS_state_combined;
    EPHYS_EB_xcorr_time_100            = EPHYS_time_100;
    BEHAVE_EB_xcorr_state_combined_100 = BEHAVE_state_combined( abs(sample_diff):end);
    BEHAVE_EB_xcorr_time_100           = BEHAVE_time_100(           abs(sample_diff):end);
end
% state_combined: make the vectors the same size
if length(BEHAVE_EB_xcorr_state_combined_100) ~= length(EPHYS_EB_xcorr_state_combined_100)
    min_length = min([ length(BEHAVE_EB_xcorr_state_combined_100),  length(EPHYS_EB_xcorr_state_combined_100)]);
    EPHYS_EB_xcorr_state_combined_100  = EPHYS_EB_xcorr_state_combined_100(  1:min_length);
    EPHYS_EB_xcorr_time_100            = EPHYS_EB_xcorr_time_100(            1:min_length);
    BEHAVE_EB_xcorr_state_combined_100 = BEHAVE_EB_xcorr_state_combined_100( 1:min_length);
    BEHAVE_EB_xcorr_time_100           = BEHAVE_EB_xcorr_time_100(           1:min_length);
end

% %% state_combined: find the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
sampling_freq = 100.0;
cutoff_freq = 25.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
EPHYS_EB_xcorr_state_combined_100_filt  = filtfilt(b_butter,a_butter,EPHYS_EB_xcorr_state_combined_100);
BEHAVE_EB_xcorr_state_combined_100_filt = filtfilt(b_butter,a_butter,BEHAVE_EB_xcorr_state_combined_100);
% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(EPHYS_EB_xcorr_time_100 ) / 500);
ind_edges = round(linspace(1, length(EPHYS_EB_xcorr_time_100 ), ind_edge_width));
ind_edges(1) = 0;
% init and loop over chunks
EPHYS_EB_inds_DTW  = cell((length(ind_edges)-1), 1);
BEHAVE_EB_inds_DTW = cell((length(ind_edges)-1), 1);
for counter_chunk = 1 : 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    EPHYS_EB_state_combined_chunk  = EPHYS_EB_xcorr_state_combined_100_filt(inds_chunk);
    BEHAVE_EB_state_combined_chunk = BEHAVE_EB_xcorr_state_combined_100_filt(inds_chunk);
    [~,ix,iy] = dtw(EPHYS_EB_state_combined_chunk,BEHAVE_EB_state_combined_chunk, 15, 'absolute');  % allow upto 15ms warp
    EPHYS_EB_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    BEHAVE_EB_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end
EPHYS_EB_inds_DTW  = cell2mat(EPHYS_EB_inds_DTW);
BEHAVE_EB_inds_DTW = cell2mat(BEHAVE_EB_inds_DTW);
EPHYS_EB_inds      = ( 1 : 1 : length(EPHYS_EB_xcorr_time_100 ) )';
BEHAVE_EB_inds     = ( 1 : 1 : length(BEHAVE_EB_xcorr_time_100) )';
% dtw works by replicating the inds to match the two signals, here we
% reverse the replicated inds to generate two matched signals but with the
% size of original signals.
EB_ind_convert_from_EPHYS_to_BEHAVE = nan(size(EPHYS_EB_inds));
EB_ind_convert_from_BEHAVE_to_EPHYS = nan(size(BEHAVE_EB_inds));
for counter_ind = 1 : 1 : length(EPHYS_EB_inds_DTW)
    ind_EPHYS_EB_DTW  = EPHYS_EB_inds_DTW(counter_ind);
    ind_BEHAVE_EB_DTW = BEHAVE_EB_inds_DTW(counter_ind);
    EB_ind_convert_from_EPHYS_to_BEHAVE(ind_BEHAVE_EB_DTW) = ind_EPHYS_EB_DTW;
    EB_ind_convert_from_BEHAVE_to_EPHYS(ind_EPHYS_EB_DTW)  = ind_BEHAVE_EB_DTW;
end

time_reference      = EPHYS_time_30K(:);
length_time         = length(time_reference);
time_EPHYS_EB_xcorr_100 = EPHYS_EB_xcorr_time_100;
time_EPHYS_EB_xcorr_100(end+1) = max([time_reference(end), time_EPHYS_EB_xcorr_100(end)])+1;
event_EPHYS_EB_xcorr_30K       = nan(length(EPHYS_EB_xcorr_time_100), 1);
counter_EPHYS_EB_xcorr     = find(time_EPHYS_EB_xcorr_100 >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_EB_xcorr_100(  counter_EPHYS_EB_xcorr)
        event_EPHYS_EB_xcorr_30K(    counter_EPHYS_EB_xcorr) = counter_time_point;
        counter_EPHYS_EB_xcorr   = counter_EPHYS_EB_xcorr   + 1;
    end
end

time_reference      = EPHYS_time_100(:);
length_time         = length(time_reference);
time_EPHYS_EB_xcorr_100 = EPHYS_EB_xcorr_time_100;
time_EPHYS_EB_xcorr_100(end+1) = max([time_reference(end), time_EPHYS_EB_xcorr_100(end)])+1;
event_EPHYS_EB_xcorr_100       = nan(length(EPHYS_EB_xcorr_time_100), 1);
counter_EPHYS_EB_xcorr     = find(time_EPHYS_EB_xcorr_100 >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_EB_xcorr_100(  counter_EPHYS_EB_xcorr)
        event_EPHYS_EB_xcorr_100(    counter_EPHYS_EB_xcorr) = counter_time_point;
        counter_EPHYS_EB_xcorr   = counter_EPHYS_EB_xcorr   + 1;
    end
end

time_reference      = BEHAVE_time_100(:);
length_time         = length(time_reference);
time_BEHAVE_EB_xcorr_100 = BEHAVE_EB_xcorr_time_100;
time_BEHAVE_EB_xcorr_100(end+1) = max([time_reference(end), time_BEHAVE_EB_xcorr_100(end)])+1;
event_BEHAVE_EB_xcorr_100       = nan(length(BEHAVE_EB_xcorr_time_100), 1);
counter_BEHAVE_EB_xcorr     = find(time_BEHAVE_EB_xcorr_100 >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_BEHAVE_EB_xcorr_100(  counter_BEHAVE_EB_xcorr)
        event_BEHAVE_EB_xcorr_100(    counter_BEHAVE_EB_xcorr) = counter_time_point;
        counter_BEHAVE_EB_xcorr   = counter_BEHAVE_EB_xcorr   + 1;
    end
end

EPHYS_EB_xcorr_ind_30K   = event_EPHYS_EB_xcorr_30K;
EPHYS_EB_xcorr_ind_100    = event_EPHYS_EB_xcorr_100;
BEHAVE_EB_xcorr_ind_100   = event_BEHAVE_EB_xcorr_100;
EPHYS_EB_aligned_ind_30K = EPHYS_EB_xcorr_ind_30K(EB_ind_convert_from_EPHYS_to_BEHAVE);
EPHYS_EB_aligned_ind_100  = EPHYS_EB_xcorr_ind_100( EB_ind_convert_from_EPHYS_to_BEHAVE);
BEHAVE_EB_aligned_ind_100 = BEHAVE_EB_xcorr_ind_100(EB_ind_convert_from_BEHAVE_to_EPHYS);

EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_30K          = EPHYS_EB_aligned_ind_30K;
EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_100           = EPHYS_EB_aligned_ind_100;
EPHYS.CH_EVE.align_states.BEHAVE_EB_aligned_ind_100          = BEHAVE_EB_aligned_ind_100;
EPHYS.CH_EVE.align_states.EB_ind_convert_from_BEHAVE_to_EPHYS = EB_ind_convert_from_BEHAVE_to_EPHYS;
EPHYS.CH_EVE.align_states.EB_ind_convert_from_EPHYS_to_BEHAVE = EB_ind_convert_from_EPHYS_to_BEHAVE;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_time_100            = EPHYS_EB_xcorr_time_100;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_time_100           = BEHAVE_EB_xcorr_time_100;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_ind_30K            = EPHYS_EB_xcorr_ind_30K;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_ind_100             = EPHYS_EB_xcorr_ind_100;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_100            = BEHAVE_EB_xcorr_ind_100;
EPHYS.CH_EVE.align_states.EPHYS_EB_xcorr_state_combined_100  = EPHYS_EB_xcorr_state_combined_100;
EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_state_combined_100 = BEHAVE_EB_xcorr_state_combined_100;
% EPHYS.CH_EVE.align_states.state_description            = BEHAVE.Alignment.state_description;
fprintf(' --> Completed. \n');

if EPHYS.debug_figures
clf(figure(1));subplot(2,1,1);plot(EPHYS_state_combined);subplot(2,1,2);plot(BEHAVE_state_combined);
clf(figure(2));subplot(2,1,1);plot(EPHYS_EB_xcorr_state_combined_100);subplot(2,1,2);plot(BEHAVE_EB_xcorr_state_combined_100);
clf(figure(3));subplot(2,1,1);plot(EPHYS_EB_xcorr_state_combined_100);subplot(2,1,2);plot(BEHAVE_EB_xcorr_state_combined_100(EB_ind_convert_from_BEHAVE_to_EPHYS));
end
disp(['length EPHYS: ' num2str( EPHYS.Alignment.time_100(end)-EPHYS.Alignment.time_100(1) )])
disp(['length BEHAVE: ' num2str( BEHAVE.Alignment.time_100(end)-BEHAVE.Alignment.time_100(1) )])
disp(['xcorr diff: ' num2str( sample_diff )])
    %% ALIGN EPHYS and BEHAVE photodiode_combined through xcorr and dtw
clearvars -except EPHYS BEHAVE
fprintf(['Aligning EPHYS and BEHAVE photodiode_combined', ' ... ']);
EPHYS_time_100              = EPHYS.Alignment.time_100;
EPHYS_time_30K             = EPHYS.time_30K;
BEHAVE_time_100             = BEHAVE.Alignment.time_100;
EPHYS_photodiode_combined  = EPHYS.Alignment.event_photodiode_combined;
BEHAVE_photodiode_combined = BEHAVE.Alignment.event_photodiode_combined;

% photodiode_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(EPHYS_photodiode_combined+1, BEHAVE_photodiode_combined+1); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);
if sample_diff > 0
    EPHYS_PD_xcorr_photodiode_combined_100  = EPHYS_photodiode_combined(  abs(sample_diff):end);
    EPHYS_PD_xcorr_time_100                 = EPHYS_time_100(                 abs(sample_diff):end);
    BEHAVE_PD_xcorr_photodiode_combined_100 = BEHAVE_photodiode_combined;
    BEHAVE_PD_xcorr_time_100                = BEHAVE_time_100;
elseif sample_diff < 0
    EPHYS_PD_xcorr_photodiode_combined_100  = EPHYS_photodiode_combined;
    EPHYS_PD_xcorr_time_100                 = EPHYS_time_100;
    BEHAVE_PD_xcorr_photodiode_combined_100 = BEHAVE_photodiode_combined( abs(sample_diff):end);
    BEHAVE_PD_xcorr_time_100                = BEHAVE_time_100(                abs(sample_diff):end);
end
% photodiode_combined: make the vectors the same size
if length(BEHAVE_PD_xcorr_photodiode_combined_100) ~= length(EPHYS_PD_xcorr_photodiode_combined_100)
    min_length = min([ length(BEHAVE_PD_xcorr_photodiode_combined_100),  length(EPHYS_PD_xcorr_photodiode_combined_100)]);
    EPHYS_PD_xcorr_photodiode_combined_100  = EPHYS_PD_xcorr_photodiode_combined_100(  1:min_length);
    EPHYS_PD_xcorr_time_100                 = EPHYS_PD_xcorr_time_100(                 1:min_length);
    BEHAVE_PD_xcorr_photodiode_combined_100 = BEHAVE_PD_xcorr_photodiode_combined_100( 1:min_length);
    BEHAVE_PD_xcorr_time_100                = BEHAVE_PD_xcorr_time_100(                1:min_length);
end

% photodiode_combined: find the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
sampling_freq = 100.0;
cutoff_freq = 25.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
EPHYS_PD_xcorr_photodiode_combined_100_filt  = filtfilt(b_butter,a_butter,EPHYS_PD_xcorr_photodiode_combined_100);
BEHAVE_PD_xcorr_photodiode_combined_100_filt = filtfilt(b_butter,a_butter,BEHAVE_PD_xcorr_photodiode_combined_100);
% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(EPHYS_PD_xcorr_time_100 ) / 500);
ind_edges = round(linspace(1, length(EPHYS_PD_xcorr_time_100 ), ind_edge_width));
ind_edges(1) = 0;
% init and loop over chunks
EPHYS_PD_inds_DTW  = cell((length(ind_edges)-1), 1);
BEHAVE_PD_inds_DTW = cell((length(ind_edges)-1), 1);
for counter_chunk = 1 : 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    EPHYS_PD_photodiode_combined_chunk  = EPHYS_PD_xcorr_photodiode_combined_100_filt(inds_chunk);
    BEHAVE_PD_photodiode_combined_chunk = BEHAVE_PD_xcorr_photodiode_combined_100_filt(inds_chunk);
    [~,ix,iy] = dtw(EPHYS_PD_photodiode_combined_chunk,BEHAVE_PD_photodiode_combined_chunk, 50, 'absolute'); % allow upto 50ms warp
    EPHYS_PD_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    BEHAVE_PD_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end
EPHYS_PD_inds_DTW  = cell2mat(EPHYS_PD_inds_DTW);
BEHAVE_PD_inds_DTW = cell2mat(BEHAVE_PD_inds_DTW);
EPHYS_PD_inds      = ( 1 : 1 : length(EPHYS_PD_xcorr_time_100 ) )';
BEHAVE_PD_inds     = ( 1 : 1 : length(BEHAVE_PD_xcorr_time_100) )';
% dtw works by replicating the inds to match the two signals, here we
% reverse the replicated inds to generate two matched signals but with the
% size of original signals.
PD_ind_convert_from_EPHYS_to_BEHAVE = nan(size(EPHYS_PD_inds));
PD_ind_convert_from_BEHAVE_to_EPHYS = nan(size(BEHAVE_PD_inds));
for counter_ind = 1 : 1 : length(EPHYS_PD_inds_DTW)
    ind_EPHYS_PD_DTW  = EPHYS_PD_inds_DTW(counter_ind);
    ind_BEHAVE_PD_DTW = BEHAVE_PD_inds_DTW(counter_ind);
    PD_ind_convert_from_EPHYS_to_BEHAVE(ind_BEHAVE_PD_DTW) = ind_EPHYS_PD_DTW;
    PD_ind_convert_from_BEHAVE_to_EPHYS(ind_EPHYS_PD_DTW)  = ind_BEHAVE_PD_DTW;
end

time_reference      = EPHYS_time_30K(:);
length_time         = length(time_reference);
time_EPHYS_PD_xcorr_100 = EPHYS_PD_xcorr_time_100;
time_EPHYS_PD_xcorr_100(end+1) = max([time_reference(end), time_EPHYS_PD_xcorr_100(end)])+1;
event_EPHYS_PD_xcorr_30K       = nan(length(EPHYS_PD_xcorr_time_100), 1);
counter_EPHYS_PD_xcorr     = find(time_EPHYS_PD_xcorr_100 >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_PD_xcorr_100(  counter_EPHYS_PD_xcorr)
        event_EPHYS_PD_xcorr_30K(    counter_EPHYS_PD_xcorr) = counter_time_point;
        counter_EPHYS_PD_xcorr   = counter_EPHYS_PD_xcorr   + 1;
    end
end

time_reference      = EPHYS_time_100(:);
length_time         = length(time_reference);
time_EPHYS_PD_xcorr_100 = EPHYS_PD_xcorr_time_100;
time_EPHYS_PD_xcorr_100(end+1) = max([time_reference(end), time_EPHYS_PD_xcorr_100(end)])+1;
event_EPHYS_PD_xcorr_100       = nan(length(EPHYS_PD_xcorr_time_100), 1);
counter_EPHYS_PD_xcorr     = find(time_EPHYS_PD_xcorr_100 >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_EPHYS_PD_xcorr_100(  counter_EPHYS_PD_xcorr)
        event_EPHYS_PD_xcorr_100(    counter_EPHYS_PD_xcorr) = counter_time_point;
        counter_EPHYS_PD_xcorr   = counter_EPHYS_PD_xcorr   + 1;
    end
end

time_reference      = BEHAVE_time_100(:);
length_time         = length(time_reference);
time_BEHAVE_PD_xcorr_100 = BEHAVE_PD_xcorr_time_100;
time_BEHAVE_PD_xcorr_100(end+1) = max([time_reference(end), time_BEHAVE_PD_xcorr_100(end)])+1;
event_BEHAVE_PD_xcorr_100       = nan(length(BEHAVE_PD_xcorr_time_100), 1);
counter_BEHAVE_PD_xcorr     = find(time_BEHAVE_PD_xcorr_100 >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_ponit_     = time_reference(counter_time_point);
    if time_ponit_ >= time_BEHAVE_PD_xcorr_100(  counter_BEHAVE_PD_xcorr)
        event_BEHAVE_PD_xcorr_100(    counter_BEHAVE_PD_xcorr) = counter_time_point;
        counter_BEHAVE_PD_xcorr   = counter_BEHAVE_PD_xcorr   + 1;
    end
end

EPHYS_PD_xcorr_ind_30K   = event_EPHYS_PD_xcorr_30K;
EPHYS_PD_xcorr_ind_100    = event_EPHYS_PD_xcorr_100;
BEHAVE_PD_xcorr_ind_100   = event_BEHAVE_PD_xcorr_100;
EPHYS_PD_aligned_ind_30K = EPHYS_PD_xcorr_ind_30K(PD_ind_convert_from_EPHYS_to_BEHAVE);
EPHYS_PD_aligned_ind_100  = EPHYS_PD_xcorr_ind_100( PD_ind_convert_from_EPHYS_to_BEHAVE);
BEHAVE_PD_aligned_ind_100 = BEHAVE_PD_xcorr_ind_100(PD_ind_convert_from_BEHAVE_to_EPHYS);

EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_30K               = EPHYS_PD_aligned_ind_30K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_100                = EPHYS_PD_aligned_ind_100;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_aligned_ind_100               = BEHAVE_PD_aligned_ind_100;
EPHYS.CH_EVE.align_photodiode.PD_ind_convert_from_BEHAVE_to_EPHYS = PD_ind_convert_from_BEHAVE_to_EPHYS;
EPHYS.CH_EVE.align_photodiode.PD_ind_convert_from_EPHYS_to_BEHAVE = PD_ind_convert_from_EPHYS_to_BEHAVE;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_time_100                 = EPHYS_PD_xcorr_time_100;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_time_100                = BEHAVE_PD_xcorr_time_100;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_ind_30K                 = EPHYS_PD_xcorr_ind_30K;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_ind_100                  = EPHYS_PD_xcorr_ind_100;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_ind_100                 = BEHAVE_PD_xcorr_ind_100;
EPHYS.CH_EVE.align_photodiode.EPHYS_PD_xcorr_photodiode_combined_100  = EPHYS_PD_xcorr_photodiode_combined_100;
EPHYS.CH_EVE.align_photodiode.BEHAVE_PD_xcorr_photodiode_combined_100 = BEHAVE_PD_xcorr_photodiode_combined_100;
EPHYS.CH_EVE.align_photodiode.photodiode_description              = BEHAVE.Alignment.photodiode_description;
fprintf(' --> Completed. \n');
    %% Save EPHYS EVENT DATA
clearvars -except EPHYS BEHAVE
EPHYS_time_100    = EPHYS.Alignment.time_100;
EPHYS_time_30K   = EPHYS.time_30K;
BEHAVE_time_1K   = BEHAVE.Alignment.time_1K;
BEHAVE_time_100   = BEHAVE.Alignment.time_100;
align_photodiode = EPHYS.CH_EVE.align_photodiode;
align_states     = EPHYS.CH_EVE.align_states;

parts = strsplit(EPHYS.file_path_CH_EVE, ["/", "\"]);
name = erase(parts(length(parts) - 2), '-');
file_name = char(extractBetween(name, 3, 15));
file_path = EPHYS.file_path_CH_EVE;
file_name = [file_name '_EVE1_ESN_aligned.mat'];
clearvars EPHYS BEHAVE
fprintf([file_name ': Saving EPHYS Event Data ...'])
save([file_path '../analyzed_data/' file_name], '-v7.3');
fprintf(' --> Completed. \n')
end