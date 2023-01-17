function PGH_event_alignment_LED_ephys(path_to_raw, path_to_save, rec_name)
%% load EPHYS EVENT DATA
file_name = 'all_channels.events';
fprintf(['Loading: ', file_name, ' ... ']);
if ~strcmp(path_to_raw(end), filesep);path_to_raw = [path_to_raw filesep];end
[ch_data, ch_time, ch_info] = load_open_ephys_data([path_to_raw file_name]);
EPHYS.CH_EVE.ch_data = ch_data;
EPHYS.CH_EVE.ch_time = ch_time;
EPHYS.CH_EVE.ch_info = ch_info;
EPHYS.file_name_CH_EVE = file_name;
EPHYS.file_path_CH_EVE = path_to_raw;
EPHYS.path_to_save = path_to_save;
EPHYS.rec_name = rec_name;

fprintf(' --> Completed. \n')

EPHYS.debug_figures = false;

%% load continuos file
file_name_ = dir([path_to_raw '*_CH*.continuous']);
file_name = file_name_(1).name;
fprintf(['Loading: ', file_name, ' ... \n']);
if ~strcmp(path_to_raw(end), filesep);path_to_raw = [path_to_raw filesep];end
[~, ch_time, ~] = load_open_ephys_data([path_to_raw file_name]);
EPHYS.ch_time = ch_time;
EPHYS.time_30k   = double(EPHYS.ch_time);
EPHYS.time_1K      = (EPHYS.ch_time(1) : 0.001 :EPHYS.ch_time(end))';
% EPHYS.time_100      = (EPHYS.ch_time(1) : 0.01 :EPHYS.ch_time(end))';

%% LED Data
clearvars -except EPHYS VID path_to_raw
fprintf('Loading LED signal ...');

path_to_tongue = [path_to_raw '..' filesep 'analyzed_data' ...
    filesep 'behavior_data' filesep 'tongue' filesep];

vid_file_name = dir([path_to_tongue '*_video.mat']);
load([path_to_tongue vid_file_name(1).name],'FPS','LED_FPS');

VID.FPS = FPS;
VID.LED_FPS = LED_FPS;

fprintf(' --> Completed. \n')

%% Alignment
fprintf('Aligning data ...')
clearvars -except EPHYS BEHAVE VID
%% extract event data in ephys
clearvars -except EPHYS BEHAVE VID
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

%% Build EPHYS Alignment events for 1K
clearvars -except EPHYS BEHAVE VID
fprintf(['Building EPHYS Alignment events for VID', ' ... ']);

time_reference      = EPHYS.time_1K;
length_time         = length(time_reference);
time_state_str_fixation   = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 3) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1);
%time_state_sac_detect_off = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1);

% variable_list = {'_state_str_fixation','_state_sac_detect_off'};
variable_list = {'_state_str_fixation'};

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
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_state_str_fixation(  counter_state_str_fixation)
        event_state_str_fixation(    counter_time_point) = true;
        counter_state_str_fixation   = counter_state_str_fixation   + 1;
    end
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'EPHYS.Alignment.event' variable_name ' = ' 'event' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'event_temp_' ' = ' 'event' variable_name ';']);
    length_time = length(event_temp_);

    inds_span_ = (0) : 1 : (200);

    ind_event_temp_  = find(event_temp_);
    inds_event_temp_ = repmat( ind_event_temp_(:), 1, length(inds_span_)) + repmat(inds_span_(:)', length(ind_event_temp_), 1);
    inds_event_temp_( inds_event_temp_ < 1 ) = 1;
    inds_event_temp_( inds_event_temp_ > length_time ) = length_time;
    event_temp_(inds_event_temp_(:)) = true;
    eval([ 'event'   variable_name  ' = ' 'event_temp_'   ';']);
end
event_LED_combined = double(event_state_str_fixation)   .* 1  ;
EPHYS.Alignment.event_LED_combined_1K = event_LED_combined;
fprintf(' --> Completed. \n');

%% Build VID Alignment events - 1K
clearvars -except EPHYS BEHAVE VID
fprintf(['Building VID Alignment events', ' ... ']);

time_vid(:,1) = ((1/VID.FPS) : (1/VID.FPS) : length(VID.LED_FPS(:,4))/VID.FPS);

% interpolation from FPS to 1K
time_1K(:,1) = time_vid(1) : 0.001 : time_vid(end);
for counter_fields = 1 : size(VID.LED_FPS,2)
    VID.LED_1K(:,counter_fields) = interp1(time_vid,VID.LED_FPS(:,counter_fields),time_1K);
end

time_reference = time_1K;
length_time = length(time_reference);
time_LED_rise = time_reference([diff(VID.LED_1K(:,4)) ; 0] > 0);

variable_list = { '_LED_rise'};

VID.Alignment.time_vid = time_vid;
VID.Alignment.time_1K= time_reference;

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'VID.Alignment.time' variable_name ' = ' 'time' variable_name ';']);
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
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_LED_rise(     counter_LED_rise)
        event_LED_rise(       counter_time_point) = true;
        counter_LED_rise      = counter_LED_rise      + 1;
    end

end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'VID.Alignment.event' variable_name ' = ' 'event' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'event_temp_' ' = ' 'event' variable_name ';']);
    length_time = length(event_temp_);

    inds_span_ = (0) : 1 : (200);

    ind_event_temp_  = find(event_temp_);
    inds_event_temp_ = repmat( ind_event_temp_(:), 1, length(inds_span_)) + repmat(inds_span_(:)', length(ind_event_temp_), 1);
    inds_event_temp_( inds_event_temp_ < 1 ) = 1;
    inds_event_temp_( inds_event_temp_ > length_time ) = length_time;
    event_temp_(inds_event_temp_(:)) = true;
    eval([ 'event'   variable_name  ' = ' 'event_temp_'   ';']);
end


event_LED_combined = double(event_LED_rise)  .* 1;

VID.Alignment.event_LED_combined_1K = event_LED_combined;


fprintf(' --> Completed. \n');


%% ALIGN EPHYS and VID through xcorr and dtw - time_1K
clearvars -except EPHYS BEHAVE VID
fprintf(['Aligning VID and EPHYS LED_combined', ' ... ']);

EPHYS_time_1K            = EPHYS.time_1K;
EPHYS_time_30k             = EPHYS.time_30k;
VID_time_1K = VID.Alignment.time_1K;

VID_LED_combined_1K = VID.Alignment.event_LED_combined_1K;
EPHYS_LED_combined_1K = EPHYS.Alignment.event_LED_combined_1K;

% LED_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(VID_LED_combined_1K+1, EPHYS_LED_combined_1K+1); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);

if sample_diff > 0
    VID_LED_xcorr_combined_1K  = VID_LED_combined_1K(abs(sample_diff):end);
    VID_LED_xcorr_time_1K                 = VID_time_1K(abs(sample_diff):end);
    EPHYS_LED_xcorr_combined_1K = EPHYS_LED_combined_1K;
    EPHYS_LED_xcorr_time_1K                = EPHYS_time_1K;
elseif sample_diff < 0
    VID_LED_xcorr_combined_1K  = VID_LED_combined_1K;
    VID_LED_xcorr_time_1K                 = VID_time_1K;
    EPHYS_LED_xcorr_combined_1K = EPHYS_LED_combined_1K(abs(sample_diff):end);
    EPHYS_LED_xcorr_time_1K                = EPHYS_time_1K(abs(sample_diff):end);
end

% LED_combined: make the vectors the same size
if length(EPHYS_LED_xcorr_combined_1K) ~= length(VID_LED_xcorr_combined_1K)
    min_length = min([ length(EPHYS_LED_xcorr_combined_1K),  length(VID_LED_xcorr_combined_1K)]);
    VID_LED_xcorr_combined_1K  = VID_LED_xcorr_combined_1K(  1:min_length);
    VID_LED_xcorr_time_1K                 = VID_LED_xcorr_time_1K(                 1:min_length);
    EPHYS_LED_xcorr_combined_1K = EPHYS_LED_xcorr_combined_1K( 1:min_length);
    EPHYS_LED_xcorr_time_1K                = EPHYS_LED_xcorr_time_1K(                1:min_length);
end


% LED_combined: find the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
sampling_freq = 1000;
cutoff_freq = 25.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
VID_LED_xcorr_combined_1K_filt  = filtfilt(b_butter,a_butter,VID_LED_xcorr_combined_1K);
EPHYS_LED_xcorr_combined_1K_filt = filtfilt(b_butter,a_butter,EPHYS_LED_xcorr_combined_1K);
% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(VID_LED_xcorr_time_1K ) / 500);
ind_edges = round(linspace(1, length(VID_LED_xcorr_time_1K ), ind_edge_width));
ind_edges(1) = 0;

% init and loop over chunks
VID_LED_inds_DTW  = cell((length(ind_edges)-1), 1);
EPHYS_LED_inds_DTW = cell((length(ind_edges)-1), 1);
for counter_chunk = 1 : 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    VID_LED_LED_combined_chunk  = VID_LED_xcorr_combined_1K_filt(inds_chunk);
    EPHYS_LED_LED_combined_chunk = EPHYS_LED_xcorr_combined_1K_filt(inds_chunk);
    [~,ix,iy] = dtw(VID_LED_LED_combined_chunk,EPHYS_LED_LED_combined_chunk, 15, 'absolute'); % allow upto 15ms warp
    VID_LED_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    EPHYS_LED_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end

VID_LED_inds_DTW  = cell2mat(VID_LED_inds_DTW);
EPHYS_LED_inds_DTW = cell2mat(EPHYS_LED_inds_DTW);
VID_LED_inds      = ( 1 : 1 : length(VID_LED_xcorr_time_1K ) )';
EPHYS_LED_inds     = ( 1 : 1 : length(EPHYS_LED_xcorr_time_1K) )';


% dtw works by replicating the inds to match the two signals, here we
% reverse the replicated inds to generate two matched signals but with the
% size of original signals.
LED_ind_convert_from_VID_to_EPHYS = nan(size(VID_LED_inds));
LED_ind_convert_from_EPHYS_to_VID = nan(size(EPHYS_LED_inds));
for counter_ind = 1 : 1 : length(VID_LED_inds_DTW)
    ind_VID_LED_DTW  = VID_LED_inds_DTW(counter_ind);
    ind_EPHYS_LED_DTW = EPHYS_LED_inds_DTW(counter_ind);
    LED_ind_convert_from_VID_to_EPHYS(ind_EPHYS_LED_DTW) = ind_VID_LED_DTW;
    LED_ind_convert_from_EPHYS_to_VID(ind_VID_LED_DTW)  = ind_EPHYS_LED_DTW;
end

time_reference      = EPHYS_time_30k(:);
length_time         = length(time_reference);
time_EPHYS_LED_xcorr_1K = EPHYS_LED_xcorr_time_1K;
time_EPHYS_LED_xcorr_1K(end+1) = max([time_reference(end), time_EPHYS_LED_xcorr_1K(end)])+1;
event_EPHYS_LED_xcorr_30k       = nan(length(EPHYS_LED_xcorr_time_1K), 1);
counter_EPHYS_LED_xcorr     = find(time_EPHYS_LED_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time

    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_EPHYS_LED_xcorr_1K(  counter_EPHYS_LED_xcorr)
        event_EPHYS_LED_xcorr_30k(    counter_EPHYS_LED_xcorr) = counter_time_point;
        counter_EPHYS_LED_xcorr   = counter_EPHYS_LED_xcorr   + 1;
    end
end

time_reference      = EPHYS_time_1K(:);
length_time         = length(time_reference);
time_EPHYS_LED_xcorr_1K = EPHYS_LED_xcorr_time_1K;
time_EPHYS_LED_xcorr_1K(end+1) = max([time_reference(end), time_EPHYS_LED_xcorr_1K(end)])+1;
event_EPHYS_LED_xcorr_1K       = nan(length(EPHYS_LED_xcorr_time_1K), 1);
counter_EPHYS_LED_xcorr     = find(time_EPHYS_LED_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_EPHYS_LED_xcorr_1K(  counter_EPHYS_LED_xcorr)
        event_EPHYS_LED_xcorr_1K(    counter_EPHYS_LED_xcorr) = counter_time_point;
        counter_EPHYS_LED_xcorr   = counter_EPHYS_LED_xcorr   + 1;
    end
end

time_reference      = VID_time_1K(:);
length_time         = length(time_reference);
time_VID_LED_xcorr_1K = VID_LED_xcorr_time_1K;
time_VID_LED_xcorr_1K(end+1) = max([time_reference(end), time_VID_LED_xcorr_1K(end)])+1;
event_VID_LED_xcorr_1K       = nan(length(VID_LED_xcorr_time_1K), 1);
counter_VID_LED_xcorr     = find(time_VID_LED_xcorr_1K >= time_reference(1), 1, 'first');

for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_VID_LED_xcorr_1K(  counter_VID_LED_xcorr)
        event_VID_LED_xcorr_1K(    counter_VID_LED_xcorr) = counter_time_point;
        counter_VID_LED_xcorr   = counter_VID_LED_xcorr   + 1;
    end
end

EPHYS_LED_xcorr_ind_30k     = event_EPHYS_LED_xcorr_30k;
EPHYS_LED_xcorr_ind_1K     = event_EPHYS_LED_xcorr_1K;
VID_LED_xcorr_ind_1K        = event_VID_LED_xcorr_1K;
EPHYS_LED_aligned_ind_30k   = EPHYS_LED_xcorr_ind_30k( LED_ind_convert_from_VID_to_EPHYS);
EPHYS_LED_aligned_ind_1K   = EPHYS_LED_xcorr_ind_1K( LED_ind_convert_from_VID_to_EPHYS);
VID_LED_aligned_ind_1K      = VID_LED_xcorr_ind_1K(LED_ind_convert_from_EPHYS_to_VID);

EPHYS.CH_EVE.align_LED.sample_diff = sample_diff;
EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_30k            = EPHYS_LED_aligned_ind_30k;
EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_1K            = EPHYS_LED_aligned_ind_1K;
EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_1K           = EPHYS_LED_aligned_ind_1K;
EPHYS.CH_EVE.align_LED.VID_LED_aligned_ind_1K              = VID_LED_aligned_ind_1K;
EPHYS.CH_EVE.align_LED.LED_ind_convert_from_EPHYS_to_VID = LED_ind_convert_from_EPHYS_to_VID;
EPHYS.CH_EVE.align_LED.LED_ind_convert_from_VID_to_EPHYS = LED_ind_convert_from_VID_to_EPHYS;
EPHYS.CH_EVE.align_LED.EPHYS_LED_xcorr_time_1K            = EPHYS_LED_xcorr_time_1K;
EPHYS.CH_EVE.align_LED.VID_LED_xcorr_time_1K               = VID_LED_xcorr_time_1K;
EPHYS.CH_EVE.align_LED.EPHYS_LED_xcorr_ind_30k             = EPHYS_LED_xcorr_ind_30k;
EPHYS.CH_EVE.align_LED.EPHYS_LED_xcorr_ind_1K             = EPHYS_LED_xcorr_ind_1K;
EPHYS.CH_EVE.align_LED.VID_LED_xcorr_ind_1K                = VID_LED_xcorr_ind_1K;
EPHYS.CH_EVE.align_LED.VID_LED_xcorr_combined_1K       = VID_LED_xcorr_combined_1K;
EPHYS.CH_EVE.align_LED.EPHYS_LED_xcorr_combined_1K    = EPHYS_LED_xcorr_combined_1K;

if EPHYS.debug_figures == 1
    figure
    subplot(3,1,1)
    hold on
    %ind
    plot(VID_LED_xcorr_combined_1K ,'*-');
    plot(EPHYS_LED_xcorr_combined_1K ,'o-')
    xlabel('Ind')
    ylim([0 1.1])
    title(['Aligned EPHYS/VID 1K | xcorr diff: ' num2str(sample_diff)  ' / ' num2str(sample_diff/1000) 's' ])
    legend('VID', 'EPHYS')
    subplot(3,1,2)
    hold on
    %time
    plot(EPHYS_LED_xcorr_time_1K , VID_LED_xcorr_combined_1K ,'*-');
    plot(EPHYS_LED_xcorr_time_1K  ,EPHYS_LED_xcorr_combined_1K ,'o-')
    legend('VID', 'EPHYS')
    xlabel('Time (s)')
    ylim([0 1.1])
    subplot(3,1,3)
    hold on
    %dtw
    plot(VID_LED_xcorr_combined_1K(LED_ind_convert_from_VID_to_EPHYS) ,'*-');
    plot(EPHYS_LED_xcorr_combined_1K ,'o-')
    xlabel('Ind')
    ylim([0 1.1])
    legend('VID', 'EPHYS')

    disp(['length EPHYS: ' num2str( EPHYS.Alignment.time_1K(end)-EPHYS.Alignment.time_1K(1) )])
    disp(['length VID: ' num2str( VID.Alignment.time_1K(end)-VID.Alignment.time_1K(1) )])
    disp(['xcorr diff: ' num2str( sample_diff )])
    ESN_Beautify_Plot(gcf, [20 10], 7)
end
fprintf(' --> Completed. \n');


%% Save EPHYS EVENT DATA
clearvars -except EPHYS VID
EPHYS_time_1K    = EPHYS.Alignment.time_1K;
EPHYS_time_30k   = EPHYS.time_30k;
VID_time_vid   = VID.Alignment.time_vid;
VID_time_1K   = VID.Alignment.time_1K;
align_LED     = EPHYS.CH_EVE.align_LED;
path_to_save = EPHYS.path_to_save;
rec_name = EPHYS.rec_name;
file_name = [rec_name '_EVE1_aligned.mat'];
clearvars EPHYS VID
fprintf([file_name ': Saving EPHYS Event Data ...'])
save([path_to_save file_name], '-v7.3');
fprintf(' --> Completed. \n')
end