function PGH_event_alignment_LED_behave(path_to_raw, path_to_save, rec_name)
%% load LED Data
fprintf('Loading LED signal ...');

path_to_tongue = [path_to_raw '..' filesep 'analyzed_data' ...
    filesep 'behavior_data' filesep 'tongue' filesep];

vid_file_name = dir([path_to_tongue '*_video.mat']);
load([path_to_tongue vid_file_name(1).name],'FPS','LED_FPS');

VID.FPS = FPS;
VID.LED_FPS = LED_FPS;
VID.rec_name = rec_name;
VID.path_to_save = path_to_save;
VID.debug_figures = 0;

fprintf(' --> Completed. \n')
%% load BEHAVE DATA
file_name = dir([path_to_raw '*.fhd']);
file_name = [file_name(1).name(1:end-4) '.mat'];
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(path_to_raw(end), filesep);path_to_raw = [path_to_raw filesep];end
data = load([path_to_raw file_name],'data');
BEHAVE = data.data;
fprintf(' --> Completed. \n')

%% Build BEHAVE Alignment events
clearvars -except BEHAVE VID
fprintf(['Building BEHAVE Alignment events', ' ... ']);

min_length = min([ ...
    length(BEHAVE.eyelink_time),...
    length(BEHAVE.t)]);
% timeseries
inds_invalid   = false(min_length, 1);
time_eyelink   = double(BEHAVE.eyelink_time(1:min_length)');          inds_invalid = isnan(time_eyelink) | inds_invalid;
time_tgt       = double(BEHAVE.t(1:min_length)');                     inds_invalid = isnan(time_tgt)     | inds_invalid;
% correct for the bias between time_eyelink and time_tgt
time_eyelink   = time_eyelink .* (time_tgt(end)-time_tgt(1)) ./ (time_eyelink(end)-time_eyelink(1));
time_eyelink   = time_eyelink - time_eyelink(1) + time_tgt(1);
time_1K        = (time_eyelink(1) : 0.001 : time_eyelink(end))';

BEHAVE.time_1K  = time_1K;

time_reference      = BEHAVE.time_1K;
BEHAVE.Alignment.time_1K = time_reference;

length_time         = length(time_reference);

time_state_str_fixation = [];
%time_state_sac_detect_off = [];

for counter_trial = 1 : length(BEHAVE.trials)-1
    time_state_str_fixation = [time_state_str_fixation; BEHAVE.trials{counter_trial}.state_start_time_str_target_fixation(:)];
    %time_state_sac_detect_off = [time_state_sac_detect_off; BEHAVE.trials{counter_trial}.state_start_time_detect_sac_end(:)];
end

% variable_list = {'_state_str_fixation','_state_sac_detect_off'};
variable_list = {'_state_str_fixation'};

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'BEHAVE.time' variable_name ' = ' 'time' variable_name ';']);
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
    %     if time_point_ >= time_state_sac_detect_off(counter_state_sac_detect_off)
    %         event_state_sac_detect_off(  counter_time_point) = true;
    %         counter_state_sac_detect_off = counter_state_sac_detect_off + 1;
    %     end
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'BEHAVE.event' variable_name ' = ' 'event' variable_name ';']);
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

event_photodiode_combined = double(event_state_str_fixation)   .* 1  ;

BEHAVE.Alignment.event_photodiode_combined_1K = event_photodiode_combined;

% photodiode_description = [
% 'STR_TARGET_FIXATION: 1 , ', ...
% 'DETECT_SACCADE_END: 1 , ', ...
% ];
photodiode_description = 'STR_TARGET_FIXATION: 1';

% BEHAVE.Alignment.state_description      = state_description;
BEHAVE.Alignment.photodiode_description = photodiode_description;

fprintf(' --> Completed. \n');

%% Build VID Alignment events - 1K
clearvars -except BEHAVE BEHAVE VID
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

%% Alignment
fprintf('Aligning data ...')
clearvars -except BEHAVE BEHAVE VID
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

%% ALIGN BEHAVE and VID through xcorr and dtw - time_1K
clearvars -except BEHAVE VID
fprintf(['Aligning VID and BEHAVE LED_combined', ' ... ']);

BEHAVE_time_1K            = BEHAVE.time_1K;
VID_time_1K = VID.Alignment.time_1K;

VID_LED_combined_1K = VID.Alignment.event_LED_combined_1K;
BEHAVE_PD_combined_1K = BEHAVE.Alignment.event_photodiode_combined_1K;

% LED_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(VID_LED_combined_1K+1, BEHAVE_PD_combined_1K+1); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);

if sample_diff > 0
    VID_LED_xcorr_combined_1K  = VID_LED_combined_1K(abs(sample_diff):end);
    VID_LED_xcorr_time_1K                 = VID_time_1K(abs(sample_diff):end);
    BEHAVE_PD_xcorr_combined_1K = BEHAVE_PD_combined_1K;
    BEHAVE_PD_xcorr_time_1K                = BEHAVE_time_1K;
elseif sample_diff < 0
    VID_LED_xcorr_combined_1K  = VID_LED_combined_1K;
    VID_LED_xcorr_time_1K                 = VID_time_1K;
    BEHAVE_PD_xcorr_combined_1K = BEHAVE_PD_combined_1K(abs(sample_diff):end);
    BEHAVE_PD_xcorr_time_1K                = BEHAVE_time_1K(abs(sample_diff):end);
end

% LED_combined: make the vectors the same size
if length(BEHAVE_PD_xcorr_combined_1K) ~= length(VID_LED_xcorr_combined_1K)
    min_length = min([ length(BEHAVE_PD_xcorr_combined_1K),  length(VID_LED_xcorr_combined_1K)]);
    VID_LED_xcorr_combined_1K  = VID_LED_xcorr_combined_1K(  1:min_length);
    VID_LED_xcorr_time_1K                 = VID_LED_xcorr_time_1K(                 1:min_length);
    BEHAVE_PD_xcorr_combined_1K = BEHAVE_PD_xcorr_combined_1K( 1:min_length);
    BEHAVE_PD_xcorr_time_1K                = BEHAVE_PD_xcorr_time_1K(                1:min_length);
end


% LED_combined: find the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
sampling_freq = 1000;
cutoff_freq = 25.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
VID_LED_xcorr_combined_1K_filt  = filtfilt(b_butter,a_butter,VID_LED_xcorr_combined_1K);
BEHAVE_PD_xcorr_combined_1K_filt = filtfilt(b_butter,a_butter,BEHAVE_PD_xcorr_combined_1K);
% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(VID_LED_xcorr_time_1K ) / 500);
ind_edges = round(linspace(1, length(VID_LED_xcorr_time_1K ), ind_edge_width));
ind_edges(1) = 0;

% init and loop over chunks
VID_LED_inds_DTW  = cell((length(ind_edges)-1), 1);
BEHAVE_PD_inds_DTW = cell((length(ind_edges)-1), 1);
for counter_chunk = 1 : 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    VID_LED_LED_combined_chunk  = VID_LED_xcorr_combined_1K_filt(inds_chunk);
    BEHAVE_PD_LED_combined_chunk = BEHAVE_PD_xcorr_combined_1K_filt(inds_chunk);
    [~,ix,iy] = dtw(VID_LED_LED_combined_chunk,BEHAVE_PD_LED_combined_chunk, 15, 'absolute'); % allow upto 15ms warp
    VID_LED_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    BEHAVE_PD_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end

VID_LED_inds_DTW  = cell2mat(VID_LED_inds_DTW);
BEHAVE_PD_inds_DTW = cell2mat(BEHAVE_PD_inds_DTW);
VID_LED_inds      = ( 1 : 1 : length(VID_LED_xcorr_time_1K ) )';
BEHAVE_PD_inds     = ( 1 : 1 : length(BEHAVE_PD_xcorr_time_1K) )';

% dtw works by replicating the inds to match the two signals, here we
% reverse the replicated inds to generate two matched signals but with the
% size of original signals.
LED_ind_convert_from_VID_to_BEHAVE = nan(size(VID_LED_inds));
LED_ind_convert_from_BEHAVE_to_VID = nan(size(BEHAVE_PD_inds));
for counter_ind = 1 : 1 : length(VID_LED_inds_DTW)
    ind_VID_LED_DTW  = VID_LED_inds_DTW(counter_ind);
    ind_BEHAVE_PD_DTW = BEHAVE_PD_inds_DTW(counter_ind);
    LED_ind_convert_from_VID_to_BEHAVE(ind_BEHAVE_PD_DTW) = ind_VID_LED_DTW;
    LED_ind_convert_from_BEHAVE_to_VID(ind_VID_LED_DTW)  = ind_BEHAVE_PD_DTW;
end

time_reference      = BEHAVE_time_1K(:);
length_time         = length(time_reference);
time_BEHAVE_PD_xcorr_1K = BEHAVE_PD_xcorr_time_1K;
time_BEHAVE_PD_xcorr_1K(end+1) = max([time_reference(end), time_BEHAVE_PD_xcorr_1K(end)])+1;
event_BEHAVE_PD_xcorr_1K       = nan(length(BEHAVE_PD_xcorr_time_1K), 1);
counter_BEHAVE_PD_xcorr     = find(time_BEHAVE_PD_xcorr_1K >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_BEHAVE_PD_xcorr_1K(  counter_BEHAVE_PD_xcorr)
        event_BEHAVE_PD_xcorr_1K(    counter_BEHAVE_PD_xcorr) = counter_time_point;
        counter_BEHAVE_PD_xcorr   = counter_BEHAVE_PD_xcorr   + 1;
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

BEHAVE_PD_xcorr_ind_1K     = event_BEHAVE_PD_xcorr_1K;
VID_LED_xcorr_ind_1K        = event_VID_LED_xcorr_1K;
BEHAVE_PD_aligned_ind_1K   = BEHAVE_PD_xcorr_ind_1K( LED_ind_convert_from_VID_to_BEHAVE);
VID_LED_aligned_ind_1K      = VID_LED_xcorr_ind_1K(LED_ind_convert_from_BEHAVE_to_VID);

VID.CH_EVE.align_PD.sample_diff = sample_diff;
VID.CH_EVE.align_PD.BEHAVE_PD_aligned_ind_1K            = BEHAVE_PD_aligned_ind_1K;
VID.CH_EVE.align_PD.BEHAVE_PD_aligned_ind_1K           = BEHAVE_PD_aligned_ind_1K;
VID.CH_EVE.align_PD.VID_LED_aligned_ind_1K              = VID_LED_aligned_ind_1K;
VID.CH_EVE.align_PD.LED_ind_convert_from_BEHAVE_to_VID = LED_ind_convert_from_BEHAVE_to_VID;
VID.CH_EVE.align_PD.LED_ind_convert_from_VID_to_BEHAVE = LED_ind_convert_from_VID_to_BEHAVE;
VID.CH_EVE.align_PD.BEHAVE_PD_xcorr_time_1K            = BEHAVE_PD_xcorr_time_1K;
VID.CH_EVE.align_PD.VID_LED_xcorr_time_1K               = VID_LED_xcorr_time_1K;
VID.CH_EVE.align_PD.BEHAVE_PD_xcorr_ind_1K             = BEHAVE_PD_xcorr_ind_1K;
VID.CH_EVE.align_PD.VID_LED_xcorr_ind_1K                = VID_LED_xcorr_ind_1K;
VID.CH_EVE.align_PD.VID_LED_xcorr_combined_1K       = VID_LED_xcorr_combined_1K;
VID.CH_EVE.align_PD.BEHAVE_PD_xcorr_combined_1K    = BEHAVE_PD_xcorr_combined_1K;

if VID.debug_figures == 1
    figure
    subplot(3,1,1)
    hold on
    %ind
    plot(VID_LED_xcorr_combined_1K ,'*-');
    plot(BEHAVE_PD_xcorr_combined_1K ,'o-')
    xlabel('Ind')
    ylim([0 1.1])
    title(['Aligned BEHAVE/VID 1K | xcorr diff: ' num2str(sample_diff)  ' / ' num2str(sample_diff/1000) 's' ])
    legend('VID', 'BEHAVE')
    subplot(3,1,2)
    hold on
    %time
    plot(BEHAVE_PD_xcorr_time_1K , VID_LED_xcorr_combined_1K ,'*-');
    plot(BEHAVE_PD_xcorr_time_1K  ,BEHAVE_PD_xcorr_combined_1K ,'o-')
    legend('VID', 'BEHAVE')
    xlabel('Time (s)')
    ylim([0 1.1])
    subplot(3,1,3)
    hold on
    %dtw
    plot(VID_LED_xcorr_combined_1K(LED_ind_convert_from_VID_to_BEHAVE) ,'*-');
    plot(BEHAVE_PD_xcorr_combined_1K ,'o-')
    xlabel('Ind')
    ylim([0 1.1])
    legend('VID', 'BEHAVE')

    disp(['length BEHAVE: ' num2str( BEHAVE.Alignment.time_1K(end)-BEHAVE.Alignment.time_1K(1) )])
    disp(['length VID: ' num2str( VID.Alignment.time_1K(end)-VID.Alignment.time_1K(1) )])
    disp(['xcorr diff: ' num2str( sample_diff )])
    ESN_Beautify_Plot(gcf, [20 10], 7)
end
fprintf(' --> Completed. \n');


%% Save BEHAVE EVENT DATA
clearvars -except BEHAVE VID
align_PD     = VID.CH_EVE.align_PD;

path_to_save = VID.path_to_save;
rec_name = VID.rec_name;
file_name = [rec_name '_EVE1_aligned.mat'];
clearvars BEHAVE VID
fprintf([file_name ': Saving BEHAVE Event Data ...'])
save([path_to_save file_name], '-append');
fprintf(' --> Completed. \n')
end