function PH_LED_FPS_alignment_COMBINED
%% clear
clear;
close all;

%% load EPHYS EVENT DATA
clearvars -except EPHYS BEHAVE VID

file_path = pwd;
[file_name,file_path] = uigetfile([file_path filesep 'all_channels.events'], 'Select all_channels.events file');
fprintf(['Loading ', file_name, ' ... ']);
[ch_data, ch_time, ch_info] = load_open_ephys_data([file_path filesep file_name]);
EPHYS.CH_EVE.ch_data = ch_data;
EPHYS.CH_EVE.ch_time = ch_time;
EPHYS.CH_EVE.ch_info = ch_info;
EPHYS.file_name_CH_EVE = file_name;
EPHYS.file_path_CH_EVE = file_path;
fprintf(' --> Completed. \n')

%% load continuos file
clearvars -except EPHYS BEHAVE VID

[file_name,file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep '*.continuous'], 'Select a continuous file');
fprintf(['Loading ', file_name, ' ... ']);
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
[~, ch_time, ~] = load_open_ephys_data([file_path file_name]);
EPHYS.ch_time = ch_time;
fprintf(' --> Completed. \n')

EPHYS.time_30K   = EPHYS.ch_time;
EPHYS.time_1K      = (EPHYS.ch_time(1) : 0.001 :EPHYS.ch_time(end))';
EPHYS.time_100      = (EPHYS.ch_time(1) : 0.01 :EPHYS.ch_time(end))';

%% load BEHAVE DATA
clearvars -except EPHYS BEHAVE VID

[file_name,file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep '*.mat'], 'Select BEHAVE file');
fprintf(['Loading ', file_name, ' ... ']);
data = load([file_path filesep file_name],'data');
BEHAVE = data.data;

time_eyelink   = double(BEHAVE.eyelink_time(1:length(BEHAVE.eyelink_time))');          
BEHAVE.time_1K       = (time_eyelink(1) : 0.001 : time_eyelink(end))';

fprintf(' --> Completed. \n')

%% load VID DATA
clearvars -except EPHYS BEHAVE VID

[file_name,file_path] = uigetfile([EPHYS.file_path_CH_EVE filesep '*.mp4'], 'Select video file');
VID.vid = VideoReader([file_path filesep file_name]);
VID.file_name = file_name;
VID.file_path = file_path;

%% LED
clearvars -except EPHYS BEHAVE VID
    %% Main program
vid = VID.vid;
imshow(read(VID.vid, 1));
rect = getrect;
close;
try
    num_frames = vid.NumFrames;
    num_frames_written= 0;
   
    for counter_frame = 1 : num_frames
        current_frame = read(vid, counter_frame);
        
        x = round(rect(1));
        y = round(rect(2));
        w = round(rect(3));
        h = round(rect(4));
        
        mean_r(counter_frame) = round(mean(mean(current_frame(y:y+h, x:x+w, 1))));
        mean_g(counter_frame) = round(mean(mean(current_frame(y:y+h, x:x+w, 2))));
        mean_b(counter_frame) = round(mean(mean(current_frame(y:y+h, x:x+w, 3))));
        
        progress_message = (sprintf('Processed frame %4d of %d.', counter_frame, num_frames));
        disp(progress_message)
        % Increment frame count (should eventually = numberOfFrames)
        num_frames_written = num_frames_written + 1;
        
    end
    
catch ME
    error_message = sprintf('Error extracting movie frames from:\n\n%s\n\nError: %s\n\n)', file_name, ME.message);
    uiwait(msgbox(error_message));
end
close;

% Cluster pixel values and select threshold based on the midpoint of the mean of each cluster
mean_r = mean_r';
clust_mean_r = kmeans(mean_r, 2);

mean_1 = mean(mean_r(clust_mean_r == 1));
mean_2 = mean(mean_r(clust_mean_r == 2));
threshold = round((mean_1 + mean_2)/2);

figure;
hold on
histogram(mean_r);
xline(mean_1, 'r', 'LineWidth', 2);
xline(mean_2, 'b', 'LineWidth', 2);
xline(threshold, 'k', 'LineWidth', 2);
title(['Red pixel value distribution | Selected threshold: ' num2str(threshold)]);
xlabel('Pixel value');
ylabel('Count');

% User manually specifies a threshold
% threshold = input('Specify threshold: ');

for counter_frame = 1:length(mean_r)
    if (mean_r(counter_frame) < threshold)
        LED_trace(counter_frame) = 0;
    else
        LED_trace(counter_frame) = 1;
    end
end

LED_trace_rising = [diff(LED_trace) 0] > 0;
LED_trace_falling = [diff(LED_trace) 0] < 0;

% figure;
% hold on;
% plot(LED_trace,'k');
% title(['LED trace | Pulses: ' num2str(sum(LED_trace_rising)) ]);
% xlabel('Frame #');
% ylabel('ON:1/OFF:0');
% ylim([0 1.1])

LED(:,1) = mean_r;
LED(:,2) = mean_g;
LED(:,3) = mean_b;
LED(:,4) = LED_trace;
LED(:,5) = LED_trace_rising;
LED(:,6) = LED_trace_falling;
    %% Save _LED file 
    parts = strsplit(VID.file_path, ["/", "\"]);
    name = erase(parts(length(parts) - 2), '-');
    file_name = char(extractBetween(name, 3, 15));
    file_path = VID.file_path;     
    save([file_path '../analyzed_data/' file_name '_LED.mat'], 'LED',  '-v7.3');
    VID.LED = LED;
    fprintf(' --> Completed. \n')
   
%% FPS
clearvars -except EPHYS BEHAVE VID
    %% Estimate FPS
fprintf(['Estimating FPS', ' ... \n']);

% fps_LB = input('Enter fps lower bound: ');
% fps_UB = input('Enter fps upper bound: ');
% fps_step = input('Enter fps step size: ');

% fps_LB = 64;
% fps_UB = 65;
fps_LB = 98;
fps_UB = 101;
fps_step = 0.001;

fps = fps_LB:fps_step:fps_UB;
    for counter_FPS = 1:length(fps)
    %% Build BEHAVE Alignment events - FPS
    time_reference      = (BEHAVE.time_1K(1) : (1/fps(counter_FPS)) : BEHAVE.time_1K(end))';
    length_time         = length(time_reference);
    
    time_state_str_fixation = [];
    %time_state_sac_detect_off = [];
    for counter_trial = 1 : length(BEHAVE.trials)-1
        time_state_str_fixation = [time_state_str_fixation; BEHAVE.trials{counter_trial}.state_start_time_str_target_fixation(:)];
        %time_state_sac_detect_off = [time_state_sac_detect_off; BEHAVE.trials{counter_trial}.state_start_time_detect_sac_end(:)];
    end  
    
    % variable_list = {'_state_str_fixation','_state_sac_detect_off'};
    variable_list = {'_state_str_fixation'};
    
    BEHAVE.time_vid = time_reference;
    
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
        
        inds_span_ = (0) : 1 : (20);
        
        ind_event_temp_  = find(event_temp_);
        inds_event_temp_ = repmat( ind_event_temp_(:), 1, length(inds_span_)) + repmat(inds_span_(:)', length(ind_event_temp_), 1);
        inds_event_temp_( inds_event_temp_ < 1 ) = 1;
        inds_event_temp_( inds_event_temp_ > length_time ) = length_time;
        event_temp_(inds_event_temp_(:)) = true;
        eval([ 'event'   variable_name  ' = ' 'event_temp_'   ';']);
    end
    
    BEHAVE_LED_combined = double(event_state_str_fixation)   .* 1  ;           
    %% Build VID Alignment events - FPS
    time_reference = ((1/fps(counter_FPS)) : (1/fps(counter_FPS)) : length(VID.LED(:,4))/fps(counter_FPS))';
    length_time = length(time_reference);
    time_LED_rise = time_reference([diff(VID.LED(:,4)) ; 0] > 0);
    
    variable_list = { '_LED_rise'};
    
    VID.time_vid = time_reference;
    
    for counter_variable = 1 : 1 : length(variable_list)
        variable_name = variable_list{counter_variable};
        eval([ 'VID.time' variable_name ' = ' 'time' variable_name ';']);
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
        if time_ponit_ >= time_LED_rise(     counter_LED_rise)
            event_LED_rise(       counter_time_point) = true;
            counter_LED_rise      = counter_LED_rise      + 1;
        end
        
    end
    
    for counter_variable = 1 : 1 : length(variable_list)
        variable_name = variable_list{counter_variable};
        eval([ 'VID.event' variable_name ' = ' 'event' variable_name ';']);
    end
    
    for counter_variable = 1 : 1 : length(variable_list)
        variable_name = variable_list{counter_variable};
        eval([ 'event_temp_' ' = ' 'event' variable_name ';']);
        length_time = length(event_temp_);
        inds_span_ = (0) : 1 : (20);
        ind_event_temp_  = find(event_temp_);
        inds_event_temp_ = repmat( ind_event_temp_(:), 1, length(inds_span_)) + repmat(inds_span_(:)', length(ind_event_temp_), 1);
        inds_event_temp_( inds_event_temp_ < 1 ) = 1;
        inds_event_temp_( inds_event_temp_ > length_time ) = length_time;
        event_temp_(inds_event_temp_(:)) = true;
        eval([ 'event'   variable_name  ' = ' 'event_temp_'   ';']);
    end
    
   VID_LED_combined = double(event_LED_rise)  .* 1;
    %% Align   
    [xcorr_value,xcorr_lag] = xcorr(VID_LED_combined+1, BEHAVE_LED_combined+1); % cross-correlate signals with each other
    [max_xcross,ind_max_xcross] = max(abs(xcorr_value));
    sample_diff = xcorr_lag(ind_max_xcross);
    
    ind_max_xcross_all(1, counter_FPS) = ind_max_xcross;
    max_xcross_all(1, counter_FPS) = max_xcross;
    end
    %% FPS estimate report
[xcross_FPS, ind_xcross_FPS] = max(max_xcross_all);
estimated_FPS = fps(ind_xcross_FPS);
fprintf(['Estimated FPS (BEHAVE): ' num2str(estimated_FPS)])

figure;
hold on;
plot(fps, max_xcross_all)
scatter(estimated_FPS, xcross_FPS,  'o')
xlabel('fps')
ylabel('xcorr value')
title(['Estimated FPS (BEHAVE): ' num2str(estimated_FPS)])

VID.fps = estimated_FPS;

fprintf(' --> Completed. \n')
    %% Save Estimated FPS and figure
    
    parts = strsplit(VID.file_path, ["/", "\"]);
    name = erase(parts(length(parts) - 2), '-');
    file_name = char(extractBetween(name, 3, 15));
    file_path = VID.file_path;
    file_name_pdf = [file_name '_FPS(BEHAVE)_' num2str(estimated_FPS) '.pdf' ];
    fprintf('Saving Estimated FPS (BEHAVE) ...')
    saveas(gcf,[file_path '../analyzed_data/' file_name_pdf])
    fprintf(' --> Completed. \n')

%% Align
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
    %% Build EPHYS Alignment events for VID - 100
clearvars -except EPHYS BEHAVE VID
fprintf(['Building EPHYS Alignment events for VID', ' ... ']);

time_reference      = EPHYS.time_100;
length_time         = length(time_reference);
time_state_str_fixation   = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 3) & ((EPHYS.CH_EVE.data(:,3) == 1)) , 1);
%time_state_sac_detect_off = EPHYS.CH_EVE.data( (EPHYS.CH_EVE.data(:,2) == 4) & ((EPHYS.CH_EVE.data(:,3) == 0)) , 1);

% variable_list = {'_state_str_fixation','_state_sac_detect_off'};
variable_list = {'_state_str_fixation'};

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
    eval([ 'EPHYS.Alignment.event' variable_name ' = ' 'event' variable_name ';']);
end

for counter_variable = 1 : 1 : length(variable_list)
    variable_name = variable_list{counter_variable};
    eval([ 'event_temp_' ' = ' 'event' variable_name ';']);
    length_time = length(event_temp_);
  
        inds_span_ = (0) : 1 : (20);
   
    ind_event_temp_  = find(event_temp_);
    inds_event_temp_ = repmat( ind_event_temp_(:), 1, length(inds_span_)) + repmat(inds_span_(:)', length(ind_event_temp_), 1);
    inds_event_temp_( inds_event_temp_ < 1 ) = 1;
    inds_event_temp_( inds_event_temp_ > length_time ) = length_time;
    event_temp_(inds_event_temp_(:)) = true;
    eval([ 'event'   variable_name  ' = ' 'event_temp_'   ';']);
end

event_LED_combined = double(event_state_str_fixation)   .* 1  ;

EPHYS.Alignment.event_LED_combined_100 = event_LED_combined;



fprintf(' --> Completed. \n');
    %% Build VID Alignment events - 100
clearvars -except EPHYS BEHAVE VID
fprintf(['Building VID Alignment events', ' ... ']);

time_reference = ((1/VID.fps) : (1/VID.fps) : length(VID.LED(:,4))/VID.fps)';
length_time = length(time_reference);
time_LED_rise = time_reference([diff(VID.LED(:,4)) ; 0] > 0);

variable_list = { '_LED_rise'};

VID.Alignment.time_vid = time_reference;

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
  
        inds_span_ = (0) : 1 : (20);
    
    ind_event_temp_  = find(event_temp_);
    inds_event_temp_ = repmat( ind_event_temp_(:), 1, length(inds_span_)) + repmat(inds_span_(:)', length(ind_event_temp_), 1);
    inds_event_temp_( inds_event_temp_ < 1 ) = 1;
    inds_event_temp_( inds_event_temp_ > length_time ) = length_time;
    event_temp_(inds_event_temp_(:)) = true;
    eval([ 'event'   variable_name  ' = ' 'event_temp_'   ';']);
end


event_LED_combined = double(event_LED_rise)  .* 1;

VID.Alignment.event_LED_combined_vid = event_LED_combined;


% Upsample video from time_vid to time_100 - 100 Hz
VID_time_100 = VID.Alignment.time_vid(1) : 0.01 : VID.Alignment.time_vid(end);
VID_LED_combined_100 = round(interp1(VID.Alignment.time_vid,  VID.Alignment.event_LED_combined_vid , VID_time_100));

VID.Alignment.event_LED_combined_100 = VID_LED_combined_100;
VID.Alignment.time_100 = VID_time_100;


fprintf(' --> Completed. \n');
    %% ALIGN EPHYS and VID through xcorr and dtw - time_100
clearvars -except EPHYS BEHAVE VID
fprintf(['Aligning VID and EPHYS LED_combined', ' ... ']);

EPHYS_time_100            = EPHYS.time_100;
EPHYS_time_30K             = EPHYS.time_30K;
VID_time_100 = VID.Alignment.time_100;

VID_LED_combined_100 = VID.Alignment.event_LED_combined_100;
EPHYS_LED_combined_100 = EPHYS.Alignment.event_LED_combined_100;


% LED_combined: find the bias between 2 signals
[xcorr_value,xcorr_lag] = xcorr(VID_LED_combined_100+1, EPHYS_LED_combined_100+1); % cross-correlate signals with each other
[~,ind_max_xcross] = max(abs(xcorr_value));
sample_diff = xcorr_lag(ind_max_xcross);

if sample_diff > 0
    VID_LED_xcorr_LED_combined_100  = VID_LED_combined_100(abs(sample_diff):end);
    VID_LED_xcorr_time_100                 = VID_time_100(abs(sample_diff):end);
    EPHYS_LED_xcorr_LED_combined_100 = EPHYS_LED_combined_100;
    EPHYS_LED_xcorr_time_100                = EPHYS_time_100;
elseif sample_diff < 0
    VID_LED_xcorr_LED_combined_100  = VID_LED_combined_100;
    VID_LED_xcorr_time_100                 = VID_time_100;
    EPHYS_LED_xcorr_LED_combined_100 = EPHYS_LED_combined_100(abs(sample_diff):end);
    EPHYS_LED_xcorr_time_100                = EPHYS_time_100(abs(sample_diff):end);
end

% LED_combined: make the vectors the same size
if length(EPHYS_LED_xcorr_LED_combined_100) ~= length(VID_LED_xcorr_LED_combined_100)
    min_length = min([ length(EPHYS_LED_xcorr_LED_combined_100),  length(VID_LED_xcorr_LED_combined_100)]);
    VID_LED_xcorr_LED_combined_100  = VID_LED_xcorr_LED_combined_100(  1:min_length);
    VID_LED_xcorr_time_100                 = VID_LED_xcorr_time_100(                 1:min_length);
    EPHYS_LED_xcorr_LED_combined_100 = EPHYS_LED_xcorr_LED_combined_100( 1:min_length);
    EPHYS_LED_xcorr_time_100                = EPHYS_LED_xcorr_time_100(                1:min_length);
end


% LED_combined: find the Dynamic Time Warp (DTW) between 2 time series
% low pass filter the signals to generate a sinusoid around rises and falls
sampling_freq = 100;
cutoff_freq = 25.0;
[b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
VID_LED_xcorr_LED_combined_100_filt  = filtfilt(b_butter,a_butter,VID_LED_xcorr_LED_combined_100);
EPHYS_LED_xcorr_LED_combined_100_filt = filtfilt(b_butter,a_butter,EPHYS_LED_xcorr_LED_combined_100);
% break the dtw analayises to smaller chunks, dtw does not work with large vectors
ind_edge_width = ceil(length(VID_LED_xcorr_time_100 ) / 500);
ind_edges = round(linspace(1, length(VID_LED_xcorr_time_100 ), ind_edge_width));
ind_edges(1) = 0;

% init and loop over chunks
VID_LED_inds_DTW  = cell((length(ind_edges)-1), 1);
EPHYS_LED_inds_DTW = cell((length(ind_edges)-1), 1);
for counter_chunk = 1 : 1 : (length(ind_edges)-1)
    inds_chunk = ( (ind_edges(counter_chunk)+1) : 1 : (ind_edges(counter_chunk+1)) )';
    VID_LED_LED_combined_chunk  = VID_LED_xcorr_LED_combined_100_filt(inds_chunk);
    EPHYS_LED_LED_combined_chunk = EPHYS_LED_xcorr_LED_combined_100_filt(inds_chunk);
    [~,ix,iy] = dtw(VID_LED_LED_combined_chunk,EPHYS_LED_LED_combined_chunk, 50, 'absolute'); % allow upto 50ms warp
    VID_LED_inds_DTW{counter_chunk}  = ix(:) + inds_chunk(1) - 1;
    EPHYS_LED_inds_DTW{counter_chunk} = iy(:) + inds_chunk(1) - 1;
end

VID_LED_inds_DTW  = cell2mat(VID_LED_inds_DTW);
EPHYS_LED_inds_DTW = cell2mat(EPHYS_LED_inds_DTW);
VID_LED_inds      = ( 1 : 1 : length(VID_LED_xcorr_time_100 ) )';
EPHYS_LED_inds     = ( 1 : 1 : length(EPHYS_LED_xcorr_time_100) )';


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

time_reference      = EPHYS_time_30K(:);
length_time         = length(time_reference);
time_EPHYS_LED_xcorr_100 = EPHYS_LED_xcorr_time_100;
time_EPHYS_LED_xcorr_100(end+1) = max([time_reference(end), time_EPHYS_LED_xcorr_100(end)])+1;
event_EPHYS_LED_xcorr_30K       = nan(length(EPHYS_LED_xcorr_time_100), 1);
counter_EPHYS_LED_xcorr     = find(time_EPHYS_LED_xcorr_100 >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_EPHYS_LED_xcorr_100(  counter_EPHYS_LED_xcorr)
        event_EPHYS_LED_xcorr_30K(    counter_EPHYS_LED_xcorr) = counter_time_point;
        counter_EPHYS_LED_xcorr   = counter_EPHYS_LED_xcorr   + 1;
    end
end

time_reference      = EPHYS_time_100(:);
length_time         = length(time_reference);
time_EPHYS_LED_xcorr_100 = EPHYS_LED_xcorr_time_100;
time_EPHYS_LED_xcorr_100(end+1) = max([time_reference(end), time_EPHYS_LED_xcorr_100(end)])+1;
event_EPHYS_LED_xcorr_100       = nan(length(EPHYS_LED_xcorr_time_100), 1);
counter_EPHYS_LED_xcorr     = find(time_EPHYS_LED_xcorr_100 >= time_reference(1), 1, 'first');
for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_EPHYS_LED_xcorr_100(  counter_EPHYS_LED_xcorr)
        event_EPHYS_LED_xcorr_100(    counter_EPHYS_LED_xcorr) = counter_time_point;
        counter_EPHYS_LED_xcorr   = counter_EPHYS_LED_xcorr   + 1;
    end
end

time_reference      = VID_time_100(:);
length_time         = length(time_reference);
time_VID_LED_xcorr_100 = VID_LED_xcorr_time_100;
time_VID_LED_xcorr_100(end+1) = max([time_reference(end), time_VID_LED_xcorr_100(end)])+1;
event_VID_LED_xcorr_100       = nan(length(VID_LED_xcorr_time_100), 1);
counter_VID_LED_xcorr     = find(time_VID_LED_xcorr_100 >= time_reference(1), 1, 'first');

for counter_time_point = 1 : length_time
    time_point_     = time_reference(counter_time_point);
    if time_point_ >= time_VID_LED_xcorr_100(  counter_VID_LED_xcorr)
        event_VID_LED_xcorr_100(    counter_VID_LED_xcorr) = counter_time_point;
        counter_VID_LED_xcorr   = counter_VID_LED_xcorr   + 1;
    end
end

EPHYS_LED_xcorr_ind_30K     = event_EPHYS_LED_xcorr_30K;
EPHYS_LED_xcorr_ind_100     = event_EPHYS_LED_xcorr_100;
VID_LED_xcorr_ind_100        = event_VID_LED_xcorr_100;
EPHYS_LED_aligned_ind_30K   = EPHYS_LED_xcorr_ind_30K( LED_ind_convert_from_VID_to_EPHYS);
EPHYS_LED_aligned_ind_100   = EPHYS_LED_xcorr_ind_100( LED_ind_convert_from_VID_to_EPHYS);
VID_LED_aligned_ind_100      = VID_LED_xcorr_ind_100(LED_ind_convert_from_EPHYS_to_VID);

EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_30K            = EPHYS_LED_aligned_ind_30K;
EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_100            = EPHYS_LED_aligned_ind_100;
EPHYS.CH_EVE.align_LED.EPHYS_LED_aligned_ind_100           = EPHYS_LED_aligned_ind_100;
EPHYS.CH_EVE.align_LED.VID_LED_aligned_ind_100              = VID_LED_aligned_ind_100;
EPHYS.CH_EVE.align_LED.LED_ind_convert_from_EPHYS_to_VID = LED_ind_convert_from_EPHYS_to_VID;
EPHYS.CH_EVE.align_LED.LED_ind_convert_from_VID_to_EPHYS = LED_ind_convert_from_VID_to_EPHYS;
EPHYS.CH_EVE.align_LED.EPHYS_LED_xcorr_time_100            = EPHYS_LED_xcorr_time_100;
EPHYS.CH_EVE.align_LED.VID_LED_xcorr_time_100               = VID_LED_xcorr_time_100;
EPHYS.CH_EVE.align_LED.EPHYS_LED_xcorr_ind_30K             = EPHYS_LED_xcorr_ind_30K;
EPHYS.CH_EVE.align_LED.EPHYS_LED_xcorr_ind_100             = EPHYS_LED_xcorr_ind_100;
EPHYS.CH_EVE.align_LED.VID_LED_xcorr_ind_100                = VID_LED_xcorr_ind_100;
EPHYS.CH_EVE.align_LED.VID_LED_xcorr_LED_combined_100       = VID_LED_xcorr_LED_combined_100;
EPHYS.CH_EVE.align_LED.EPHYS_LED_xcorr_LED_combined_100    = EPHYS_LED_xcorr_LED_combined_100;

fprintf(' --> Completed. \n');


%PLOT RAW
% figure(1); 
% hold on;
% %plot(VID.LED(:,4), 'k');
% plot(VID.Alignment.event_LED_combined_100, 'r');
% plot(EPHYS.Alignment.event_LED_combined_100, 'b');
% title(['LED | Raw: ' num2str(sum(VID.LED(:,5) > 0)) ', Built: '  num2str(sum(diff(VID.Alignment.event_LED_combined_100)> 0)) ...
%     '. EPHYS | Raw: ' num2str(sum(diff(VID.Alignment.event_LED_combined_100) > 0))]);
% legend( 'LED', 'EPHYS')
% ylim([0 1.1])

%PLOT ALIGNED
figure
subplot(3,1,1)
hold on
%ind
plot(VID_LED_xcorr_LED_combined_100 ,'*-');
plot(EPHYS_LED_xcorr_LED_combined_100 ,'o-')
xlabel('Ind')
ylim([0 1.1])
title(['Aligned EPHYS/VID 100 | xcorr diff: ' num2str(sample_diff)  ' / ' num2str(sample_diff/100) 's' ])
legend('VID', 'EPHYS')
subplot(3,1,2)
hold on
%time
plot(EPHYS_LED_xcorr_time_100 , VID_LED_xcorr_LED_combined_100 ,'*-');
plot(EPHYS_LED_xcorr_time_100  ,EPHYS_LED_xcorr_LED_combined_100 ,'o-')
legend('VID', 'EPHYS')
xlabel('Time (s)')
ylim([0 1.1])
subplot(3,1,3)
hold on
%dtw
plot(VID_LED_xcorr_LED_combined_100(LED_ind_convert_from_VID_to_EPHYS) ,'*-');
plot(EPHYS_LED_xcorr_LED_combined_100 ,'o-')
xlabel('Ind')
ylim([0 1.1])
legend('VID', 'EPHYS')

disp(['length EPHYS: ' num2str( EPHYS.Alignment.time_100(end)-EPHYS.Alignment.time_100(1) )])
disp(['length VID: ' num2str( VID.Alignment.time_100(end)-VID.Alignment.time_100(1) )])
disp(['xcorr diff: ' num2str( sample_diff )])


fprintf(' --> Completed. \n');
    %% Save EPHYS EVENT DATA
clearvars -except EPHYS BEHAVE VID
% EPHYS_time_vid    = EPHYS.Alignment.time_vid;
EPHYS_time_100    = EPHYS.Alignment.time_100;
EPHYS_time_1K   = EPHYS.time_1K;
EPHYS_time_30K   = EPHYS.time_30K;

VID_time_vid   = VID.Alignment.time_vid;
VID_time_100   = VID.Alignment.time_100;

align_LED     = EPHYS.CH_EVE.align_LED;

parts = strsplit(EPHYS.file_path_CH_EVE, ["/", "\"]);
name = erase(parts(length(parts) - 2), '-');
file_name = char(extractBetween(name, 3, 15));
file_path = EPHYS.file_path_CH_EVE;
file_name = [file_name '_EVE1_PH_aligned.mat'];
clearvars EPHYS BEHAVE VID
fprintf([file_name ': Saving EPHYS Event Data ...'])
save([file_path '../analyzed_data/' file_name], '-v7.3');
fprintf(' --> Completed. \n')


end