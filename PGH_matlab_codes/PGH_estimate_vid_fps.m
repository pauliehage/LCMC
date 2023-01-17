function [LED_FPS, FPS, height, width, duration, num_frames] = PGH_estimate_vid_fps(path_to_raw, flag_figure, fps_LB, fps_UB, fps_step)
%% load BEHAVE DATA
if ~strcmp(path_to_raw(end), filesep);path_to_raw = [path_to_raw filesep];end
file_name_ = dir([path_to_raw '*.fhd']);
file_name = [file_name_(1).name(1:end-4) '.mat'];
fprintf(['Loading: ', file_name, ' ... ']);
data = load([path_to_raw file_name],'data');
BEHAVE = data.data;
time_eyelink   = double(BEHAVE.eyelink_time(1:length(BEHAVE.eyelink_time))');
BEHAVE.time_1K       = (time_eyelink(1) : 0.001 : time_eyelink(end))';
VID.flag_figure = flag_figure;

fprintf(' --> Completed. \n')
%% load VID DATA
if ~strcmp(path_to_raw(end), filesep);path_to_raw = [path_to_raw filesep];end
path_to_analyzed_tongue = [path_to_raw, '..' filesep 'analyzed_data', filesep, ...
    'behavior_data', filesep, 'tongue', filesep];

path_to_analyzed_figs_tongue = [path_to_raw, '..' filesep 'analyzed_figs', filesep, ...
    'behavior_data', filesep, 'tongue', filesep];
try
    file_name_ = dir([path_to_analyzed_tongue '*_DLC.mp4']);
    file_name = [file_name_(1).name];
    file_path = [path_to_analyzed_tongue filesep];
    VID.vid = VideoReader([file_path file_name]);
    fprintf(['Loading: ', file_name, ' ... ']);
    num_frames = 0;
    while hasFrame(VID.vid)
        num_frames = num_frames + 1;
        frames{num_frames} = readFrame(VID.vid);
    end
    VID.num_frames = num_frames;
catch
    try
        file_name_ = dir([path_to_raw '*.mp4']);
        file_name = [file_name_(1).name];
        file_path = [path_to_raw filesep];
        VID.vid = VideoReader([file_path file_name]);
        fprintf(['Loading: ', file_name, ' ... ']);
        num_frames = 0;
        while hasFrame(VID.vid)
            num_frames = num_frames + 1;
            frames{num_frames} = readFrame(VID.vid);
        end
        VID.num_frames = num_frames;
    catch
        file_name_ = dir([path_to_raw '*.avi']);
        file_name = [file_name_(1).name];
        file_path = [path_to_raw filesep];
        VID.vid = VideoReader([file_path file_name]);
        fprintf(['Loading: ', file_name, ' ... ']);
        num_frames = 0;
        while hasFrame(VID.vid)
            num_frames = num_frames + 1;
            frames{num_frames} = readFrame(VID.vid);
        end
        VID.num_frames = num_frames;
    end
end

VID.frames = frames;
VID.debug_figures = true;
VID.file_path = path_to_raw;
VID.file_name = file_name;
VID.fps_LB = fps_LB;
VID.fps_UB = fps_UB;
VID.fps_step = fps_step;
VID.path_to_analyzed_tongue = path_to_analyzed_tongue;
VID.path_to_analyzed_figs_tongue = path_to_analyzed_figs_tongue;
VID.height = VID.vid.Height;
VID.width = VID.vid.Width;
VID.duration = VID.vid.Duration;

fprintf(' --> Completed. \n')

%% LED Data
clearvars -except BEHAVE VID
fprintf('Extracting LED signal ...');

auto_find_rect = 1;
if auto_find_rect == 1
    rect_ = findrect(VID);
    rect_1 = rect_(1); rect_2 = rect_(2); rect_3 = rect_(3) - rect_(1); rect_4 = rect_(4) - rect_(2);
    rect = [rect_1 rect_2 rect_3 rect_4]; % [xmin, ymin, width, height]
else
    fig = figure;
    fig.WindowState = 'maximized';
    imshow(VID.frames{1});
    rect = getrect;
    close;
end
try
    num_frames = VID.num_frames;
    for counter_frame = 1 : num_frames
        current_frame = VID.frames{counter_frame};
        x = round(rect(1));y = round(rect(2)); w = round(rect(3));h = round(rect(4));
        mean_r(counter_frame) = round(mean(mean(current_frame(y:y+h, x:x+w, 1))));
        mean_g(counter_frame) = round(mean(mean(current_frame(y:y+h, x:x+w, 2))));
        mean_b(counter_frame) = round(mean(mean(current_frame(y:y+h, x:x+w, 3))));
    end
catch
    fprintf(' --> Error in extracting movie frames. \n');
end
close;

% Cluster pixel values and select threshold based on the midpoint of the mean of each cluster
mean_r = mean_r';
clust_mean_r = kmeans(mean_r, 2);
mean_1 = mean(mean_r(clust_mean_r == 1));
mean_2 = mean(mean_r(clust_mean_r == 2));
threshold = round((mean_1 + mean_2)/2);

if VID.debug_figures == 1
    figure;
    hold on
    histogram(mean_r);
    xline(mean_1, 'r', 'LineWidth', 2);
    xline(mean_2, 'b', 'LineWidth', 2);
    xline(threshold, 'k', 'LineWidth', 2);
    title(['Red pixel value distribution | Selected threshold: ' num2str(threshold)]);
    xlabel('Pixel value');
    ylabel('Count');
    ESN_Beautify_Plot(gcf, [20 10], 7)
end

for counter_frame = 1:length(mean_r)
    if (mean_r(counter_frame) < threshold)
        LED_trace(counter_frame) = 0;
    else
        LED_trace(counter_frame) = 1;
    end
end

LED_trace_rising = [diff(LED_trace) 0] > 0;
LED_trace_falling = [diff(LED_trace) 0] < 0;

LED(:,1) = mean_r;
LED(:,2) = mean_g;
LED(:,3) = mean_b;
LED(:,4) = LED_trace;
LED(:,5) = LED_trace_rising;
LED(:,6) = LED_trace_falling;

VID.LED_FPS = LED;
fprintf(' --> Completed. \n')


%% FPS
clearvars -except BEHAVE VID
%% Estimate FPS
fprintf(['Estimating FPS', ' ... \n']);

% fps_LB = input('Enter fps lower bound: ');
% fps_UB = input('Enter fps upper bound: ');
% fps_step = input('Enter fps step size: ');


fps = VID.fps_LB:VID.fps_step:VID.fps_UB;
%% CALCULATE FPS COARSE
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
    time_reference = ((1/fps(counter_FPS)) : (1/fps(counter_FPS)) : length(VID.LED_FPS(:,4))/fps(counter_FPS))';
    length_time = length(time_reference);
    time_LED_rise = time_reference([diff(VID.LED_FPS(:,4)) ; 0] > 0);

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
[~, ind_xcross_FPS_] = max(max_xcross_all);
estimated_FPS_ = fps(ind_xcross_FPS_);

if estimated_FPS_ > 103 || estimated_FPS_ < 90
    estimated_FPS_ = 100;
end
%% CALCULATE FPS FINE
%estimated_FPS_ = 40;
fps = estimated_FPS_ - 3 : 0.001 : estimated_FPS_ +3;
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
    time_reference = ((1/fps(counter_FPS)) : (1/fps(counter_FPS)) : length(VID.LED_FPS(:,4))/fps(counter_FPS))';
    length_time = length(time_reference);
    time_LED_rise = time_reference([diff(VID.LED_FPS(:,4)) ; 0] > 0);

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
[xcross_FPS, ind_xcross_FPS] = max(max_xcross_all);
estimated_FPS = fps(ind_xcross_FPS);

fprintf(['Estimated FPS (BEHAVE): ' num2str(estimated_FPS)])
fprintf(' --> Completed. \n')

%% PLOT FPS estimate report
if VID.flag_figure
    figure;
    hold on;
    plot(fps, max_xcross_all,'k')
    scatter(estimated_FPS, xcross_FPS,  'ro')
    xlabel('fps')
    ylabel('xcorr value')
    title(['Estimated FPS (BEHAVE): ' num2str(estimated_FPS)])
    ESN_Beautify_Plot(gcf, [20 10], 8)
    VID.FPS = estimated_FPS;
end
FPS = VID.FPS;
LED_FPS = VID.LED_FPS;
height = VID.height;
width = VID.width;
duration = VID.duration;
num_frames = VID.num_frames;


end

%% function findrect -  output: [x1, y1, x2, y2, frame of detection]
function rect = findrect(VID)
for i = 1 : VID.num_frames - 1

    frame = VID.frames{i};
    frame_next = VID.frames{i+1};

    rows = [];
    cols = [];

    for r = (size(frame, 1) - 50):size(frame, 1)
        for c = (size(frame, 2) - 100):size(frame, 2)
            red = frame(r, c, 1);
            red_next = frame_next(r, c, 1);

            green = frame(r, c, 2);
            green_next = frame_next(r, c, 2);

            blue = frame(r, c, 3);
            blue_next = frame_next(r, c, 3);

            if (green < 150) && (green_next > 200) && ((green_next - green) > 100)
                rows(end + 1) = r;
                cols(end + 1) = c;
            end
        end
    end

    if (~isempty(rows)) && (~isempty(cols))
        if ((rows(end) - rows(1)) >= 3) && ((cols(end) - cols(1)) >= 3)
            rect = [ceil(mean(cols)-1.5), ceil(mean(rows)-1.5), floor(mean(cols)+1.5), floor(mean(rows)+1.5), i];
            break
        end
    end
end
end
