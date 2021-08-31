%% Authors: PH, VL
%% Last modified: 04/28/2021
%
%%  DLC_behavior.m 
%   Analyzes the tracking data from the DLC network to
%   produce a .mat package and a .png report
%
%%  Inputs (optional):
%   var_name1 :  variable name ('file_path' / 'session_type' / 'with_nose_markers')
%   arg1      :	 value for the var_name1
%   var_name2 :  variable name ('file_path' / 'session_type' / 'with_nose_markers')
%   arg2      :	 value for the var_name2
%   var_name3 :  variable name ('file_path' / 'session_type' / 'with_nose_markers')
%   arg3      :	 value for the var_name3
%
%%  Note:
%   -   If 'file_path' is not specified, users will be prompted to select
%       a .csv file from a directory
%   -   The default value for 'session_type' is 2 if it is not specified
%   -   The default value for 'with_nose_markers' is 1 if it is not
%       specified
%
%%  Example use:
%
%   file_path = C:\Users\Vivian Looi\Johns Hopkins\Paul Hage - Mirza\2019-09-23_sorted_analyzed(eyeB&tongueRB&RP_1)_standard\2019-09-23_14-37-25\raw_data
%   DLC_behavior('file_path', file_path);
%

function DLC_behavior(var_name1, arg1, var_name2, arg2, var_name3, arg3)
    close all;
    warning('off','all');
    rng(1);
    
    if exist('var_name1', 'var') && exist('arg1', 'var')
        if strcmp(var_name1, 'file_path')
            file_path = char(arg1);
        elseif strcmp(var_name1, 'session_type')
            session_type = arg1;
        elseif strcmp(var_name1, 'with_nose_markers')
            with_nose_markers = arg1;
        end
    end
    
    if exist('var_name2', 'var') && exist('arg2', 'var')
        if strcmp(var_name2, 'file_path')
            file_path = char(arg2);
        elseif strcmp(var_name2, 'session_type')
            session_type = arg2;
        elseif strcmp(var_name2, 'with_nose_markers')
            with_nose_markers = arg2;
        end
    end
    
    if  exist('var_name3', 'var') && exist('arg3', 'var')
        if strcmp(var_name3, 'file_path')
            file_path = char(arg3);
        elseif strcmp(var_name3, 'session_type')
            session_type = arg3;
        elseif strcmp(var_name3, 'with_nose_markers')
            with_nose_markers = arg3;
        end
    end
    
    if ~exist('file_path','var')
        file_path = 0;
    elseif file_path(length(file_path)) ~= '\' || file_path(length(file_path)) ~= '/'
        file_path = [file_path '/'];
    end
    
    if ~exist('session_type','var')
        session_type = 2;
    end
    
    if ~exist('with_nose_markers','var')
        with_nose_markers = 1;
    end    
    
    % FUNCTION CALLS %
    DLC = load_extract_shift(file_path, session_type, with_nose_markers);
    DLC = evaluate(DLC);
    %DLC = quality_assurance(DLC);
    save_DLC(DLC);
end 

%% Load, extract, and shift DeepLabCut data
function DLC = load_extract_shift(file_path, session_type, with_nose_markers)
    fprintf('Loading, extracting, and shifting DLC data: ');
    
    % load csv file
    if length(file_path) == 1
        % file path not provided
        [file_name, file_path] = uigetfile([pwd filesep '*.csv'], 'Select CSV file');
    else
        % get .csv file name from the file path given
        dir_info = dir(file_path);
        for i = 1:length(dir_info)
            if regexp(dir_info(i).name, regexptranslate('wildcard', '*.csv'))
                file_name = dir_info(i).name;
            end
        end
    end    
    fprintf([file_name '... ']);
    DLC_data = xlsread([file_path filesep file_name]);

    %webcam
    %scal = 7/70;
    scal = 7/47;
    
    % get FPS
    dir_info = dir(strcat(file_path, '/../analyzed_data'));
    for i = 1:length(dir_info)
        if regexp(dir_info(i).name, regexptranslate('wildcard', '*.pdf'))
            name = erase(dir_info(i).name, '.pdf');
            num_str = regexprep(regexprep(name, '[^_.0-9(,)/]',''), ' \D* ',' ');
            parts = strsplit(num_str, '_');
            FPS = str2double(char(parts(length(parts))));
            break;
        end
    end  
    if ~exist('FPS','var')
        error("The FPS of the video hasn't been analyzed yet.");
    end

    for counter_frame = 1:length(DLC_data)
        time_vid(counter_frame) =  length(DLC_data(1:counter_frame))/FPS;
    end

    tip_tongue_x = scal*DLC_data(:,2);
    tip_tongue_y = scal*DLC_data(:,3);
    r_tongue_x = scal*DLC_data(:,5);
    r_tongue_y = scal*DLC_data(:,6);
    l_tongue_x = scal*DLC_data(:,8);
    l_tongue_y = scal*DLC_data(:,9);
    mid_tongue_x = scal*DLC_data(:,11);
    mid_tongue_y = scal*DLC_data(:,12);

    if session_type == 2 && with_nose_markers
        r_nose_x = scal*DLC_data(:,14);
        r_nose_y = scal*DLC_data(:,15);
        l_nose_x = scal*DLC_data(:,17);
        l_nose_y = scal*DLC_data(:,18);
        r_food_x = scal*DLC_data(:,20);
        r_food_y = scal*DLC_data(:,21);
        l_food_x = scal*DLC_data(:,23);
        l_food_y = scal*DLC_data(:,24);
        r_tube_r_x = scal*DLC_data(:,26);
        r_tube_r_y = scal*DLC_data(:,27);
        r_tube_l_x = scal*DLC_data(:,29);
        r_tube_l_y = scal*DLC_data(:,30);
        l_tube_r_x = scal*DLC_data(:,32);
        l_tube_r_y = scal*DLC_data(:,33);
        l_tube_l_x = scal*DLC_data(:,35);
        l_tube_l_y = scal*DLC_data(:,36);
    elseif session_type == 2 && with_nose_markers == 0
        r_food_x = scal*DLC_data(:,14);
        r_food_y = scal*DLC_data(:,15);
        l_food_x = scal*DLC_data(:,17);
        l_food_y = scal*DLC_data(:,18);
        r_tube_r_x = scal*DLC_data(:,20);
        r_tube_r_y = scal*DLC_data(:,21);
        r_tube_l_x = scal*DLC_data(:,23);
        r_tube_l_y = scal*DLC_data(:,24);
        l_tube_r_x = scal*DLC_data(:,26);
        l_tube_r_y = scal*DLC_data(:,27);
        l_tube_l_x = scal*DLC_data(:,29);
        l_tube_l_y = scal*DLC_data(:,30);
    elseif session_type == 1 && with_nose_markers
        r_nose_x = scal*DLC_data(:,14);
        r_nose_y = scal*DLC_data(:,15);
        l_nose_x = scal*DLC_data(:,17);
        l_nose_y = scal*DLC_data(:,18);
        r_food_x = scal*DLC_data(:,20);
        r_food_y = scal*DLC_data(:,21);
        r_tube_r_x = scal*DLC_data(:,23);
        r_tube_r_y = scal*DLC_data(:,24);
        r_tube_l_x = scal*DLC_data(:,26);
        r_tube_l_y = scal*DLC_data(:,27);
    elseif session_type == 1 && with_nose_markers == 0
        r_food_x = scal*DLC_data(:,14);
        r_food_y = scal*DLC_data(:,15);
        r_tube_r_x = scal*DLC_data(:,17);
        r_tube_r_y = scal*DLC_data(:,18);
        r_tube_l_x = scal*DLC_data(:,20);
        r_tube_l_y = scal*DLC_data(:,21);
    end

    [ind_cluster_x, cent_cluster_x] = kmeans(tip_tongue_x, 3);
    [ind_cluster_y, cent_cluster_y] = kmeans(tip_tongue_y, 3);
    x0 = min(cent_cluster_x);

    midpoint = (mean(r_nose_y) + mean(l_nose_y)) / 2;
    if abs(cent_cluster_y(1) - midpoint) < abs(cent_cluster_y(2) - midpoint) && abs(cent_cluster_y(1) - midpoint) < abs(cent_cluster_y(3) - midpoint)
        y0 = cent_cluster_y(1);
    elseif abs(cent_cluster_y(2) - midpoint) < abs(cent_cluster_y(1) - midpoint) && abs(cent_cluster_y(2) - midpoint) < abs(cent_cluster_y(3) - midpoint)
        y0 = cent_cluster_y(2);
    else
        y0 = cent_cluster_y(3);
    end

    tip_tongue_x = tip_tongue_x - x0;
    tip_tongue_y = tip_tongue_y - y0;
    r_tongue_x = r_tongue_x - x0;
    r_tongue_y = r_tongue_y - y0;
    l_tongue_x = l_tongue_x - x0;
    l_tongue_y = l_tongue_y - y0;
    mid_tongue_x = mid_tongue_x - x0;
    mid_tongue_y = mid_tongue_y - y0;
    if with_nose_markers
        r_nose_x = r_nose_x - x0;
        r_nose_y = r_nose_y - y0;
        l_nose_x = l_nose_x - x0;
        l_nose_y = l_nose_y - y0;
    end
    if session_type == 2
        l_food_x = l_food_x - x0;
        l_food_y = l_food_y - y0;
        l_tube_r_x = l_tube_r_x - x0;
        l_tube_r_y = l_tube_r_y - y0;
        l_tube_l_x = l_tube_l_x - x0;
        l_tube_l_y = l_tube_l_y - y0;
    end
    r_food_x = r_food_x - x0;
    r_food_y = r_food_y - y0;
    r_tube_r_x = r_tube_r_x - x0;
    r_tube_r_y = r_tube_r_y - y0;
    r_tube_l_x = r_tube_l_x - x0;
    r_tube_l_y = r_tube_l_y - y0;

    DLC.FILE.file_name = file_name;
    DLC.FILE.file_path = file_path;

    DLC.TIME.time_vid = time_vid;

    DLC.POINTS.tip_tongue_x = tip_tongue_x;
    DLC.POINTS.tip_tongue_y = tip_tongue_y;
    DLC.POINTS.r_tongue_x = r_tongue_x;
    DLC.POINTS.r_tongue_y = r_tongue_y;
    DLC.POINTS.l_tongue_x = l_tongue_x;
    DLC.POINTS.l_tongue_y = l_tongue_y;
    DLC.POINTS.mid_tongue_x = mid_tongue_x;
    DLC.POINTS.mid_tongue_y = mid_tongue_y;
    if with_nose_markers
        DLC.POINTS.r_nose_x = r_nose_x;
        DLC.POINTS.r_nose_y = r_nose_y;
        DLC.POINTS.l_nose_x = l_nose_x;
        DLC.POINTS.l_nose_y = l_nose_y;
    end
    if session_type == 2
        DLC.POINTS.l_food_x = l_food_x;
        DLC.POINTS.l_food_y = l_food_y;
        DLC.POINTS.l_tube_r_x = l_tube_r_x;
        DLC.POINTS.l_tube_r_y = l_tube_r_y;
        DLC.POINTS.l_tube_l_x = l_tube_l_x;
        DLC.POINTS.l_tube_l_y = l_tube_l_y;
    end
    DLC.POINTS.r_food_x = r_food_x;
    DLC.POINTS.r_food_y = r_food_y;
    DLC.POINTS.r_tube_r_x = r_tube_r_x;
    DLC.POINTS.r_tube_r_y = r_tube_r_y;
    DLC.POINTS.r_tube_l_x = r_tube_l_x;
    DLC.POINTS.r_tube_l_y = r_tube_l_y;
    DLC.POINTS.origin_x = x0;
    DLC.POINTS.origin_y = y0;

    DLC.FILE.session_type = session_type;
    DLC.FILE.with_nose_markers = with_nose_markers;
    DLC.FILE.FPS = FPS;
    clearvars -except DLC
    fprintf(' --> Completed. \n');
end

%% Suggest threshold
function [DLC, threshold_forward, threshold_backward] = suggest_threshold(DLC)
    session_type = DLC.FILE.session_type;
    with_nose_markers = DLC.FILE.with_nose_markers;
    
    % forward licks
    tip_tongue_x_forward = DLC.POINTS.tip_tongue_x(find(DLC.POINTS.tip_tongue_x>=0));
    tip_tongue_y_forward = DLC.POINTS.tip_tongue_y(find(DLC.POINTS.tip_tongue_x>=0));            
    d_tip_forward = sqrt(tip_tongue_x_forward.^2 + tip_tongue_y_forward.^2).';
    distribution = [];
            
   %find 2mm interval
    for i = 0:0.5:3.5
        distribution = [distribution; [i length(find(d_tip_forward < i))]];
        d_tip_forward(find(d_tip_forward < i)) = [];
    end

    for i = 1:length(distribution)
        distribution(i, 3) = var(distribution([i:length(distribution)],2));
    end

    for i = 2:length(distribution)
        if distribution(i, 3) == 0
            distribution(i, 4) = 0;
        else       
            distribution(i, 4) = distribution(i - 1, 3)/distribution(i, 3);
        end
    end

    [~, max_ind] = max(distribution(:,4));

    %filter data
    lower_bound = distribution(max_ind, 1);
    upper_bound = lower_bound + 1;

    %obtain distribution
    d_tip_forward = sqrt(tip_tongue_x_forward.^2 + tip_tongue_y_forward.^2).';
    d_tip_forward(find(d_tip_forward < lower_bound)) = [];
    distribution = [];
    for i = lower_bound:0.1:upper_bound
        distribution = [distribution; [i length(find(d_tip_forward < i))]];
        d_tip_forward(find(d_tip_forward < i)) = [];
    end

    for i = 1:length(distribution)
        distribution(i, 3) = var(distribution([i:length(distribution)],2));
    end

    for i = 2:length(distribution)
        if distribution(i, 3) == 0
            distribution(i, 4) = 0;
        else       
            distribution(i, 4) = distribution(i - 1, 3)/distribution(i, 3);
        end
    end

    [~,idx] = min(abs(distribution(:,4) - 1));
    threshold_forward = distribution(idx, 1) - 0.1;


    % backward licks
    tip_tongue_x_backward = DLC.POINTS.tip_tongue_x(find(DLC.POINTS.tip_tongue_x<0));
    tip_tongue_y_backward = DLC.POINTS.tip_tongue_y(find(DLC.POINTS.tip_tongue_x<0));            
    d_tip_backward = sqrt(tip_tongue_x_backward.^2 + tip_tongue_y_backward.^2).';
    distribution = [];

   %find 2mm interval
    for i = 0:0.5:3.5
        distribution = [distribution; [i length(find(d_tip_backward < i))]];
        d_tip_backward(find(d_tip_backward < i)) = [];
    end

    for i = 1:length(distribution)
        distribution(i, 3) = var(distribution([i:length(distribution)],2));
    end

    for i = 2:length(distribution)
        if distribution(i, 3) == 0
            distribution(i, 4) = 0;
        else       
            distribution(i, 4) = distribution(i - 1, 3)/distribution(i, 3);
        end
    end

    [~, max_ind] = max(distribution(:,4));

    %filter data
    lower_bound = distribution(max_ind, 1);
    upper_bound = lower_bound + 1;

    %obtain distribution
    d_tip_backward = sqrt(tip_tongue_x_backward.^2 + tip_tongue_y_backward.^2).';
    d_tip_backward(find(d_tip_backward < lower_bound)) = [];
    distribution = [];
    for i = lower_bound:0.1:upper_bound
        distribution = [distribution; [i length(find(d_tip_backward < i))]];
        d_tip_backward(find(d_tip_backward < i)) = [];
    end

    for i = 1:length(distribution)
        distribution(i, 3) = var(distribution([i:length(distribution)],2));
    end

    for i = 2:length(distribution)
        if distribution(i, 3) == 0
            distribution(i, 4) = 0;
        else       
            distribution(i, 4) = distribution(i - 1, 3)/distribution(i, 3);
        end
    end

    [~,idx] = min(abs(distribution(:,4) - 1));
    if distribution(idx, 1) == 0
        threshold_backward = 0;
    else
        threshold_backward = distribution(idx, 1) - 0.1;
    end

    d_tip = sqrt(DLC.POINTS.tip_tongue_x.^2 + DLC.POINTS.tip_tongue_y.^2).';
    idx = find(d_tip <= threshold_forward);


    %Nan for frames with tongue out before the first lick onset
    if length(idx) ~=0 && idx(1) ~= 1
        DLC.POINTS.tip_tongue_x(1:idx(1)-1) = nan;
        DLC.POINTS.tip_tongue_y(1:idx(1)-1) = nan;
        DLC.POINTS.r_tongue_x(1:idx(1)-1) = nan;
        DLC.POINTS.r_tongue_y(1:idx(1)-1) = nan;
        DLC.POINTS.l_tongue_x(1:idx(1)-1) = nan;
        DLC.POINTS.l_tongue_y(1:idx(1)-1) = nan;
        DLC.POINTS.mid_tongue_x(1:idx(1)-1) = nan;
        DLC.POINTS.mid_tongue_y(1:idx(1)-1) = nan;
        if with_nose_markers
            DLC.POINTS.r_nose_x(1:idx(1)-1) = nan;
            DLC.POINTS.r_nose_y(1:idx(1)-1) = nan;
            DLC.POINTS.l_nose_x(1:idx(1)-1) = nan;
            DLC.POINTS.l_nose_y(1:idx(1)-1) = nan;
        end
        DLC.POINTS.r_food_x(1:idx(1)-1) = nan;
        DLC.POINTS.r_food_y(1:idx(1)-1) = nan;
        DLC.POINTS.r_tube_r_x(1:idx(1)-1) = nan;
        DLC.POINTS.r_tube_r_y(1:idx(1)-1) = nan;
        DLC.POINTS.r_tube_l_x(1:idx(1)-1) = nan;
        DLC.POINTS.r_tube_l_y(1:idx(1)-1) = nan;
        if session_type == 2
            DLC.POINTS.l_food_x(1:idx(1)-1) = nan;
            DLC.POINTS.l_food_y(1:idx(1)-1) = nan;
            DLC.POINTS.l_tube_r_x(1:idx(1)-1) = nan;
            DLC.POINTS.l_tube_r_y(1:idx(1)-1) = nan;
            DLC.POINTS.l_tube_l_x(1:idx(1)-1) = nan;
            DLC.POINTS.l_tube_l_y(1:idx(1)-1) = nan;
        end
    end
end

%% Detect licks and bouts
function DLC = detect_licks_and_bouts(DLC, threshold_forward, threshold_backward)
    fprintf('Detecting licks and bouts ... ');
    
    tip_tongue_x = DLC.POINTS.tip_tongue_x;
    tip_tongue_y = DLC.POINTS.tip_tongue_y;
    
    d_tip = sqrt(tip_tongue_x.^2 + tip_tongue_y.^2);
    d_tip_forward = zeros(size(tip_tongue_x));
    d_tip_backward = zeros(size(tip_tongue_x));
    d_tip_forward(find(tip_tongue_x > 0)) = sqrt(tip_tongue_x(find(tip_tongue_x > 0)).^2 + tip_tongue_y(find(tip_tongue_x > 0)).^2);
    d_tip_backward(find(tip_tongue_x < 0)) = sqrt(tip_tongue_x(find(tip_tongue_x < 0)).^2 + tip_tongue_y(find(tip_tongue_x < 0)).^2);
    time_vid = DLC.TIME.time_vid';

    tongue_out = d_tip_forward > threshold_forward;
    backward_ind = find(d_tip_backward > threshold_backward);
    tongue_out(backward_ind) = 1;
    if ~isempty(backward_ind)
        for i = 1:length(backward_ind)
            j = backward_ind(i) - 1;
            while j ~= 0 && tongue_out(j) == 0 
                tongue_out(j) = 1;
                j = j - 1;
            end
        end
    end

    ind_lick_onset = find(diff(tongue_out) == 1);
    ind_lick_offset = find(diff(tongue_out) == -1) +1;
    if length(ind_lick_offset) < length(ind_lick_onset)
        ind_lick_offset = [ind_lick_offset; length(d_tip)];
    end
    
    for i = 1:length(ind_lick_onset)
        ind = ind_lick_onset(i);
        while ind ~= 1 && d_tip(ind - 1) < d_tip(ind)
            ind = ind - 1;
        end
        ind_lick_onset(i) = ind;
    end
    
    for i = 1:length(ind_lick_offset)
        ind = ind_lick_offset(i);
        while ind ~= length(d_tip) && d_tip(ind + 1) < d_tip(ind)
            ind = ind + 1;
        end
        ind_lick_offset(i) = ind;
    end
    
    num_lick = length(ind_lick_onset);

    time_lick_onset = time_vid(ind_lick_onset);
    time_lick_offset =  time_vid(ind_lick_offset);
    time_lick_duration = time_lick_offset - time_lick_onset;

    ind_lick_onset_str_bout_ = [1; find(diff(ind_lick_onset) >100) + 1];
    ind_lick_onset_str_bout = ind_lick_onset(ind_lick_onset_str_bout_);
    ind_lick_onset_end_bout_ = [find(diff(ind_lick_onset) >100); length(ind_lick_onset)];
    ind_lick_onset_end_bout = ind_lick_onset(ind_lick_onset_end_bout_);

    time_lick_onset_str_bout = time_vid(ind_lick_onset_str_bout);
    time_lick_onset_end_bout = time_vid(ind_lick_onset_end_bout);
    time_bout_duration = time_lick_onset_end_bout - time_lick_onset_str_bout;

    num_lick_bout = (ind_lick_onset_end_bout_ - ind_lick_onset_str_bout_) + 1;

    num_bout = length(ind_lick_onset_str_bout);

    DLC.KINEMATIC.d_tip = d_tip;
    DLC.IND.num_lick = num_lick;
    DLC.IND.num_bout = num_bout;
    DLC.IND.num_lick_bout = num_lick_bout;
    DLC.IND.ind_lick_onset = ind_lick_onset;
    DLC.IND.ind_lick_offset = ind_lick_offset;
    DLC.IND.ind_lick_onset_str_bout = ind_lick_onset_str_bout;
    DLC.IND.ind_lick_onset_end_bout = ind_lick_onset_end_bout;
    DLC.TIME.time_lick_onset_str_bout = time_lick_onset_str_bout;
    DLC.TIME.time_lick_onset_end_bout = time_lick_onset_end_bout;
    DLC.TIME.time_lick_onset = time_lick_onset;
    DLC.TIME.time_lick_offset = time_lick_offset;
    DLC.TIME.time_lick_duration = time_lick_duration;
    DLC.TIME.time_bout_duration = time_bout_duration;

    figure
    hold on;
    plot(time_vid,d_tip,'k');
    plot(time_lick_onset_str_bout,21,'*g');
    plot(time_lick_onset_end_bout,21,'*r');
    plot(time_lick_onset, 20, '.k');
    xlabel('Time (s)');
    ylabel('Distance (mm)');
    % ESN_Beautify_Plot

    fprintf(' --> Completed. \n')
end

%% Calculate lick kinematics: d, v, a, angle, ILI(inter lick interval), ILR(instantaneous lick rate = 1/ILI)
function DLC = calculate_lick_kinematics(DLC)
    fprintf('Calculating lick kinematics ...');

    session_type = DLC.FILE.session_type;
    with_nose_markers = DLC.FILE.with_nose_markers;
    d_tip = DLC.KINEMATIC.d_tip;
    time_lick_duration = DLC.TIME.time_lick_duration;
    time_bout_duration = DLC.TIME.time_bout_duration;

    num_lick = DLC.IND.num_lick;
    num_bout = DLC.IND.num_bout;
    num_lick_bout = DLC.IND.num_lick_bout;
    time_vid = DLC.TIME.time_vid';
    ind_lick_onset = DLC.IND.ind_lick_onset;
    ind_lick_offset = DLC.IND.ind_lick_offset;
    ind_lick_onset_str_bout = DLC.IND.ind_lick_onset_str_bout;
    ind_lick_onset_end_bout = DLC.IND.ind_lick_onset_end_bout;
    tip_tongue_x = DLC.POINTS.tip_tongue_x;
    tip_tongue_y = DLC.POINTS.tip_tongue_y;
    mid_tongue_x = DLC.POINTS.mid_tongue_x;
    mid_tongue_y = DLC.POINTS.mid_tongue_y;
    midtip_y = tip_tongue_y - mid_tongue_y;
    midtip_x = tip_tongue_x - mid_tongue_x;

    time_lick_onset = DLC.TIME.time_lick_onset;
    time_lick_offset = DLC.TIME.time_lick_offset;

    d_tip(1:ind_lick_onset(1)) = 0;
    for i = 1:length(ind_lick_onset) - 1
        d_tip(ind_lick_offset(i):ind_lick_onset(i+1)) = 0;
    end
    
    for counter_lick =  1 : 1 : num_lick
        inds_ = ind_lick_onset(counter_lick):ind_lick_offset(counter_lick);
        for counter_inds_ = 1 : 1: length(inds_)
            d_lick_(counter_inds_) = d_tip(inds_(counter_inds_));
            d_lick(counter_lick, counter_inds_) = d_lick_(counter_inds_);
        end
        [d_lick_max_local, ind_d_lick_max_local] = max(d_tip(inds_));
        ind_d_lick_max(counter_lick,1) = ind_lick_onset(counter_lick) + ind_d_lick_max_local - 1;
        d_lick_max(counter_lick,1) =  d_lick_max_local;
    end

    for counter_lick =  1 : 1 : num_lick
        inds_ = ind_lick_onset(counter_lick):ind_lick_offset(counter_lick);
        v_lick_ = [0 (diff(d_tip(inds_))./ time_vid(1))'];
        v_lick(counter_lick, 1:length(inds_)) = v_lick_;
        [v_lick_max_local, ind_v_lick_max_local] = max(v_lick(counter_lick,:));
        ind_v_lick_max(counter_lick,1) = ind_lick_onset(counter_lick) + ind_v_lick_max_local - 1;
        v_lick_max(counter_lick,1) =  v_lick_max_local;

        [v_lick_min_local, ind_v_lick_min_local] = min(v_lick(counter_lick,:));
        ind_v_lick_min(counter_lick,1) = ind_lick_onset(counter_lick) + ind_v_lick_min_local - 1;
        v_lick_min(counter_lick,1) =  v_lick_min_local;
    end

    for counter_lick =  1 : 1 : num_lick
        inds_ = ind_lick_onset(counter_lick):ind_lick_offset(counter_lick) ;
        a_lick_ = [0 0 (diff(d_tip(inds_),2)./ time_vid(1))'];
        a_lick(counter_lick, 1:length(inds_)) = a_lick_;
        [a_lick_max_local, ind_a_lick_max_local] = max(a_lick(counter_lick,:));
        ind_a_lick_max(counter_lick,1) = ind_lick_onset(counter_lick) + ind_a_lick_max_local - 1;
        a_lick_max(counter_lick,1) =  a_lick_max_local;
    end

    for counter_lick = 1 : 1 : num_lick
        inds_ = ind_lick_onset(counter_lick):ind_lick_offset(counter_lick);
        for counter_inds_ = 1 : 1: length(inds_)
            angle_lick_(counter_inds_) = rad2deg(atan2(midtip_y(inds_(counter_inds_)),...
                midtip_x(inds_(counter_inds_))));
            angle_lick(counter_lick, counter_inds_) = angle_lick_(counter_inds_);
        end
    end

    for counter_bout = 1 : 1 : num_bout
        ILI_bout(counter_bout,1) = mean(diff(time_lick_onset(find(ind_lick_onset >= ind_lick_onset_str_bout(counter_bout) & ind_lick_onset <= ind_lick_onset_end_bout(counter_bout)))));
    end

    v_tip = ([0 (diff(d_tip)./ time_vid(1))'])';
    a_tip = ([0 0 (diff(d_tip, 2)./ time_vid(1))'])';
    angle_midtip = rad2deg(atan2(midtip_y, midtip_x));
    angle_lick_max = angle_midtip(ind_d_lick_max);

    ILR_bout = 1./ILI_bout;
    ILI_bout(find(ILR_bout==Inf)) = nan;
    ILR_bout(find(ILR_bout==Inf)) = nan;
    
    ILI_lick = diff(time_lick_onset);
    ILI_lick = [ILI_lick; nan];
    ILR_lick = 1./ILI_lick;

    time_d_lick_max_abs = time_vid(ind_d_lick_max);
    time_d_lick_max_rel = (time_d_lick_max_abs - time_lick_onset);
    time_v_lick_max_abs = time_vid(ind_v_lick_max);
    time_v_lick_max_rel = (time_v_lick_max_abs - time_lick_onset);
    time_v_lick_min_abs = time_vid(ind_v_lick_min);
    time_v_lick_min_rel = (time_v_lick_min_abs - time_lick_onset);

    DLC.KINEMATIC.d_tip = d_tip;
    DLC.KINEMATIC.d_lick = d_lick;
    DLC.KINEMATIC.d_lick_max = d_lick_max;
    DLC.IND.ind_d_lick_max = ind_d_lick_max;
    DLC.KINEMATIC.v_tip = v_tip;
    DLC.KINEMATIC.a_tip = a_tip;
    DLC.KINEMATIC.angle_lick = angle_lick;
    DLC.KINEMATIC.angle_midtip = angle_midtip;
    DLC.KINEMATIC.angle_lick_max = angle_lick_max;
    DLC.KINEMATIC.v_lick = v_lick;
    DLC.KINEMATIC.v_lick_max = v_lick_max;
    DLC.IND.ind_v_lick_max = ind_v_lick_max;
    DLC.KINEMATIC.v_lick_min = v_lick_min;
    DLC.IND.ind_v_lick_min = ind_v_lick_min;
    DLC.KINEMATIC.a_lick = a_lick;
    DLC.KINEMATIC.a_lick_max = a_lick_max;
    DLC.IND.ind_a_lick_max = ind_a_lick_max;
    DLC.TIME.time_d_lick_max_abs = time_d_lick_max_abs;
    DLC.TIME.time_d_lick_max_rel = time_d_lick_max_rel;
    DLC.TIME.time_v_lick_max_abs = time_v_lick_max_abs;
    DLC.TIME.time_v_lick_max_rel = time_v_lick_max_rel;
    DLC.TIME.time_v_lick_min_abs = time_v_lick_min_abs;
    DLC.TIME.time_v_lick_min_rel = time_v_lick_min_rel;
    DLC.KINEMATIC.ILI_bout = ILI_bout;
    DLC.KINEMATIC.ILR_bout = ILR_bout;
    DLC.KINEMATIC.ILI_lick = ILI_lick;
    DLC.KINEMATIC.ILR_lick = ILR_lick;

    figure;
    hold on;
    subplot(1,6,1)
    boxplot(d_lick_max, {'D max'});
    ylabel('mm');
    ylim([0 25])
    subplot(1,6,2)
    boxplot(v_lick_max, {'V max'});
    ylabel('mm/s');
    ylim([0 650])
    subplot(1,6,3)
    boxplot(v_lick_min, {'V min'});
    ylabel('mm/s');
    ylim([-650 0])
    subplot(1,6,4)
    boxplot(time_lick_duration*1000, {'Duration'});
    ylabel('ms');
    ylim([0 650])
    subplot(1,6,5)
    boxplot(ILI_bout*1000, {'ILI'});
    ylabel('ms');
    ylim([0 650])
    subplot(1,6,6)
    boxplot(ILR_bout, {'ILR'});
    ylabel('hz');
    ylim([0 6])
    % ESN_Beautify_Plot

    fprintf(' --> Completed. \n')
end

%% Geometrization
function DLC = geometrization(DLC)
    fprintf('Geometrizing frames ...');

    session_type = DLC.FILE.session_type;
    with_nose_markers = DLC.FILE.with_nose_markers;
    tip_tongue_x = DLC.POINTS.tip_tongue_x;
    tip_tongue_y = DLC.POINTS.tip_tongue_y;
    r_tongue_x = DLC.POINTS.r_tongue_x;
    r_tongue_y = DLC.POINTS.r_tongue_y;
    l_tongue_x = DLC.POINTS.l_tongue_x;
    l_tongue_y = DLC.POINTS.l_tongue_y;
    mid_tongue_x = DLC.POINTS.mid_tongue_x;
    mid_tongue_y = DLC.POINTS.mid_tongue_y;
    if session_type == 2
        l_tube_r_x = DLC.POINTS.l_tube_r_x;
        l_tube_r_y = DLC.POINTS.l_tube_r_y;
        l_tube_l_x = DLC.POINTS.l_tube_l_x;
        l_tube_l_y = DLC.POINTS.l_tube_l_y;
        l_food_x = DLC.POINTS.l_food_x;
        l_food_y = DLC.POINTS.l_food_y;
    end
    r_tube_r_x = DLC.POINTS.r_tube_r_x;
    r_tube_r_y = DLC.POINTS.r_tube_r_y;
    r_tube_l_x = DLC.POINTS.r_tube_l_x;
    r_tube_l_y = DLC.POINTS.r_tube_l_y;
    r_food_x = DLC.POINTS.r_food_x;
    r_food_y = DLC.POINTS.r_food_y;
    num_lick = DLC.IND.num_lick;
    ind_d_lick_max = DLC.IND.ind_d_lick_max;


    for counter_lick = 1:1:num_lick
        counter_frame = ind_d_lick_max(counter_lick);

        geo_tongue = polyshape([tip_tongue_x(counter_frame),l_tongue_x(counter_frame),...
            mid_tongue_x(counter_frame),r_tongue_x(counter_frame)],[tip_tongue_y(counter_frame),...
            l_tongue_y(counter_frame),mid_tongue_y(counter_frame),r_tongue_y(counter_frame)]);

        geo_r_tube_empty = polyshape([r_food_x(counter_frame),r_tube_r_x(counter_frame),...
            r_tube_r_x(counter_frame),r_tube_l_x(counter_frame), r_tube_l_x(counter_frame)],[r_food_y(counter_frame),...
            r_food_y(counter_frame), r_tube_r_y(counter_frame),r_tube_l_y(counter_frame), r_food_y(counter_frame)]);

        geo_r_tube_full = polyshape([r_tube_r_x(counter_frame), r_tube_r_x(counter_frame),...
            r_tube_l_x(counter_frame), r_tube_l_x(counter_frame), r_food_x(counter_frame)],...
            [r_food_y(counter_frame), r_tube_r_y(counter_frame) + abs(max(r_food_y)-min(r_food_y)), ...
            r_tube_l_y(counter_frame) + abs(max(r_food_y)-min(r_food_y)), r_food_y(counter_frame), r_food_y(counter_frame)]);

        geo_inter_tongue_r_tube_empty = intersect(geo_tongue, geo_r_tube_empty);

        geo_inter_tongue_r_tube_full = intersect(geo_tongue, geo_r_tube_full);

        geo_all = [geo_tongue geo_r_tube_empty geo_r_tube_full ...
            geo_inter_tongue_r_tube_empty geo_inter_tongue_r_tube_full];

        if session_type == 2
            geo_l_tube_empty = polyshape([l_food_x(counter_frame),l_tube_r_x(counter_frame),...
                l_tube_r_x(counter_frame),l_tube_l_x(counter_frame), l_tube_l_x(counter_frame)],[l_food_y(counter_frame),...
                l_food_y(counter_frame), l_tube_r_y(counter_frame),l_tube_l_y(counter_frame), l_food_y(counter_frame)]);

            geo_l_tube_full = polyshape([l_tube_r_x(counter_frame), l_tube_r_x(counter_frame),...
                l_tube_l_x(counter_frame), l_tube_l_x(counter_frame), l_food_x(counter_frame)],...
                [l_food_y(counter_frame), l_tube_r_y(counter_frame) - abs(max(l_food_y)-min(l_food_y)), ...
                l_tube_l_y(counter_frame) - abs(max(l_food_y)-min(l_food_y)), l_food_y(counter_frame), l_food_y(counter_frame)]);

            geo_inter_tongue_l_tube_empty = intersect(geo_tongue, geo_l_tube_empty);

            geo_inter_tongue_l_tube_full = intersect(geo_tongue, geo_l_tube_full);

            geo_all = [geo_tongue geo_r_tube_empty geo_r_tube_full geo_l_tube_empty geo_l_tube_full...
                geo_inter_tongue_r_tube_empty geo_inter_tongue_r_tube_full geo_inter_tongue_l_tube_empty...
                geo_inter_tongue_l_tube_full ];
        end

        [cent_tongue_r_tube_empty_x(counter_lick, 1), cent_tongue_r_tube_empty_y(counter_lick, 1) ] = centroid(geo_inter_tongue_r_tube_empty);
        [cent_tongue_r_tube_full_x(counter_lick, 1), cent_tongue_r_tube_full_y(counter_lick, 1)] = centroid(geo_inter_tongue_r_tube_full);
        bool_overlaps_all = overlaps(geo_all);
        bool_tongue_r_tube_empty(counter_lick, 1) = bool_overlaps_all(1,2);
        bool_tongue_r_tube_full(counter_lick, 1) = bool_overlaps_all(1,3);
        area_tongue(counter_lick, 1) = area(geo_tongue);
        area_r_tube_empty(counter_lick, 1) = area(geo_r_tube_empty);
        area_r_tube_full(counter_lick, 1) = area(geo_r_tube_full);
        area_inter_tongue_r_tube_empty(counter_lick, 1) = area(geo_inter_tongue_r_tube_empty);
        area_inter_tongue_r_tube_full(counter_lick, 1) = area(geo_inter_tongue_r_tube_full);

        if session_type == 2
            [cent_tongue_l_tube_empty_x(counter_lick, 1),cent_tongue_l_tube_empty_y(counter_lick, 1)]  = centroid(geo_inter_tongue_l_tube_empty);
            [cent_tongue_l_tube_full_x(counter_lick, 1),cent_tongue_l_tube_full_y(counter_lick, 1)] = centroid(geo_inter_tongue_l_tube_full);
            bool_tongue_l_tube_empty(counter_lick, 1)= bool_overlaps_all(1,4);
            bool_tongue_l_tube_full(counter_lick, 1) = bool_overlaps_all(1,5);
            area_l_tube_empty(counter_lick, 1) = area(geo_l_tube_empty);
            area_l_tube_full(counter_lick, 1) = area(geo_l_tube_full);
            area_inter_tongue_l_tube_empty(counter_lick, 1) = area(geo_inter_tongue_l_tube_empty);
            area_inter_tongue_l_tube_full(counter_lick, 1) = area(geo_inter_tongue_l_tube_full);
        end
    end

    DLC.GEO.area_tongue = area_tongue;
    DLC.GEO.area_r_tube_empty = area_r_tube_empty;
    DLC.GEO.area_r_tube_full = area_r_tube_full;
    DLC.GEO.area_inter_tongue_r_tube_empty = area_inter_tongue_r_tube_empty;
    DLC.GEO.area_inter_tongue_r_tube_full = area_inter_tongue_r_tube_full;
    DLC.GEO.bool_overlaps_all = bool_overlaps_all;
    DLC.GEO.bool_tongue_r_tube_empty = bool_tongue_r_tube_empty;
    DLC.GEO.bool_tongue_r_tube_full = bool_tongue_r_tube_full;
    DLC.GEO.cent_tongue_r_tube_empty_x = cent_tongue_r_tube_empty_x;
    DLC.GEO.cent_tongue_r_tube_empty_y = cent_tongue_r_tube_empty_y;
    DLC.GEO.cent_tongue_r_tube_full_x = cent_tongue_r_tube_full_x;
    DLC.GEO.cent_tongue_r_tube_full_y = cent_tongue_r_tube_full_y;

    if session_type == 2
        DLC.GEO.area_l_tube_empty = area_l_tube_empty;
        DLC.GEO.area_l_tube_full = area_l_tube_full;
        DLC.GEO.area_inter_tongue_l_tube_empty = area_inter_tongue_l_tube_empty;
        DLC.GEO.area_inter_tongue_l_tube_full = area_inter_tongue_l_tube_full;
        DLC.GEO.bool_tongue_l_tube_empty = bool_tongue_l_tube_empty;
        DLC.GEO.bool_tongue_l_tube_full = bool_tongue_l_tube_full;
        DLC.GEO.cent_tongue_l_tube_empty_x = cent_tongue_l_tube_empty_x;
        DLC.GEO.cent_tongue_l_tube_empty_y = cent_tongue_l_tube_empty_y;
        DLC.GEO.cent_tongue_l_tube_full_x = cent_tongue_l_tube_full_x;
        DLC.GEO.cent_tongue_l_tube_full_y = cent_tongue_l_tube_full_y;
    end

    fprintf(' --> Completed. \n')
end

%% Sort licks
function DLC = sort_licks(DLC)
    fprintf('Sorting licks ...');

    session_type = DLC.FILE.session_type;
    with_nose_markers = DLC.FILE.with_nose_markers;

    if (session_type == 2)

        bool_tongue_r_tube_empty = DLC.GEO.bool_tongue_r_tube_empty;
        bool_tongue_r_tube_full = DLC.GEO.bool_tongue_r_tube_full;
        bool_tongue_l_tube_empty = DLC.GEO.bool_tongue_l_tube_empty;
        bool_tongue_l_tube_full = DLC.GEO.bool_tongue_l_tube_full;

        tip_tongue_x = DLC.POINTS.tip_tongue_x;
        tip_tongue_y = DLC.POINTS.tip_tongue_y;

        if with_nose_markers
            r_nose_x =  DLC.POINTS.r_nose_x;
            r_nose_y = DLC.POINTS.r_nose_y;
            l_nose_x =  DLC.POINTS.l_nose_x;
            l_nose_y = DLC.POINTS.l_nose_y;
        end

        r_tube_r_x = DLC.POINTS.r_tube_r_x;
        r_tube_r_y = DLC.POINTS.r_tube_r_y;
        r_tube_l_x = DLC.POINTS.r_tube_l_x;
        r_tube_l_y = DLC.POINTS.r_tube_l_y;
        r_food_y = DLC.POINTS.r_food_y;

        l_tube_r_x = DLC.POINTS.l_tube_r_x;
        l_tube_r_y = DLC.POINTS.l_tube_r_y;
        l_tube_l_x = DLC.POINTS.l_tube_l_x;
        l_tube_l_y = DLC.POINTS.l_tube_l_y;
        l_food_y = DLC.POINTS.l_food_y;

        ind_d_lick_max = DLC.IND.ind_d_lick_max;
        ind_lick_onset = DLC.IND.ind_lick_onset;
        ind_lick_offset = DLC.IND.ind_lick_offset;
        time_vid = DLC.TIME.time_vid';

        ind_grooming_lick = find(bool_tongue_r_tube_empty == 0 & bool_tongue_l_tube_empty == 0 ...
            & bool_tongue_r_tube_full == 0 & bool_tongue_l_tube_full == 0 &...
            tip_tongue_y(ind_d_lick_max) < (r_tube_r_y(ind_d_lick_max)) & ...
            tip_tongue_y(ind_d_lick_max) > (l_tube_r_y(ind_d_lick_max)));
        is_grooming_lick = false(size(ind_lick_onset));
        is_grooming_lick(ind_grooming_lick) = true;
        
        ind_r_reward_inner_tube_lick = find(bool_tongue_r_tube_full == 1 & ...
            tip_tongue_x(ind_d_lick_max) < r_tube_r_x(ind_d_lick_max) & tip_tongue_x(ind_d_lick_max) > r_tube_l_x(ind_d_lick_max));
        is_r_reward_inner_tube_lick = false(size(ind_lick_onset));
        is_r_reward_inner_tube_lick(ind_r_reward_inner_tube_lick) = true;
        
        ind_r_reward_outer_tube_lick = find(bool_tongue_r_tube_full == 1 &...
            (tip_tongue_x(ind_d_lick_max) > r_tube_r_x(ind_d_lick_max) | tip_tongue_x(ind_d_lick_max) < r_tube_l_x(ind_d_lick_max)));
        is_r_reward_outer_tube_lick = false(size(ind_lick_onset));
        is_r_reward_outer_tube_lick(ind_r_reward_outer_tube_lick) = true;

        ind_r_noreward_inner_tube_lick = find(bool_tongue_r_tube_empty == 1 & bool_tongue_r_tube_full == 0 &...
            tip_tongue_x(ind_d_lick_max) < r_tube_r_x(ind_d_lick_max) & tip_tongue_x(ind_d_lick_max) > r_tube_l_x(ind_d_lick_max));
        is_r_noreward_inner_tube_lick = false(size(ind_lick_onset));
        is_r_noreward_inner_tube_lick(ind_r_noreward_inner_tube_lick) = true;

        ind_r_noreward_outer_tube_lick = find(bool_tongue_r_tube_full == 0 & tip_tongue_y(ind_d_lick_max) > r_tube_r_y(ind_d_lick_max) & ...
            (tip_tongue_x(ind_d_lick_max) > r_tube_r_x(ind_d_lick_max) | tip_tongue_x(ind_d_lick_max) < r_tube_l_x(ind_d_lick_max)));
        is_r_noreward_outer_tube_lick = false(size(ind_lick_onset));
        is_r_noreward_outer_tube_lick(ind_r_noreward_outer_tube_lick) = true;

        ind_lick_onset_r_reward = ind_lick_onset(sort([ind_r_reward_inner_tube_lick;ind_r_reward_outer_tube_lick;ind_r_noreward_inner_tube_lick;ind_r_noreward_outer_tube_lick]));
        time_lick_onset_r_reward = time_vid(ind_lick_onset_r_reward);

        ind_lick_offset_r_reward = ind_lick_offset(sort([ind_r_reward_inner_tube_lick;ind_r_reward_outer_tube_lick;ind_r_noreward_inner_tube_lick;ind_r_noreward_outer_tube_lick]));
        time_lick_offset_r_reward = time_vid(ind_lick_offset_r_reward);

        ind_l_reward_inner_tube_lick = find(bool_tongue_l_tube_full == 1 & ...
            tip_tongue_x(ind_d_lick_max) < l_tube_r_x(ind_d_lick_max) & tip_tongue_x(ind_d_lick_max) > l_tube_l_x(ind_d_lick_max));
        is_l_reward_inner_tube_lick = false(size(ind_lick_onset));
        is_l_reward_inner_tube_lick(ind_l_reward_inner_tube_lick) = true;

        ind_l_reward_outer_tube_lick = find(bool_tongue_l_tube_full == 1 & ...
            (tip_tongue_x(ind_d_lick_max) > l_tube_r_x(ind_d_lick_max) | tip_tongue_x(ind_d_lick_max) < l_tube_l_x(ind_d_lick_max)));
        is_l_reward_outer_tube_lick = false(size(ind_lick_onset));
        is_l_reward_outer_tube_lick(ind_l_reward_outer_tube_lick) = true;

        ind_l_noreward_inner_tube_lick = find(bool_tongue_l_tube_empty == 1 & bool_tongue_l_tube_full == 0 & ...
            tip_tongue_x(ind_d_lick_max) < l_tube_r_x(ind_d_lick_max) & tip_tongue_x(ind_d_lick_max) > l_tube_l_x(ind_d_lick_max));
        is_l_noreward_inner_tube_lick = false(size(ind_lick_onset));
        is_l_noreward_inner_tube_lick(ind_l_noreward_inner_tube_lick) = true;

        ind_l_noreward_outer_tube_lick = find(bool_tongue_l_tube_full == 0 & tip_tongue_y(ind_d_lick_max) < l_tube_r_y(ind_d_lick_max) &...
            (tip_tongue_x(ind_d_lick_max) > l_tube_r_x(ind_d_lick_max) | tip_tongue_x(ind_d_lick_max) < l_tube_l_x(ind_d_lick_max)));
        is_l_noreward_outer_tube_lick = false(size(ind_lick_onset));
        is_l_noreward_outer_tube_lick(ind_l_noreward_outer_tube_lick) = true;

        
        %     ind_l_noreward_inner_tube_lick(ismember(ind_l_noreward_inner_tube_lick, ind_l_reward_inner_tube_lick)) = [];
        %     ind_r_noreward_inner_tube_lick(ismember(ind_r_noreward_inner_tube_lick, ind_r_reward_inner_tube_lick)) = [];

        ind_lick_onset_l_reward = ind_lick_onset(sort([ind_l_reward_inner_tube_lick;ind_l_reward_outer_tube_lick;ind_l_noreward_inner_tube_lick;ind_l_noreward_outer_tube_lick]));
        time_lick_onset_l_reward = time_vid(ind_lick_onset_l_reward);

        ind_lick_offset_l_reward = ind_lick_offset(sort([ind_l_reward_inner_tube_lick;ind_l_reward_outer_tube_lick;ind_l_noreward_inner_tube_lick;ind_l_noreward_outer_tube_lick]));
        time_lick_offset_l_reward = time_vid(ind_lick_offset_l_reward);

        DLC.CLASS.is_grooming_lick = is_grooming_lick;
        DLC.CLASS.is_r_reward_inner_tube_lick = is_r_reward_inner_tube_lick;
        DLC.CLASS.is_r_reward_outer_tube_lick = is_r_reward_outer_tube_lick;
        DLC.CLASS.is_r_noreward_inner_tube_lick = is_r_noreward_inner_tube_lick;
        DLC.CLASS.is_r_noreward_outer_tube_lick = is_r_noreward_outer_tube_lick;
        DLC.CLASS.ind_lick_onset_r_reward = ind_lick_onset_r_reward;
        DLC.TIME.time_lick_onset_r_reward = time_lick_onset_r_reward;
        DLC.CLASS.ind_lick_offset_r_reward = ind_lick_offset_r_reward;
        DLC.TIME.time_lick_offset_r_reward = time_lick_offset_r_reward;

        DLC.CLASS.is_l_reward_inner_tube_lick = is_l_reward_inner_tube_lick;
        DLC.CLASS.is_l_reward_outer_tube_lick = is_l_reward_outer_tube_lick;
        DLC.CLASS.is_l_noreward_inner_tube_lick = is_l_noreward_inner_tube_lick;
        DLC.CLASS.is_l_noreward_outer_tube_lick = is_l_noreward_outer_tube_lick;
        DLC.CLASS.ind_lick_onset_l_reward = ind_lick_onset_l_reward;
        DLC.TIME.time_lick_onset_l_reward = time_lick_onset_l_reward;
        DLC.CLASS.ind_lick_offset_l_reward = ind_lick_offset_l_reward;
        DLC.TIME.time_lick_offset_l_reward = time_lick_offset_l_reward;

        fig = figure;
        hold on;
        % plot(tip_tongue_x(ind_d_lick_max),tip_tongue_y(ind_d_lick_max), '.k');
        plot(tip_tongue_x(ind_d_lick_max(ind_grooming_lick)),tip_tongue_y(ind_d_lick_max(ind_grooming_lick)), '.b');
        plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_inner_tube_lick)), '.g');
        plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_inner_tube_lick)), '.r');
        plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_outer_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_outer_tube_lick)), '.c');
        plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_outer_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_outer_tube_lick)), '.m');
        if with_nose_markers
            plot(r_nose_x, r_nose_y,'ok');
            plot(l_nose_x, l_nose_y,'ok');
        end
        plot(r_tube_r_x(ind_d_lick_max),r_tube_r_y(ind_d_lick_max),'sk');
        plot(r_tube_l_x(ind_d_lick_max),r_tube_l_y(ind_d_lick_max),'sk');

        plot(tip_tongue_x(ind_d_lick_max(ind_l_noreward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_noreward_inner_tube_lick)), '.g');
        plot(tip_tongue_x(ind_d_lick_max(ind_l_reward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_reward_inner_tube_lick)), '.r');
        plot(tip_tongue_x(ind_d_lick_max(ind_l_noreward_outer_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_noreward_outer_tube_lick)), '.c');

        plot(tip_tongue_x(ind_d_lick_max(ind_l_reward_outer_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_reward_outer_tube_lick)), '.m');
        plot(l_tube_r_x(ind_d_lick_max), l_tube_r_y(ind_d_lick_max),'sk');
        plot(l_tube_l_x(ind_d_lick_max), l_tube_l_y(ind_d_lick_max),'sk');
        xlabel('x position (mm)');
        ylabel('y position (mm)');
        set(gca, 'YDir','reverse')
        xlim([-10 15]);
        ylim([-25 25]);
        title('Max lick positions');
        
        parts = strsplit(DLC.FILE.file_path, ["/", "\"]);
        name = erase(parts(length(parts) - 2), '-');
        file_name = char(extractBetween(name, 3, 15));
        file_path = DLC.FILE.file_path;
        saveas(fig,[file_path '../analyzed_figs/' file_name '_DLC_Report.png'])
        %     ESN_Beautify_Plot

        
    end

    if (session_type == 1)
        bool_tongue_r_tube_empty = DLC.GEO.bool_tongue_r_tube_empty;
        bool_tongue_r_tube_full = DLC.GEO.bool_tongue_r_tube_full;

        tip_tongue_x = DLC.POINTS.tip_tongue_x;
        tip_tongue_y = DLC.POINTS.tip_tongue_y;
        if with_nose_markers
            r_nose_x =  DLC.POINTS.r_nose_x;
            r_nose_y = DLC.POINTS.r_nose_y;
            l_nose_x =  DLC.POINTS.l_nose_x;
            l_nose_y = DLC.POINTS.l_nose_y;
        end
        r_tube_r_x = DLC.POINTS.r_tube_r_x;
        r_tube_r_y = DLC.POINTS.r_tube_r_y;
        r_tube_l_x = DLC.POINTS.r_tube_l_x;
        r_tube_l_y = DLC.POINTS.r_tube_l_y;
        r_food_y = DLC.POINTS.r_food_y;

        ind_d_lick_max = DLC.IND.ind_d_lick_max;
        ind_lick_onset = DLC.IND.ind_lick_onset;
        ind_lick_offset = DLC.IND.ind_lick_offset;
        time_vid = DLC.TIME.time_vid';


        ind_grooming_lick = find(bool_tongue_r_tube_empty == 0 & bool_tongue_r_tube_full == 0 & ...
            tip_tongue_y(ind_d_lick_max) < (r_tube_r_y(ind_d_lick_max)));
        is_grooming_lick = false(size(ind_lick_onset));
        is_grooming_lick(ind_grooming_lick) = true;
        
        ind_r_reward_inner_tube_lick = find(bool_tongue_r_tube_full == 1 & ...
            tip_tongue_x(ind_d_lick_max) < r_tube_r_x(ind_d_lick_max) & tip_tongue_x(ind_d_lick_max) > r_tube_l_x(ind_d_lick_max));
        is_r_reward_inner_tube_lick = false(size(ind_lick_onset));
        is_r_reward_inner_tube_lick(ind_r_reward_inner_tube_lick) = true;
        
        ind_r_reward_outer_tube_lick = find(bool_tongue_r_tube_full == 1 &...
            (tip_tongue_x(ind_d_lick_max) > r_tube_r_x(ind_d_lick_max) | tip_tongue_x(ind_d_lick_max) < r_tube_l_x(ind_d_lick_max)));
        is_r_reward_outer_tube_lick = false(size(ind_lick_onset));
        is_r_reward_outer_tube_lick(ind_r_reward_outer_tube_lick) = true;

        ind_r_noreward_inner_tube_lick = find(bool_tongue_r_tube_empty == 1 & bool_tongue_r_tube_full == 0 &...
            tip_tongue_x(ind_d_lick_max) < r_tube_r_x(ind_d_lick_max) & tip_tongue_x(ind_d_lick_max) > r_tube_l_x(ind_d_lick_max));
        is_r_noreward_inner_tube_lick = false(size(ind_lick_onset));
        is_r_noreward_inner_tube_lick(ind_r_noreward_inner_tube_lick) = true;

        ind_r_noreward_outer_tube_lick = find(bool_tongue_r_tube_full == 0 & tip_tongue_y(ind_d_lick_max) > r_tube_r_y(ind_d_lick_max) & ...
            (tip_tongue_x(ind_d_lick_max) > r_tube_r_x(ind_d_lick_max) | tip_tongue_x(ind_d_lick_max) < r_tube_l_x(ind_d_lick_max)));
        is_r_noreward_outer_tube_lick = false(size(ind_lick_onset));
        ind_r_noreward_outer_tube_lick(is_r_noreward_outer_tube_lick) = true;

        ind_lick_onset_r_reward = ind_lick_onset(sort([ind_r_reward_inner_tube_lick;ind_r_reward_outer_tube_lick;ind_r_noreward_inner_tube_lick;ind_r_noreward_outer_tube_lick]));
        time_lick_onset_r_reward = time_vid(ind_lick_onset_r_reward);
        
        ind_lick_offset_r_reward = ind_lick_offset(sort([ind_r_reward_inner_tube_lick;ind_r_reward_outer_tube_lick;ind_r_noreward_inner_tube_lick;ind_r_noreward_outer_tube_lick]));
        time_lick_offset_r_reward = time_vid(ind_lick_offset_r_reward);
        
        
        DLC.CLASS.is_grooming_lick = is_grooming_lick;
        DLC.CLASS.is_r_reward_inner_tube_lick = is_r_reward_inner_tube_lick;
        DLC.CLASS.is_r_reward_outer_tube_lick = is_r_reward_outer_tube_lick;
        DLC.CLASS.is_r_noreward_inner_tube_lick = is_r_noreward_inner_tube_lick;
        DLC.CLASS.is_r_noreward_outer_tube_lick = is_r_noreward_outer_tube_lick;
        DLC.CLASS.ind_lick_onset_r_reward = ind_lick_onset_r_reward;
        DLC.TIME.time_lick_onset_r_reward = time_lick_onset_r_reward;
        DLC.CLASS.ind_lick_offset_r_reward = ind_lick_offset_r_reward;
        DLC.TIME.time_lick_offset_r_reward = time_lick_offset_r_reward;

        fig = figure;
        hold on;
        % plot(tip_tongue_x(ind_d_lick_max),tip_tongue_y(ind_d_lick_max), '.k');
        plot(tip_tongue_x(ind_d_lick_max(ind_grooming_lick)),tip_tongue_y(ind_d_lick_max(ind_grooming_lick)), '.b');
        plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_inner_tube_lick)), '.g');
        plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_inner_tube_lick)), '.r');
        plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_outer_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_outer_tube_lick)), '.c');
        plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_outer_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_outer_tube_lick)), '.m');
        if with_nose_markers
            plot(r_nose_x,r_nose_y,'ok');
            plot(l_nose_x,l_nose_y,'ok');
        end
        plot(r_tube_r_x(ind_d_lick_max),r_tube_r_y(ind_d_lick_max),'sk');
        plot(r_tube_l_x(ind_d_lick_max),r_tube_l_y(ind_d_lick_max),'sk');

        xlabel('x position (mm)');
        ylabel('y position (mm)');
        set(gca, 'YDir','reverse')
        xlim([-10 15]);
        ylim([-25 25]);
        title('Max lick positions');
        
        [~,file_name,~] = fileparts(DLC.FILE.file_name);
        file_path = DLC.FILE.file_path;
        saveas(fig,[file_path '../analyzed_figs/' file_name '_Report.png'])
        %     ESN_Beautify_Plot
    end

    fprintf(' --> Completed. \n')
end

%% Quantify food in tube
function DLC = quantify_food(DLC)
    fprintf('Quantifying food in tube ...');

    session_type = DLC.FILE.session_type;
    with_nose_markers = DLC.FILE.with_nose_markers;
    time_vid = DLC.TIME.time_vid;
    r_tube_r_y = DLC.POINTS.r_tube_r_y;
    r_food_x = DLC.POINTS.r_food_x;
    r_food_y = DLC.POINTS.r_food_y;
    r_food_x(r_food_y < 0) = nan;
    r_food_y(r_food_y < 0) = nan;
    r_food_max = max(r_food_y);
    r_food_min = min(r_food_y);
    r_tube_food = (r_food_y - r_food_max)./(r_tube_r_y - r_food_max);
    DLC.FOOD.r_tube_food = r_tube_food;
    DLC.FOOD.r_food_x = r_food_x;
    DLC.FOOD.r_food_y = r_food_y;

    if session_type == 2
        l_tube_r_y = DLC.POINTS.l_tube_r_y;
        l_food_x = DLC.POINTS.l_food_x;
        l_food_y = DLC.POINTS.l_food_y;
        l_food_x(l_food_y > 0) = nan;
        l_food_y(l_food_y > 0) = nan;
        l_food_max = max(l_food_y);
        l_food_min = min(l_food_y);
        l_tube_food = (l_food_y - l_food_min)./(l_tube_r_y - l_food_min);
        DLC.FOOD.l_tube_food = l_tube_food;
        DLC.FOOD.l_food_x = l_food_x;
        DLC.FOOD.l_food_y = l_food_y;
    end

    % figure;
    % hold on;
    % plot(time_vid, r_tube_food);
    % if session_type == 2
    %     plot(time_vid, l_tube_food);
    % end
    % xlabel('Time (s)');
    % ylabel('Food in tube (empty:0 | full:1)');
    % % ESN_Beautify_Plot

    fprintf(' --> Completed. \n')
end

%% Detect harvest str, end, num, times, and duration
function DLC = detect_harvest(DLC)
    fprintf('Detecting harvest ...');

    session_type = DLC.FILE.session_type;
    with_nose_markers = DLC.FILE.with_nose_markers;
    d_tip = DLC.KINEMATIC.d_tip;
    time_vid = DLC.TIME.time_vid';
    ind_lick_onset = DLC.IND.ind_lick_onset;
    is_grooming_lick = DLC.CLASS.is_grooming_lick;
    ind_lick_onset_r_reward = DLC.CLASS.ind_lick_onset_r_reward;
    time_lick_onset_r_reward = DLC.TIME.time_lick_onset_r_reward;
    r_tube_food = DLC.FOOD.r_tube_food;
    if session_type == 2
        ind_lick_onset_l_reward = DLC.CLASS.ind_lick_onset_l_reward;
        time_lick_onset_l_reward = DLC.TIME.time_lick_onset_l_reward;
        l_tube_food = DLC.FOOD.l_tube_food;
    end
    ind_lick_onset_str_bout = DLC.IND.ind_lick_onset_str_bout;
    ind_lick_onset_end_bout = DLC.IND.ind_lick_onset_end_bout;

    ind_lick_onset = ind_lick_onset(is_grooming_lick);
   
    ind_lick_onset_str_harvest_ = [1; find(diff(ind_lick_onset) >100) + 1];
    ind_lick_onset_str_harvest = ind_lick_onset(ind_lick_onset_str_harvest_);
    ind_lick_onset_end_harvest_ = [find(diff(ind_lick_onset) >100); length(ind_lick_onset)];
    ind_lick_onset_end_harvest = ind_lick_onset(ind_lick_onset_end_harvest_);

    num_lick_harvest = (ind_lick_onset_end_harvest_ - ind_lick_onset_str_harvest_) + 1;

    time_lick_onset_str_harvest = time_vid(ind_lick_onset_str_harvest);
    time_lick_onset_end_harvest = time_vid(ind_lick_onset_end_harvest);
    time_harvest_duration = time_vid(ind_lick_onset_end_harvest) - ...
        time_vid(ind_lick_onset_str_harvest);

    inds_del_harvest = num_lick_harvest<3;
    num_lick_harvest(inds_del_harvest) = [];
    time_harvest_duration(inds_del_harvest) = [];
    ind_lick_onset_str_harvest(inds_del_harvest) = [];
    ind_lick_onset_end_harvest(inds_del_harvest) = [];
    time_lick_onset_str_harvest(inds_del_harvest) = [];
    time_lick_onset_end_harvest(inds_del_harvest) = [];
    
    % Determine direction of bouts and harvest
    num_r = [];
    if session_type == 1
        num_l = zeros(size(ind_lick_onset_str_bout));
    else
        num_l = [];
    end
    for i = 1:length(ind_lick_onset_str_bout)
        inds_lick_onset = ind_lick_onset_str_bout(i):ind_lick_onset_end_bout(i);
        num_r = [num_r; length(find(ismember(inds_lick_onset, ind_lick_onset_r_reward)))];
        if session_type == 2
            num_l = [num_l; length(find(ismember(inds_lick_onset, ind_lick_onset_l_reward)))];
        end   
    end
    is_bout_r = num_r > num_l;
    if session_type == 2
        is_bout_l = num_l > num_r;
    end
    
    num_r = [];
    if session_type == 1
        num_l = zeros(size(ind_lick_onset_str_harvest));
    else
        num_l = [];
    end
    for i = 1:length(ind_lick_onset_str_harvest)
        inds_lick_onset = ind_lick_onset_str_harvest(i):ind_lick_onset_end_harvest(i);
        num_r = [num_r; length(find(ismember(inds_lick_onset, ind_lick_onset_r_reward)))];
        if session_type == 2
            num_l = [num_l; length(find(ismember(inds_lick_onset, ind_lick_onset_l_reward)))];
        end
    end
    is_harvest_r = num_r > num_l;
    if session_type == 2
        is_harvest_l = num_l > num_r;
    end
      
    % time_r_harvest_duration = time_harvest_duration(ismember(ind_lick_onset_str_harvest,ind_lick_onset_r_reward));
    % time_l_harvest_duration = time_harvest_duration(ismember(ind_lick_onset_str_harvest,ind_lick_onset_l_reward));
    %
    % num_lick_r_harvest = num_lick_harvest(ismember(ind_lick_onset_str_harvest,ind_lick_onset_r_reward));
    % num_lick_l_harvest = num_lick_harvest(ismember(ind_lick_onset_str_harvest,ind_lick_onset_l_reward));

    % r_tube_food_str_harvest = r_tube_food(ind_lick_onset_str_harvest(ismember(ind_lick_onset_str_harvest,ind_lick_onset_r_reward)));
    % r_tube_food_end_harvest = r_tube_food(ind_lick_onset_end_harvest(ismember(ind_lick_onset_end_harvest,ind_lick_onset_r_reward)));
    % l_tube_food_str_harvest = l_tube_food(ind_lick_onset_str_harvest(ismember(ind_lick_onset_str_harvest,ind_lick_onset_l_reward)));
    % l_tube_food_end_harvest = l_tube_food(ind_lick_onset_end_harvest(ismember(ind_lick_onset_end_harvest,ind_lick_onset_l_reward)));

    % r_tube_food_consumed_harvest = r_tube_food_str_harvest - [r_tube_food_end_harvest;0;0];
    % l_tube_food_consumed_harvest = l_tube_food_str_harvest - l_tube_food_end_harvest;
    
    DLC.IND.ind_lick_onset_str_harvest = ind_lick_onset_str_harvest;
    DLC.IND.ind_lick_onset_end_harvest = ind_lick_onset_end_harvest;
    DLC.IND.num_lick_harvest = num_lick_harvest;
    % DLC.IND.num_lick_r_harvest = num_lick_r_harvest;
    % DLC.IND.num_lick_l_harvest = num_lick_l_harvest;
    DLC.TIME.time_lick_onset_str_harvest = time_lick_onset_str_harvest;
    DLC.TIME.time_lick_onset_end_harvest = time_lick_onset_end_harvest;
    DLC.TIME.time_harvest_duration = time_harvest_duration;
    % DLC.TIME.time_r_harvest_duration = time_r_harvest_duration;
    % DLC.TIME.time_l_harvest_duration = time_l_harvest_duration;
    % DLC.FOOD.r_tube_food_str_harvest = r_tube_food_str_harvest;
    % DLC.FOOD.r_tube_food_end_harvest = r_tube_food_end_harvest;
    % DLC.FOOD.l_tube_food_str_harvest = l_tube_food_str_harvest;
    % DLC.FOOD.l_tube_food_end_harvest = l_tube_food_end_harvest;
    DLC.CLASS.is_bout_r = is_bout_r;
    DLC.CLASS.is_harvest_r = is_harvest_r;
    if session_type == 2
        DLC.CLASS.is_bout_l = is_bout_l;
        DLC.CLASS.is_harvest_l = is_harvest_l;
    end
    
    time_lick_onset_str_bout =  DLC.TIME.time_lick_onset_str_bout;
    time_lick_onset_end_bout= DLC.TIME.time_lick_onset_end_bout;

    figure;
    hold on;
    plot(time_vid,d_tip,'.-k');
    xlabel('Time (s)');
    ylabel('Displacement (mm)');

    yyaxis right
    plot(time_lick_onset_str_bout,2,'*g');
    plot(time_lick_onset_end_bout,2,'*r');
    if ~isempty(time_lick_onset_r_reward)
        plot(time_lick_onset_r_reward,1.9,'.r');
    end
    plot(time_vid,r_tube_food,'.-r');

    if session_type == 2
        if ~isempty(time_lick_onset_l_reward)
        	plot(time_lick_onset_l_reward,1.9,'.b');
        end
        plot(time_vid,l_tube_food,'.-b');
    end
    ylabel('Reward capacity (0:Empty | 1:Full)')
    % ESN_Beautify_Plot

    fprintf(' --> Completed. \n')
end

%% This function calls all other functions that analyze the tracking data to product DLC.mat
function DLC = evaluate(DLC)
    [DLC, threshold_forward, threshold_backward] = suggest_threshold(DLC);
    DLC = detect_licks_and_bouts(DLC, threshold_forward, threshold_backward);
    DLC = calculate_lick_kinematics(DLC);
    DLC = geometrization(DLC);
    DLC = sort_licks(DLC);
    DLC = quantify_food(DLC);
    DLC = detect_harvest(DLC);
end

%% Save preprocessed DLC data
function save_DLC(DLC)
    fprintf('Saving: ');
    parts = strsplit(DLC.FILE.file_path, ["/", "\"]);
    name = erase(parts(length(parts) - 2), '-');
    file_name = char(extractBetween(name, 3, 15));
    file_path = DLC.FILE.file_path;
    fprintf([file_name '.mat ...'])
    save([ file_path  '../analyzed_data/' file_name  '_DLC.mat'], 'DLC', '-v7.3');
    fprintf(' --> Completed. \n')   
end

%% Quality Assurance
function DLC = quality_assurance(DLC)
    fprintf("Performing Quality Assurance...");

    % read from video
    [vidfile_name, vidfile_path] = uigetfile([pwd filesep '*.mp4;*.avi'], 'Select a video file');
    vid = VideoReader(fullfile(vidfile_path, vidfile_name));
    scal = 7/47;
    origin_x = DLC.POINTS.origin_x;
    origin_y = DLC.POINTS.origin_y;

    % initialization
    ind_inaccurate_tip_tongue = [];
    ind_inaccurate_l_tongue = [];
    ind_inaccurate_r_tongue = [];
    ind_inaccurate_mid_tongue = [];
    ind_inaccurate_l_food = [];
    jumping_l_food = [];
    ind_inaccurate_r_food = [];
    jumping_r_food = [];

    session_type = DLC.FILE.session_type;
    with_nose_markers = DLC.FILE.with_nose_markers;
    ind_lick_onset = DLC.IND.ind_lick_onset;
    ind_lick_offset = DLC.IND.ind_lick_offset;

    x_upper_lim = vid.Width - 10;
    y_lower_lim = 10;
    y_upper_lim = vid.Height - 10;
    tip_tongue_x = (DLC.POINTS.tip_tongue_x + origin_x) / scal;
    tip_tongue_y = (DLC.POINTS.tip_tongue_y + origin_y) / scal;
    r_tongue_x = (DLC.POINTS.r_tongue_x + origin_x) / scal;
    r_tongue_y = (DLC.POINTS.r_tongue_y + origin_y) / scal;
    l_tongue_x = (DLC.POINTS.l_tongue_x + origin_x) / scal;
    l_tongue_y = (DLC.POINTS.l_tongue_y + origin_y) / scal;
    mid_tongue_x = (DLC.POINTS.mid_tongue_x + origin_x) / scal;
    mid_tongue_y = (DLC.POINTS.mid_tongue_y + origin_y) / scal;           
    r_food_x = (DLC.POINTS.r_food_x + origin_x) / scal;
    r_food_y = (DLC.POINTS.r_food_y + origin_y) / scal;
    r_tube_l_x = (DLC.POINTS.r_tube_l_x + origin_x) / scal;
    r_tube_l_y = (DLC.POINTS.r_tube_l_y + origin_y) / scal;
    r_tube_r_x = (DLC.POINTS.r_tube_r_x + origin_x) / scal;
    r_tube_r_y = (DLC.POINTS.r_tube_r_y + origin_y) / scal;
    if with_nose_markers
        r_nose_x = (DLC.POINTS.r_nose_x + origin_x) / scal;
        r_nose_y = (DLC.POINTS.r_nose_y + origin_y) / scal;
        l_nose_x = (DLC.POINTS.l_nose_x + origin_x) / scal;
        l_nose_y = (DLC.POINTS.l_nose_y + origin_y) / scal;
        mean_nose_x = (mean(r_nose_x) + mean(l_nose_x)) / 2;
    end
    if session_type == 2
        l_food_x = (DLC.POINTS.l_food_x + origin_x) / scal;
        l_food_y = (DLC.POINTS.l_food_y + origin_y) / scal;
        l_tube_l_x = (DLC.POINTS.l_tube_l_x + origin_x) / scal;
        l_tube_l_y = (DLC.POINTS.l_tube_l_y + origin_y) / scal;
        l_tube_r_x = (DLC.POINTS.l_tube_r_x + origin_x) / scal;
        l_tube_r_y = (DLC.POINTS.l_tube_r_y + origin_y) / scal;
    end
    y0 = origin_y / scal;

    mean_r_tube_l_x = mean(r_tube_l_x);
    mean_r_tube_r_x = mean(r_tube_r_x);
    mean_r_tube_l_y = mean(r_tube_l_y);
    if session_type == 2
        mean_l_tube_l_x = mean(l_tube_l_x);
        mean_l_tube_r_x = mean(l_tube_r_x);
        mean_l_tube_l_y = mean(l_tube_l_y);
    end

    % Correct stationary markers
    if with_nose_markers
        [DLC.POINTS.l_nose_x, DLC.POINTS.l_nose_y] = correct_stationary_markers(l_nose_x, l_nose_y, DLC.POINTS.l_nose_x, DLC.POINTS.l_nose_y);
        [DLC.POINTS.r_nose_x, DLC.POINTS.r_nose_y] = correct_stationary_markers(r_nose_x, r_nose_y, DLC.POINTS.r_nose_x, DLC.POINTS.r_nose_y);        
    end
    [DLC.POINTS.r_tube_l_x, DLC.POINTS.r_tube_l_y] = correct_stationary_markers(r_tube_l_x, r_tube_l_y, DLC.POINTS.r_tube_l_x, DLC.POINTS.r_tube_l_y);
    [DLC.POINTS.r_tube_r_x, DLC.POINTS.r_tube_r_y] = correct_stationary_markers(r_tube_r_x, r_tube_r_y, DLC.POINTS.r_tube_r_x, DLC.POINTS.r_tube_r_y);
    if session_type == 2
        [DLC.POINTS.l_tube_l_x, DLC.POINTS.l_tube_l_y] = correct_stationary_markers(l_tube_l_x, l_tube_l_y, DLC.POINTS.l_tube_l_x, DLC.POINTS.l_tube_l_y);
        [DLC.POINTS.l_tube_r_x, DLC.POINTS.l_tube_r_y] = correct_stationary_markers(l_tube_r_x, l_tube_r_y, DLC.POINTS.l_tube_r_x, DLC.POINTS.l_tube_r_y);
    end

    % Find inaccurate tongue markers
    midtip_y = tip_tongue_y - mid_tongue_y;
    midtip_x = tip_tongue_x - mid_tongue_x;
    angle = zeros(size(midtip_y));

    for i = 1:length(midtip_y)
        angle(i) = rad2deg(atan2(midtip_y(i), midtip_x(i)));
    end

     for i = 2:length(tip_tongue_y)
       if with_nose_markers
            if tip_tongue_y(i) <= y_lower_lim || tip_tongue_x(i) < mean_nose_x || tip_tongue_y(i) >= y_upper_lim || abs(angle(i) - angle(i-1)) > 80
                ind_inaccurate_tip_tongue = [ind_inaccurate_tip_tongue; i];
            end

            if l_tongue_y(i) <= y_lower_lim || l_tongue_x(i) < mean_nose_x || l_tongue_y(i) >= y_upper_lim
                ind_inaccurate_l_tongue = [ind_inaccurate_l_tongue; i];
            end

            if r_tongue_y(i) <= y_lower_lim || r_tongue_x(i) < mean_nose_x || r_tongue_y(i) >= y_upper_lim
                ind_inaccurate_r_tongue = [ind_inaccurate_r_tongue; i];
            end

            if mid_tongue_y(i) <= y_lower_lim || mid_tongue_x(i) < mean_nose_x || mid_tongue_x(i) >= x_upper_lim || mid_tongue_y(i) >= y_upper_lim
                ind_inaccurate_mid_tongue = [ind_inaccurate_mid_tongue; i];
            end
       else
           if tip_tongue_y(i) <= y_lower_lim || tip_tongue_y(i) >= y_upper_lim || abs(angle(i) - angle(i-1)) > 80
                ind_inaccurate_tip_tongue = [ind_inaccurate_tip_tongue; i];
            end

            if l_tongue_y(i) <= y_lower_lim || l_tongue_y(i) >= y_upper_lim
                ind_inaccurate_l_tongue = [ind_inaccurate_l_tongue; i];
            end

            if r_tongue_y(i) <= y_lower_lim || r_tongue_y(i) >= y_upper_lim
                ind_inaccurate_r_tongue = [ind_inaccurate_r_tongue; i];
            end

            if mid_tongue_y(i) <= y_lower_lim || mid_tongue_x(i) >= x_upper_lim || mid_tongue_y(i) >= y_upper_lim
                ind_inaccurate_mid_tongue = [ind_inaccurate_mid_tongue; i];
            end
       end
        % Find inaccurate right food markers
        if r_food_y(i) < y0 || r_food_x(i) < mean_r_tube_l_x ||  r_food_x(i) > mean_r_tube_r_x || r_food_y(i) >= y_upper_lim
            ind_inaccurate_r_food = [ind_inaccurate_r_food; i];
        end
        if abs(r_food_y(i) - r_food_y(i - 1)) > (vid.Height - mean_r_tube_l_y)/2
            jumping_r_food = [jumping_r_food; i];
        end

        if session_type == 2
            % Find inaccurate left food markers
            if l_food_y(i) > y0 || l_food_x(i) < mean_l_tube_l_x ||  l_food_x(i) > mean_l_tube_r_x || l_food_y(i) <= y_lower_lim
                ind_inaccurate_l_food = [ind_inaccurate_l_food; i];
            end
            
            if abs(l_food_y(i) - l_food_y(i - 1)) > mean_l_tube_l_y/2
                jumping_l_food = [jumping_l_food; i];
            end

            if r_food_y(i) - l_food_y(i) < mean_r_tube_l_y - mean_l_tube_l_y && ~((r_food_y(i) < mean_r_tube_l_y && r_food_y(i) > y0) || (l_food_y(i) > mean_l_tube_l_y && l_food_y(i) < y0))
                ind_inaccurate_r_food = [ind_inaccurate_r_food; i];
                ind_inaccurate_l_food = [ind_inaccurate_l_food; i];
            end
        end
    end

%     ind_inaccurate_tip_tongue = sortrows(ind_inaccurate_tip_tongue);
%     ind_inaccurate_l_tongue = sortrows(ind_inaccurate_l_tongue);
%     ind_inaccurate_r_tongue = sortrows(ind_inaccurate_r_tongue);
%     ind_inaccurate_mid_tongue = sortrows(ind_inaccurate_mid_tongue);
    
    if mod(length(jumping_r_food), 2) == 1
        fprintf('Please check the r_food marker at frame %d\n', jumping_r_food(length(jumping_r_food)));
        jumping_r_food(length(jumping_r_food)) = [];
    end   
    for i = 1:2:length(jumping_r_food)
        to_add = jumping_r_food(i):jumping_r_food(i+1);
        ind_inaccurate_r_food = [ind_inaccurate_r_food; to_add(:)];
    end
    ind_inaccurate_r_food = unique(sortrows(ind_inaccurate_r_food));
    
    if session_type == 2
        if mod(length(jumping_l_food), 2) == 1
            fprintf('Please check the l_food marker at frame %d\n', jumping_l_food(length(jumping_l_food)));
            jumping_l_food(length(jumping_l_food)) = [];
        end
        for i = 1:2:length(jumping_l_food)
            to_add = jumping_l_food(i):jumping_l_food(i+1);
            ind_inaccurate_l_food = [ind_inaccurate_l_food; to_add(:)];
        end
        ind_inaccurate_l_food = unique(sortrows(ind_inaccurate_l_food));
    end

    % Correct tongue markers manually
    vid = VideoReader(fullfile(vidfile_path, vidfile_name));
    for i = 1:length(ind_inaccurate_tip_tongue)
        counter_frame = ind_inaccurate_tip_tongue(i);
        frame = read(vid, counter_frame);  
        fig = figure;
        set(fig, 'ToolBar', 'none');
        imshow(frame)
        title('Indicate the position of tongue tip marker')
        [x, y] = getpts;
        close(fig)

        DLC.POINTS.tip_tongue_x(counter_frame) = scal * x - origin_x;
        DLC.POINTS.tip_tongue_y(counter_frame) = scal * y - origin_y;
    end

    for i = 1:length(ind_inaccurate_l_tongue)
        counter_frame = ind_inaccurate_l_tongue(i);
        frame = read(vid, counter_frame);  
        fig = figure;
        set(fig, 'ToolBar', 'none');
        imshow(frame)
        title('Indicate the position of tongue left marker')
        [x, y] = getpts;
        close(fig)

        DLC.POINTS.l_tongue_x(counter_frame) = scal * x - origin_x;
        DLC.POINTS.l_tongue_y(counter_frame) = scal * y - origin_y;
    end

    for i = 1:length(ind_inaccurate_r_tongue)
        counter_frame = ind_inaccurate_r_tongue(i);
        frame = read(vid, counter_frame);  
        fig = figure;
        set(fig, 'ToolBar', 'none');
        imshow(frame)
        title('Indicate the position of tongue right marker')
        [x, y] = getpts;
        close(fig)

        DLC.POINTS.r_tongue_x(counter_frame) = scal * x - origin_x;
        DLC.POINTS.r_tongue_y(counter_frame) = scal * y - origin_y;
    end

    for i = 1:length(ind_inaccurate_mid_tongue)
        counter_frame = ind_inaccurate_mid_tongue(i);
        frame = read(vid, counter_frame);  
        fig = figure;
        set(fig, 'ToolBar', 'none');
        imshow(frame)
        title('Indicate the position of tongue mid marker')
        [x, y] = getpts;
        close(fig)

        DLC.POINTS.mid_tongue_x(counter_frame) = scal * x - origin_x;
        DLC.POINTS.mid_tongue_y(counter_frame) = scal * y - origin_y;
    end

    if ~isempty(ind_inaccurate_tip_tongue) || ~isempty(ind_inaccurate_l_tongue) || ~isempty(ind_inaccurate_r_tongue) || ~isempty(ind_inaccurate_mid_tongue)
        DLC = evaluate(DLC); 
    end

    % Correct food markers
    rightward = [DLC.CLASS.ind_r_reward_inner_tube_lick; DLC.CLASS.ind_r_noreward_inner_tube_lick; DLC.CLASS.ind_r_reward_outer_tube_lick; DLC.CLASS.ind_r_noreward_outer_tube_lick];
    if session_type == 2
        leftward = [DLC.CLASS.ind_l_reward_inner_tube_lick; DLC.CLASS.ind_l_noreward_inner_tube_lick; DLC.CLASS.ind_l_reward_outer_tube_lick; DLC.CLASS.ind_l_noreward_outer_tube_lick]; 
    end
    
    if ~isempty(ind_inaccurate_r_food)
        i = 1;
        while i <= size(ind_inaccurate_r_food, 1)
            ind = ind_inaccurate_r_food(i);
            n = ind + 1;
            while ismember(n, ind_inaccurate_r_food)
                n = n + 1;
                i = i + 1;
            end

            previous_onset_ind = find(ind_lick_onset > ind - 1, 1);
            previous_offset_ind = find(ind_lick_offset > ind - 1, 1);
            previous_in_lick = 1;
            if previous_onset_ind ~= 1 && (ind - 1) < ind_lick_onset(previous_onset_ind - 1) && (ind - 1) > ind_lick_offset(previous_offset_ind - 1)
                previous_in_lick = 0;
            end
            next_onset_ind = find(ind_lick_onset > n, 1);
            next_offset_ind = find(ind_lick_offset > n, 1);
            next_in_lick = 1;
            if n < ind_lick_onset(next_onset_ind - 1) & n > ind_lick_offset(next_offset_ind - 1)
                next_in_lick = 0;
            end

            if previous_onset_ind == 1 || (previous_in_lick && next_in_lick && ismember(ind_lick_onset(previous_onset_ind - 1), rightward) && ismember(ind_lick_onset(next_onset_ind - 1), rightward))
                frame = read(vid, ind);  
                fig = figure;
                imshow(frame)
                set(fig, 'ToolBar', 'none');
                title('Indicate the position of right food marker')
                [x, y] = getpts;
                close(fig)

                DLC.POINTS.r_food_x(ind) = scal * x - origin_x;
                DLC.POINTS.r_food_y(ind) = scal * y - origin_y;
            else
                for j = ind : n - 1
                    DLC.POINTS.r_food_x(j) = DLC.POINTS.r_food_x(j-1) + ((DLC.POINTS.r_food_x(n) - DLC.POINTS.r_food_x(ind - 1)) / (n - ind + 1));
                    DLC.POINTS.r_food_y(j) = DLC.POINTS.r_food_y(j-1) + ((DLC.POINTS.r_food_y(n) - DLC.POINTS.r_food_y(ind - 1)) / (n - ind + 1));
                end
            end
            i = i + 1;
        end
    end

    if session_type == 2
        if ~isempty(ind_inaccurate_l_food)
            i = 1;
            while i <= size(ind_inaccurate_l_food, 1)
                ind = ind_inaccurate_l_food(i);
                n = ind + 1;
                while ismember(n, ind_inaccurate_l_food)
                    n = n + 1;
                    i = i + 1;
                end

                if n == ind + 1
                    i = i + 1;
                end

                previous_onset_ind = find(ind_lick_onset > ind - 1, 1);
                previous_offset_ind = find(ind_lick_offset > ind - 1, 1);
                previous_in_lick = 1;
                if previous_onset_ind ~= 1 && (ind - 1) < ind_lick_onset(previous_onset_ind - 1) && (ind - 1) > ind_lick_offset(previous_offset_ind - 1)
                    previous_in_lick = 0;
                end
                next_onset_ind = find(ind_lick_onset > n, 1);
                next_offset_ind = find(ind_lick_offset > n, 1);
                next_in_lick = 1;
                if n < ind_lick_onset(next_onset_ind - 1) & n > ind_lick_offset(next_offset_ind - 1)
                    next_in_lick = 0;
                end

                if previous_onset_ind == 1 || (previous_in_lick && next_in_lick && ismember(ind_lick_onset(previous_onset_ind - 1), leftward) && ismember(ind_lick_onset(next_onset_ind - 1), leftward))
                    frame = read(vid, ind);  
                    fig = figure;
                    imshow(frame)
                    set(fig, 'ToolBar', 'none');
                    title('Indicate the position of right food marker')
                    [x, y] = getpts;
                    close(fig)

                    DLC.POINTS.l_food_x(ind) = scal * x - origin_x;
                    DLC.POINTS.l_food_y(ind) = scal * y - origin_y;

                else
                    for j = ind : n - 1
                        DLC.POINTS.l_food_x(j) = DLC.POINTS.l_food_x(j-1) + ((DLC.POINTS.l_food_x(n) - DLC.POINTS.l_food_x(ind - 1)) / (n - ind + 1));
                        DLC.POINTS.l_food_y(j) = DLC.POINTS.l_food_y(j-1) + ((DLC.POINTS.l_food_y(n) - DLC.POINTS.l_food_y(ind - 1)) / (n - ind + 1));
                    end
                end
                i = i + 1;
            end
        end  
    end
    
    % Save reevaluated DLC data
    file_path = DLC.FILE.file_path;
    file_name = DLC.FILE.file_name;
    session_type = DLC.FILE.session_type;
    with_nose_markers = DLC.FILE.with_nose_markers;
    scal = 7/47;
    origin_x = DLC.POINTS.origin_x;
    origin_y = DLC.POINTS.origin_y;
        
    CSV(:,1) = [1:length(DLC.POINTS.tip_tongue_x)];
    CSV(:,2) = (DLC.POINTS.tip_tongue_x + origin_x) / scal;
    CSV(:,3) = (DLC.POINTS.tip_tongue_y + origin_y) / scal;
    CSV(:,5) = (DLC.POINTS.r_tongue_x + origin_x) / scal;
    CSV(:,6) = (DLC.POINTS.r_tongue_y + origin_y) / scal;
    CSV(:,8) = (DLC.POINTS.l_tongue_x + origin_x) / scal;
    CSV(:,9) = (DLC.POINTS.l_tongue_y + origin_y) / scal;
    CSV(:,11) = (DLC.POINTS.mid_tongue_x + origin_x) / scal;
    CSV(:,12) = (DLC.POINTS.mid_tongue_y + origin_y) / scal;
    if session_type == 2 && with_nose_markers
        CSV(:,14) = (DLC.POINTS.r_nose_x + origin_x) / scal;
        CSV(:,15) = (DLC.POINTS.r_nose_y + origin_y) / scal;
        CSV(:,17) = (DLC.POINTS.l_nose_x + origin_x) / scal;
        CSV(:,18) = (DLC.POINTS.l_nose_y + origin_y) / scal;
        CSV(:,20) = (DLC.POINTS.r_food_x + origin_x) / scal;
        CSV(:,21) = (DLC.POINTS.r_food_y + origin_y) / scal;
        CSV(:,23) = (DLC.POINTS.l_food_x + origin_x) / scal;
        CSV(:,24) = (DLC.POINTS.l_food_y + origin_y) / scal;
        CSV(:,26) = (DLC.POINTS.r_tube_r_x + origin_x) / scal;
        CSV(:,27) = (DLC.POINTS.r_tube_r_y + origin_y) / scal;
        CSV(:,29) = (DLC.POINTS.r_tube_l_x + origin_x) / scal;
        CSV(:,30) = (DLC.POINTS.r_tube_l_y + origin_y) / scal;            
        CSV(:,32) = (DLC.POINTS.l_tube_r_x + origin_x) / scal;
        CSV(:,33) = (DLC.POINTS.l_tube_r_y + origin_y) / scal;
        CSV(:,35) = (DLC.POINTS.l_tube_l_x + origin_x) / scal;
        CSV(:,36) = (DLC.POINTS.l_tube_l_y + origin_y) / scal;
    elseif session_type == 2 && with_nose_markers == 0
        CSV(:,14) = (DLC.POINTS.r_food_x + origin_x) / scal;
        CSV(:,15) = (DLC.POINTS.r_food_y + origin_y) / scal;
        CSV(:,17) = (DLC.POINTS.l_food_x + origin_x) / scal;
        CSV(:,18) = (DLC.POINTS.l_food_y + origin_y) / scal;
        CSV(:,20) = (DLC.POINTS.r_tube_r_x + origin_x) / scal;
        CSV(:,21) = (DLC.POINTS.r_tube_r_y + origin_y) / scal;
        CSV(:,23) = (DLC.POINTS.r_tube_l_x + origin_x) / scal;
        CSV(:,24) = (DLC.POINTS.r_tube_l_y + origin_y) / scal;            
        CSV(:,26) = (DLC.POINTS.l_tube_r_x + origin_x) / scal;
        CSV(:,27) = (DLC.POINTS.l_tube_r_y + origin_y) / scal;
        CSV(:,29) = (DLC.POINTS.l_tube_l_x + origin_x) / scal;
        CSV(:,30) = (DLC.POINTS.l_tube_l_y + origin_y) / scal;
    elseif session_type == 1 && with_nose_markers 
        CSV(:,14) = (DLC.POINTS.r_nose_x + origin_x) / scal;
        CSV(:,15) = (DLC.POINTS.r_nose_y + origin_y) / scal;
        CSV(:,17) = (DLC.POINTS.l_nose_x + origin_x) / scal;
        CSV(:,18) = (DLC.POINTS.l_nose_y + origin_y) / scal;
        CSV(:,20) = (DLC.POINTS.r_food_x + origin_x) / scal;
        CSV(:,21) = (DLC.POINTS.r_food_y + origin_y) / scal;
        CSV(:,23) = (DLC.POINTS.r_tube_r_x + origin_x) / scal;
        CSV(:,24) = (DLC.POINTS.r_tube_r_y + origin_y) / scal;
        CSV(:,26) = (DLC.POINTS.r_tube_l_x + origin_x) / scal;
        CSV(:,27) = (DLC.POINTS.r_tube_l_y + origin_y) / scal; 
    elseif session_type == 1 && with_nose_markers == 0
        CSV(:,14) = (DLC.POINTS.r_food_x + origin_x) / scal;
        CSV(:,15) = (DLC.POINTS.r_food_y + origin_y) / scal;
        CSV(:,17) = (DLC.POINTS.r_tube_r_x + origin_x) / scal;
        CSV(:,18) = (DLC.POINTS.r_tube_r_y + origin_y) / scal;
        CSV(:,20) = (DLC.POINTS.r_tube_l_x + origin_x) / scal;
        CSV(:,21) = (DLC.POINTS.r_tube_l_y + origin_y) / scal; 
    end
    writematrix(CSV, [file_path file_name(1:length(file_name) - 4) '_Revaluated' '.csv']);
    
    if session_type == 2
        if ~isempty(ind_inaccurate_r_food) || ~isempty(ind_inaccurate_l_food) 
            DLC = evaluate(DLC); 
        end 
    else
        if ~isempty(ind_inaccurate_r_food)
            DLC = evaluate(DLC); 
        end 
    end

    fprintf(' --> Completed. \n')
    
    function [points_x, points_y] = correct_stationary_markers(x, y, points_x, points_y)
        rect = [mean(x) - 5 mean(y) - 5 10 10];
        ind_inaccurate = find(x < rect(1) | x > rect(1) + rect(3) | y < rect(2) | y > rect(2) + rect(4));
        
        if ~isempty(ind_inaccurate)
            i = 1;
            while i <= size(ind_inaccurate, 1)
                ind = ind_inaccurate(i);
                n = ind + 1;
                while ismember(n, ind_inaccurate)
                    n = n + 1;
                    i = i + 1;
                end
                if n == ind + 1
                    i = i + 1;
                end
                
                if ind == 1 || n > length(points_x)
                    points_x(ind) = mean(points_x);
                    points_y(ind) =  mean(points_y);
                    i = i + 1;
                else
                    for j = ind : n - 1
                        points_x(j) = points_x(j-1) + ((points_x(n) - points_x(ind - 1)) / (n - ind + 1));
                        points_y(j) = points_y(j-1) + ((points_y(n) - points_y(ind - 1)) / (n - ind + 1));
                    end
                end
            end
        end
    end
end