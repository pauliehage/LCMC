%% Combines ESN_monkey_behavior with JSP_Sac_Sorter
function [SACS_ALL_DATA, TRIALS_DATA, EXPERIMENT_PARAMS] = MAF_monkey_behavior_sac(mat_file_address, flag_figure, params, funcs)
%% Get file_name and file_path
[file_path,file_name,file_ext] = fileparts(mat_file_address);
file_name = [file_name file_ext];
% add filesep ('/' or '\') to the end of file_path
if ~strcmp(file_path(end), filesep)
    file_path = [file_path filesep];
end

EXPERIMENT_PARAMS.mat_FileName = file_name;
EXPERIMENT_PARAMS.mat_PathName = file_path;
[~, EXPERIMENT_PARAMS.file_name, ~] = fileparts(file_name);

%% Load Data
fprintf('Loading ...\n')
clearvars -except EXPERIMENT_PARAMS TRIALS_DATA_ALL flag_figure funcs;
filename = EXPERIMENT_PARAMS.mat_FileName;
pathname = EXPERIMENT_PARAMS.mat_PathName;
load([pathname filename], 'data');
fprintf('%s: Loading complete.\n', EXPERIMENT_PARAMS.file_name);
[~, foldername, ~] = fileparts(pathname(1:end-1)); % (1:end-1) to delete the '/' character
EXPERIMENT_PARAMS.folder_name = foldername;
EXPERIMENT_PARAMS.num_trials  = length(data.trials);

%% Analyze trials
clearvars -except EXPERIMENT_PARAMS TRIALS_DATA_ALL data flag_figure funcs;
num_trials = length(data.trials);
fprintf([EXPERIMENT_PARAMS.file_name ': Analyzing TRIALS ...'])
for counter_trial = 1 : 1 : num_trials-1
    %% Extract Trial Varibales
    clearvars -except EXPERIMENT_PARAMS ...
        TRIALS_DATA_ALL TRIALS data counter_trial flag_figure funcs;
    if ~exist('counter_trial','var'); counter_trial = 1; end
    % get trial struct
    trial_struct = data.trials{1, counter_trial};
    % trial start/end time
    TRIAL.time_start = double(trial_struct.trial_start_time);
    TRIAL.time_end   = double(trial_struct.trial_end_time);
    % trial variables
    TRIAL.start_x                    = ESN_Round(trial_struct.start_x, 0.01);
    TRIAL.start_y                    = ESN_Round(trial_struct.start_y, 0.01);
    TRIAL.cue_x                      = ESN_Round(trial_struct.cue_x, 0.01);
    TRIAL.cue_y                      = ESN_Round(trial_struct.cue_y, 0.01);
    TRIAL.end_x                      = ESN_Round(trial_struct.end_x, 0.01);
    TRIAL.end_y                      = ESN_Round(trial_struct.end_y, 0.01);
    TRIAL.iss_x                      = ESN_Round(trial_struct.iss_x, 0.01);
    TRIAL.iss_y                      = ESN_Round(trial_struct.iss_y, 0.01);
    if isfield(TRIAL, 'pursuit_x')
        TRIAL.pursuit_x                  = ESN_Round(trial_struct.pursuit_x, 0.01);
        TRIAL.pursuit_y                  = ESN_Round(trial_struct.pursuit_y, 0.01);
        TRIAL.reward_area                = ESN_Round(trial_struct.reward_area, 0.01);
        TRIAL.time_pursuit               = double(trial_struct.pursuit_duration);
    else
        TRIAL.pursuit_x                  = nan;
        TRIAL.pursuit_y                  = nan;
        TRIAL.reward_area                = nan;
        TRIAL.time_pursuit               = nan;
    end
    TRIAL.time_state_str_pursuit     = double(trial_struct.state_start_time_str_target_pursuit);
    TRIAL.time_state_str_present     = double(trial_struct.state_start_time_str_target_present);
    TRIAL.time_state_str_fixation    = double(trial_struct.state_start_time_str_target_fixation);
    TRIAL.time_state_cue_present     = double(trial_struct.state_start_time_cue_target_present);
    TRIAL.time_state_sac_detect_on   = double(trial_struct.state_start_time_detect_sac_start);
    TRIAL.time_state_sac_onset       = double(trial_struct.state_start_time_saccade);
    TRIAL.time_state_sac_detect_off  = double(trial_struct.state_start_time_detect_sac_end);
    TRIAL.time_state_reward          = double(trial_struct.state_start_time_deliver_reward);
    TRIAL.time_state_end_fixation    = double(trial_struct.state_start_time_end_target_fixation);
    TRIAL.time_state_iti             = double(trial_struct.state_start_time_iti);
    TRIAL.time_state_next_trial      = double(trial_struct.state_start_time_next_trial);
    TRIAL.time_iti                   = double(trial_struct.iti);
    TRIAL.time_punishment            = double(trial_struct.punishment_time);
    TRIAL.time_fixation              = double(trial_struct.fixation_time);
    % trial start/end indices
    length_data = min([length(data.eyelink_time) length(data.t)]);
    time_array = double(data.t(1: length_data));
    TRIAL.ind_trial_str   = find(time_array>TRIAL.time_start, 1, 'first');
    TRIAL.ind_trial_end   = find(time_array<TRIAL.time_end, 1, 'last');
    TRIAL.inds_trial          = TRIAL.ind_trial_str:TRIAL.ind_trial_end;
    % trial timeseries
    TRIAL.inds_invalid   = false(1, length(TRIAL.inds_trial));
    TRIAL.time_eyelink   = double(data.eyelink_time(1, TRIAL.inds_trial));                           TRIAL.inds_invalid = isnan(TRIAL.time_eyelink) | TRIAL.inds_invalid;
    TRIAL.eye_r_px       = double(data.right_horizontal_eye(1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_r_px)     | TRIAL.inds_invalid;
    TRIAL.eye_r_py       = double(data.right_vertical_eye(  1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_r_py)     | TRIAL.inds_invalid;
    TRIAL.eye_l_px       = double(data.left_horizontal_eye( 1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_l_px)     | TRIAL.inds_invalid;
    TRIAL.eye_l_py       = double(data.left_vertical_eye(   1, TRIAL.inds_trial));                   TRIAL.inds_invalid = isnan(TRIAL.eye_l_py)     | TRIAL.inds_invalid;
    TRIAL.eye_r_vx       = double(data.right_horizontal_eye_velocity_filtered(1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_r_vx)     | TRIAL.inds_invalid;
    TRIAL.eye_r_vy       = double(data.right_vertical_eye_velocity_filtered(  1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_r_vy)     | TRIAL.inds_invalid;
    TRIAL.eye_l_vx       = double(data.left_horizontal_eye_velocity_filtered( 1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_l_vx)     | TRIAL.inds_invalid;
    TRIAL.eye_l_vy       = double(data.left_vertical_eye_velocity_filtered(   1, TRIAL.inds_trial)); TRIAL.inds_invalid = isnan(TRIAL.eye_l_vy)     | TRIAL.inds_invalid;
    TRIAL.eye_r_vm       = sqrt(TRIAL.eye_r_vx.^2 + TRIAL.eye_r_vy.^2);
    TRIAL.eye_l_vm       = sqrt(TRIAL.eye_l_vx.^2 + TRIAL.eye_l_vy.^2);
    TRIAL.time_tgt       = double(data.t(1, TRIAL.inds_trial));                                      TRIAL.inds_invalid = isnan(TRIAL.time_tgt)     | TRIAL.inds_invalid;
    TRIAL.tgt_px         = double(data.target_x(1, TRIAL.inds_trial));                               TRIAL.inds_invalid = isnan(TRIAL.tgt_px)       | TRIAL.inds_invalid;
    TRIAL.tgt_py         = double(data.target_y(1, TRIAL.inds_trial));                               TRIAL.inds_invalid = isnan(TRIAL.tgt_py)       | TRIAL.inds_invalid;
    TRIAL.reward         = double(data.reward(1, TRIAL.inds_trial));                                 TRIAL.inds_invalid = isnan(TRIAL.reward)       | TRIAL.inds_invalid;
    TRIAL.target_visible = logical(double(data.target_visible(1, TRIAL.inds_trial)));
    % correct for the bias between time_eyelink and time_tgt
    TRIAL.time_eyelink   = TRIAL.time_eyelink .* (TRIAL.time_tgt(end)-TRIAL.time_tgt(1)) ./ (TRIAL.time_eyelink(end)-TRIAL.time_eyelink(1));
    TRIAL.time_eyelink   = TRIAL.time_eyelink - TRIAL.time_eyelink(1) + TRIAL.time_tgt(1);
    TRIAL.time_1K        = TRIAL.time_eyelink(1) : 0.001 : TRIAL.time_eyelink(end);
    % make non unique points of eye traces invalid
    TRIAL.inds_invalid = ([false (diff(TRIAL.time_eyelink)==0)])       | TRIAL.inds_invalid;
    % remove invalid values
    TRIAL.time_eyelink(TRIAL.inds_invalid) = [];
    TRIAL.eye_r_px(    TRIAL.inds_invalid) = [];
    TRIAL.eye_r_py(    TRIAL.inds_invalid) = [];
    TRIAL.eye_l_px(    TRIAL.inds_invalid) = [];
    TRIAL.eye_l_py(    TRIAL.inds_invalid) = [];
    % reconstruct eye_r data
    TRIAL.eye_r_px = interp1(TRIAL.time_eyelink, TRIAL.eye_r_px, TRIAL.time_1K, 'linear', 'extrap');
    TRIAL.eye_r_py = interp1(TRIAL.time_eyelink, TRIAL.eye_r_py, TRIAL.time_1K, 'linear', 'extrap');
    TRIAL.eye_r_vx = diff(TRIAL.eye_r_px)./diff(TRIAL.time_1K); TRIAL.eye_r_vx=[TRIAL.eye_r_vx(1) TRIAL.eye_r_vx];
    TRIAL.eye_r_vy = diff(TRIAL.eye_r_py)./diff(TRIAL.time_1K); TRIAL.eye_r_vy=[TRIAL.eye_r_vy(1) TRIAL.eye_r_vy];
    TRIAL.eye_r_vm = sqrt(TRIAL.eye_r_vx.^2 + TRIAL.eye_r_vy.^2);
    % reconstruct eye_l data
    TRIAL.eye_l_px = interp1(TRIAL.time_eyelink, TRIAL.eye_l_px, TRIAL.time_1K, 'linear', 'extrap');
    TRIAL.eye_l_py = interp1(TRIAL.time_eyelink, TRIAL.eye_l_py, TRIAL.time_1K, 'linear', 'extrap');
    TRIAL.eye_l_vx = diff(TRIAL.eye_l_px)./diff(TRIAL.time_1K); TRIAL.eye_l_vx=[TRIAL.eye_l_vx(1) TRIAL.eye_l_vx];
    TRIAL.eye_l_vy = diff(TRIAL.eye_l_py)./diff(TRIAL.time_1K); TRIAL.eye_l_vy=[TRIAL.eye_l_vy(1) TRIAL.eye_l_vy];
    TRIAL.eye_l_vm = sqrt(TRIAL.eye_l_vx.^2 + TRIAL.eye_l_vy.^2);
    TRIAL.time     = TRIAL.time_1K;
    % filter params
    sampling_freq = 1000.0;
    cutoff_freq = 100.0;
    [b_butter,a_butter] = butter(3,(cutoff_freq/(sampling_freq/2)), 'low');
    % filter eye_r data
    TRIAL.eye_r_px_filt = filtfilt(b_butter,a_butter,TRIAL.eye_r_px);
    TRIAL.eye_r_py_filt = filtfilt(b_butter,a_butter,TRIAL.eye_r_py);
    TRIAL.eye_r_vx_filt = diff(TRIAL.eye_r_px_filt)./diff(TRIAL.time_1K); TRIAL.eye_r_vx_filt=[TRIAL.eye_r_vx_filt(1) TRIAL.eye_r_vx_filt];
    TRIAL.eye_r_vy_filt = diff(TRIAL.eye_r_py_filt)./diff(TRIAL.time_1K); TRIAL.eye_r_vy_filt=[TRIAL.eye_r_vy_filt(1) TRIAL.eye_r_vy_filt];
    TRIAL.eye_r_vm_filt = sqrt(TRIAL.eye_r_vx_filt.^2 + TRIAL.eye_r_vy_filt.^2);
    % filter eye_l data
    TRIAL.eye_l_px_filt = filtfilt(b_butter,a_butter,TRIAL.eye_l_px);
    TRIAL.eye_l_py_filt = filtfilt(b_butter,a_butter,TRIAL.eye_l_py);
    TRIAL.eye_l_vx_filt = diff(TRIAL.eye_l_px_filt)./diff(TRIAL.time_1K); TRIAL.eye_l_vx_filt=[TRIAL.eye_l_vx_filt(1) TRIAL.eye_l_vx_filt];
    TRIAL.eye_l_vy_filt = diff(TRIAL.eye_l_py_filt)./diff(TRIAL.time_1K); TRIAL.eye_l_vy_filt=[TRIAL.eye_l_vy_filt(1) TRIAL.eye_l_vy_filt];
    TRIAL.eye_l_vm_filt = sqrt(TRIAL.eye_l_vx_filt.^2 + TRIAL.eye_l_vy_filt.^2);
    % trial state indices
    TRIAL.ind_state_str_pursuit    = find(TRIAL.time>TRIAL.time_state_str_pursuit(end), 1, 'first');
    TRIAL.ind_state_str_present    = find(TRIAL.time>TRIAL.time_state_str_present(end), 1, 'first');
    TRIAL.ind_state_str_fixation   = find(TRIAL.time>TRIAL.time_state_str_fixation(end), 1, 'first');
    TRIAL.ind_state_cue_present    = find(TRIAL.time>TRIAL.time_state_cue_present(end), 1, 'first');
    TRIAL.ind_state_sac_detect_on  = find(TRIAL.time>TRIAL.time_state_sac_detect_on(end), 1, 'first');
    TRIAL.ind_state_sac_onset      = find(TRIAL.time>TRIAL.time_state_sac_onset(end), 1, 'first');
    TRIAL.ind_state_sac_detect_off = find(TRIAL.time>TRIAL.time_state_sac_detect_off(end), 1, 'first');
    TRIAL.ind_state_reward         = find(TRIAL.time>TRIAL.time_state_reward(end), 1, 'first');
    TRIAL.ind_state_end_fixation   = find(TRIAL.time>TRIAL.time_state_end_fixation(end), 1, 'first');
    TRIAL.ind_state_iti            = find(TRIAL.time>TRIAL.time_state_iti(end), 1, 'first');
    TRIAL.ind_state_next_trial     = find(TRIAL.time>TRIAL.time_state_next_trial(end), 1, 'first');
    
    %% Compute the start position Bias
    inds_start_fixation = TRIAL.ind_state_cue_present - round(TRIAL.time_fixation * 1000) : 1 : TRIAL.ind_state_cue_present;
    state_start_eye_r_px_start_fixation = TRIAL.eye_l_px_filt(inds_start_fixation);
    state_start_eye_r_py_start_fixation = TRIAL.eye_l_py_filt(inds_start_fixation);
    tgt_start_x = TRIAL.start_x;
    tgt_start_y = TRIAL.start_y;
    start_x_bias = mean(state_start_eye_r_px_start_fixation) - tgt_start_x;
    start_y_bias = mean(state_start_eye_r_py_start_fixation) - tgt_start_y;
    TRIAL.start_x_bias = start_x_bias;
    TRIAL.start_y_bias = start_y_bias;
    
    %% Build TRIALS
    TRIALS(counter_trial) = TRIAL;
    
    % print a dot every 20 trials
    if rem(counter_trial, 20) == 0
        fprintf('.');
    end
    
end
fprintf(' --> Completed. \n')

%% Arrange 'TRIALS_DATA'
clearvars -except EXPERIMENT_PARAMS TRIALS_DATA_ALL TRIALS_DATA TRIALS flag_figure funcs
fprintf([EXPERIMENT_PARAMS.file_name ': Arranging TRIALS_DATA ...'])
clearvars('TRIALS_DATA'); TRIALS_DATA = struct;
field_names_TRIALS = fieldnames(TRIALS);
for counter_fields = 1 : 1 : length(field_names_TRIALS)
    for counter_trials = 1 : 1 : length(TRIALS)
        variable_TRIALS_ = TRIALS(counter_trials).(field_names_TRIALS{counter_fields});
        % handling an error which the variable_TRIALS_ was []
        if isempty(variable_TRIALS_)
            variable_TRIALS_ = nan;
        end
        variable_TRIALS_ = variable_TRIALS_(:);
        if max(size(variable_TRIALS_)) > 1
            variable_TRIALS_ = mat2cell(variable_TRIALS_, size(variable_TRIALS_,1), size(variable_TRIALS_,2));
        end
        % the field does not exist in TRIALS_DATA
        if ~isfield(TRIALS_DATA, field_names_TRIALS{counter_fields})
            TRIALS_DATA.(field_names_TRIALS{counter_fields}) = [];
        end
        variable_TRIALS_DATA_ = TRIALS_DATA.(field_names_TRIALS{counter_fields});
        % variable_TRIALS_ is cell array
        if iscell(variable_TRIALS_)
            % variable_TRIALS_DATA_ is cell array
            % variables are compatible (both are cell), add new data
            if iscell(variable_TRIALS_DATA_)
                variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
            end
            % variable_TRIALS_DATA_ is matrix array
            % convert variable_TRIALS_DATA_ to cell, and add new data
            if isnumeric(variable_TRIALS_DATA_)
                variable_TRIALS_DATA_ = num2cell(variable_TRIALS_DATA_);
                variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
            end
        end
        % variable_TRIALS_ is matrix array
        if isnumeric(variable_TRIALS_)
            % variable_TRIALS_DATA_ is cell array
            % convert variable_TRIALS_ to cell, and add new data
            if iscell(variable_TRIALS_DATA_)
                variable_TRIALS_ = mat2cell(variable_TRIALS_, size(variable_TRIALS_,1), size(variable_TRIALS_,2));
                variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
            end
            % variable_TRIALS_DATA_ is matrix array
            % variables are compatible (both are matrix), add new data
            if isnumeric(variable_TRIALS_DATA_)
                variable_TRIALS_DATA_(1, counter_trials) = variable_TRIALS_;
            end
        end
        TRIALS_DATA.(field_names_TRIALS{counter_fields}) = variable_TRIALS_DATA_;
    end
    fprintf('.');
end
fprintf(' --> Completed. \n')

%% Build TRIALS_DATA_ALL
TRIALS_DATA_ALL    = TRIALS_DATA;

%% Arrange 'TRIALS_DATA'
clearvars -except EXPERIMENT_PARAMS TRIALS_DATA_ALL TRIALS_DATA flag_figure funcs
fprintf('Arranging TRIALS_DATA ...')
clearvars('TRIALS_DATA'); TRIALS_DATA = struct;
field_names_TRIALS_DATA_ALL = fieldnames(TRIALS_DATA_ALL);
for counter_fields = 1 : 1 : length(field_names_TRIALS_DATA_ALL)
    for counter_files = 1 : 1 : length(TRIALS_DATA_ALL)
        variable_TRIALS_DATA_ALL_ = TRIALS_DATA_ALL(counter_files).(field_names_TRIALS_DATA_ALL{counter_fields});
        % the field does not exist in TRIALS_DATA
        if ~isfield(TRIALS_DATA, field_names_TRIALS_DATA_ALL{counter_fields})
            TRIALS_DATA.(field_names_TRIALS_DATA_ALL{counter_fields}) = [];
        end
        variable_TRIALS_DATA_ = TRIALS_DATA.(field_names_TRIALS_DATA_ALL{counter_fields});
        variable_TRIALS_DATA_ = horzcat(variable_TRIALS_DATA_, variable_TRIALS_DATA_ALL_);
        TRIALS_DATA.(field_names_TRIALS_DATA_ALL{counter_fields}) = variable_TRIALS_DATA_;
    end
    fprintf('.');
end
fprintf(' --> Completed. \n')

[SACS_ALL_DATA, TRIALS_DATA, EXPERIMENT_PARAMS] = funcs.Sac_Sorter(TRIALS_DATA, EXPERIMENT_PARAMS, flag_figure);
