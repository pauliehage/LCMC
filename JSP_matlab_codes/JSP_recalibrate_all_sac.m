function [SACS_ALL_DATA_M_RECAL,TRIALS_DATA_M_RECAL,EXPERIMENT_PARAMS_M, CAL_MATRIX] = JSP_recalibrate_all_sac(SACS_ALL_DATA_M,TRIALS_DATA_M,EXPERIMENT_PARAMS_M, flag_use_avg,min_num_sample, params, funcs, pre_cal_matrix)
%JSP_RECALIBRATE_ALL_SAC
%
%ARGS:
%   BEHAVE_DATA_SORTED - array of behavior data where each element is a structure array for each behavior session,
%                        with fields consisting of exp. parameters, primary and corrective saccade, and
%                        other trial data (struct)
%                      - data in each field is either a vector or matrix where each column is data
%                        for each trial
%                      - output from 'backstep_exp_preprocess_behave_data.m'
%                      - data are sorted for the presence of corrective saccade
%                      - 'BEHAVE_DATA_SORTED.TRIALS_DATA.is_corrSac' is a vector where 
%                        each element indicates the presence of corrective saccade
%                        0: no corrSac; 1: corrSac; 2: ambiguous or invalid
%                           - use this info. as sole determinant of validity of trials
%   use_avg - indicates whether to use the mean pos. of corr. sac. to recalibrate (bool)
%   min_num_sample - minimum number of samples of corr. sac. per unique end target position to count as it is a reliable position to be used for calibration
%                  - mitigates few outlier corr. sac. offset position to bias the calibration
%   pre_cal_matrix - 3 x 3 calibration matrix (optional)
%OUTPUTS:
%   BEHAVE_DATA_RECAL - same as above, except the following behavior data for RIGHT EYE are recalibrated:
%                       - primSacOnset, primSacOffset, primSac, corrSacOnset, corrSacOffset, corrSac 
%NOTE: 
%   Each of recalibrated behavior session data is saved in its respective folder as specified in 
%   'EXPERIMENT_PARAMS.mat_PathName' folder, and its file name contains '_ANALYZED_corrSac_sorted_RECAL' tag at the end.
%   Data are recalibrated primarily using corrective saccade offset positions
%   as the references with respect to end target positions.
%   If different offsets or calibration matrix used during experiment across sessions,
%   different session needs to be separately recalibrated.
%   The number of trials per end target used to compute calibration matrix needs to be roughly equal.
%   Otherwise, calibration matrix is biased to the data with more number of examples 
%% PREPROCESS DATA
num_recording = length(SACS_ALL_DATA_M);
field_names_SACS_ALL_DATA = fieldnames(SACS_ALL_DATA_M(1));
for counter_field = 1 : 1 : length(field_names_SACS_ALL_DATA)
    field_name_SACS_ALL_DATA = field_names_SACS_ALL_DATA{counter_field};
    SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = [];
    for counter_recording = 1 : 1 : num_recording
        SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = horzcat(...
            SACS_ALL_DATA.(field_name_SACS_ALL_DATA), ...
            SACS_ALL_DATA_M(counter_recording).(field_name_SACS_ALL_DATA));
    end
end
%% CENTER POSITIONS
%% GET VARIABLES
sac_cal_tag = 4; % successful corr. sac.
sac_cal_ind = (SACS_ALL_DATA.tag == sac_cal_tag);% | (SACS_ALL_DATA.tag == 1);
fprintf('Num. samples for calibration (before): %d\n',length(find(SACS_ALL_DATA.tag == sac_cal_tag)));
tgt_px = SACS_ALL_DATA.end_x(sac_cal_ind);
tgt_py = SACS_ALL_DATA.end_y(sac_cal_ind);
tgt_p = [tgt_px', tgt_py']; 
eye_px = SACS_ALL_DATA.eye_r_px_offset(sac_cal_ind);
eye_py = SACS_ALL_DATA.eye_r_py_offset(sac_cal_ind); 
%% CHECK TO MAKE SURE CALIBRATION WAS NOT CHANGED DURING EXPERIMENT

%% FIND THE MEAN OF EYE POS. FOR EACH UNIQUE TGT. POS. (OPTIONAL)
if flag_use_avg
    [unique_tgt_pos,ia,unique_tgt_grp] =  unique(round(tgt_p,2),'rows');
    num_unique_tgt_pos = size(unique_tgt_pos,1);
    for p_idx = 1:num_unique_tgt_pos
        unique_tgt_ind = (unique_tgt_grp' == p_idx); % indices of each of unique tgt. pos.
        unique_eye_px = eye_px(unique_tgt_ind);
        unique_eye_py = eye_py(unique_tgt_ind);
        % If there are less than min. required number of samples for a position, do not compute the mean
        if length(unique_eye_px) < min_num_sample
            mean_eye_px(1,p_idx) = NaN;
            mean_eye_py(1,p_idx) = NaN;
        else
            mean_eye_px(1,p_idx) = mean(unique_eye_px); % mean of unique pos. w/o any data will be NaN
            mean_eye_py(1,p_idx) = mean(unique_eye_py);
        end
    end
    non_nan_idx = find(~isnan(mean_eye_px));
    eye_px = mean_eye_px(non_nan_idx);
    eye_py = mean_eye_py(non_nan_idx);
    tgt_px = unique_tgt_pos(non_nan_idx,1);
    tgt_px = tgt_px';
    tgt_py = unique_tgt_pos(non_nan_idx,2);
    tgt_py = tgt_py';
    % Plot the mean pos.
    avg_fig = figure('Name','Average of eye pos.','NumberTitle','off');
    avg_axes = axes(avg_fig);
    hold(avg_axes,'on');
    title(avg_axes, 'Average position of eye pos.');
    axis(avg_axes,'equal');
    xlim(avg_axes,[-12,12]); ylim(avg_axes,[-12,12]);
    xlabel(avg_axes, 'Horizontal Eye Position (deg)'); ylabel(avg_axes, 'Vertical Eye Position (deg)');
    plot(avg_axes, tgt_px,tgt_py,'dk');
    plot(avg_axes, eye_px, eye_py, '.r', 'MarkerSize',14);
end
%% EYE TRAJECTORY BEFORE RECALIBRATION
figure('WindowState','maximized');
pre_ax = subplot(1,2,1);
hold(pre_ax,'on');
title(pre_ax, 'Before recalibration');
axis(pre_ax,'equal');
xlim(pre_ax,[-12,12]); ylim(pre_ax,[-12,12]);
xlabel(pre_ax, 'Horizontal Eye Position (deg)'); ylabel(pre_ax, 'Vertical Eye Position (deg)');
% Experiment parameters
plot(pre_ax, SACS_ALL_DATA.start_x,SACS_ALL_DATA.start_y,'dk','MarkerFaceColor','k');
plot(pre_ax, SACS_ALL_DATA.cue_x,SACS_ALL_DATA.cue_y,'dk','MarkerFaceColor','k');
plot(pre_ax, SACS_ALL_DATA.end_x,SACS_ALL_DATA.end_y,'dk');
% Prim. sac. and corr. sac. finish pos.
prim_sac_tag = 1;
prim_sac_ind = SACS_ALL_DATA.tag == prim_sac_tag;
corr_sac_tag = 4;
corr_sac_ind = SACS_ALL_DATA.tag == corr_sac_tag;
plot(pre_ax, SACS_ALL_DATA.eye_r_px_offset(prim_sac_ind), SACS_ALL_DATA.eye_r_py_offset(prim_sac_ind),'.b','MarkerSize',14);
plot(pre_ax, SACS_ALL_DATA.eye_r_px_offset(corr_sac_ind), SACS_ALL_DATA.eye_r_py_offset(corr_sac_ind),'.r','MarkerSize',14);
%% COMPUTE CALIBRATION MATRIX
eye_p = [eye_px', eye_py', ones(size(eye_px'))];
tgt_p = [tgt_px', tgt_py', ones(size(tgt_px'))];
cal_matrix = ((eye_p' * eye_p)^(-1)) * (eye_p') * tgt_p;
%% OPTIONAL INPUT CALIBRATION MATRIX
if nargin == 8
    cal_matrix = pre_cal_matrix;
end
%% RECALIBRATE 
variable_list = { ...
                'eye_r_px', 'eye_r_py'; ...
                'eye_r_px_filt', 'eye_r_py_filt'; ...
                'eye_r_vx','eye_r_vy';...
                'eye_r_vx_filt','eye_r_vy_filt'};
for r_idx = 1:num_recording
    TRIALS_DATA = TRIALS_DATA_M(r_idx);
    num_trial = length(TRIALS_DATA.start_x);
    for t_idx = 1:num_trial
        for v_idx = 1:size(variable_list,1)
            x_to_cal = TRIALS_DATA.(variable_list{v_idx,1}){t_idx};
            y_to_cal = TRIALS_DATA.(variable_list{v_idx,2}){t_idx};
            data_to_cal = [x_to_cal, y_to_cal, ones(size(x_to_cal))];
            data_cal = data_to_cal * cal_matrix;
            x_cal = data_cal(:,1);
            y_cal = data_cal(:,2);
            TRIALS_DATA.(variable_list{v_idx,1}){t_idx} = x_cal;
            TRIALS_DATA.(variable_list{v_idx,2}){t_idx} = y_cal;
            if contains(variable_list{v_idx,1},'eye_r_v')
                if contains(variable_list{v_idx,1} ,'filt')
                    variable_name = [variable_list{v_idx,1}(1:7),'m_filt'];
                else
                    variable_name = [variable_list{v_idx,1}(1:7),'m'];
                end
                TRIALS_DATA.(variable_name){t_idx} = sqrt(x_cal.^2 + y_cal.^2);
            end
        end
    end
    TRIALS_DATA_M(r_idx) = TRIALS_DATA;
end
%% RE-RUN SAC SORTER
fprintf('Recalibrating...\n');
for r_idx = 1:num_recording
    [SACS_ALL_DATA, TRIALS_DATA, EXPERIMENT_PARAMS] = funcs.Sac_Sorter(TRIALS_DATA_M(r_idx),EXPERIMENT_PARAMS_M(r_idx),false);
    SACS_ALL_DATA_M(r_idx) = SACS_ALL_DATA;
    TRIALS_DATA_M(r_idx) = TRIALS_DATA;
    EXPERIMENT_PARAMS_M(r_idx) = EXPERIMENT_PARAMS;
    
end
%% GET VARIABLES AGAIN
clearvars SACS_ALL_DATA TRIALS_DATA EXPERIMENT_PARAMS
field_names_SACS_ALL_DATA = fieldnames(SACS_ALL_DATA_M(1));
for counter_field = 1 : 1 : length(field_names_SACS_ALL_DATA)
    field_name_SACS_ALL_DATA = field_names_SACS_ALL_DATA{counter_field};
    SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = [];
    for counter_recording = 1 : 1 : num_recording
        SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = horzcat(...
            SACS_ALL_DATA.(field_name_SACS_ALL_DATA), ...
            SACS_ALL_DATA_M(counter_recording).(field_name_SACS_ALL_DATA));
    end
end
sac_cal_ind = SACS_ALL_DATA.tag == sac_cal_tag;
fprintf('Num. samples for calibration (after): %d\n',length(find(sac_cal_ind)));
tgt_px = SACS_ALL_DATA.tgt_px_offset(sac_cal_ind);
tgt_py = SACS_ALL_DATA.tgt_py_offset(sac_cal_ind);
tgt_p = [tgt_px', tgt_py']; 
eye_px = SACS_ALL_DATA.eye_r_px_offset(sac_cal_ind);
eye_py = SACS_ALL_DATA.eye_r_py_offset(sac_cal_ind); 
%% EYE TRAJECTORY AFTER RECALIBRATION
post_ax = subplot(1,2,2);
hold(post_ax,'on');
title(post_ax, 'After recalibration');
axis(post_ax,'equal');
xlim(post_ax,[-12,12]); ylim(post_ax,[-12,12]);
xlabel(post_ax, 'Horizontal Eye Position (deg)'); ylabel(post_ax, 'Vertical Eye Position (deg)');
% Experiment parameters
plot(post_ax, SACS_ALL_DATA.start_x,SACS_ALL_DATA.start_y,'dk','MarkerFaceColor','k');
plot(post_ax, SACS_ALL_DATA.cue_x,SACS_ALL_DATA.cue_y,'dk','MarkerFaceColor','k');
plot(post_ax, SACS_ALL_DATA.end_x,SACS_ALL_DATA.end_y,'dk');
% Prim. sac. and corr. sac. finish pos.
prim_sac_tag = 1;
prim_sac_ind = SACS_ALL_DATA.tag == prim_sac_tag;
corr_sac_tag = 4;
corr_sac_ind = SACS_ALL_DATA.tag == corr_sac_tag;
plot(post_ax, SACS_ALL_DATA.eye_r_px_offset(prim_sac_ind),SACS_ALL_DATA.eye_r_py_offset(prim_sac_ind),'.b','MarkerSize',14);
plot(post_ax, SACS_ALL_DATA.eye_r_px_offset(corr_sac_ind),SACS_ALL_DATA.eye_r_py_offset(corr_sac_ind),'.r','MarkerSize',14);
%% RETURN DATA
SACS_ALL_DATA_M_RECAL = SACS_ALL_DATA_M;
TRIALS_DATA_M_RECAL = TRIALS_DATA_M;
CAL_MATRIX = cal_matrix;
end
