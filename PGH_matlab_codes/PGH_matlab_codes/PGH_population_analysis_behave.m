%% MASTER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function PGH_population_analysis
function PGH_population_analysis_behave
clc; clear; close all;
tic
path_data_monkey_sorted = 'data_125d';
%% params
params.sac.length_trace = 500;
params.sac.inds_span        = ((-(params.sac.length_trace/2)+1) : 1 : (params.sac.length_trace/2))';
params.sac.ang_step         = 45;
params.sac.ang_edges        = (0 - (params.sac.ang_step/2)) : params.sac.ang_step : (360 + (params.sac.ang_step/2));
params.sac.ang_values       = (0) : params.sac.ang_step : (360 - params.sac.ang_step);
params.sac.tags_CS_ang_avg  = [1,4,6];
params.sac.tag_name_list    = { ...
    'prim_success', ... % tag 1
    'prim_attempt', ... % tag 2
    'prim_fail', ... % tag 3
    'corr_success', ... % tag 4
    'corr_fail', ... % tag 5
    'back_center_success', ... % tag 6
    'back_center_prim', ... % tag 7
    'back_center_irrelev', ... % tag 8
    'target_irrelev', ... % tag 9
    'other_irrelev', ... % tag 10
    'prim_no_corr',... % tag 11; prim. sac. that is not followed by corr. sac.
    'db_corr_success',... % tag 12; sac. that follows first corr. sac., back to 2nd jumped cue
    'corr_no_db_corr',... % tag 13; corr. sac. that is not followed by another corr. sac.
    'other_irrelev_visual',... % tag 14; like tag 10, but visual ang. based on the pursuit target present on the screen
    'back_center_irrelev_visual',... % tag 15; like tag 10, but visual start position based on offset of previous saccade
    ... % Add more tags here, do not reorder or change the tags defined above.
    };

params.lick.length_trace    = 2000;
params.lick.inds_span       = ((-(params.lick.length_trace/2)+1) : 1 : (params.lick.length_trace/2))';
params.lick.ang_step        = 45;
params.lick.ang_edges       = -90 - params.lick.ang_step /2 : params.lick.ang_step  : 90 + params.lick.ang_step /2 ;
params.lick.ang_values      = -90 : params.lick.ang_step : 90;
params.lick.amp_step    = 5;
params.lick.amp_edges    = 0 : params.lick.amp_step : 25;
params.lick.vel_step    = 150;
params.lick.vel_edges    = 0 : params.lick.vel_step : 750;
params.lick.dur_step = 50;
params.lick.dur_edges = 0 : params.lick.dur_step : 500;
params.lick.tags_CS_ang_avg  = [1:5];
params.lick.tag_name_list = {  ...
    'groom', ... % tag 1
    'inner_tube_success', ... % tag 2
    'inner_tube_fail', ... % tag 3
    'outer_tube_success' ..., % tag 4
    'outer_tube_fail', ... % tag 5
    };
params.lick.tag_bout_name_list = {  ...
    'bout_start', ... % tag 1
    'bout_end', ... % tag 2
    };
params.lick.tag_harvest_name_list = {  ...
    'harvest_start',  ...% tag 1
    'harvest_end', ... % tag 2
    };
params.lick.line_colors = [0,0,0; pink(round(1.5*length(params.lick.amp_edges)))];



%% Build functions
%extract_population_data(path_data_monkey_sorted);
%build_population_data(path_data_monkey_sorted)
%generate_population_meta_data(path_data_monkey_sorted)
%% Plot Functions
%plot_session_analysis(params,path_data_monkey_sorted)
plot_harvest_analysis(params,path_data_monkey_sorted)

toc
end

%% UTIL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function extract_population_data
function extract_population_data(path_data_monkey_sorted)
out_path = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted filesep 'session_data' filesep];
mkdir(out_path);

current_monkey_path_ = path_data_monkey_sorted;
fprintf(['### ' 'Analyzing subject ', current_monkey_path_ ' ### \n']);
if ~strcmp(current_monkey_path_(end), filesep)
    current_monkey_path = [current_monkey_path_ filesep];
end
session_list = dir([current_monkey_path filesep '20*' filesep '20*']);
session_list([session_list.isdir] == 0) = [];
sess_list = {session_list.name};

num_sess = length(sess_list);
num_trial_sess = nan(num_sess,1);
num_unit_sess = nan(num_sess,1);
probe_type_sess = nan(num_sess,1);
% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    clearvars EXPERIMENT_DATA LICKS_ALL_DATA SACS_ALL_DATA data_recording
    current_sess = sess_list{counter_sess};
    fprintf(['  ### ' 'sess ', current_sess, ' num. ',...
        num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [current_monkey_path, current_sess(1:7), filesep, current_sess, filesep];

    sess_name = regexprep(current_sess(3:end),'-','');
    sess_meta_data = readtable([sess_path sess_name '.xls']);

    rec_list = sess_meta_data.folder_name(logical...
        (sess_meta_data.ephys .* sess_meta_data.eye));

    n_recs = length(rec_list);

    % EXPERIMENT_DATA
    path_to_units = [sess_path 'units' filesep];

    % Load session trial count and electrode type from session meta dat
    session_meta_data_name_ = dir([sess_path '*.xls']);
    session_meta_data_name = session_meta_data_name_.name;
    session_meta_data = readtable([sess_path session_meta_data_name]);

    %EXPERIMENT_DATA.num_trial_sess_rec = (session_meta_data.num_trial);
    EXPERIMENT_DATA.num_trial_sess = sum(session_meta_data.num_trial);
    EXPERIMENT_DATA.probe_type_sess = session_meta_data.elec(1);
    EXPERIMENT_DATA.experiment_type_sess = session_meta_data.exp(1);

    % Count number of units from rec unit summary
    session_rec_unit_summary_name_ = dir([path_to_units '*.mat']);
    session_rec_unit_summary_name = session_rec_unit_summary_name_.name;
    session_rec_unit_summary = load([path_to_units session_rec_unit_summary_name]);

    EXPERIMENT_DATA.num_unit_sess = length(session_rec_unit_summary.cell_ids);
    EXPERIMENT_DATA.sess_date = string(current_sess(3:end));
    EXPERIMENT_DATA.id =  string([num2str(current_sess(3:4)) num2str(current_sess(6:7)) ...
        num2str(current_sess(9:10)) '_combine_' num2str(n_recs)]);
    for counter_rec = 1 : n_recs
        EXPERIMENT_DATA.rec_time(counter_rec, 1) = timeofday(datetime(rec_list{counter_rec}, 'InputFormat', 'yyyy-MM-dd_HH-mm-ss'));

        path_to_rec_eye = [sess_path rec_list{counter_rec} filesep 'analyzed_data' filesep 'behavior_data' filesep 'eye' filesep];
        current_rec_data_eye = dir([path_to_rec_eye '*_ANALYZED_RECAL.mat']);
        current_rec_data_eye_name = current_rec_data_eye.name;
        load([path_to_rec_eye current_rec_data_eye_name], 'TRIALS_DATA');

        EXPERIMENT_DATA.rec_duration(counter_rec, 1) = duration(seconds(TRIALS_DATA.time_end(end) - TRIALS_DATA.time_start(1)));

        path_to_rec_tongue = [sess_path rec_list{counter_rec} filesep 'analyzed_data' filesep 'behavior_data' filesep 'tongue' filesep];
        current_rec_data_tongue_align = dir([path_to_rec_tongue '*_aligned.mat']);
        current_rec_data_tongue_align_name = current_rec_data_tongue_align.name;
        load([path_to_rec_tongue current_rec_data_tongue_align_name], 'align_PD');
        EXPERIMENT_DATA.sample_diff(counter_rec, 1) = align_PD.sample_diff;
        EXPERIMENT_DATA.exp_start_time(counter_rec, 1) = align_PD.BEHAVE_PD_xcorr_time_1K(1,1);

        current_rec_data_tongue = dir([path_to_rec_tongue '*_ANALYZED.mat']);
        current_rec_data_tongue_name = current_rec_data_tongue.name;

        %LICKS_ALL_DATA
        load([path_to_rec_tongue current_rec_data_tongue_name], 'LICKS_ALL_DATA');

        data_recording(counter_rec).LICKS_ALL_DATA.tag = LICKS_ALL_DATA.tag;
        data_recording(counter_rec).LICKS_ALL_DATA.tag_bout = LICKS_ALL_DATA.tag_bout;
        data_recording(counter_rec).LICKS_ALL_DATA.tag_harvest = LICKS_ALL_DATA.tag_harvest;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_tip_px = LICKS_ALL_DATA.tongue_tip_px;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_tip_py = LICKS_ALL_DATA.tongue_tip_py;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_mid_px = LICKS_ALL_DATA.tongue_mid_px;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_mid_py = LICKS_ALL_DATA.tongue_mid_py;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_r_px = LICKS_ALL_DATA.tongue_r_px;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_r_py = LICKS_ALL_DATA.tongue_r_py;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_l_px = LICKS_ALL_DATA.tongue_l_px;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_l_py = LICKS_ALL_DATA.tongue_l_py;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_tip_px_dmax = LICKS_ALL_DATA.tongue_tip_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_tip_py_dmax = LICKS_ALL_DATA.tongue_tip_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_mid_px_dmax = LICKS_ALL_DATA.tongue_mid_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_mid_py_dmax = LICKS_ALL_DATA.tongue_mid_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_r_px_dmax = LICKS_ALL_DATA.tongue_r_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_r_py_dmax = LICKS_ALL_DATA.tongue_r_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_l_px_dmax = LICKS_ALL_DATA.tongue_l_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_l_py_dmax = LICKS_ALL_DATA.tongue_l_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.nose_r_px_dmax = LICKS_ALL_DATA.nose_r_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.nose_r_py_dmax = LICKS_ALL_DATA.nose_r_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.nose_l_px_dmax = LICKS_ALL_DATA.nose_l_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.nose_l_py_dmax = LICKS_ALL_DATA.nose_l_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_r_px_dmax = LICKS_ALL_DATA.rtube_r_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_r_py_dmax = LICKS_ALL_DATA.rtube_r_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_l_px_dmax = LICKS_ALL_DATA.rtube_l_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_l_py_dmax = LICKS_ALL_DATA.rtube_l_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_r_px_dmax = LICKS_ALL_DATA.ltube_r_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_r_py_dmax = LICKS_ALL_DATA.ltube_r_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_l_px_dmax = LICKS_ALL_DATA.ltube_l_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_l_py_dmax = LICKS_ALL_DATA.ltube_l_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_r_px_dmax = LICKS_ALL_DATA.rew_r_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_r_py_dmax = LICKS_ALL_DATA.rew_r_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_l_px_dmax = LICKS_ALL_DATA.rew_l_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_l_py_dmax = LICKS_ALL_DATA.rew_l_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.time_onset = LICKS_ALL_DATA.time_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.time_vmax = LICKS_ALL_DATA.time_vmax;
        data_recording(counter_rec).LICKS_ALL_DATA.time_dmax = LICKS_ALL_DATA.time_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.time_vmin = LICKS_ALL_DATA.time_vmin;
        data_recording(counter_rec).LICKS_ALL_DATA.time_offset = LICKS_ALL_DATA.time_offset;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_dm = LICKS_ALL_DATA.tongue_dm;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_vm = LICKS_ALL_DATA.tongue_vm;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_ang = LICKS_ALL_DATA.tongue_ang;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_dm_max = LICKS_ALL_DATA.tongue_dm_max;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_vm_max = LICKS_ALL_DATA.tongue_vm_max;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_vm_min = LICKS_ALL_DATA.tongue_vm_min;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_ang_max = LICKS_ALL_DATA.tongue_ang_max;
        data_recording(counter_rec).LICKS_ALL_DATA.duration_lick = LICKS_ALL_DATA.duration_lick;
        data_recording(counter_rec).LICKS_ALL_DATA.duration_harvest= LICKS_ALL_DATA.duration_harvest;
        data_recording(counter_rec).LICKS_ALL_DATA.duration_harvest= LICKS_ALL_DATA.duration_harvest;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_capacity_r_lick_onset= LICKS_ALL_DATA.rew_capacity_r_lick_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_capacity_r_lick_offset= LICKS_ALL_DATA.rew_capacity_r_lick_offset;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_capacity_l_lick_onset= LICKS_ALL_DATA.rew_capacity_l_lick_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_capacity_l_lick_offset= LICKS_ALL_DATA.rew_capacity_l_lick_offset;

        EXPERIMENT_DATA.num_licks_sess_rec(counter_rec, 1) = length(LICKS_ALL_DATA.tag);
        EXPERIMENT_DATA.num_harvests_sess_rec(counter_rec, 1) = sum(LICKS_ALL_DATA.tag_harvest == 1);

        %SACS_ALL_DATA
        load([path_to_rec_eye current_rec_data_eye_name], 'SACS_ALL_DATA');
        data_recording(counter_rec).SACS_ALL_DATA.tag = SACS_ALL_DATA.tag  ;
        data_recording(counter_rec).SACS_ALL_DATA.trial_num = SACS_ALL_DATA.trial_num  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_start = TRIALS_DATA.time_start  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_end = TRIALS_DATA.time_end  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_onset = SACS_ALL_DATA.time_onset  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_vmax = SACS_ALL_DATA.time_vmax  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_offset = SACS_ALL_DATA.time_offset  ;
        data_recording(counter_rec).SACS_ALL_DATA.reaction = SACS_ALL_DATA.reaction  ;
        data_recording(counter_rec).SACS_ALL_DATA.duration = SACS_ALL_DATA.duration  ;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_vm = SACS_ALL_DATA.eye_r_vm  ;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_amp_m = SACS_ALL_DATA.eye_r_amp_m  ;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_ang = SACS_ALL_DATA.eye_r_ang  ;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_px_onset = SACS_ALL_DATA.eye_r_px_onset;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_py_onset = SACS_ALL_DATA.eye_r_py_onset;
        data_recording(counter_rec).SACS_ALL_DATA.tgt_px_onset = SACS_ALL_DATA.visual_px_onset  ;
        data_recording(counter_rec).SACS_ALL_DATA.tgt_py_onset = SACS_ALL_DATA.visual_py_onset  ;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_px_offset = SACS_ALL_DATA.eye_r_px_offset;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_py_offset = SACS_ALL_DATA.eye_r_py_offset;
        data_recording(counter_rec).SACS_ALL_DATA.tgt_px_offset = SACS_ALL_DATA.visual_px_offset  ;
        data_recording(counter_rec).SACS_ALL_DATA.tgt_py_offset = SACS_ALL_DATA.visual_py_offset  ;
        data_recording(counter_rec).SACS_ALL_DATA.validity = SACS_ALL_DATA.validity  ;

        EXPERIMENT_DATA.num_sacs_sess_rec(counter_rec, 1) = length(SACS_ALL_DATA.tag);
        EXPERIMENT_DATA.num_trial_sess_rec(counter_rec, 1) = length(TRIALS_DATA.time_start);

        path_to_raw = [path_to_rec_tongue '..' filesep '..' filesep '..' filesep 'raw_data' filesep];

        % add vid frame to EXPERIMENT_DATA
        % Save data
        if ~isempty(dir([path_to_rec_tongue '*_DLC.mp4']))
            vid = dir([path_to_rec_tongue '*_DLC.mp4']);
        elseif ~isempty(dir([path_to_raw '*.mp4']))
            vid = dir([path_to_raw '*.mp4']);
        else
            vid = dir([path_to_raw '*.avi']);
        end

        vid_name = vid.name;
        vid_path = vid.folder;

        EXPERIMENT_DATA.vid_frame = readFrame(VideoReader([vid_path filesep vid_name]));
    end

    % compute break duration -> time between recordings
    for counter_rec = 1 : n_recs - 1
        EXPERIMENT_DATA.break_duration(counter_rec,1) = EXPERIMENT_DATA.rec_time(counter_rec + 1) - ...
            EXPERIMENT_DATA.rec_time(counter_rec) - EXPERIMENT_DATA.rec_duration(counter_rec);
    end
    EXPERIMENT_DATA.break_duration(n_recs,1) = 0;

    % concatenate LICKS_ALL_DATA and SACS_ALL_DATA
    [LICKS_ALL_DATA,SACS_ALL_DATA] = PGH_combine_recs(data_recording);

    % add num_licks and num_sacs to EXPERIMENT_DATA
    EXPERIMENT_DATA.num_licks_sess = length(LICKS_ALL_DATA.tag);
    EXPERIMENT_DATA.num_sacs_sess = length(SACS_ALL_DATA.tag);

    save([out_path  char(EXPERIMENT_DATA.id) '.mat'], 'LICKS_ALL_DATA','SACS_ALL_DATA','EXPERIMENT_DATA','-v7.3')

end
fprintf('### ALL DONE. ###\n')
end
%% function build_population_data
function build_population_data(path_data_monkey_sorted)
%% Load data
path_BEHAVE_session_data = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted filesep 'session_data'];
out_path = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted filesep];

if ~strcmp(path_BEHAVE_session_data(end), filesep);path_BEHAVE_session_data = [path_BEHAVE_session_data filesep];end
session_list = dir([path_BEHAVE_session_data '*mat']);
num_session = numel(session_list);

%% Loop over sessions
for counter_session = 1 : num_session
    fprintf(['### ' 'Building population data session num. ' num2str(counter_session), ' / ' num2str(num_session) ' ### \n']);
    load([session_list(counter_session).folder filesep session_list(counter_session).name]);

    % population_experiment_data
    fieldnames_EXPERIMENT_DATA = fieldnames(EXPERIMENT_DATA);
    for counter_fieldnames = 1 : length(fieldnames_EXPERIMENT_DATA)
        population_experiment_data.(fieldnames_EXPERIMENT_DATA{counter_fieldnames}){counter_session,1} = EXPERIMENT_DATA.(fieldnames_EXPERIMENT_DATA{counter_fieldnames});
    end

    % population_lick_data
    fieldnames_LICK_DATA = fieldnames(LICKS_ALL_DATA);
    for counter_fieldnames = 1 : length(fieldnames_LICK_DATA)
        population_lick_data.(fieldnames_LICK_DATA{counter_fieldnames}){counter_session,1} = LICKS_ALL_DATA.(fieldnames_LICK_DATA{counter_fieldnames});
    end

    % population_sac_data
    fieldnames_SAC_DATA = fieldnames(SACS_ALL_DATA);
    for counter_fieldnames = 1 : length(fieldnames_SAC_DATA)
        population_sac_data.(fieldnames_SAC_DATA{counter_fieldnames}){counter_session,1} = SACS_ALL_DATA.(fieldnames_SAC_DATA{counter_fieldnames});
    end
end
%% Save data
fprintf(['Saving .mat file' ' ...'])
save([out_path 'population_experiment_data' '.mat'], 'population_experiment_data', '-v7.3');
save([out_path 'population_lick_data' '.mat'], 'population_lick_data', '-v7.3');
save([out_path 'population_sac_data' '.mat'], 'population_sac_data', '-v7.3');

fprintf('### ALL DONE. ###\n')

end
%% function generate_population_meta_data
function generate_population_meta_data(path_data_monkey_sorted)
%% Load data
path_BEHAVE_population_data = ['Z:\video_10TB\Paul\BEHAVE\' path_data_monkey_sorted ];
if ~strcmp(path_BEHAVE_population_data(end), filesep);path_BEHAVE_population_data = [path_BEHAVE_population_data filesep];end
out_path = [path_BEHAVE_population_data];
load([path_BEHAVE_population_data 'population_experiment_data.mat']);
%% Build inclusion table
num_session = length(population_experiment_data.sess_date);
include_tongue = ones(num_session,1);
include_eye = ones(num_session,1);

session_date = population_experiment_data.sess_date;
num_trial = population_experiment_data.num_trial_sess;
num_lick = population_experiment_data.num_licks_sess;
num_sac = population_experiment_data.num_sacs_sess;
exp_type = population_experiment_data.experiment_type_sess;

meta_data = table(session_date, exp_type, num_trial, num_lick, num_sac, include_tongue, include_eye);

%% Save
writetable(meta_data, [out_path 'meta_data.xls'])
end
%% function combine_recs
function [LICKS_ALL_DATA,SACS_ALL_DATA] = PGH_combine_recs(data_recording)

num_recording = numel(data_recording);

% LICKS_ALL_DATA
field_names_LICKS_ALL_DATA = fieldnames(data_recording(1).LICKS_ALL_DATA);
for counter_field = 1 : 1 : length(field_names_LICKS_ALL_DATA)
    field_name_LICKS_ALL_DATA = field_names_LICKS_ALL_DATA{counter_field};
    data_cell.LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA) = [];
    for counter_recording = 1 : 1 : num_recording
        data_cell.LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA) = horzcat(...
            data_cell.LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA), ...
            data_recording(counter_recording).LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA));
    end
end
LICKS_ALL_DATA = data_cell.LICKS_ALL_DATA;

% SACS_ALL_DATA
field_names_SACS_ALL_DATA = fieldnames(data_recording(1).SACS_ALL_DATA);
for counter_field = 1 : 1 : length(field_names_SACS_ALL_DATA)
    field_name_SACS_ALL_DATA = field_names_SACS_ALL_DATA{counter_field};
    data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = [];
    for counter_recording = 1 : 1 : num_recording
        data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = horzcat(...
            data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA), ...
            data_recording(counter_recording).SACS_ALL_DATA.(field_name_SACS_ALL_DATA));
    end
end
SACS_ALL_DATA = data_cell.SACS_ALL_DATA;

end

%% function build_variables
function [count,lick,sac] = build_variables(population_experiment_data, population_lick_data,population_sac_data, path)
count.num_session =  length(population_lick_data.tag);
count.num_type = 1; % determine 1 = between harvest, 2 = within harvest trial/sac
count.num_tag_lick = 5; % 5 = consider tags 1 : 5
count.num_tag_sac = 10;
%%% compute vigor fit %%%
[FIT] = compute_vigor(population_lick_data, population_sac_data, count, path, 0);

for counter_session = 1 : count.num_session
    %%% build lick data %%%
    ind_harvest_str = find(population_lick_data.tag_harvest{counter_session, 1} == 1);
    ind_harvest_end = find(population_lick_data.tag_harvest{counter_session, 1} == 2);
    count.num_harvest(counter_session,1) = length(ind_harvest_str);
    sac.duration_work{counter_session,1} = num2cell([population_lick_data.time_onset{counter_session, 1}(ind_harvest_str(2:end)) - population_lick_data.time_onset{counter_session, 1}(ind_harvest_end(1:end-1)) nan]'); ...
        sac.duration_work{counter_session,1}(cell2mat(sac.duration_work{counter_session,1}) < 0,1) = {nan};sac.duration_work{counter_session,1}(cell2mat(sac.duration_work{counter_session,1}) > 60,1) = {nan};
    for counter_harvest = 1 : count.num_harvest(counter_session,1)
        ind_harvest_span = ind_harvest_str(counter_harvest):ind_harvest_end(counter_harvest); % span of inds from harvest str to harvest end
        lick.tag_lick{counter_session,1}{counter_harvest,:} = population_lick_data.tag{counter_session, 1}(ind_harvest_span);
        lick.tag_harvest{counter_session,1}{counter_harvest,:} = population_lick_data.tag_harvest{counter_session, 1}(ind_harvest_span);
        lick.time_onset_lick{counter_session,1}{counter_harvest,:} = population_lick_data.time_onset{counter_session, 1}(ind_harvest_span);
        lick.time_offset_lick{counter_session,1}{counter_harvest,:} = population_lick_data.time_offset{counter_session, 1}(ind_harvest_span);
        lick.tongue_dm_max{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_dm_max{counter_session, 1}(ind_harvest_span);
        lick.tongue_vm_max{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_vm_max{counter_session, 1}(ind_harvest_span);
        lick.tongue_vm_min{counter_session,1}{counter_harvest,:} = abs(population_lick_data.tongue_vm_min{counter_session, 1}(ind_harvest_span));
        lick.tongue_ang_max{counter_session,1}{counter_harvest,:} = abs(population_lick_data.tongue_ang_max{counter_session, 1}(ind_harvest_span));
        lick.ILI{counter_session,1}{counter_harvest,:} = [diff(lick.time_onset_lick{counter_session,1}{counter_harvest,:}) nan]; ...
            lick.ILI{counter_session,1}{counter_harvest,1}(1,lick.ILI{counter_session,1}{counter_harvest,1}(1,:) < 0.1) = nan; lick.ILI{counter_session,1}{counter_harvest,1}(1,lick.ILI{counter_session,1}{counter_harvest,1}(1,:) > 0.6) = nan;
        lick.ILR{counter_session,1}{counter_harvest,:} = 1./lick.ILI{counter_session,1}{counter_harvest,:};
        count.num_lick_harvest{counter_session,1}{counter_harvest,1} = length(lick.time_onset_lick{counter_session,1}{counter_harvest,:});
        lick.num_lick_harvest{counter_session,1}{counter_harvest,1} =  count.num_lick_harvest{counter_session,1}{counter_harvest,1}; % copy for sake of figures
        lick.duration_harvest{counter_session,1}{counter_harvest,1} =  lick.time_offset_lick{counter_session,1}{counter_harvest,1}(1,end) -  lick.time_onset_lick{counter_session,1}{counter_harvest,1}(1,1);
        lick.tongue_dm{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_dm{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_vm{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_vm{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_ang{counter_session,1}{counter_harvest,1} = abs(population_lick_data.tongue_ang{counter_session,1}(:,ind_harvest_span))';
        lick.tongue_duration{counter_session,1}{counter_harvest,1} = lick.time_offset_lick{counter_session,1}{counter_harvest,:} - lick.time_onset_lick{counter_session,1}{counter_harvest,:};
        lick.tongue_tip_px{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_tip_px{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_tip_py{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_tip_py{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_mid_px{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_mid_px{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_mid_py{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_mid_py{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_r_px{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_r_px{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_r_py{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_r_py{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_l_px{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_l_px{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_l_py{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_l_py{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_tip_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_tip_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_tip_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_tip_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_mid_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_mid_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_mid_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_mid_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_r_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_r_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_r_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_r_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_l_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_l_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_l_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_l_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.nose_r_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.nose_r_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.nose_r_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.nose_r_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.nose_l_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.nose_l_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.nose_l_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.nose_l_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.rtube_r_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_r_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.rtube_r_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_r_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.rtube_l_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_l_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.rtube_l_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_l_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.ltube_r_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_r_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.ltube_r_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_r_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.ltube_l_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_l_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.ltube_l_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_l_py_dmax{counter_session, 1}(ind_harvest_span);

        %find if a bout is made to the right or left tube using sum of angles, then nan out the opposite reward data. if equal, completley nan out
        if sum(lick.tongue_ang_max{counter_session,1}{counter_harvest,:} > 0) > sum(lick.tongue_ang_max{counter_session,1}{counter_harvest,:} < 0)
            rew_r{counter_session,1}{counter_harvest,:} = population_lick_data.rew_capacity_r_lick_onset{counter_session, 1}(ind_harvest_span);
            rew_l{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            lick.rew_r_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rew_r_px_dmax{counter_session, 1}(ind_harvest_span);
            lick.rew_r_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rew_r_py_dmax{counter_session, 1}(ind_harvest_span);
            lick.rew_l_px_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            lick.rew_l_py_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            %filter out bad rew value is_nan_r = rew_r{counter_session,1}{counter_harvest,:}>1.5;
            %rew_r{counter_session,1}{counter_harvest,:}(is_nan_r) = nan;
        elseif sum(lick.tongue_ang_max{counter_session,1}{counter_harvest,:} > 0) < sum(lick.tongue_ang_max{counter_session,1}{counter_harvest,:} < 0)
            rew_r{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            rew_l{counter_session,1}{counter_harvest,:} = population_lick_data.rew_capacity_l_lick_onset{counter_session, 1}(ind_harvest_span);
            lick.rew_r_px_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            lick.rew_r_py_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            lick.rew_l_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rew_l_px_dmax{counter_session, 1}(ind_harvest_span);
            lick.rew_l_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rew_l_py_dmax{counter_session, 1}(ind_harvest_span);
            %filter out bad rew value is_nan_l = rew_r{counter_session,1}{counter_harvest,:}>1.5;
            %rew_l{counter_session,1}{counter_harvest,:}(is_nan_l) = nan;
        else
            rew_r{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            rew_l{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            lick.rew_r_px_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            lick.rew_r_py_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            lick.rew_l_px_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
            lick.rew_l_py_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}{counter_harvest,1});
        end
        rew_r_harvest_str{counter_session,1}{counter_harvest,1} = rew_r{counter_session,1}{counter_harvest,:}(1); % as harvest series
        rew_r_harvest_end{counter_session,1}{counter_harvest,1} = rew_r{counter_session,1}{counter_harvest,:}(end); % as harvest series
        rew_l_harvest_str{counter_session,1}{counter_harvest,1} = rew_l{counter_session,1}{counter_harvest,:}(1); % as harvest series
        rew_l_harvest_end{counter_session,1}{counter_harvest,1} = rew_l{counter_session,1}{counter_harvest,:}(end); % as harvest series
        rew_r_harvest_consumed{counter_session,1}{counter_harvest,1} = rew_r_harvest_str{counter_session,1}{counter_harvest,1} - rew_r_harvest_end{counter_session,1}{counter_harvest,1}; % as harvest series
        rew_l_harvest_consumed{counter_session,1}{counter_harvest,1} = rew_l_harvest_str{counter_session,1}{counter_harvest,1} - rew_l_harvest_end{counter_session,1}{counter_harvest,1}; % as harvest series
        rew_r_consumed{counter_session,1}{counter_harvest,:} = [rew_r_harvest_str{counter_session,1}{counter_harvest,1} -  rew_r{counter_session,1}{counter_harvest,:}(1,2:end) nan]; % as lick series
        rew_l_consumed{counter_session,1}{counter_harvest,:} = [rew_l_harvest_str{counter_session,1}{counter_harvest,1} -  rew_l{counter_session,1}{counter_harvest,:}(1,2:end) nan]; % as lick series
        % combine food tube data
        if ~isnan(rew_r_harvest_str{counter_session,1}{counter_harvest,1})
            lick.rew{counter_session,1}{counter_harvest,:} = rew_r{counter_session,1}{counter_harvest,:};
            lick.rew_harvest_str{counter_session,1}{counter_harvest,1} = rew_r_harvest_str{counter_session,1}{counter_harvest,1};
            lick.rew_harvest_end{counter_session,1}{counter_harvest,1} = rew_r_harvest_end{counter_session,1}{counter_harvest,1};
            lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} = rew_r_harvest_consumed{counter_session,1}{counter_harvest,1};
            lick.rew_consumed{counter_session,1}{counter_harvest,:} = rew_r_consumed{counter_session,1}{counter_harvest,:};
        elseif ~isnan(rew_l_harvest_str{counter_session,1}{counter_harvest,1})
            lick.rew{counter_session,1}{counter_harvest,:} = rew_l{counter_session,1}{counter_harvest,:};
            lick.rew_harvest_str{counter_session,1}{counter_harvest,1} = rew_l_harvest_str{counter_session,1}{counter_harvest,1};
            lick.rew_harvest_end{counter_session,1}{counter_harvest,1} = rew_l_harvest_end{counter_session,1}{counter_harvest,1};
            lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} = rew_l_harvest_consumed{counter_session,1}{counter_harvest,1};
            lick.rew_consumed{counter_session,1}{counter_harvest,:} = rew_l_consumed{counter_session,1}{counter_harvest,:};
        else % use r since it will be nan anyway
            lick.rew{counter_session,1}{counter_harvest,:} = rew_r{counter_session,1}{counter_harvest,:};
            lick.rew_harvest_str{counter_session,1}{counter_harvest,1} = rew_r_harvest_str{counter_session,1}{counter_harvest,1};
            lick.rew_harvest_end{counter_session,1}{counter_harvest,1} = rew_r_harvest_end{counter_session,1}{counter_harvest,1};
            lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} = rew_r_harvest_consumed{counter_session,1}{counter_harvest,1};
            lick.rew_consumed{counter_session,1}{counter_harvest,:} = rew_r_consumed{counter_session,1}{counter_harvest,:};
        end

        % load and compute vigor of pro and ret
        load([path.out_path 'VIGOR' filesep path.path_data_monkey_sorted '_FIT.mat' ]);
        for counter_lick = 1 : count.num_lick_harvest{counter_session,1}{counter_harvest,1}
            if (lick.tag_lick{counter_session,1}{counter_harvest,:}(1,counter_lick) == 1)
                fit_pro = FIT.fit_lick_groom_pro ;
                fit_ret = FIT.fit_lick_groom_ret ;
            else
                fit_pro = FIT.fit_lick_rew_pro ;
                fit_ret = FIT.fit_lick_rew_ret ;
            end
            lick.tongue_vigor_pro{counter_session,1}{counter_harvest,:}(1,counter_lick) = lick.tongue_vm_max{counter_session,1}{counter_harvest,:}(1,counter_lick) ...
                /fit_pro(lick.tongue_dm_max{counter_session,1}{counter_harvest,:}(1,counter_lick));
            lick.tongue_vigor_ret{counter_session,1}{counter_harvest,:}(1,counter_lick) = lick.tongue_vm_min{counter_session,1}{counter_harvest,:}(1,counter_lick) ...
                /fit_ret(lick.tongue_dm_max{counter_session,1}{counter_harvest,:}(1,counter_lick));
        end

        lick.tongue_vigor_pro{counter_session,1}{counter_harvest,1}(1, lick.tongue_vigor_pro{counter_session,1}{counter_harvest,:}(1,:) > 1.5 | lick.tongue_vigor_pro{counter_session,1}{counter_harvest,:}(1,:) < 0.5) = nan;
        lick.tongue_vigor_ret{counter_session,1}{counter_harvest,1}(1, lick.tongue_vigor_ret{counter_session,1}{counter_harvest,:}(1,:) > 1.5 | lick.tongue_vigor_ret{counter_session,1}{counter_harvest,:}(1,:) < 0.5) = nan;

        % compute reward gained during work
        if counter_harvest == 1
            sac.rew_gained_work{counter_session,1}{counter_harvest,1} = nan ; % as harvest series
        else
            if ~isnan(lick.rew_harvest_str{counter_session,1}{counter_harvest,1}) && ~isnan(lick.rew_harvest_end{counter_session,1}{counter_harvest-1,1})
                sac.rew_gained_work{counter_session,1}{counter_harvest,1} = lick.rew_harvest_str{counter_session,1}{counter_harvest,1} - lick.rew_harvest_end{counter_session,1}{counter_harvest-1,1};
                if  sac.rew_gained_work{counter_session,1}{counter_harvest,1} < 0 || sac.rew_gained_work{counter_session,1}{counter_harvest,1} > 1.5
                    sac.rew_gained_work{counter_session,1}{counter_harvest,1} = nan;
                end
            else
                sac.rew_gained_work{counter_session,1}{counter_harvest,1} = nan;
            end
        end
        % filter rew
        lick.rew_consumed{counter_session,1}{counter_harvest,1}(1, lick.rew_consumed{counter_session,1}{counter_harvest,:}(1,:) > 1.5 | lick.rew_consumed{counter_session,1}{counter_harvest,:}(1,:) < -1) = nan;
        lick.rew{counter_session,1}{counter_harvest,1}(1, lick.rew{counter_session,1}{counter_harvest,:}(1,:) > 1.5 | lick.rew{counter_session,1}{counter_harvest,:}(1,:) < -1) = nan;

        if (lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} > 1.5 || lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} < - 1)
            lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} = nan;
        end
        if (lick.rew_harvest_str{counter_session,1}{counter_harvest,1} > 1.5 || lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} < - 1 )
            lick.rew_harvest_str{counter_session,1}{counter_harvest,1} = nan;
        end
        if (lick.rew_harvest_end{counter_session,1}{counter_harvest,1} > 1.5 || lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} < - 1 )
            lick.rew_harvest_end{counter_session,1}{counter_harvest,1} = nan;
        end
    end

    %%% build eye data %%%
    % shift eye data to align eyelink data with video data using sample_diff
    sample_diff = population_experiment_data.sample_diff{counter_session, 1};
    num_sac_rec{counter_session, 1} = population_experiment_data.num_sacs_sess_rec{counter_session, 1};
    exp_start_time   = population_experiment_data.exp_start_time{counter_session, 1};
    num_trial_rec{counter_session, 1} = population_experiment_data.num_trial_sess_rec{counter_session, 1};
    num_harvest_rec{counter_session,1} = population_experiment_data.num_harvests_sess_rec{counter_session, 1};
    time_end_trial_{counter_session,1} = population_sac_data.time_end{counter_session,1};
    tag_sac_{counter_session,1} = population_sac_data.tag{counter_session,1};
    time_onset_sac_{counter_session,1} = population_sac_data.time_onset{counter_session,1};
    eye_vm_{counter_session,1} = population_sac_data.eye_r_vm{counter_session,1};
    eye_ang_{counter_session,1} = population_sac_data.eye_r_ang{counter_session,1};
    eye_vm_max_{counter_session,1} = max(population_sac_data.eye_r_vm{counter_session,1});
    eye_dm_max_{counter_session,1} = population_sac_data.eye_r_amp_m{counter_session,1};
    trial_num_{counter_session,1} = population_sac_data.trial_num{counter_session,1};
    reaction_{counter_session,1} = population_sac_data.reaction{counter_session,1};
    eye_px_onset_{counter_session,1} = population_sac_data.eye_r_px_onset{counter_session,1};
    eye_py_onset_{counter_session,1} = population_sac_data.eye_r_py_onset{counter_session,1};
    tgt_px_onset_{counter_session,1} = population_sac_data.tgt_px_onset{counter_session,1};
    tgt_py_onset_{counter_session,1} = population_sac_data.tgt_py_onset{counter_session,1};
    eye_px_offset_{counter_session,1} = population_sac_data.eye_r_px_offset{counter_session,1};
    eye_py_offset_{counter_session,1} = population_sac_data.eye_r_py_offset{counter_session,1};
    tgt_px_offset_{counter_session,1} = population_sac_data.tgt_px_offset{counter_session,1};
    tgt_py_offset_{counter_session,1} = population_sac_data.tgt_py_offset{counter_session,1};
    validity_{counter_session,1} = population_sac_data.validity{counter_session,1};

    shift_sac = 1;
    shift_trial = 1;

    for counter_rec = 1 :  length(num_trial_rec{counter_session})
        time_end_trial__{counter_session,1}{1,counter_rec} = (time_end_trial_{counter_session,1}(1,shift_trial : num_trial_rec{counter_session, 1}(counter_rec) + shift_trial - 1 ) ...
            - exp_start_time(counter_rec) + sample_diff(counter_rec)/1000);
        time_onset_sac__{counter_session,1}{1,counter_rec} = (time_onset_sac_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ) ...
            - exp_start_time(counter_rec) + sample_diff(counter_rec)/1000);
        tag_sac__{counter_session,1}{1,counter_rec} = (tag_sac_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_vm__{counter_session,1}{1,counter_rec} = (eye_vm_{counter_session,1}(:,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_ang__{counter_session,1}{1,counter_rec} = (eye_ang_{counter_session,1}(:,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_vm_max__{counter_session,1}{1,counter_rec} = (eye_vm_max_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_dm_max__{counter_session,1}{1,counter_rec} = (eye_dm_max_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        trial_num__{counter_session,1}{1,counter_rec} = (trial_num_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        reaction__{counter_session,1}{1,counter_rec} = (reaction_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_px_onset__{counter_session,1}{1,counter_rec} = (eye_px_onset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_py_onset__{counter_session,1}{1,counter_rec} = (eye_py_onset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        tgt_px_onset__{counter_session,1}{1,counter_rec} = (tgt_px_onset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        tgt_py_onset__{counter_session,1}{1,counter_rec} = (tgt_py_onset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_px_offset__{counter_session,1}{1,counter_rec} = (eye_px_offset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_py_offset__{counter_session,1}{1,counter_rec} = (eye_py_offset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        tgt_px_offset__{counter_session,1}{1,counter_rec} = (tgt_px_offset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        tgt_py_offset__{counter_session,1}{1,counter_rec} = (tgt_py_offset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        validity__{counter_session,1}{1,counter_rec} = (validity_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));

        shift_sac = shift_sac + num_sac_rec{counter_session, 1}(counter_rec);
        shift_trial = shift_trial + num_trial_rec{counter_session, 1}(counter_rec);

        for counter_harvest_rec = 1 : num_harvest_rec{counter_session,1}(counter_rec)
            shift = 0;
            if counter_rec > 1
                counter_harvest_rec = counter_harvest_rec + sum(num_harvest_rec{counter_session,1}(1 : counter_rec - 1));
                shift = sum(num_harvest_rec{counter_session,1}(1 : counter_rec - 1));
            end

            % detect trials and saccs occuring either pre or within bout
            if counter_harvest_rec == 1
                is_pre_bout_trial = time_end_trial__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1);
                is_pre_bout_sac = time_onset_sac__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1);
            elseif counter_harvest_rec > 1
                is_pre_bout_trial = time_end_trial__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                    & time_end_trial__{counter_session,1}{counter_rec} > lick.time_onset_lick{counter_session,1}{counter_harvest_rec-1,:}(1,end);
                is_pre_bout_sac = time_onset_sac__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                    & time_onset_sac__{counter_session,1}{counter_rec} > lick.time_onset_lick{counter_session,1}{counter_harvest_rec-1,:}(1,end);
            end
            is_in_bout_trial = time_end_trial__{counter_session,1}{counter_rec} > lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                & time_end_trial__{counter_session,1}{counter_rec} < lick.time_offset_lick{counter_session,1}{counter_harvest_rec,:}(1,end);
            is_in_bout_sac = time_onset_sac__{counter_session,1}{counter_rec} > lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                & time_onset_sac__{counter_session,1}{counter_rec} < lick.time_offset_lick{counter_session,1}{counter_harvest_rec,:}(1,end);

            sac.time_end_trial{counter_session,1}{counter_harvest_rec,:} = time_end_trial__{counter_session,1}{counter_rec}(1,is_pre_bout_trial);
            count.num_trial_work{counter_session,1}{counter_harvest_rec,1} = length(sac.time_end_trial{counter_session,1}{counter_harvest_rec,:});
            sac.num_trial_work{counter_session,1}{counter_harvest_rec,1} = count.num_trial_work{counter_session,1}{counter_harvest_rec,1};

            %        sac.time_end_trial{counter_session,2}{counter_harvest_rec ,:} = time_end_trial__{counter_session,1}{counter_rec}(1,is_in_bout_trial);
            %          count.num_trial_work{counter_session,2}{counter_harvest_rec,1} = length(sac.time_end_trial{counter_session,2}{counter_harvest_rec,:});

            % handle situations where no trials occur
            if isempty(sac.time_end_trial{counter_session,1}{counter_harvest_rec,:})
                sac.ITI{counter_session,1}{counter_harvest_rec,:} = [];
                sac.ITR{counter_session,1}{counter_harvest_rec,:} = [];
            else
                sac.ITI{counter_session,1}{counter_harvest_rec,:} = [diff(sac.time_end_trial{counter_session,1}{counter_harvest_rec,:}) nan];
                %                     sac.ITI{counter_session,1}{counter_harvest_rec,:}(sac.ITI{counter_session,1}{counter_harvest_rec,:}>20) = nan; % filter out high values
                sac.ITR{counter_session,1}{counter_harvest_rec,:} = 1./sac.ITI{counter_session,1}{counter_harvest_rec,:};
            end
            %             if isempty(sac.time_end_trial{counter_session,2}{counter_harvest_rec,:})
            %                 sac.ITI{counter_session,2}{counter_harvest_rec,:} = [];
            %                 sac.ITR{counter_session,2}{counter_harvest_rec,:} = [];
            %             else
            %                 sac.ITI{counter_session,2}{counter_harvest_rec,:} = [diff(sac.time_end_trial{counter_session,2}{counter_harvest_rec,:}) nan];
            %                 %                     sac.ITI{counter_session,2}{counter_harvest_rec,:}(sac.ITI{counter_session,2}{counter_harvest_rec,:}>20) = nan; % filter out high values
            %                 sac.ITR{counter_session,2}{counter_harvest_rec,:} = 1./sac.ITI{counter_session,2}{counter_harvest_rec,:};
            %             end

            % specify saccade type to analyze
            tag = [1];
            is_tag = tag_sac__{counter_session,1}{counter_rec} == tag;
            is_valid = validity__{counter_session,1}{counter_rec} == 1;

            sac.tag_sac{counter_session,1}{counter_harvest_rec,:} = tag_sac__{counter_session,1}{counter_rec}(1,is_pre_bout_sac & is_tag & is_valid);
            %             sac.tag_sac{counter_session,2}{counter_harvest_rec ,:} = tag_sac__{counter_session,1}{counter_rec}(1,is_in_bout_sac & is_tag & is_valid);

            sac.trial_num{counter_session,1}{counter_harvest_rec,:} = trial_num__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.trial_num{counter_session,2}{counter_harvest_rec ,:} = trial_num__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.time_onset_sac{counter_session,1}{counter_harvest_rec,:} = time_onset_sac__{counter_session,1}{counter_rec}(1,is_pre_bout_sac & is_tag & is_valid);
            %             sac.time_onset_sac{counter_session,2}{counter_harvest_rec ,:} = time_onset_sac__{counter_session,1}{counter_rec}(1,is_in_bout_sac & is_tag & is_valid);

            count.num_sac_work{counter_session,1}{counter_harvest_rec,1} = length(sac.time_onset_sac{counter_session,1}{counter_harvest_rec,:});
            %             count.num_sac_work{counter_session,2}{counter_harvest_rec,1} = length(sac.time_onset_sac{counter_session,2}{counter_harvest_rec,:});
            sac.num_sac_work{counter_session,1}{counter_harvest_rec,1} = count.num_sac_work{counter_session,1}{counter_harvest_rec,1};
            %             sac.num_sac_work{counter_session,2}{counter_harvest_rec,1} = count.num_sac_work{counter_session,2}{counter_harvest_rec,1};

            sac.eye_vm_max{counter_session,1}{counter_harvest_rec,:} = eye_vm_max__{counter_session,1}{counter_rec}(1,is_pre_bout_sac & is_tag & is_valid);
            %             sac.eye_vm_max{counter_session,2}{counter_harvest_rec ,:} = eye_vm_max__{counter_session,1}{counter_rec}(1,is_in_bout_sac & is_tag & is_valid);

            sac.eye_ang{counter_session,1}{counter_harvest_rec,:} = eye_ang__{counter_session,1}{counter_rec}(1,is_pre_bout_sac & is_tag & is_valid);
            %             sac.eye_ang{counter_session,2}{counter_harvest_rec ,:} = eye_ang__{counter_session,1}{counter_rec}(1,is_in_bout_sac & is_tag & is_valid);

            sac.eye_dm_max{counter_session,1}{counter_harvest_rec,:} = eye_dm_max__{counter_session,1}{counter_rec}(1,is_pre_bout_sac & is_tag & is_valid);
            %             sac.eye_dm_max{counter_session,2}{counter_harvest_rec ,:} = eye_dm_max__{counter_session,1}{counter_rec}(1,is_in_bout_sac & is_tag & is_valid);

            sac.eye_vm{counter_session,1}{counter_harvest_rec,:} = eye_vm__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.eye_vm{counter_session,2}{counter_harvest_rec ,:} = eye_vm__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.reaction{counter_session,1}{counter_harvest_rec,:} = reaction__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.reaction{counter_session,2}{counter_harvest_rec ,:} = reaction__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.eye_px_onset{counter_session,1}{counter_harvest_rec,:} = eye_px_onset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.eye_px_onset{counter_session,2}{counter_harvest_rec ,:} = eye_px_onset__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.eye_py_onset{counter_session,1}{counter_harvest_rec,:} = eye_py_onset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.eye_py_onset{counter_session,2}{counter_harvest_rec ,:} = eye_py_onset__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.tgt_px_onset{counter_session,1}{counter_harvest_rec,:} = tgt_px_onset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.tgt_px_onset{counter_session,2}{counter_harvest_rec ,:} = tgt_px_onset__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.tgt_py_onset{counter_session,1}{counter_harvest_rec,:} = tgt_py_onset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.tgt_py_onset{counter_session,2}{counter_harvest_rec ,:} = tgt_py_onset__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.eye_px_offset{counter_session,1}{counter_harvest_rec,:} = eye_px_offset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.eye_px_offset{counter_session,2}{counter_harvest_rec ,:} = eye_px_offset__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.eye_py_offset{counter_session,1}{counter_harvest_rec,:} = eye_py_offset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.eye_py_offset{counter_session,2}{counter_harvest_rec ,:} = eye_py_offset__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.tgt_px_offset{counter_session,1}{counter_harvest_rec,:} = tgt_px_offset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.tgt_px_offset{counter_session,2}{counter_harvest_rec ,:} = tgt_px_offset__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.tgt_py_offset{counter_session,1}{counter_harvest_rec,:} = tgt_py_offset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag & is_valid);
            %             sac.tgt_py_offset{counter_session,2}{counter_harvest_rec ,:} = tgt_py_offset__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag & is_valid);

            sac.eye_vigor{counter_session,1}{counter_harvest_rec,:} = sac.eye_vm_max{counter_session,1}{counter_harvest_rec,:} ./ ...
                FIT.fit_sac_rel(sac.eye_dm_max{counter_session,1}{counter_harvest_rec,:})';
            %             sac.eye_vigor{counter_session,2}{counter_harvest_rec,:} = sac.eye_vm_max{counter_session,2}{counter_harvest_rec,:} ./ ...
            %                 FIT.fit_sac_rel(sac.eye_dm_max{counter_session,2}{counter_harvest_rec,:})';

            %             if isempty(sac.time_end_trial{counter_session,1}{counter_harvest_rec,:} | isempty(sac.time_onset_sac{counter_session,1}{counter_harvest_rec,:}))
            %                 sac.duration_work{counter_session,1}{counter_harvest_rec,:} = {nan};
            %             else
            %                 sac.duration_work{counter_session,1} = sac.time_end_trial{counter_session,1}{counter_harvest_rec,:}(end) - sac.time_onset_sac{counter_session,1}{counter_harvest_rec,:}(1);
            %             end

            for counter_type = 1 : count.num_type
                % handle situations where no saccades occur
                if isempty(sac.time_onset_sac{counter_session,counter_type}{counter_harvest_rec,:})
                    sac.ISI{counter_session,counter_type}{counter_harvest_rec,:} = [];
                    sac.ISR{counter_session,counter_type}{counter_harvest_rec,:} = [];
                else
                    sac.ISI{counter_session,counter_type}{counter_harvest_rec,:} = [diff(sac.time_onset_sac{counter_session,counter_type}{counter_harvest_rec,:}) nan];
                    sac.ISR{counter_session,counter_type}{counter_harvest_rec,:} = 1./sac.ISI{counter_session,counter_type}{counter_harvest_rec,:};
                end
            end

        end
    end
    % can be used for padding
    %     max_num_lick_harvest_sess(counter_session,1) = max(cell2mat(count.num_lick_harvest{counter_session,1}));
    %     max_num_trial_work_sess(counter_session,1) = max(cell2mat(count.num_trial_work{counter_session,1}));
    %     max_num_trial_work_sess(counter_session,2) = max(cell2mat(count.num_trial_work{counter_session,2}));
    %     max_num_sac_work_sess(counter_session,1) = max(cell2mat(count.num_sac_work{counter_session,1}));
    %     max_num_sac_work_sess(counter_session,2) = max(cell2mat(count.num_sac_work{counter_session,2}));
end

%%% set max number for later padding %%%
count.max_num_harvest = 250;
count.max_num_harvest_early = 250/2;
count.max_num_harvest_late = 250/2;

count.max_num_lick_harvest= 20;

count.max_num_trial_work(1,1) = 7;
count.max_num_sac_work(1,1) = 7;

if count.num_type == 2
    count.max_num_trial_work(1,2) = 3;
    count.max_num_sac_work(1,2) = 3;
end
end
%% function endpoint_error_mag
function [data] = endpoint_error_mag(data, count)
fieldname = fieldnames(data);
for counter_session = 1 : count.num_session
    for counter_type = 1 : count.num_type
        if counter_type == 1
            for counter_fieldname = 1 : length(fieldname)
                % check if lick or sac related data
                if sum(contains(string(fieldname),'_lick')) > 0
                    for counter_harvest = 1 : count.num_harvest(counter_session)
                        for counter_lick = 1 : count.num_lick_harvest{counter_session,1}{counter_harvest,1}
                            % compute x and y of error vector after rotation
                            % check which tongue point is closest

                            % tip
                            err_tongue_tip_rew_r_px_ = (data.rew_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_tip_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_tip_rew_r_py_ = (data.rew_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_tip_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_tip_rew_l_px_ = ( data.rew_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_tip_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_tip_rew_l_py_ = (data.rew_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_tip_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            dist_err_tongue_tip_rew_ = sqrt(nansum([err_tongue_tip_rew_r_px_; err_tongue_tip_rew_l_px_]).^2 ...
                                +  nansum([err_tongue_tip_rew_r_py_; err_tongue_tip_rew_l_py_]).^2);
                            % mid
                            err_tongue_mid_rew_r_px_ = (data.rew_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_mid_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_mid_rew_r_py_ = (data.rew_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_mid_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_mid_rew_l_px_ = ( data.rew_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_mid_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_mid_rew_l_py_ = (data.rew_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_mid_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            dist_err_tongue_mid_rew_ = sqrt(nansum([err_tongue_mid_rew_r_px_; err_tongue_mid_rew_l_px_]).^2 ...
                                +  nansum([err_tongue_mid_rew_r_py_; err_tongue_mid_rew_l_py_]).^2);
                            % r
                            err_tongue_r_rew_r_px_ = (data.rew_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_r_rew_r_py_ = (data.rew_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_r_rew_l_px_ = ( data.rew_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_r_rew_l_py_ = (data.rew_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            dist_err_tongue_r_rew_ = sqrt(nansum([err_tongue_r_rew_r_px_; err_tongue_r_rew_l_px_]).^2 ...
                                +  nansum([err_tongue_r_rew_r_py_; err_tongue_r_rew_l_py_]).^2);
                            % l
                            err_tongue_l_rew_r_px_ = (data.rew_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_l_rew_r_py_ = (data.rew_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_l_rew_l_px_ = ( data.rew_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_l_rew_l_py_ = (data.rew_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            dist_err_tongue_l_rew_ = sqrt(nansum([err_tongue_l_rew_r_px_; err_tongue_l_rew_l_px_]).^2 ...
                                +  nansum([err_tongue_l_rew_r_py_; err_tongue_l_rew_l_py_]).^2);

                            dist_err_tongue_rew_ = [dist_err_tongue_tip_rew_; dist_err_tongue_mid_rew_; dist_err_tongue_r_rew_; dist_err_tongue_l_rew_];

                            if (min(dist_err_tongue_rew_) == dist_err_tongue_tip_rew_)
                                err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_r_px_;
                                err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_r_py_;
                                err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_l_px_;
                                err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_l_py_;
                            elseif (min(dist_err_tongue_rew_) == dist_err_tongue_mid_rew_)
                                err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_mid_rew_r_px_;
                                err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_mid_rew_r_py_;
                                err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_mid_rew_l_px_;
                                err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_mid_rew_l_py_;
                            elseif (min(dist_err_tongue_rew_) == dist_err_tongue_r_rew_)
                                err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_r_rew_r_px_;
                                err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_r_rew_r_py_;
                                err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_r_rew_l_px_;
                                err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_r_rew_l_py_;
                            elseif (min(dist_err_tongue_rew_) == dist_err_tongue_l_rew_)
                                err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_l_rew_r_px_;
                                err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_l_rew_r_py_;
                                err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_l_rew_l_px_;
                                err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_l_rew_l_py_;
                            else
                                err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_r_px_;
                                err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_r_py_;
                                err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_l_px_;
                                err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_l_py_;
                            end

                            theta{counter_session,1}{counter_harvest,1}(1,counter_lick) = data.tongue_ang_max{counter_session,1}{counter_harvest,1}(1,counter_lick);

                            R = [cosd(theta{counter_session,1}{counter_harvest,1}(1,counter_lick)) -sind(theta{counter_session,1}{counter_harvest,1}(1,counter_lick)) ; ...
                                sind(theta{counter_session,1}{counter_harvest,1}(1,counter_lick)) cosd(theta{counter_session,1}{counter_harvest,1}(1,counter_lick))];

                            err_rot = R * [err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                ; err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick)];

                            data.err_tongue_rew_px_rot{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_rot(1);
                            data.err_tongue_rew_py_rot{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_rot(2);

                            data.dist_err_tongue_rew{counter_session,1}{counter_harvest,1}(1,counter_lick) = sqrt(data.err_tongue_rew_px_rot{counter_session,1}{counter_harvest,1}(1,counter_lick).^2 +  data.err_tongue_rew_py_rot{counter_session,1}{counter_harvest,1}(1,counter_lick).^2);
                            %  data. det_var_err_tongue_rew_{counter_session,1}{counter_harvest,1}(1,counter_lick) = det(nancov(data.err_tongue_rew_px_rot_{counter_session,1}{counter_harvest,1}(1,counter_lick), data.err_tongue_rew_py_rot_{counter_session,1}{counter_harvest,1}(1,counter_lick)));
                        end
                    end
                elseif sum(contains(string(fieldname),'_sac')) > 0
                    for counter_harvest = 1 : count.num_harvest(counter_session)
                        % initialize to handle situations where harvest has
                        % no saccades
                        %                         data.err_eye_tgt_px_rot{counter_session,1}(counter_harvest,1 = nan(size(count.num_sac_work{counter_harvest,1},1),1);
                        %                         data.err_eye_tgt_py_rot{counter_session,1}{counter_harvest,1} =nan(1,count.num_sac_work(counter_harvest));
                        %                         data.dist_err_eye_tgt{counter_session,1}{counter_harvest,1} = nan(1,count.num_sac_work(counter_harvest));
                        for counter_sac = 1 : count.num_sac_work{counter_session,1}{counter_harvest,1}
                            % compute x and y of error vector after rotation
                            err_eye_tgt_px{counter_session,1}{counter_harvest,1}(1,counter_sac) = (data.tgt_px_offset{counter_session,1}{counter_harvest,1}(1,counter_sac) ...
                                - data.eye_px_offset{counter_session,1}{counter_harvest,1}(1,counter_sac)) ;
                            err_eye_tgt_py{counter_session,1}{counter_harvest,1}(1,counter_sac) = (data.tgt_py_offset{counter_session,1}{counter_harvest,1}(1,counter_sac) ...
                                - data.eye_py_offset{counter_session,1}{counter_harvest,1}(1,counter_sac)) ;

                            theta{counter_session,1}{counter_harvest,1}(1,counter_sac) = data.eye_ang{counter_session,1}{counter_harvest,1}(1,counter_sac);

                            R = [cosd(theta{counter_session,1}{counter_harvest,1}(1,counter_sac)) -sind(theta{counter_session,1}{counter_harvest,1}(1,counter_sac)) ; ...
                                sind(theta{counter_session,1}{counter_harvest,1}(1,counter_sac)) cosd(theta{counter_session,1}{counter_harvest,1}(1,counter_sac))];

                            err_rot = R * [err_eye_tgt_px{counter_session,1}{counter_harvest,1}(1,counter_sac); err_eye_tgt_py{counter_session,1}{counter_harvest,1}(1,counter_sac)];

                            data.err_eye_tgt_px_rot{counter_session,1}{counter_harvest,1}(1,counter_sac) = err_rot(1);
                            data.err_eye_tgt_py_rot{counter_session,1}{counter_harvest,1}(1,counter_sac) = err_rot(2);
                            data.dist_err_eye_tgt{counter_session,1}{counter_harvest,1}(1,counter_sac) = sqrt(data.err_eye_tgt_px_rot{counter_session,1}{counter_harvest,1}(1,counter_sac)^2 + data.err_eye_tgt_py_rot{counter_session,1}{counter_harvest,1}(1,counter_sac)^2);
                            %  data.det_var_err_eye_tgt_{counter_session,1}{counter_harvest,1}(1,counter_sac) = det(nancov(data.err_eye_tgt_px_rot_{counter_session,1}{counter_harvest,1}(1,counter_sac), data.err_eye_tgt_py_rot_{counter_session,1}{counter_harvest,1}(1,counter_sac)));
                        end
                        if isempty(counter_sac)
                            data.err_eye_tgt_px_rot{counter_session,1}{counter_harvest,1} = [];
                            data.err_eye_tgt_py_rot{counter_session,1}{counter_harvest,1} = [];
                            data.dist_err_eye_tgt{counter_session,1}{counter_harvest,1} = [];
                        end
                    end
                end
            end
        end
    end
end
end

%% function endpoint_error_var
function [data] = endpoint_error_var(data, count)
fieldname = fieldnames(data);
fieldname(~contains(fieldname, 'err') | contains(fieldname, 'dist')) = [];

fieldname_x = fieldname(contains(fieldname,'px'));
fieldname_y = fieldname(contains(fieldname,'py'));

if sum(contains(string(fieldname),'tongue')) > 0
    count_event = count.max_num_lick_harvest;
else
    count_event = count.max_num_sac_work;
end

for counter_session = 1 : count.num_session
    for counter_type = 1 : count.num_type
        if counter_type == 1
            for counter_fieldname = 1 : length(fieldname_x)
                fieldname_output = strcat(fieldname_x{counter_fieldname}(1:strfind(fieldname_x{1}, '_p')-1),fieldname_x{counter_fieldname}(strfind(fieldname_x{1}, 'rot')+3:end));
                for counter_event = 1 : count_event
                    data.(strcat('det_var_',fieldname_output)){counter_session,1}(:,counter_event) = det(nancov(data.(string(fieldname_x(counter_fieldname))){counter_session,1}(:,counter_event), data.(string(fieldname_y(counter_fieldname))){counter_session,1}(:,counter_event)));
                end
            end
        end
    end
end
end


%% function split_early_late
function [data,count] = split_early_late(data, count)
fieldname = fieldnames(data);
for counter_session = 1 : count.num_session
    count.num_harvest_early(counter_session,1) = floor(count.num_harvest(counter_session,1)/2);
    count.num_harvest_late(counter_session,1) = ceil(count.num_harvest(counter_session,1)/2);

    flag_same = 0;

    % detect when same due to even number
    if (count.num_harvest_early(counter_session,1) == count.num_harvest_late(counter_session,1))
        flag_same = 1;
    end

    for counter_fieldname = 1 : length(fieldname)
        %%% split early data %%%
        for counter_early = 1 : count.num_harvest_early(counter_session,1)
            data.(strcat(string(fieldname(counter_fieldname)), '_early')){counter_session,1}{counter_early,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_early,:};
        end

        %%% split late data %%%
        for counter_late = count.num_harvest_late(counter_session,1) : count.num_harvest(counter_session,1)
            data.(strcat(string(fieldname(counter_fieldname)), '_late')){counter_session,1}{counter_late - count.num_harvest_late(counter_session,1) + 1,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_late,:};
        end

        % save number of licks for early/late
        count.num_lick_harvest_early{counter_session,1}(1:count.num_harvest_early(counter_session,1),1) = count.num_lick_harvest{counter_session,1}(1:count.num_harvest_early(counter_session,1),1);
        count.num_trial_work_early{counter_session,1}(1:count.num_harvest_early(counter_session,1),1) = count.num_trial_work{counter_session,1}(1:count.num_harvest_early(counter_session,1),1);
        count.num_sac_work_early{counter_session,1}(1:count.num_harvest_early(counter_session,1),1) = count.num_sac_work{counter_session,1}(1:count.num_harvest_early(counter_session,1),1);
        if flag_same == 1
            % cut of redundant bout due to even number
            data.(strcat(string(fieldname(counter_fieldname)), '_late')){counter_session,1}(1,:) = [];
            count.num_lick_harvest_late{counter_session,1}(1:count.num_harvest_late(counter_session,1),1) = count.num_lick_harvest{counter_session,1}(count.num_harvest_late(counter_session,1) + 1 : count.num_harvest(counter_session,1),1);
            count.num_trial_work_late{counter_session,1}(1:count.num_harvest_late(counter_session,1),1) = count.num_trial_work{counter_session,1}(count.num_harvest_late(counter_session,1) + 1 : count.num_harvest(counter_session,1),1);
            count.num_sac_work_late{counter_session,1}(1:count.num_harvest_late(counter_session,1),1) = count.num_sac_work{counter_session,1}(count.num_harvest_late(counter_session,1)+ 1 : count.num_harvest(counter_session,1),1);

        else
            count.num_lick_harvest_late{counter_session,1}(1:count.num_harvest_late(counter_session,1),1) = count.num_lick_harvest{counter_session,1}(count.num_harvest_late(counter_session,1) : count.num_harvest(counter_session,1),1);
            count.num_trial_work_late{counter_session,1}(1:count.num_harvest_late(counter_session,1),1) = count.num_trial_work{counter_session,1}(count.num_harvest_late(counter_session,1) : count.num_harvest(counter_session,1),1);
            count.num_sac_work_late{counter_session,1}(1:count.num_harvest_late(counter_session,1),1) = count.num_sac_work{counter_session,1}(count.num_harvest_late(counter_session,1) : count.num_harvest(counter_session,1),1);

        end
    end
end
end
%% function pad_data
function [data] = pad_data(data, count)
fieldname = fieldnames(data);
for counter_session = 1 : count.num_session
    for counter_type = 1 : count.num_type
        if counter_type == 1
            for counter_fieldname = 1 : length(fieldname)
                % check if lick or sac related data
                if sum(contains(string(fieldname),'_lick')) > 0
                    % check if combined or early/late data
                    if(~contains(string(fieldname(counter_fieldname)),'_early') && ~contains(string(fieldname(counter_fieldname)),'_late'))
                        % check if bout related data or lick related data
                        if ~contains(string(fieldname(counter_fieldname)),[ "_harvest"]) && size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},2) < 500
                            for counter_harvest = 1 : count.num_harvest(counter_session,1)
                                pad_val_within = count.max_num_lick_harvest - count.num_lick_harvest{counter_session,1}{counter_harvest,1};
                                if pad_val_within >= 0
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1},[0 pad_val_within], nan, 'post');
                                else
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1}(:,1:1:count.max_num_lick_harvest);
                                end
                            end
                            data.(string(fieldname(counter_fieldname))){counter_session,1} =cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1});
                        elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},2) == 500
                            for counter_harvest = 1 : count.num_harvest(counter_session,1)
                                pad_val_within = count.max_num_lick_harvest - count.num_lick_harvest{counter_session,1}{counter_harvest,1};
                                if pad_val_within >= 0
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1},[pad_val_within 0], nan, 'post');
                                else
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1}(1:1:count.max_num_lick_harvest,:);
                                end
                            end
                        elseif contains(string(fieldname(counter_fieldname)),["_harvest"]) && ~contains(string(fieldname(counter_fieldname)),["tag"])
                            pad_val_accross = count.max_num_harvest - count.num_harvest(counter_session,1);
                            if pad_val_accross >= 0
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = padarray(cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}),[pad_val_accross 0], nan, 'post')';
                            else
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}(1:1:count.max_num_harvest,1))';
                            end
                        end
                    elseif(contains(string(fieldname(counter_fieldname)),'_early'))
                        if ~contains(string(fieldname(counter_fieldname)),["_harvest"]) && size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},2) < 500
                            for counter_harvest_early = 1 : count.num_harvest_early(counter_session,1)
                                pad_val_within = count.max_num_lick_harvest - count.num_lick_harvest_early{counter_session,1}{counter_harvest_early,1};
                                if pad_val_within >= 0
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,1},[0 pad_val_within], nan, 'post');
                                else
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,1}(:,1:1:count.max_num_lick_harvest);
                                end
                            end
                            data.(string(fieldname(counter_fieldname))){counter_session,1} =horzcat(cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}));
                        elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},2) == 500
                            for counter_harvest_early = 1 : count.num_harvest_early(counter_session,1)
                                pad_val_within = count.max_num_lick_harvest - count.num_lick_harvest_early{counter_session,1}{counter_harvest_early,1};
                                if pad_val_within >= 0
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,1},[pad_val_within 0], nan, 'post');
                                else
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,1}(1:1:count.max_num_lick_harvest,:);
                                end
                            end
                        elseif contains(string(fieldname(counter_fieldname)),["_harvest"])  && ~contains(string(fieldname(counter_fieldname)),["tag"])
                            pad_val_accross = count.max_num_harvest_early - count.num_harvest_early(counter_session,1);
                            if pad_val_accross >= 0
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = padarray(cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}),[pad_val_accross 0], nan, 'post')';
                            else
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}(1:1:count.max_num_harvest_early,1))';
                            end
                        end
                    elseif(contains(string(fieldname(counter_fieldname)),'_late'))
                        if ~contains(string(fieldname(counter_fieldname)),["_harvest"]) && size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},2) < 500
                            for counter_harvest_late = 1 : count.num_harvest_late(counter_session,1)
                                pad_val_within = count.max_num_lick_harvest - count.num_lick_harvest_late{counter_session,1}{counter_harvest_late,1};
                                if pad_val_within >= 0
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,1},[0 pad_val_within], nan, 'post');
                                else
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,1}(:,1:1:count.max_num_lick_harvest);
                                end
                            end
                            data.(string(fieldname(counter_fieldname))){counter_session,1} =horzcat(cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}));
                        elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},2) == 500
                            for counter_harvest_late = 1 : count.num_harvest_late(counter_session,1)
                                pad_val_within = count.max_num_lick_harvest - count.num_lick_harvest_late{counter_session,1}{counter_harvest_late,1};
                                if pad_val_within >= 0
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,1},[pad_val_within 0], nan, 'post');
                                else
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,1}(1:1:count.max_num_lick_harvest,:);
                                end
                            end
                        elseif contains(string(fieldname(counter_fieldname)),["_harvest"])  && ~contains(string(fieldname(counter_fieldname)),["tag"])
                            pad_val_accross = count.max_num_harvest_late - count.num_harvest_late(counter_session,1);
                            if pad_val_accross >= 0
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = padarray(cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}),[pad_val_accross 0], nan, 'post')';
                            else
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}(1:1:count.max_num_harvest_late,1))';
                            end
                        end
                    end
                elseif sum(contains(string(fieldname),'_sac')) > 0
                    % check if combined or early/late data
                    if(~contains(string(fieldname(counter_fieldname)),'_early') && ~contains(string(fieldname(counter_fieldname)),'_late'))
                        % check if bout related data or lick related data
                        if ~contains(string(fieldname(counter_fieldname)),["duration_work", "num_sac_work", "rew_gained_work", "num_trial_work"]) && size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},1) < 150
                            for counter_harvest = 1 : count.num_harvest(counter_session,1)
                                if(contains(string(fieldname(counter_fieldname)),["time_end_trial", "ITI","ITR"]))
                                    pad_val_within = count.max_num_trial_work - count.num_trial_work{counter_session,1}{counter_harvest,1};
                                else
                                    pad_val_within = count.max_num_sac_work - count.num_sac_work{counter_session,1}{counter_harvest,1};
                                end
                                if isempty(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1})
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = nan(1, pad_val_within);
                                else
                                    if pad_val_within >= 0
                                        data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1},[0 pad_val_within], nan, 'post');
                                    else
                                        data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1}(:,1:1:count.max_num_sac_work);
                                    end
                                end
                            end
                            data.(string(fieldname(counter_fieldname))){counter_session,1} =cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1});
                        elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},1) == 150
                            for counter_harvest = 1 : count.num_harvest(counter_session,1)
                                pad_val_within = count.max_num_sac_work - count.num_sac_work{counter_session,1}{counter_harvest,1};
                                if pad_val_within >= 0
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1}',[pad_val_within 0], nan, 'post');
                                else
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1}(:,1:1:count.max_num_sac_work)';
                                end
                            end
                        elseif contains(string(fieldname(counter_fieldname)),["duration_work", "num_sac_work", "rew_gained_work", "num_trial_work"])
                            pad_val_accross = count.max_num_harvest - count.num_harvest(counter_session,1);
                            if pad_val_accross >= 0
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = padarray(cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}),[pad_val_accross 0], nan, 'post')';
                            else
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}(1:1:count.max_num_harvest,1))';
                            end
                        end
                    elseif(contains(string(fieldname(counter_fieldname)),'_early'))
                        % check if bout related data or lick related data
                        if ~contains(string(fieldname(counter_fieldname)),["duration_work", "num_sac_work", "rew_gained_work", "num_trial_work"]) && size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},1) < 150
                            for counter_harvest_early = 1 : count.num_harvest_early(counter_session,1)
                                if(contains(string(fieldname(counter_fieldname)),["time_end_trial", "ITI","ITR"]))
                                    pad_val_within = count.max_num_trial_work - count.num_trial_work_early{counter_session,1}{counter_harvest_early,1};
                                else
                                    pad_val_within = count.max_num_sac_work - count.num_sac_work_early{counter_session,1}{counter_harvest_early,1};
                                end
                                if isempty(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,1})
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,:} = nan(1, pad_val_within);
                                else
                                    if pad_val_within >= 0
                                        data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,1},[0 pad_val_within], nan, 'post');
                                    else
                                        data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,1}(:,1:1:count.max_num_sac_work);
                                    end
                                end
                            end
                            data.(string(fieldname(counter_fieldname))){counter_session,1} =cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1});
                        elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},1) == 150
                            for counter_harvest_early = 1 : count.num_harvest_early(counter_session,1)
                                pad_val_within = count.max_num_sac_work - count.num_sac_work_early{counter_session,1}{counter_harvest_early,1};
                                if pad_val_within >= 0
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,1}',[pad_val_within 0], nan, 'post');
                                else
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_early,1}(:,1:1:count.max_num_sac_work)';
                                end
                            end
                        elseif contains(string(fieldname(counter_fieldname)),["duration_work", "num_sac_work", "rew_gained_work", "num_trial_work"])
                            pad_val_accross = count.max_num_harvest_early - count.num_harvest_early(counter_session,1);
                            if pad_val_accross >= 0
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = padarray(cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}),[pad_val_accross 0], nan, 'post')';
                            else
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}(1:1:count.max_num_harvest_early,1))';
                            end
                        end
                    elseif(contains(string(fieldname(counter_fieldname)),'_late'))
                        % check if bout related data or lick related data
                        if ~contains(string(fieldname(counter_fieldname)),["duration_work", "num_sac_work", "rew_gained_work", "num_trial_work"]) && size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},1) < 150
                            for counter_harvest_late = 1 : count.num_harvest_late(counter_session,1)
                                if(contains(string(fieldname(counter_fieldname)),["time_end_trial", "ITI","ITR"]))
                                    pad_val_within = count.max_num_trial_work - count.num_trial_work_late{counter_session,1}{counter_harvest_late,1};
                                else
                                    pad_val_within = count.max_num_sac_work - count.num_sac_work_late{counter_session,1}{counter_harvest_late,1};
                                end
                                if isempty(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,1})
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,:} = nan(1, pad_val_within);
                                else
                                    if pad_val_within >= 0
                                        data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,1},[0 pad_val_within], nan, 'post');
                                    else
                                        data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,1}(:,1:1:count.max_num_sac_work);
                                    end
                                end
                            end
                            data.(string(fieldname(counter_fieldname))){counter_session,1} =cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1});
                        elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},1) == 150
                            for counter_harvest_late = 1 : count.num_harvest_late(counter_session,1)
                                pad_val_within = count.max_num_sac_work - count.num_sac_work_late{counter_session,1}{counter_harvest_late,1};
                                if pad_val_within >= 0
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,1}',[pad_val_within 0], nan, 'post');
                                else
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest_late,1}(:,1:1:count.max_num_sac_work)';
                                end
                            end
                        elseif contains(string(fieldname(counter_fieldname)),["duration_work", "num_sac_work", "rew_gained_work", "num_trial_work"])
                            pad_val_accross = count.max_num_harvest_late - count.num_harvest_late(counter_session,1);
                            if pad_val_accross >= 0
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = padarray(cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}),[pad_val_accross 0], nan, 'post')';
                            else
                                data.(string(fieldname(counter_fieldname))){counter_session,1} = cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}(1:1:count.max_num_harvest_late,1))';
                            end
                        end
                    end
                end
            end
        end
    end
end
end
%% function compute_vigor
function [FIT] = compute_vigor(data_lick, data_sac, count, path, run_vigor)
if run_vigor == 1
    % build and concatenate amp and vel data
    for counter_session = 1 : count.num_session
        tag_lick_{1,counter_session} = data_lick.tag{counter_session, 1};
        tongue_dm_max_{1,counter_session} = data_lick.tongue_dm_max{counter_session, 1};
        tongue_vm_max_{1,counter_session} = data_lick.tongue_vm_max{counter_session, 1};
        tongue_vm_min_{1,counter_session} = abs(data_lick.tongue_vm_min{counter_session, 1});

        tag_sac_{1,counter_session} = data_sac.tag{counter_session, 1};
        eye_dm_max_{1,counter_session} = data_sac.eye_r_amp_m{counter_session, 1};
        eye_vm_max_{1,counter_session} =  max(data_sac.eye_r_vm{counter_session, 1});
    end
    tongue_dm_max = cell2mat(tongue_dm_max_);
    tongue_vm_max = cell2mat(tongue_vm_max_);
    tongue_vm_min = cell2mat(tongue_vm_min_);
    tag_lick = cell2mat(tag_lick_);

    eye_dm_max = cell2mat(eye_dm_max_);
    eye_vm_max = cell2mat(eye_vm_max_);
    tag_sac = cell2mat(tag_sac_);

    % separate based on tag
    tag_lick_use = 1; % grooming
    tongue_dm_max_groom = tongue_dm_max(ismember(tag_lick,tag_lick_use));
    tongue_vm_max_groom  = tongue_vm_max(ismember(tag_lick,tag_lick_use));
    tongue_vm_min_groom  = tongue_vm_min(ismember(tag_lick,tag_lick_use));

    tag_lick_use = 2:5; % reward
    tongue_dm_max_rew = tongue_dm_max(ismember(tag_lick,tag_lick_use));
    tongue_vm_max_rew  = tongue_vm_max(ismember(tag_lick,tag_lick_use));
    tongue_vm_min_rew  = tongue_vm_min(ismember(tag_lick,tag_lick_use));

    tag_sac_use = 10; % irrelevant
    eye_dm_max_irr = eye_dm_max(ismember(tag_sac,tag_sac_use));
    eye_vm_max_irr  = eye_vm_max(ismember(tag_sac,tag_sac_use));

    tag_sac_use = [1 4 6]; % relevant
    eye_dm_max_rel= eye_dm_max(ismember(tag_sac,tag_sac_use));
    eye_vm_max_rel = eye_vm_max(ismember(tag_sac,tag_sac_use));

    % amplitude based bins for vigor plot
    tongue_bins = [2.5 7.5 12.5 17.5 22.5 ];
    tongue_edges = [0:5:25];

    tongue_bins = [2.5:5:22.5];
    tongue_edges = [0:5:25];

    bin_tongue_groom = discretize(tongue_dm_max_groom,tongue_edges);
    bin_tongue_rew = discretize(tongue_dm_max_rew,tongue_edges);

    for counter_bin = 1:length(tongue_bins)
        mean_tongue_vm_max_groom_bin(counter_bin) = nanmean(tongue_vm_max_groom(bin_tongue_groom == counter_bin)) ;
        sem_tongue_vm_max_groom_bin(counter_bin) = nanstd(tongue_vm_max_groom(bin_tongue_groom == counter_bin))/sqrt(sum(bin_tongue_groom == counter_bin)) ;
        mean_tongue_vm_min_groom_bin(counter_bin) = nanmean(tongue_vm_min_groom(bin_tongue_groom == counter_bin)) ;
        sem_tongue_vm_min_groom_bin(counter_bin) = nanstd(tongue_vm_min_groom(bin_tongue_groom == counter_bin))/sqrt(sum(bin_tongue_groom == counter_bin)) ;
        mean_tongue_diff_groom_bin(counter_bin) = nanmean(tongue_vm_min_groom(bin_tongue_groom == counter_bin) - tongue_vm_max_groom(bin_tongue_groom == counter_bin)) ;
        sem_tongue_diff_groom_bin(counter_bin) = nanstd(tongue_vm_min_groom(bin_tongue_groom == counter_bin) - tongue_vm_max_groom(bin_tongue_groom == counter_bin))/sqrt(sum(bin_tongue_groom == counter_bin)) ;

        mean_tongue_vm_max_rew_bin(counter_bin) = nanmean(tongue_vm_max_rew(bin_tongue_rew == counter_bin)) ;
        sem_tongue_vm_max_rew_bin(counter_bin) = nanstd(tongue_vm_max_rew(bin_tongue_rew == counter_bin))/sqrt(sum(bin_tongue_rew == counter_bin)) ;
        mean_tongue_vm_min_rew_bin(counter_bin) = nanmean(tongue_vm_min_rew(bin_tongue_rew == counter_bin)) ;
        sem_tongue_vm_min_rew_bin(counter_bin) = nanstd(tongue_vm_min_rew(bin_tongue_rew == counter_bin))/sqrt(sum(bin_tongue_rew == counter_bin)) ;
        mean_tongue_diff_rew_bin(counter_bin) = nanmean(tongue_vm_min_rew(bin_tongue_rew == counter_bin) - tongue_vm_max_rew(bin_tongue_rew == counter_bin)) ;
        sem_tongue_diff_rew_bin(counter_bin) = nanstd(tongue_vm_min_rew(bin_tongue_rew == counter_bin) - tongue_vm_max_rew(bin_tongue_rew == counter_bin))/sqrt(sum(bin_tongue_rew == counter_bin)) ;

    end

    eye_bins = [1:20];
    eye_edges = [0:20];

    bin_eye_irr = discretize(eye_dm_max_irr,eye_edges);
    bin_eye_rel = discretize(eye_dm_max_rel,eye_edges);

    for counter_bin = 1:length(eye_bins)
        mean_eye_vm_max_irr_bin(counter_bin) = nanmean(eye_vm_max_irr(bin_eye_irr == counter_bin)) ;
        sem_eye_vm_max_irr_bin(counter_bin) = nanstd(eye_vm_max_irr(bin_eye_irr == counter_bin))/sqrt(sum(bin_eye_irr == counter_bin)) ;

        mean_eye_vm_max_rel_bin(counter_bin) = nanmean(eye_vm_max_rel(bin_eye_rel == counter_bin)) ;
        sem_eye_vm_max_rel_bin(counter_bin) = nanstd(eye_vm_max_rel(bin_eye_rel == counter_bin))/sqrt(sum(eye_vm_max_rel == counter_bin)) ;

    end

    plot_vigor = 1;
    if plot_vigor == 1
        % plot
        fig_vigor = figure;
        s1 = subplot(3,3,1);
        hold on
        x = tongue_dm_max_groom ;
        y = tongue_vm_max_groom ;
        B0 = [100; 0.1];
        f = fittype('a*(1-(1./(1+b*x)))');
        FIT.fit_lick_groom_pro = fit(x',y',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
        plot(x,y,'or', 'MarkerSize', 0.1)
        plot(FIT.fit_lick_groom_pro, 'r')

        x = tongue_dm_max_rew ;
        y = tongue_vm_max_rew ;
        B0 = [100; 0.1];
        f = fittype('a*(1-(1./(1+b*x)))');
        FIT.fit_lick_rew_pro = fit(x',y',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
        plot(x,y,'ob', 'MarkerSize', 0.1)
        plot(FIT.fit_lick_rew_pro, 'b')
        xlabel('lick amp. (mm)')
        ylabel('lick pro. speed (mm/s)')
        legend('off')
        title('groom (r) vs rew (b)')
        ylim([0 800])
        yticks([0 : 100: 800])
        xlim([0 25])
        xticks([0 : 5: 25])

        s2 = subplot(3,3,2);
        hold on
        x = tongue_dm_max_groom ;
        y = tongue_vm_min_groom ;
        B0 = [100; 0.1];
        f = fittype('a*(1-(1./(1+b*x)))');
        FIT.fit_lick_groom_ret = fit(x',y',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
        plot(x,y,'or', 'MarkerSize', 0.1)
        plot(FIT.fit_lick_groom_ret, 'r')

        x = tongue_dm_max_rew ;
        y = tongue_vm_min_rew ;
        B0 = [100; 0.1];
        f = fittype('a*(1-(1./(1+b*x)))');
        FIT.fit_lick_rew_ret = fit(x',y',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
        plot(x,y,'ob', 'MarkerSize', 0.1)
        plot(FIT.fit_lick_rew_ret, 'b')
        xlabel('lick amp. (mm)')
        ylabel('lick ret. speed (mm/s)')
        legend('off')
        title('groom (r) vs rew (b)')
        ylim([0 800])
        yticks([0 : 100: 800])
        xlim([0 25])
        xticks([0 : 5: 25])

        s3 = subplot(3,3,3);
        hold on
        x = eye_dm_max_irr ;
        y = eye_vm_max_irr ;
        B0 = [100; 0.1];
        f = fittype('a*(1-(1./(1+b*x)))');
        FIT.fit_sac_irr = fit(x',y',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
        plot(x,y,'or', 'MarkerSize', 0.1)
        plot(FIT.fit_sac_irr, 'r')

        x = eye_dm_max_rel ;
        y = eye_vm_max_rel ;
        B0 = [100; 0.1];
        f = fittype('a*(1-(1./(1+b*x)))');
        FIT.fit_sac_rel = fit(x',y',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
        plot(x,y,'ob', 'MarkerSize', 0.1)
        plot(FIT.fit_sac_rel, 'b')
        xlabel('saccade amp. (deg)')
        ylabel('saccade vel. (deg/s)')
        legend('off')
        title('task irr (r) vs task rel (b)')

        s4 = subplot(3,3,4);
        hold on
        plot(tongue_bins,mean_tongue_vm_max_groom_bin,'r')
        errorbar(tongue_bins,mean_tongue_vm_max_groom_bin, sem_tongue_vm_max_groom_bin, '.r', 'LineWidth', 2, 'MarkerSize', 20);
        plot(tongue_bins,mean_tongue_vm_max_rew_bin,'b')
        errorbar(tongue_bins,mean_tongue_vm_max_rew_bin, sem_tongue_vm_max_rew_bin, '.b', 'LineWidth', 2, 'MarkerSize', 20);
        xlabel('lick amp. (mm)')
        ylabel('lick pro. speed (mm/s)')

        s5 = subplot(3,3,5);
        hold on
        plot(tongue_bins,mean_tongue_vm_min_groom_bin,'r')
        errorbar(tongue_bins,mean_tongue_vm_min_groom_bin, sem_tongue_vm_min_groom_bin, '.r', 'LineWidth', 2, 'MarkerSize', 20);
        plot(tongue_bins,mean_tongue_vm_min_rew_bin,'b')
        errorbar(tongue_bins,mean_tongue_vm_min_rew_bin, sem_tongue_vm_min_rew_bin, '.b', 'LineWidth', 2, 'MarkerSize', 20);
        xlabel('lick amp. (mm)')
        ylabel('lick ret. speed (mm/s)')
        ylim([0 800])
        yticks([0 : 100: 800])
        xlim([0 25])
        xticks([0 : 5: 25])

        s6 = subplot(3,3,6);
        hold on
        plot(eye_bins,mean_eye_vm_max_irr_bin,'r')
        errorbar(eye_bins,mean_eye_vm_max_irr_bin, sem_eye_vm_max_irr_bin, '.r', 'LineWidth', 2, 'MarkerSize', 20);
        plot(eye_bins,mean_eye_vm_max_rel_bin,'b')
        errorbar(eye_bins,mean_eye_vm_max_rel_bin, sem_eye_vm_max_rel_bin, '.b', 'LineWidth', 2, 'MarkerSize', 20);
        xlabel('lick amp. (mm)')
        xlabel('saccade amp. (deg)')
        ylabel('saccade vel. (deg/s)')

        s7 = subplot(3,3,7);
        hold on
        plot(tongue_bins,mean_tongue_diff_groom_bin,'r')
        errorbar(tongue_bins,mean_tongue_diff_groom_bin, sem_tongue_diff_groom_bin, '.r', 'LineWidth', 2, 'MarkerSize', 20);
        plot(tongue_bins,mean_tongue_diff_rew_bin,'b')
        errorbar(tongue_bins,mean_tongue_diff_rew_bin, sem_tongue_diff_rew_bin, '.b', 'LineWidth', 2, 'MarkerSize', 20);
        xlabel('lick amp. (mm)')
        ylabel('diff. ret.-pro. (mm/s)')
        ylim([-100 300])
        yticks([-100 : 50: 300])
        xlim([0 25])
        xticks([0 : 5: 25])


        linkaxes([s1 s2 s4 s5], ['x','y']);
        linkaxes([s3 s6], ['x' ,'y']);

        sgtitle([path.path_data_monkey_sorted(6:end) ' | sessions: ' num2str(count.num_session) ' | licks (groom, rew): ' num2str(length(bin_tongue_groom)), ', ' num2str(length(bin_tongue_rew)) ' | sacs (irr, rel): ' num2str(length(bin_eye_irr)), ', ' num2str(length(bin_eye_rel))], 'Interpreter', 'none')
        ESN_Beautify_Plot
        fig_vigor.WindowState = 'maximized';
        saveas(gcf, [path.out_path 'VIGOR' filesep path.path_data_monkey_sorted '_vigor' ], 'pdf');
        close(gcf)
    end

    save([path.out_path 'VIGOR' filesep path.path_data_monkey_sorted '_FIT.mat' ], 'FIT', '-v7');
else
    FIT = [];
end
end
%% function split_high_low_rew_avail
function [data] = split_high_low_rew_avail(data,data_lick, count)
fieldname = fieldnames(data);
fieldname(contains(fieldname, '_rew_avail') | contains(fieldname, '_rew_hist')) = [];
for counter_session = 1 : count.num_session
    for counter_fieldname = 1 : length(fieldname)
        if sum(contains(string(fieldname),'_lick')) > 0
            count_event = count.max_num_lick_harvest;
        else
            count_event = count.max_num_sac_work;
        end
        if(~contains(string(fieldname(counter_fieldname)),'_early') && ~contains(string(fieldname(counter_fieldname)),'_late')) && ~contains(string(fieldname(counter_fieldname)),'_px') && ~contains(string(fieldname(counter_fieldname)),'_py')
            % check if bout related data or lick related data
            if size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) > 1 && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1} = nan(count.num_harvest(counter_session), count_event);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1} = nan(count.num_harvest(counter_session), count_event);
                % build
                is_high_rew = find(data_lick.rew{counter_session, 1}(:,1) > nanmedian(data_lick.rew{counter_session, 1}(:,1)));
                is_low_rew = find(data_lick.rew{counter_session, 1}(:,1) < nanmedian(data_lick.rew{counter_session, 1}(:,1)));
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1}(is_high_rew,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_high_rew,:);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1}(is_low_rew,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_low_rew,:);
            elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) == 1  && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1} = nan(1,count.max_num_harvest);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1} = nan(1,count.max_num_harvest);

                is_high_rew = find(data_lick.rew_harvest_str{counter_session, 1} > nanmedian(data_lick.rew_harvest_str{counter_session, 1}));
                is_low_rew = find(data_lick.rew_harvest_str{counter_session, 1} < nanmedian(data_lick.rew_harvest_str{counter_session, 1}));
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1}(1,is_high_rew) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_high_rew);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1}(1,is_low_rew) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_low_rew);
            end

        elseif (contains(string(fieldname(counter_fieldname)),'_early')) && ~contains(string(fieldname(counter_fieldname)),'_px') && ~contains(string(fieldname(counter_fieldname)),'_py')
            if size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) > 1 && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1} = nan(count.num_harvest_early(counter_session), count_event);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1} = nan(count.num_harvest_early(counter_session), count_event);
                % build
                is_high_rew_early = find(data_lick.rew_early{counter_session, 1}(:,1) > nanmedian(data_lick.rew_early{counter_session, 1}(:,1)));
                is_low_rew_early = find(data_lick.rew_early{counter_session, 1}(:,1) < nanmedian(data_lick.rew_early{counter_session, 1}(:,1)));
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1}(is_high_rew_early,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_high_rew_early,:);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1}(is_low_rew_early,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_low_rew_early,:);
            elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) == 1  && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1} = nan(1,count.max_num_harvest_early);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1} = nan(1,count.max_num_harvest_early);

                is_high_rew_early = find(data_lick.rew_harvest_str_early{counter_session, 1} > nanmedian(data_lick.rew_harvest_str_early{counter_session, 1}));
                is_low_rew_early = find(data_lick.rew_harvest_str_early{counter_session, 1} < nanmedian(data_lick.rew_harvest_str_early{counter_session, 1}));
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1}(1,is_high_rew_early) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_high_rew_early);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1}(1,is_low_rew_early) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_low_rew_early);
            end

        elseif (contains(string(fieldname(counter_fieldname)),'_late')) && ~contains(string(fieldname(counter_fieldname)),'_px') && ~contains(string(fieldname(counter_fieldname)),'_py')
            if size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) > 1 && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1} = nan(count.num_harvest_late(counter_session), count_event);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1} = nan(count.num_harvest_late(counter_session), count_event);
                % build
                is_high_rew_late = find(data_lick.rew_late{counter_session, 1}(:,1) > nanmedian(data_lick.rew_late{counter_session, 1}(:,1)));
                is_low_rew_late = find(data_lick.rew_late{counter_session, 1}(:,1) < nanmedian(data_lick.rew_late{counter_session, 1}(:,1)));
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1}(is_high_rew_late,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_high_rew_late,:);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1}(is_low_rew_late,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_low_rew_late,:);
            elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) == 1  && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1} = nan(1,count.max_num_harvest_late);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1} = nan(1,count.max_num_harvest_late);

                is_high_rew_late = find(data_lick.rew_harvest_str_late{counter_session, 1} > nanmedian(data_lick.rew_harvest_str_late{counter_session, 1}));
                is_low_rew_late = find(data_lick.rew_harvest_str_late{counter_session, 1} < nanmedian(data_lick.rew_harvest_str_late{counter_session, 1}));
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_avail')){counter_session,1}(1,is_high_rew_late) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_high_rew_late);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_avail')){counter_session,1}(1,is_low_rew_late) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_low_rew_late);
            end
        end
    end
end
end
%% function split_high_low_rew_hist
function [data] = split_high_low_rew_hist(data, data_lick, count)
fieldname = fieldnames(data);
fieldname(contains(fieldname, '_rew_avail') | contains(fieldname, '_rew_hist')) = [];
for counter_session = 1 : count.num_session
    for counter_fieldname = 1 : length(fieldname)
        if sum(contains(string(fieldname),'_lick')) > 0
            count_event = count.max_num_lick_harvest;
        else
            count_event = count.max_num_sac_work;
        end
        if(~contains(string(fieldname(counter_fieldname)),'_early') && ~contains(string(fieldname(counter_fieldname)),'_late')) && ~contains(string(fieldname(counter_fieldname)),'_px') && ~contains(string(fieldname(counter_fieldname)),'_py')
            % check if bout related data or sac related data
            if size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) > 1 && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1} = nan(count.num_harvest(counter_session), count_event);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1} = nan(count.num_harvest(counter_session), count_event);
                % build
                is_high_rew_ = data_lick.rew_consumed_harvest{counter_session, 1} > nanmedian(data_lick.rew_consumed_harvest{counter_session, 1});
                is_low_rew_ = data_lick.rew_consumed_harvest{counter_session, 1} < nanmedian(data_lick.rew_consumed_harvest{counter_session, 1});
                is_high_rew = logical(zeros(count.num_harvest(counter_session,1),1));
                is_low_rew = logical(zeros(count.num_harvest(counter_session,1),1));
                is_high_rew(find(is_high_rew_ == 1) + 1) = 1; is_high_rew = find(is_high_rew(1:count.num_harvest(counter_session,1)));
                is_low_rew(find(is_low_rew_ == 1) + 1) = 1; is_low_rew = find(is_low_rew(1:count.num_harvest(counter_session,1)));

                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1}(is_high_rew,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_high_rew,:);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1}(is_low_rew,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_low_rew,:);

            elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) == 1
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1} = nan(1,count.max_num_harvest);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1} = nan(1,count.max_num_harvest);
                %build
                is_high_rew_ = data_lick.rew_consumed_harvest{counter_session, 1} > nanmedian(data_lick.rew_consumed_harvest{counter_session, 1});
                is_low_rew_ = data_lick.rew_consumed_harvest{counter_session, 1} < nanmedian(data_lick.rew_consumed_harvest{counter_session, 1});
                is_high_rew = logical(zeros(1,count.max_num_harvest));
                is_low_rew = logical(zeros(1,count.max_num_harvest));
                is_high_rew(find(is_high_rew_ == 1) + 1) = 1;  is_high_rew = find(is_high_rew(1:count.max_num_harvest));
                is_low_rew(find(is_low_rew_ == 1) + 1) = 1;  is_low_rew = find(is_low_rew(1:count.max_num_harvest));

                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1}(1,is_high_rew) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_high_rew);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1}(1,is_low_rew) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_low_rew);
            end

        elseif (contains(string(fieldname(counter_fieldname)),'_early')) && ~contains(string(fieldname(counter_fieldname)),'_px') && ~contains(string(fieldname(counter_fieldname)),'_py')
            if size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) > 1 && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1} = nan(count.num_harvest(counter_session), count_event);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1} = nan(count.num_harvest(counter_session), count_event);
                % build
                is_high_rew_early_ = data_lick.rew_consumed_harvest_early{counter_session, 1} > nanmedian(data_lick.rew_consumed_harvest_early{counter_session, 1});
                is_low_rew_early_ = data_lick.rew_consumed_harvest_early{counter_session, 1} < nanmedian(data_lick.rew_consumed_harvest_early{counter_session, 1});
                is_high_rew_early = logical(zeros(count.num_harvest_early(counter_session,1),1));
                is_low_rew_early = logical(zeros(count.num_harvest_early(counter_session,1),1));
                is_high_rew_early(find(is_high_rew_early_ == 1) + 1) = 1; is_high_rew_early = find(is_high_rew_early(1:count.num_harvest_early(counter_session,1)));
                is_low_rew_early(find(is_low_rew_early_ == 1) + 1) = 1; is_low_rew_early = find(is_low_rew_early(1:count.num_harvest_early(counter_session,1)));

                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1}(is_high_rew_early,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_high_rew_early,:);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1}(is_low_rew_early,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_low_rew_early,:);

            elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) == 1
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1} = nan(1,count.max_num_harvest_early);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1} = nan(1,count.max_num_harvest_early);
                %build
                is_high_rew_early_ = data_lick.rew_consumed_harvest_early{counter_session, 1} > nanmedian(data_lick.rew_consumed_harvest_early{counter_session, 1});
                is_low_rew_early_ = data_lick.rew_consumed_harvest_early{counter_session, 1} < nanmedian(data_lick.rew_consumed_harvest_early{counter_session, 1});
                is_high_rew_early = logical(zeros(1,count.max_num_harvest_early));
                is_low_rew_early = logical(zeros(1,count.max_num_harvest_early));
                is_high_rew_early(find(is_high_rew_early_ == 1) + 1) = 1; is_high_rew_early = find(is_high_rew_early(1:count.max_num_harvest_early));
                is_low_rew_early(find(is_low_rew_early_ == 1) + 1) = 1; is_low_rew_early = find(is_low_rew_early(1:count.max_num_harvest_early));

                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1}(1,is_high_rew_early) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_high_rew_early);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1}(1,is_low_rew_early) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_low_rew_early);
            end

        elseif (contains(string(fieldname(counter_fieldname)),'_late')) && ~contains(string(fieldname(counter_fieldname)),'_px') && ~contains(string(fieldname(counter_fieldname)),'_py')
            if size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) > 1 && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                %initialize
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1} = nan(count.num_harvest(counter_session), count_event);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1} = nan(count.num_harvest(counter_session), count_event);
                % build
                is_high_rew_late_ = data_lick.rew_consumed_harvest_late{counter_session, 1} > nanmedian(data_lick.rew_consumed_harvest_late{counter_session, 1});
                is_low_rew_late_ = data_lick.rew_consumed_harvest_late{counter_session, 1} < nanmedian(data_lick.rew_consumed_harvest_late{counter_session, 1});
                is_high_rew_late = logical(zeros(count.num_harvest_late(counter_session,1),1));
                is_low_rew_late = logical(zeros(count.num_harvest_late(counter_session,1),1));
                is_high_rew_late(find(is_high_rew_late_ == 1) + 1) = 1; is_high_rew_late = find(is_high_rew_late(1:count.num_harvest_late(counter_session,1)));
                is_low_rew_late(find(is_low_rew_late_ == 1) + 1) = 1; is_low_rew_late = find(is_low_rew_late(1:count.num_harvest_late(counter_session,1)));

                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1}(is_high_rew_late,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_high_rew_late,:);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1}(is_low_rew_late,:) = data.(string(fieldname(counter_fieldname))){counter_session,1}(is_low_rew_late,:);


            elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) == 1
                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1} = nan(1,count.max_num_harvest_late);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1} = nan(1,count.max_num_harvest_late);
                %build
                is_high_rew_late_ = data_lick.rew_consumed_harvest_late{counter_session, 1} > nanmedian(data_lick.rew_consumed_harvest_late{counter_session, 1});
                is_low_rew_late_ = data_lick.rew_consumed_harvest_late{counter_session, 1} < nanmedian(data_lick.rew_consumed_harvest_late{counter_session, 1});
                is_high_rew_late = logical(zeros(1,count.max_num_harvest_late));
                is_low_rew_late = logical(zeros(1,count.max_num_harvest_late));
                is_high_rew_late(find(is_high_rew_late_ == 1) + 1) = 1; is_high_rew_late = find(is_high_rew_late(1:count.max_num_harvest_late));
                is_low_rew_late(find(is_low_rew_late_ == 1) + 1) = 1; is_low_rew_late = find(is_low_rew_late(1:count.max_num_harvest_late));

                data.(strcat(string(fieldname(counter_fieldname)), '_high_rew_hist')){counter_session,1}(1,is_high_rew_late) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_high_rew_late);
                data.(strcat(string(fieldname(counter_fieldname)), '_low_rew_hist')){counter_session,1}(1,is_low_rew_late) = data.(string(fieldname(counter_fieldname))){counter_session,1}(1,is_low_rew_late);
            end
        end
    end
end
end

%% function split_RPE
function [data] = split_RPE(data, count)
fieldname = fieldnames(data);
% fieldname = {'ILR';'tongue_vigor_pro'; 'tongue_vigor_ret'; 'tongue_vm_max'; 'tongue_vm_min'; 'tongue_ang_max'; 'tongue_dm_max'; 'tongue_duration';'dist_err_tongue_rew'; 'rew'; 'rew_consumed'; ...
%     'ILR_early_high_rew_avail';'tongue_vigor_pro_early_high_rew_avail'; 'tongue_vigor_ret_early_high_rew_avail'; 'tongue_vm_max_early_high_rew_avail'; 'tongue_vm_min_early_high_rew_avail'; 'tongue_ang_max_early_high_rew_avail'; 'tongue_dm_max_early_high_rew_avail'; 'tongue_duration_early_high_rew_avail';'dist_err_tongue_rew_early_high_rew_avail'; 'rew_early_high_rew_avail'; 'rew_consumed_early_high_rew_avail'; ...
%     'ILR_early_low_rew_avail';'tongue_vigor_pro_early_low_rew_avail'; 'tongue_vigor_ret_early_low_rew_avail'; 'tongue_vm_max_early_low_rew_avail'; 'tongue_vm_min_early_low_rew_avail'; 'tongue_ang_max_early_low_rew_avail'; 'tongue_dm_max_early_low_rew_avail'; 'tongue_duration_early_low_rew_avail';'dist_err_tongue_rew_early_low_rew_avail'; 'rew_early_low_rew_avail'; 'rew_consumed_early_low_rew_avail'; ...
%     'ILR_late_high_rew_avail';'tongue_vigor_pro_late_high_rew_avail'; 'tongue_vigor_ret_late_high_rew_avail'; 'tongue_vm_max_late_high_rew_avail'; 'tongue_vm_min_late_high_rew_avail'; 'tongue_ang_max_late_high_rew_avail'; 'tongue_dm_max_late_high_rew_avail'; 'tongue_duration_late_high_rew_avail';'dist_err_tongue_rew_late_high_rew_avail'; 'rew_late_high_rew_avail'; 'rew_consumed_late_high_rew_avail'; ...
%     'ILR_late_low_rew_avail';'tongue_vigor_pro_late_low_rew_avail'; 'tongue_vigor_ret_late_low_rew_avail'; 'tongue_vm_max_late_low_rew_avail'; 'tongue_vm_min_late_low_rew_avail'; 'tongue_ang_max_late_low_rew_avail'; 'tongue_dm_max_late_low_rew_avail'; 'tongue_duration_late_low_rew_avail';'dist_err_tongue_rew_late_low_rew_avail'; 'rew_late_low_rew_avail'; 'rew_consumed_late_low_rew_avail'};


for counter_session = 1 : count.num_session
    for counter_fieldname = 1 : length(fieldname)
        if(~contains(string(fieldname(counter_fieldname)),'_early') && ~contains(string(fieldname(counter_fieldname)),'_late')) && ~contains(string(fieldname(counter_fieldname)),'_px') && ~contains(string(fieldname(counter_fieldname)),'_py')
            for counter_harvest = 1 : count.num_harvest(counter_session)
                % check if bout related data or lick related data
                if size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) > 1 && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                    % initialize
                    is_rpe_TT_curr = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FF_curr = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FT_curr = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_TF_curr = logical(zeros(1,count.max_num_lick_harvest));

                    is_rpe_TT_prev = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FF_prev = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FT_prev = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_TF_prev = logical(zeros(1,count.max_num_lick_harvest));

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_curr')){counter_session,1} = nan(count.num_harvest(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_curr')){counter_session,1} = nan(count.num_harvest(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_curr')){counter_session,1} = nan(count.num_harvest(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_curr')){counter_session,1} = nan(count.num_harvest(counter_session),count.max_num_lick_harvest);

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_prev')){counter_session,1} = nan(count.num_harvest(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_prev')){counter_session,1} = nan(count.num_harvest(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_prev')){counter_session,1} = nan(count.num_harvest(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_prev')){counter_session,1} = nan(count.num_harvest(counter_session),count.max_num_lick_harvest);

                    for counter_lick = 2 : count.max_num_lick_harvest
                        is_rpe_TT_curr(1,counter_lick) = (data.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 2 || data.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 4 ) && ...
                            (data.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 2 || data.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 4);
                        is_rpe_FF_curr(1,counter_lick) = (data.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 3 || data.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 5 ) && ...
                            (data.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 3 || data.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 5);
                        is_rpe_FT_curr(1,counter_lick) = (data.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 2 || data.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 4 ) && ...
                            (data.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 3 || data.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 5);
                        is_rpe_TF_curr(1,counter_lick) = (data.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 3 || data.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 5 ) && ...
                            (data.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 2 || data.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 4);
                    end

                    % is for previous lick
                    is_rpe_TT_prev(find(is_rpe_TT_curr == 1) - 1) = 1;
                    is_rpe_FF_prev(find(is_rpe_FF_curr == 1) - 1) = 1;
                    is_rpe_FT_prev(find(is_rpe_FT_curr == 1) - 1) = 1;
                    is_rpe_TF_prev(find(is_rpe_TF_curr == 1) - 1) = 1;

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_curr')){counter_session,1}(counter_harvest,1:sum(is_rpe_TT_curr)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest, is_rpe_TT_curr);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_curr')){counter_session,1}(counter_harvest,1:sum(is_rpe_FF_curr)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest, is_rpe_FF_curr);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_curr')){counter_session,1}(counter_harvest,1:sum(is_rpe_FT_curr)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest, is_rpe_FT_curr);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_curr')){counter_session,1}(counter_harvest,1:sum(is_rpe_TF_curr)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest, is_rpe_TF_curr);

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_prev')){counter_session,1}(counter_harvest,1:sum(is_rpe_TT_prev)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest, is_rpe_TT_prev);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_prev')){counter_session,1}(counter_harvest,1:sum(is_rpe_FF_prev)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest, is_rpe_FF_prev);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_prev')){counter_session,1}(counter_harvest,1:sum(is_rpe_FT_prev)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest, is_rpe_FT_prev);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_prev')){counter_session,1}(counter_harvest,1:sum(is_rpe_TF_prev)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest, is_rpe_TF_prev);
                end
            end
        elseif (contains(string(fieldname(counter_fieldname)),'_early')) && ~contains(string(fieldname(counter_fieldname)),'_px') && ~contains(string(fieldname(counter_fieldname)),'_py')
            for counter_harvest_early = 1 : count.num_harvest_early(counter_session)
                % check if bout related data or lick related data
                if size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) > 1 && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                    % initialize
                    is_rpe_TT_curr_early = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FF_curr_early = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FT_curr_early = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_TF_curr_early = logical(zeros(1,count.max_num_lick_harvest));

                    is_rpe_TT_prev_early = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FF_prev_early = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FT_prev_early = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_TF_prev_early = logical(zeros(1,count.max_num_lick_harvest));

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_curr')){counter_session,1} = nan(count.num_harvest_early(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_curr')){counter_session,1} = nan(count.num_harvest_early(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_curr')){counter_session,1} = nan(count.num_harvest_early(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_curr')){counter_session,1} = nan(count.num_harvest_early(counter_session),count.max_num_lick_harvest);

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_prev')){counter_session,1} = nan(count.num_harvest_early(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_prev')){counter_session,1} = nan(count.num_harvest_early(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_prev')){counter_session,1} = nan(count.num_harvest_early(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_prev')){counter_session,1} = nan(count.num_harvest_early(counter_session),count.max_num_lick_harvest);

                    for counter_lick = 2 : count.max_num_lick_harvest
                        is_rpe_TT_curr_early(1,counter_lick) = (data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick) == 2 || data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick) == 4 ) && ...
                            (data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick - 1) == 2 || data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick - 1) == 4);
                        is_rpe_FF_curr_early(1,counter_lick) = (data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick) == 3 || data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick) == 5 ) && ...
                            (data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick - 1) == 3 || data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick - 1) == 5);
                        is_rpe_FT_curr_early(1,counter_lick) = (data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick) == 2 || data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick) == 4 ) && ...
                            (data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick - 1) == 3 || data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick - 1) == 5);
                        is_rpe_TF_curr_early(1,counter_lick) = (data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick) == 3 || data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick) == 5 ) && ...
                            (data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick - 1) == 2 || data.tag_lick_early{counter_session, 1}(counter_harvest_early,counter_lick - 1) == 4);
                    end

                    % is for previous lick
                    is_rpe_TT_prev_early(find(is_rpe_TT_curr_early == 1) - 1) = 1;
                    is_rpe_FF_prev_early(find(is_rpe_FF_curr_early == 1) - 1) = 1;
                    is_rpe_FT_prev_early(find(is_rpe_FT_curr_early == 1) - 1) = 1;
                    is_rpe_TF_prev_early(find(is_rpe_TF_curr_early == 1) - 1) = 1;

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_curr')){counter_session,1}(counter_harvest_early,1:sum(is_rpe_TT_curr_early)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_early, is_rpe_TT_curr_early);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_curr')){counter_session,1}(counter_harvest_early,1:sum(is_rpe_FF_curr_early)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_early, is_rpe_FF_curr_early);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_curr')){counter_session,1}(counter_harvest_early,1:sum(is_rpe_FT_curr_early)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_early, is_rpe_FT_curr_early);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_curr')){counter_session,1}(counter_harvest_early,1:sum(is_rpe_TF_curr_early)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_early, is_rpe_TF_curr_early);

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_prev')){counter_session,1}(counter_harvest_early,1:sum(is_rpe_TT_prev_early)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_early, is_rpe_TT_prev_early);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_prev')){counter_session,1}(counter_harvest_early,1:sum(is_rpe_FF_prev_early)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_early, is_rpe_FF_prev_early);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_prev')){counter_session,1}(counter_harvest_early,1:sum(is_rpe_FT_prev_early)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_early, is_rpe_FT_prev_early);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_prev')){counter_session,1}(counter_harvest_early,1:sum(is_rpe_TF_prev_early)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_early, is_rpe_TF_prev_early);
                end
            end

        elseif (contains(string(fieldname(counter_fieldname)),'_late')) && ~contains(string(fieldname(counter_fieldname)),'_px') && ~contains(string(fieldname(counter_fieldname)),'_py')
            for counter_harvest_late = 1 : count.num_harvest_late(counter_session)
                % check if bout related data or lick related data
                if size(data.(string(fieldname(counter_fieldname))){counter_session,1},1) > 1 && ~iscell(data.(string(fieldname(counter_fieldname))){counter_session,1})
                    % initialize
                    is_rpe_TT_curr_late = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FF_curr_late = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FT_curr_late = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_TF_curr_late = logical(zeros(1,count.max_num_lick_harvest));

                    is_rpe_TT_prev_late = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FF_prev_late = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_FT_prev_late = logical(zeros(1,count.max_num_lick_harvest));
                    is_rpe_TF_prev_late = logical(zeros(1,count.max_num_lick_harvest));

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_curr')){counter_session,1} = nan(count.num_harvest_late(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_curr')){counter_session,1} = nan(count.num_harvest_late(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_curr')){counter_session,1} = nan(count.num_harvest_late(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_curr')){counter_session,1} = nan(count.num_harvest_late(counter_session),count.max_num_lick_harvest);

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_prev')){counter_session,1} = nan(count.num_harvest_late(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_prev')){counter_session,1} = nan(count.num_harvest_late(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_prev')){counter_session,1} = nan(count.num_harvest_late(counter_session),count.max_num_lick_harvest);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_prev')){counter_session,1} = nan(count.num_harvest_late(counter_session),count.max_num_lick_harvest);

                    for counter_lick = 2 : count.max_num_lick_harvest
                        is_rpe_TT_curr_late(1,counter_lick) = (data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick) == 2 || data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick) == 4 ) && ...
                            (data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick - 1) == 2 || data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick - 1) == 4);
                        is_rpe_FF_curr_late(1,counter_lick) = (data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick) == 3 || data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick) == 5 ) && ...
                            (data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick - 1) == 3 || data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick - 1) == 5);
                        is_rpe_FT_curr_late(1,counter_lick) = (data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick) == 2 || data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick) == 4 ) && ...
                            (data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick - 1) == 3 || data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick - 1) == 5);
                        is_rpe_TF_curr_late(1,counter_lick) = (data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick) == 3 || data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick) == 5 ) && ...
                            (data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick - 1) == 2 || data.tag_lick_late{counter_session, 1}(counter_harvest_late,counter_lick - 1) == 4);
                    end

                    % is for previous lick
                    is_rpe_TT_prev_late(find(is_rpe_TT_curr_late == 1) - 1) = 1;
                    is_rpe_FF_prev_late(find(is_rpe_FF_curr_late == 1) - 1) = 1;
                    is_rpe_FT_prev_late(find(is_rpe_FT_curr_late == 1) - 1) = 1;
                    is_rpe_TF_prev_late(find(is_rpe_TF_curr_late == 1) - 1) = 1;

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_curr')){counter_session,1}(counter_harvest_late,1:sum(is_rpe_TT_curr_late)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_late, is_rpe_TT_curr_late);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_curr')){counter_session,1}(counter_harvest_late,1:sum(is_rpe_FF_curr_late)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_late, is_rpe_FF_curr_late);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_curr')){counter_session,1}(counter_harvest_late,1:sum(is_rpe_FT_curr_late)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_late, is_rpe_FT_curr_late);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_curr')){counter_session,1}(counter_harvest_late,1:sum(is_rpe_TF_curr_late)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_late, is_rpe_TF_curr_late);

                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TT_prev')){counter_session,1}(counter_harvest_late,1:sum(is_rpe_TT_prev_late)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_late, is_rpe_TT_prev_late);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FF_prev')){counter_session,1}(counter_harvest_late,1:sum(is_rpe_FF_prev_late)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_late, is_rpe_FF_prev_late);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_FT_prev')){counter_session,1}(counter_harvest_late,1:sum(is_rpe_FT_prev_late)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_late, is_rpe_FT_prev_late);
                    data.(strcat(string(fieldname(counter_fieldname)), '_rpe_TF_prev')){counter_session,1}(counter_harvest_late,1:sum(is_rpe_TF_prev_late)) = data.(string(fieldname(counter_fieldname))){counter_session,1}(counter_harvest_late, is_rpe_TF_prev_late);

                end
            end
        end
    end
end
end


%% function move_type_prob
function [data] = move_type_prob(data, count)
fieldname = fieldnames(data);
fieldname(~contains(fieldname, 'tag_lick') | contains(fieldname, 'prev') | contains(fieldname, 'curr')) = [];

for counter_session = 1 : count.num_session
    for counter_fieldname = 1 : length(fieldname)

        count_valid_harvest = sum(~isnan(data.(string(fieldname(counter_fieldname))){counter_session,1}(:,1)));

        max_tag = max(data.(string(fieldname(counter_fieldname))){counter_session,1}, [], 'all');

        for counter_tag = 1 : max_tag
            data.(strcat('prob_', string(fieldname(counter_fieldname)),'_', num2str(counter_tag))){counter_session,1} = sum(data.(string(fieldname(counter_fieldname))){counter_session,1} == counter_tag)/count_valid_harvest;
        end

        % form rew seeking and grooming
        data.(strcat('prob_', string(fieldname(counter_fieldname)),'_', 'rewseek')){counter_session,1} = sum(data.(string(fieldname(counter_fieldname))){counter_session,1} ~= 1)/count_valid_harvest;


    end

end
end

%% PLOT FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function plot_session_analysis
function plot_session_analysis(params,path_data_monkey_sorted)
%% Load data
path_BEHAVE_population_data = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted];
if ~strcmp(path_BEHAVE_population_data(end), filesep);path_BEHAVE_population_data = [path_BEHAVE_population_data filesep];end
out_path = [path_BEHAVE_population_data filesep 'SESSION' filesep];
load([path_BEHAVE_population_data 'population_experiment_data.mat']);
%% Build data
probe_type_sess = cell2mat(population_experiment_data.probe_type_sess);
num_trial_sess = cell2mat(population_experiment_data.num_trial_sess);
num_unit_sess   = cell2mat(population_experiment_data.num_unit_sess);
num_licks = cell2mat(population_experiment_data.num_licks_sess);
num_sacs = cell2mat(population_experiment_data.num_sacs_sess);

mean_num_sacs = nanmean(num_sacs);
std_num_sacs = nanstd(num_sacs);

mean_num_licks = nanmean(num_licks);
std_num_licks = nanstd(num_licks);

mean_num_unit_sess = nanmean(num_unit_sess);
std_num_unit_sess = nanstd(num_unit_sess);

% silicon array
is_sa = probe_type_sess == 1 | probe_type_sess == 2;
num_sess_sa = sum(is_sa);

num_trial_sa = sum(num_trial_sess(is_sa));
mean_num_trial_sa = nanmean(num_trial_sess(is_sa));
std_num_trial_sa = nanstd(num_trial_sess(is_sa));
se_num_trial_sa = std_num_trial_sa/sqrt(num_sess_sa);

num_unit_sa = sum(num_unit_sess(is_sa));
mean_num_unit_sa = nanmean(num_unit_sess(is_sa));
std_num_unit_sa = nanstd(num_unit_sess(is_sa));
se_num_unit_sa = std_num_unit_sa/sqrt(num_sess_sa);

if num_unit_sa == 0
    mean_num_unit_sa = 0;
    std_num_unit_sa = 0;
    se_num_unit_sa = 0;
end

% hept/tet
is_ht = probe_type_sess == 3 | probe_type_sess == 4;
num_sess_ht = sum(is_ht);

num_trial_ht = sum(num_trial_sess(is_ht));
mean_num_trial_ht = nanmean(num_trial_sess(is_ht));
std_num_trial_ht = nanstd(num_trial_sess(is_ht));
se_num_trial_ht = std_num_trial_ht/sqrt(num_sess_ht);

num_unit_ht = sum(num_unit_sess(is_ht));
mean_num_unit_ht = nanmean(num_unit_sess(is_ht));
std_num_unit_ht = nanstd(num_unit_sess(is_ht));
se_num_unit_ht = std_num_unit_ht/sqrt(num_sess_ht);

% combinded
num_sess_all = num_sess_ht + num_sess_sa;

num_trial_all = num_trial_ht + num_trial_sa;
mean_num_trial_all = nanmean(num_trial_sess);
std_num_trial_all = nanstd(num_trial_sess);
se_num_trial_all = std_num_trial_all/sqrt(num_sess_all);

num_unit_all = num_unit_ht + num_unit_sa;

% extract num_unit_sess_type
num_unit_sess_sa = num_unit_sess;
num_unit_sess_sa(~is_sa) = nan;

num_unit_sess_ht = num_unit_sess;
num_unit_sess_ht(~is_ht) = nan;

for counter_session = 1 : num_sess_all
    rec_duration(counter_session, 1) = sum(population_experiment_data.rec_duration{counter_session, 1});
    break_duration(counter_session, 1) = sum(population_experiment_data.break_duration{counter_session, 1});
    session_duration(counter_session, 1) = rec_duration(counter_session, 1) + break_duration(counter_session, 1);
end

mean_rec_duration = mean(rec_duration);
std_rec_duration = std(rec_duration);

mean_break_duration = mean(break_duration);
std_break_duration = std(break_duration);

mean_session_duration = mean(session_duration);
std_session_duration = std(session_duration);

rec_duration.Format = 'hh:mm:ss';
break_duration.Format = 'hh:mm:ss';
session_duration.Format = 'hh:mm:ss';

%% Plot data - Fig 1
x_axis_session = 1:num_sess_all;
fig1 = figure;
fig1.WindowState = 'maximized';

subplot(3,3,1)
bar(x_axis_session, session_duration, 'FaceColor', 'black')
xlabel('recording session')
ylabel('time')
title('session duration')
xlim([0 length(num_trial_sess)+1])
ylim([duration([0 0 0]) duration([4 0 0])])
xticks(1:1:length(num_trial_sess))
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,2)
bar(x_axis_session, rec_duration, 'FaceColor', 'black')
xlabel('recording session')
ylabel('time')
title('rec duration')
xlim([0 length(num_trial_sess)+1])
ylim([duration([0 0 0]) duration([4 0 0])])
xticks(1:1:length(num_trial_sess))
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,3)
bar(x_axis_session, break_duration, 'FaceColor', 'black')
xlabel('recording session')
ylabel('time')
title('break duration')
xlim([0 length(num_trial_sess)+1])
ylim([duration([0 0 0]) duration([4 0 0])])
xticks(1:1:length(num_trial_sess))
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,4)
bar(x_axis_session, num_trial_sess, 'FaceColor', 'black')
xlabel('recording session')
ylabel('count')
title('num. correct trials')
xlim([0 length(num_trial_sess)+1])
xticks(1:1:length(num_trial_sess))
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,5)
bar(x_axis_session, num_sacs, 'FaceColor', 'black')
xlabel('recording session')
ylabel('count')
title('num. sacs')
xlim([0 length(num_trial_sess)+1])
xticks(1:1:length(num_trial_sess))
xlim([0 length(num_trial_sess)+1])
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,6)
bar(x_axis_session, num_licks, 'FaceColor', 'black')
xlabel('recording session')
ylabel('count')
title('num. licks')
xlim([0 length(num_trial_sess)+1])
xticks(1:1:length(num_trial_sess))
xlim([0 length(num_trial_sess)+1])
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,7)
bar(x_axis_session, num_unit_sess_sa, 'FaceColor', 'black')
xlabel('recording session')
ylabel('count')
title('num. neurons: silicon array')
xlim([0 length(num_trial_sess)+1])
xticks(1:1:length(num_trial_sess))
xlim([0 length(num_trial_sess)+1])
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,8)
bar(x_axis_session, num_unit_sess_ht, 'FaceColor', 'black')
xlabel('recording session')
ylabel('count')
title('num. neurons: tetrode/heptode')
xlim([0 length(num_trial_sess)+1])
xticks(1:1:length(num_trial_sess))
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,9)
pie([num_sess_sa, num_sess_ht], {['silicon array: ' num2str(num_sess_sa)], ['terode/heptode: ' num2str(num_sess_ht)]})
title('num. sessions w/ probe type')

sgtitle([path_data_monkey_sorted ' | sessions: ' num2str(num_sess_all) ], 'Interpreter', 'none')
ESN_Beautify_Plot
saveas(gcf, [out_path 'session_' num2str(get(gcf).Number)], 'pdf')

%% Plot data - Fig 2
fig2 = figure;
fig2.WindowState = 'maximized';
subplot(3,3,1)
histogram(session_duration, 'FaceColor', 'k');
xline(mean_session_duration, 'r', 'LineWidth', 1)
xline(mean_session_duration + std_session_duration, '--r', 'LineWidth', 0.5)
xline(mean_session_duration - std_session_duration, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('sess. duration')

subplot(3,3,2)
histogram(rec_duration, 'FaceColor', 'k');
xline(mean_rec_duration, 'r', 'LineWidth', 1)
xline(mean_rec_duration + std_rec_duration, '--r', 'LineWidth', 0.5)
xline(mean_rec_duration - std_rec_duration, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('rec. duration')

subplot(3,3,3)
histogram(break_duration, 'FaceColor', 'k');
xline(mean_break_duration, 'r', 'LineWidth', 1)
xline(mean_break_duration + std_break_duration, '--r', 'LineWidth', 0.5)
xline(mean_break_duration - std_break_duration, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('break duration')

subplot(3,3,4)
histogram(num_trial_sess, 'FaceColor', 'k');
xline(mean_num_trial_all, 'r', 'LineWidth', 1)
xline(mean_num_trial_all + std_num_trial_all, '--r', 'LineWidth', 0.5)
xline(mean_num_trial_all - std_num_trial_all, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('num. correct trials')

subplot(3,3,5)
histogram(num_sacs, 'FaceColor', 'k');
xline(mean_num_sacs, 'r', 'LineWidth', 1)
xline(mean_num_sacs + std_num_sacs, '--r', 'LineWidth', 0.5)
xline(mean_num_sacs - std_num_sacs, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('num. sacs')

subplot(3,3,6)
histogram(num_licks, 'FaceColor', 'k');
xline(mean_num_licks, 'r', 'LineWidth', 1)
xline(mean_num_licks + std_num_licks, '--r', 'LineWidth', 0.5)
xline(mean_num_licks - std_num_licks, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('num. licks')

subplot(3,3,7)
histogram(num_unit_sess(is_ht), 'FaceColor', 'k');
xline(mean_num_unit_ht, 'r', 'LineWidth', 1)
xline(mean_num_unit_ht + std_num_unit_ht, '--r', 'LineWidth', 0.5)
xline(mean_num_unit_ht - std_num_unit_ht, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('num. units w/ hept/tet')

subplot(3,3,8)
histogram(num_unit_sess(is_sa), 'FaceColor', 'k');
xline(mean_num_unit_sa, 'r', 'LineWidth', 1)
xline(mean_num_unit_sa + std_num_unit_sa, '--r', 'LineWidth', 0.5)
xline(mean_num_unit_sa - std_num_unit_sa, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('num. units w/ silicon array')

sgtitle([path_data_monkey_sorted ' | sessions: ' num2str(num_sess_all) ], 'Interpreter', 'none')
ESN_Beautify_Plot
saveas(gcf, [out_path 'session_' num2str(get(gcf).Number)], 'pdf')
%% Plot data - Fig 3
fig3 = figure;
fig3.WindowState = 'maximized';

subplot(3,3,1)
plot(rec_duration,num_trial_sess, '.k', 'MarkerSize', 20);
xlabel('rec. duration')
ylabel('num. correct trials')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,2)
plot(rec_duration,num_sacs, '.k', 'MarkerSize', 20);
xlabel('rec. duration')
ylabel('num. sacs')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,3)
plot(rec_duration,num_licks, '.k', 'MarkerSize', 20);
xlabel('rec. duration')
ylabel('num. licks')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,4)
plot(break_duration,num_trial_sess, '.k', 'MarkerSize', 20);
xlabel('break duration')
ylabel('num. correct trials')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,5)
plot(break_duration,num_sacs, '.k', 'MarkerSize', 20);
xlabel('break duration')
ylabel('num. sacs')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,6)
plot(break_duration,num_licks, '.k', 'MarkerSize', 20);
xlabel('break duration')
ylabel('num. licks')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,7)
plot(num_trial_sess,num_sacs, '.k', 'MarkerSize', 20);
xlabel('num. correct trials')
ylabel('num. sacs')
ylim([0 inf])
xlim([0 inf])

subplot(3,3,8)
plot(num_trial_sess,num_licks, '.k', 'MarkerSize', 20);
xlabel('num. correct trials')
ylabel('num. licks')
ylim([0 inf])
xlim([0 inf])

subplot(3,3,9)
plot(num_sacs,num_licks, '.k', 'MarkerSize', 20);
xlabel('num. sacs')
ylabel('num. licks')
ylim([0 inf])
xlim([0 inf])

sgtitle([path_data_monkey_sorted ' | sessions: ' num2str(num_sess_all) ], 'Interpreter', 'none')
ESN_Beautify_Plot
saveas(gcf, [out_path 'session_' num2str(get(gcf).Number)], 'pdf')

%% Plot data - Fig 4
fig4 = figure;
fig4.WindowState = 'maximized';

session_list = population_experiment_data.sess_date;
num_session = length(session_list);

num_col = floor(num_session/5);
num_row = ceil(num_session/floor(num_session/5));

for counter_session = 1 : num_session
    subplot(num_row, num_col, counter_session)
    imshow(population_experiment_data.vid_frame{counter_session, 1}  )
    title(session_list(counter_session,1))
end


sgtitle([path_data_monkey_sorted ' | sessions: ' num2str(num_sess_all) ], 'Interpreter', 'none')
ESN_Beautify_Plot
saveas(gcf, [out_path 'session_' num2str(get(gcf).Number)], 'pdf')
end
%% function plot_harvest_analysis
function plot_harvest_analysis(params,path_data_monkey_sorted)
%% Load data
fprintf('Loading data ... ')
path.path_data_monkey_sorted = path_data_monkey_sorted;
path_BEHAVE_population_data = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted];
if ~strcmp(path_BEHAVE_population_data(end), filesep);path_BEHAVE_population_data = [path_BEHAVE_population_data filesep];end
path.out_path = [path_BEHAVE_population_data];
load([path_BEHAVE_population_data 'population_experiment_data.mat']); % exp data
load([path_BEHAVE_population_data 'population_lick_data.mat']); % lick data
load([path_BEHAVE_population_data 'population_sac_data.mat']); % sac data
meta_data = readtable([path_BEHAVE_population_data 'meta_data.xls']); % meta data
include_tongue = find(meta_data{:,6} == 1);
include_eye = find(meta_data{:,7} == 1);

include = include_tongue(ismember(include_tongue,include_eye));

% exp data
field_name = string(fieldnames(population_experiment_data));
for counter_field_name = 1 : length(field_name)
    population_experiment_data.(field_name(counter_field_name)) = population_experiment_data.(field_name(counter_field_name))(include,1);
end

% lick data
field_name = string(fieldnames(population_lick_data));
for counter_field_name = 1 : length(field_name)
    population_lick_data.(field_name(counter_field_name)) = population_lick_data.(field_name(counter_field_name))(include,1);
end

% sac data
field_name = string(fieldnames(population_sac_data));
for counter_field_name = 1 : length(field_name)
    population_sac_data.(field_name(counter_field_name)) = population_sac_data.(field_name(counter_field_name))(include,1);
end

fprintf('---> Complete. \n')

%% Build variables
fprintf('Building variables ... ')
[count,lick,sac] = build_variables(population_experiment_data, population_lick_data,population_sac_data,path);

clearvars -except params lick sac path count
fprintf('---> Complete. \n')

%% compute endpoint error mag
fprintf('Computing endpoint error magnitude ... ')
[lick] = endpoint_error_mag(lick,count);
[sac] = endpoint_error_mag(sac,count);

clearvars -except params lick sac path count
fprintf('---> Complete. \n')
%% split bouts into early vs late
fprintf('Splitting data ... ')
[lick,count] = split_early_late(lick, count); % split lick data
[sac,count] = split_early_late(sac, count); % split sac data


clearvars -except params lick sac path count
fprintf('---> Complete. \n')

%% pad/cut variables
fprintf('Padding data ... ')
[lick] = pad_data(lick, count); % pad lick data
[sac] = pad_data(sac, count); % pad sac data

clearvars -except params lick sac path count
fprintf('---> Complete. \n')

%% split high/low rew
fprintf('Split based on reward  ... ')
[lick] = split_high_low_rew_avail(lick,lick, count); % split lick data based on reward available at harvest start
[sac] = split_high_low_rew_avail(sac,lick, count); % split sac data based on reward available at harvest start

[lick] = split_high_low_rew_hist(lick, lick, count); % split lick data based on reward consumed during previous harvest
[sac] = split_high_low_rew_hist(sac, lick, count); % split sac data based on reward consumed during previous harvest

clearvars -except params lick sac path count
fprintf('---> Complete. \n')
%% split RPE
fprintf('Split based on RPE  ... ')
[lick] = split_RPE(lick, count); % split lick data based on RPE

clearvars -except params lick sac path count
fprintf('---> Complete. \n')

%% compute end point var
fprintf('Computing endpoint error variance ... ')
% [lick] = endpoint_error_var(lick,count);
% [sac] = endpoint_error_var(sac,count);

clearvars -except params lick sac path count
fprintf('---> Complete. \n')

%% compute move type prob
fprintf('Computing move type probability... ')
[lick] = move_type_prob(lick,count);
[sac] = move_type_prob(sac,count);

clearvars -except params lick sac path count
fprintf('---> Complete. \n')
%% Figure - harvest period
variable_list = ["duration_harvest", "num_lick_harvest", "rew_consumed_harvest"];

for counter_variable = 1 : length(variable_list)
    eval(['fig' num2str(counter_variable) ' = figure;'])
    num_type = 1;
    num_category = 1;
    clearvars mean_ mean__ err_ub err_lb
    x_label = categorical({'all', 'early', 'late', 'early/high rew curr.', 'early/low rew curr.', 'late/high rew curr.', 'late/low rew curr.' , 'early/high rew hist.', 'early/low rew hist.', 'late/high rew hist.', 'late/low rew hist.' });
    x_label = reordercats(x_label,{'all', 'early', 'late', 'early/high rew curr.', 'early/low rew curr.', 'late/high rew curr.', 'late/low rew curr.', 'early/high rew hist.', 'early/low rew hist.', 'late/high rew hist.', 'late/low rew hist.' });

    mean_ = [nanmean(cell2mat(lick.(variable_list(counter_variable))),2) nanmean(cell2mat(lick.(strcat(variable_list(counter_variable), '_early'))),2) nanmean(cell2mat(lick.(strcat(variable_list(counter_variable), '_late'))),2) ...
        nanmean(cell2mat(lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail'))),2)  nanmean(cell2mat(lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail'))),2) ...
        nanmean(cell2mat(lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail'))),2)  nanmean(cell2mat(lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail'))),2)  ...
        nanmean(cell2mat(lick.(strcat(variable_list(counter_variable), '_early_high_rew_hist'))),2)  nanmean(cell2mat(lick.(strcat(variable_list(counter_variable), '_early_low_rew_hist'))),2) ...
        nanmean(cell2mat(lick.(strcat(variable_list(counter_variable), '_late_high_rew_hist'))),2)  nanmean(cell2mat(lick.(strcat(variable_list(counter_variable), '_late_low_rew_hist'))),2)];

    mean__ = nanmean(mean_);
    sem_ = nanstd(mean_)/sqrt(count.num_session);

    hold on
    plot(1 : size(mean_,2), mean_, '.r', 'MarkerSize', 10)
    errorbar(mean__, sem_, '.k', 'LineWidth', 2, 'MarkerSize', 30);
    xticklabels(x_label)
    %     xlabel('harvest type')
    ylabel(variable_list(counter_variable), 'Interpreter', 'none')
    xlim([0.5 11.5])
    xticks([1 : 11])

    sgtitle([path.path_data_monkey_sorted(6:end) ' | sessions: ' num2str(count.num_session) ], 'Interpreter', 'none')
    ESN_Beautify_Plot
    eval(['fig' num2str(counter_variable) '.WindowState = ''maximized'';'])
    saveas(gcf, [path.out_path  'HARVEST' filesep  path.path_data_monkey_sorted '_' char(variable_list(counter_variable))], 'pdf');

    is_for_paper = 1;
    if is_for_paper == 1
        % filter data based on variable and specify axes
        data_ = cell2mat(lick.(variable_list(counter_variable)));
        if contains(variable_list(counter_variable), 'duration')
            is_filter = data_ > 40;
            x_lim = [-0.5 40.5];
            x_ticks = [0 : 5 : 40];

        elseif contains(variable_list(counter_variable), 'num')
            is_filter = data_ > 80;
            x_lim = [-0.5 80.5];
            x_ticks = [0 : 5 : 80];
        else
            x_lim = [-0.5 1.2];
            x_ticks = [-0.6 : 0.2 : 1.2];
        end
        data_(is_filter ) = nan;

        figure
        % histogram
        subplot(1,2,1)
        histogram(data_,25,'Normalization','probability','FaceColor', 'black')
        xlabel(variable_list(counter_variable), 'Interpreter', 'none')
        ylabel('probability')
        xlim(x_lim)
        xticks(x_ticks)
        xline(nanmean(data_, 'all'), '-r', 'LineWidth', 1);
        xline(nanmean(data_, 'all') + nanstd(data_,0, 'all'), '--r', 'LineWidth', 0.5);
        xline(nanmean(data_, 'all') - nanstd(data_,0, 'all'), '--r', 'LineWidth', 0.5);

        % violin plot
        subplot(1,2,2)
        data__ = reshape(data_, [size(data_,1) * size(data_,2) 1]);
        violinplot(data__);
        ylim(x_lim)
        yticks(x_ticks)
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlabel('ID')


        sgtitle([path.path_data_monkey_sorted(6:end) ' | sessions: ' num2str(count.num_session) ' | harvest intervals: ' num2str(sum(count.num_harvest)) ' | mean+/-std: ' num2str(nanmean(data_, 'all')) '+/-' num2str(nanstd(data_,0, 'all')) ' | median+/-std: ' num2str(nanmedian(data_, 'all')) '+/-' num2str(nanstd(data_,0, 'all'))  ], 'Interpreter', 'none')
        ESN_Beautify_Plot(gcf, [10 8])
        saveas(gcf, [path.out_path  'HARVEST' filesep  path.path_data_monkey_sorted '_' char(variable_list(counter_variable)) '_hist'], 'pdf');






    end
end
%% Figure - work period
variable_list = ["duration_work","num_trial_work", "num_sac_work", "rew_gained_work"];

for counter_variable = 1 : length(variable_list)
    eval(['fig' num2str(counter_variable) ' = figure;'])
    num_type = 1;
    num_category = 1;
    clearvars mean_ mean__ err_ub err_lb
    x_label = categorical({'all', 'early', 'late', 'early/high rew curr.', 'early/low rew curr.', 'late/high rew curr.', 'late/low rew curr.' , 'early/high rew hist.', 'early/low rew hist.', 'late/high rew hist.', 'late/low rew hist.' });
    x_label = reordercats(x_label,{'all', 'early', 'late', 'early/high rew curr.', 'early/low rew curr.', 'late/high rew curr.', 'late/low rew curr.', 'early/high rew hist.', 'early/low rew hist.', 'late/high rew hist.', 'late/low rew hist.' });

    mean_ = [nanmean(cell2mat(sac.(variable_list(counter_variable))),2) nanmean(cell2mat(sac.(strcat(variable_list(counter_variable), '_early'))),2) nanmean(cell2mat(sac.(strcat(variable_list(counter_variable), '_late'))),2) ...
        nanmean(cell2mat(sac.(strcat(variable_list(counter_variable), '_early_high_rew_avail'))),2)  nanmean(cell2mat(sac.(strcat(variable_list(counter_variable), '_early_low_rew_avail'))),2) ...
        nanmean(cell2mat(sac.(strcat(variable_list(counter_variable), '_late_high_rew_avail'))),2)  nanmean(cell2mat(sac.(strcat(variable_list(counter_variable), '_late_low_rew_avail'))),2)  ...
        nanmean(cell2mat(sac.(strcat(variable_list(counter_variable), '_early_high_rew_hist'))),2)  nanmean(cell2mat(sac.(strcat(variable_list(counter_variable), '_early_low_rew_hist'))),2) ...
        nanmean(cell2mat(sac.(strcat(variable_list(counter_variable), '_late_high_rew_hist'))),2)  nanmean(cell2mat(sac.(strcat(variable_list(counter_variable), '_late_low_rew_hist'))),2)];

    mean__ = nanmean(mean_);
    sem_ = nanstd(mean_)/sqrt(count.num_session);

    hold on
    plot(1 : size(mean_,2), mean_, '.r', 'MarkerSize', 10)
    errorbar(mean__, sem_, '.k', 'LineWidth', 2, 'MarkerSize', 30);
    xticklabels(x_label)
    %     xlabel('work type')
    ylabel(variable_list(counter_variable), 'Interpreter', 'none')
    xlim([0.5 11.5])
    xticks([1 : 11])

    sgtitle([path.path_data_monkey_sorted(6:end) ' | sessions: ' num2str(count.num_session) ], 'Interpreter', 'none')
    ESN_Beautify_Plot
    eval(['fig' num2str(counter_variable) '.WindowState = ''maximized'';'])
    saveas(gcf, [path.out_path  'WORK' filesep  path.path_data_monkey_sorted '_' char(variable_list(counter_variable))], 'pdf');

    is_for_paper = 1;
    if is_for_paper == 1
        % filter data based on variable and specify axes
        data_ = cell2mat(sac.(variable_list(counter_variable)));
        if contains(variable_list(counter_variable), 'duration')
            is_filter = data_ > 60;
            x_lim = [-0.5 60.5];
            x_ticks = [0 : 5 : 60];

        elseif contains(variable_list(counter_variable), 'num')
            is_filter = data_ > 20 | data_ <= 0;
            x_lim = [-0.5 20.5];
            x_ticks = [0 : 5 : 20];
        else
            x_lim = [-0.5 1.2];
            x_ticks = [-0.6 : 0.2 : 1.2];
        end
        data_(is_filter ) = nan;
        figure
        % histogram
        subplot(1,2,1)
        histogram(data_,25,'Normalization','probability','FaceColor', 'black')
        xlabel(variable_list(counter_variable), 'Interpreter', 'none')
        ylabel('probability')
        xlim(x_lim)
        xticks(x_ticks)
        xline(nanmean(data_, 'all'), '-r', 'LineWidth', 1);
        xline(nanmean(data_, 'all') + nanstd(data_,0, 'all'), '--r', 'LineWidth', 0.5);
        xline(nanmean(data_, 'all') - nanstd(data_,0, 'all'), '--r', 'LineWidth', 0.5);

        % violin plot
        subplot(1,2,2)
        data__ = reshape(data_, [size(data_,1) * size(data_,2) 1]);
        violinplot(data__);
        ylim(x_lim)
        yticks(x_ticks)
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlabel('ID')
        sgtitle([path.path_data_monkey_sorted(6:end) ' | sessions: ' num2str(count.num_session) ' | work intervals: ' num2str(sum(count.num_harvest)) ' | mean+/-std: ' num2str(nanmean(data_, 'all')) '+/-' num2str(nanstd(data_,0, 'all'))  ' | median+/-std: ' num2str(nanmedian(data_, 'all')) '+/-' num2str(nanstd(data_,0, 'all')) ], 'Interpreter', 'none')
        ESN_Beautify_Plot(gcf, [10 8])
        saveas(gcf, [path.out_path  'WORK' filesep  path.path_data_monkey_sorted '_' char(variable_list(counter_variable)) '_hist'], 'pdf');
    end

end
%% Figure - lick harvest
% variable_list = ["ILR","tongue_vigor_pro", "tongue_vigor_ret","tongue_vm_max", "tongue_vm_min", "tongue_ang_max", "tongue_dm_max", "tongue_duration", "dist_err_tongue_rew", "det_var_err_tongue_rew" , "rew", "rew_consumed"];
variable_list = ["ILR","tongue_vigor_pro", "tongue_vigor_ret","tongue_vm_max", "tongue_vm_min", "tongue_ang_max", "tongue_dm_max", "tongue_duration", "dist_err_tongue_rew" , "rew", "rew_consumed", "prob_tag_lick"];
for counter_variable = 1 : length(variable_list)
    clearvars mean_ err_ub err_lb
    x_err = 1:count.max_num_lick_harvest;
    line_color = ['k';'r'; 'g'; 'b'; 'c'; 'm'; 'y';];
    eval(['fig' num2str(counter_variable) ' = figure;'])

    if ~strcmp(variable_list(counter_variable), "prob_tag_lick")
        % average accross sessions
        % all harvests
        s1 = subplot(2,3,1);
        hold on
        for counter_session = 1 : count.num_session
            mean_(counter_session,:) = nanmean(lick.(variable_list(counter_variable)){counter_session,1});
        end
        err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
        err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
        shade(x_err, err_ub, 'k', x_err, err_lb, 'k', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','k', 'FillColor', 'k');
        plot(nanmean(mean_),'k', 'LineWidth', 1);
        xlabel('lick num. in harvest')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlim([0 21])
        xticks([1 5:5:20])
        title('all harvests')

        s2 = subplot(2,3,2);
        % early vs late
        hold on
        for counter_session = 1 : count.num_session
            mean_(counter_session,:) = nanmean(lick.(strcat(variable_list(counter_variable), '_early')){counter_session,1});
        end
        err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
        err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
        shade(x_err, err_ub, 'r', x_err, err_lb, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
        plot(nanmean(mean_),'r', 'LineWidth', 1);

        for counter_session = 1 : count.num_session
            mean_(counter_session,:) = nanmean(lick.(strcat(variable_list(counter_variable), '_late')){counter_session,1});
        end
        err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
        err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
        shade(x_err, err_ub, 'b', x_err, err_lb, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
        plot(nanmean(mean_),'b', 'LineWidth', 1);
        xlabel('lick num. in harvest')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlim([0 21])
        xticks([1 5:5:20])
        title('early (r) vs late (b) harvests')

        s3 = subplot(2,3,3);
        % high vs low rew
        hold on
        for counter_session = 1 : count.num_session
            mean_(counter_session,:) = nanmean(lick.(strcat(variable_list(counter_variable), '_high_rew_avail')){counter_session,1});
        end
        err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
        err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
        shade(x_err, err_ub, 'r', x_err, err_lb, 'r', 'LineWidth',  0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
        plot(nanmean(mean_),'r', 'LineWidth', 1);

        for counter_session = 1 : count.num_session
            mean_(counter_session,:) = nanmean(lick.(strcat(variable_list(counter_variable), '_low_rew_avail')){counter_session,1});
        end
        err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
        err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
        shade(x_err, err_ub, 'b', x_err, err_lb, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
        plot(nanmean(mean_),'b', 'LineWidth', 1);
        xlabel('lick num. in harvest')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlim([0 21])
        xticks([1 5:5:20])
        title('high rew available (r) vs low rew available (b)')

        s4 = subplot(2,3,4);
        % high vs low rew during early
        hold on
        for counter_session = 1 : count.num_session
            mean_(counter_session,:) = nanmean(lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail')){counter_session,1});
        end
        err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
        err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
        shade(x_err, err_ub, 'r', x_err, err_lb, 'r', 'LineWidth',  0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
        plot(nanmean(mean_),'r', 'LineWidth', 1);

        for counter_session = 1 : count.num_session
            mean_(counter_session,:) = nanmean(lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail')){counter_session,1});
        end
        err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
        err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
        shade(x_err, err_ub, 'b', x_err, err_lb, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
        plot(nanmean(mean_),'b', 'LineWidth', 1);
        xlabel('lick num. in harvest')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlim([0 21])
        xticks([1 5:5:20])
        title('early harvests: high rew available (r) vs low rew available (b)')

        % high vs low rew during late
        s5 = subplot(2,3,5);
        hold on
        for counter_session = 1 : count.num_session
            mean_(counter_session,:) = nanmean(lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail')){counter_session,1});
        end
        err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
        err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
        shade(x_err, err_ub, 'r', x_err, err_lb, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
        plot(nanmean(mean_),'r', 'LineWidth', 1);

        for counter_session = 1 : count.num_session
            mean_(counter_session,:) = nanmean(lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail')){counter_session,1});
        end
        err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
        err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
        shade(x_err, err_ub, 'b', x_err, err_lb, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
        plot(nanmean(mean_),'b', 'LineWidth', 1);
        xlabel('lick num. in harvest')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlim([0 21])
        xticks([1 5:5:20])
        title('late harvests: high rew available (r) vs low rew available (b)')
    else
        % average accross sessions
        % all harvests
        s1 = subplot(2,3,1);
        hold on
        for counter_tag = 1 : count.num_tag_lick
            for counter_session = 1 : count.num_session
                prob_(counter_session,:) = (lick.(strcat(variable_list(counter_variable), '_', num2str(counter_tag))){counter_session,1});
            end
            err_ub = nanmean(prob_) + nanstd(prob_)/sqrt(count.num_session);
            err_lb = nanmean(prob_) - nanstd(prob_)/sqrt(count.num_session);
            shade(x_err, err_ub, line_color(counter_tag), x_err, err_lb, line_color(counter_tag), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',line_color(counter_tag), 'FillColor', line_color(counter_tag));
            plot(nanmean(prob_),line_color(counter_tag), 'LineWidth', 1);
        end
        xlabel('lick num. in harvest')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlim([0 21])
        xticks([1 5:5:20])
        title('all harvests | Gr(k), IS(r), IF(g), OS(b), OF(c)')

        s2 = subplot(2,3,2);
        % early
        hold on
        for counter_tag = 1 : count.num_tag_lick
            for counter_session = 1 : count.num_session
                prob_(counter_session,:) = (lick.(strcat(variable_list(counter_variable), '_early_', num2str(counter_tag))){counter_session,1});
            end
            err_ub = nanmean(prob_) + nanstd(prob_)/sqrt(count.num_session);
            err_lb = nanmean(prob_) - nanstd(prob_)/sqrt(count.num_session);
            shade(x_err, err_ub, line_color(counter_tag), x_err, err_lb, line_color(counter_tag), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',line_color(counter_tag), 'FillColor', line_color(counter_tag));
            plot(nanmean(prob_),line_color(counter_tag), 'LineWidth', 1);
        end
        xlabel('lick num. in harvest')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlim([0 21])
        xticks([1 5:5:20])
        title('early harvests | Gr(k), IS(r), IF(g), OS(b), OF(c)')

        s3 = subplot(2,3,3);
        % late
        hold on
        for counter_tag = 1 : count.num_tag_lick
            for counter_session = 1 : count.num_session
                prob_(counter_session,:) = (lick.(strcat(variable_list(counter_variable), '_late_', num2str(counter_tag))){counter_session,1});
            end
            err_ub = nanmean(prob_) + nanstd(prob_)/sqrt(count.num_session);
            err_lb = nanmean(prob_) - nanstd(prob_)/sqrt(count.num_session);
            shade(x_err, err_ub, line_color(counter_tag), x_err, err_lb, line_color(counter_tag), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',line_color(counter_tag), 'FillColor', line_color(counter_tag));
            plot(nanmean(prob_),line_color(counter_tag), 'LineWidth', 1);
        end
        xlabel('lick num. in harvest')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlim([0 21])
        xticks([1 5:5:20])
        title('late harvests | Gr(k), IS(r), IF(g), OS(b), OF(c)')

        s4 = subplot(2,3,4);
        % high rew available
        hold on
        for counter_tag = 1 : count.num_tag_lick
            for counter_session = 1 : count.num_session
                prob_(counter_session,:) = (lick.(strcat(variable_list(counter_variable), '_high_rew_avail_', num2str(counter_tag))){counter_session,1});
            end
            err_ub = nanmean(prob_) + nanstd(prob_)/sqrt(count.num_session);
            err_lb = nanmean(prob_) - nanstd(prob_)/sqrt(count.num_session);
            shade(x_err, err_ub, line_color(counter_tag), x_err, err_lb, line_color(counter_tag), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',line_color(counter_tag), 'FillColor', line_color(counter_tag));
            plot(nanmean(prob_),line_color(counter_tag), 'LineWidth', 1);
        end
        xlabel('lick num. in harvest')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlim([0 21])
        xticks([1 5:5:20])
        title('high rew available | Gr(k), IS(r), IF(g), OS(b), OF(c)')


        s5 = subplot(2,3,5);
        % low rew available
        hold on
        for counter_tag = 1 : count.num_tag_lick
            for counter_session = 1 : count.num_session
                prob_(counter_session,:) = (lick.(strcat(variable_list(counter_variable), '_low_rew_avail_', num2str(counter_tag))){counter_session,1});
            end
            err_ub = nanmean(prob_) + nanstd(prob_)/sqrt(count.num_session);
            err_lb = nanmean(prob_) - nanstd(prob_)/sqrt(count.num_session);
            shade(x_err, err_ub, line_color(counter_tag), x_err, err_lb, line_color(counter_tag), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',line_color(counter_tag), 'FillColor', line_color(counter_tag));
            plot(nanmean(prob_),line_color(counter_tag), 'LineWidth', 1);
        end
        xlabel('lick num. in harvest')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        xlim([0 21])
        xticks([1 5:5:20])
        title('low rew available | Gr(k), IS(r), IF(g), OS(b), OF(c)')
    end

    linkaxes([s1 s2 s3 s4 s5], 'y')
    sgtitle([path.path_data_monkey_sorted(6:end) ' | sessions: ' num2str(count.num_session) ], 'Interpreter', 'none')
    ESN_Beautify_Plot
    eval(['fig' num2str(counter_variable) '.WindowState = ''maximized'';'])
    saveas(gcf, [path.out_path  'HARVEST' filesep  path.path_data_monkey_sorted '_' char(variable_list(counter_variable))], 'pdf');
end

%% Figure - lick rpe
%variable_list = ["ILR","tongue_vigor_pro", "tongue_vigor_ret", "tongue_vm_max", "tongue_vm_min", "tongue_ang_max", "tongue_dm_max", "tongue_duration", "dist_err_tongue_rew", "det_var_err_tongue_rew", "rew", "rew_consumed"];
variable_list = ["ILR","tongue_vigor_pro", "tongue_vigor_ret", "tongue_vm_max", "tongue_vm_min", "tongue_ang_max", "tongue_dm_max", "tongue_duration", "dist_err_tongue_rew", "rew", "rew_consumed"];

for counter_variable = 1 : length(variable_list)
    x_label = categorical({'TT', 'FF', 'FT', 'TF','T', 'F'});
    x_label = reordercats(x_label,{'TT', 'FF', 'FT', 'TF', 'T' , 'F'});

    eval(['fig' num2str(counter_variable) ' = figure;'])

    % average accross sessions
    for counter_session = 1 : count.num_session
        % all harvests
        mean_prev_TT_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_rpe_TT_prev')){counter_session,1}));
        mean_prev_FF_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_rpe_FF_prev')){counter_session,1}));
        mean_prev_FT_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_rpe_FT_prev')){counter_session,1}));
        mean_prev_TF_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_rpe_TF_prev')){counter_session,1}));
        mean_prev_T_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_rpe_TF_prev')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_rpe_TT_prev')){counter_session,1}]));
        mean_prev_F_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_rpe_FT_prev')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_rpe_FF_prev')){counter_session,1}]));

        mean_curr_TT_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_rpe_TT_curr')){counter_session,1}));
        mean_curr_FF_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_rpe_FF_curr')){counter_session,1}));
        mean_curr_FT_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_rpe_FT_curr')){counter_session,1}));
        mean_curr_TF_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_rpe_TF_curr')){counter_session,1}));
        mean_curr_T_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_rpe_TF_curr')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_rpe_TT_curr')){counter_session,1}]));
        mean_curr_F_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_rpe_FT_curr')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_rpe_FF_curr')){counter_session,1}]));

        % early harvests
        mean_prev_TT_early_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_TT_prev')){counter_session,1}));
        mean_prev_FF_early_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_FF_prev')){counter_session,1}));
        mean_prev_FT_early_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_FT_prev')){counter_session,1}));
        mean_prev_TF_early_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_TF_prev')){counter_session,1}));
        mean_prev_T_early_high_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_TF_prev')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_TT_prev')){counter_session,1}]));
        mean_prev_F_early_high_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_FT_prev')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_FF_prev')){counter_session,1}]));

        mean_prev_TT_early_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_TT_prev')){counter_session,1}));
        mean_prev_FF_early_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_FF_prev')){counter_session,1}));
        mean_prev_FT_early_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_FT_prev')){counter_session,1}));
        mean_prev_TF_early_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_TF_prev')){counter_session,1}));
        mean_prev_T_early_low_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_TF_prev')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_TT_prev')){counter_session,1}]));
        mean_prev_F_early_low_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_FT_prev')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_FF_prev')){counter_session,1}]));

        mean_curr_TT_early_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_TT_curr')){counter_session,1}));
        mean_curr_FF_early_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_FF_curr')){counter_session,1}));
        mean_curr_FT_early_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_FT_curr')){counter_session,1}));
        mean_curr_TF_early_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_TF_curr')){counter_session,1}));
        mean_curr_T_early_high_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_TF_curr')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_TT_curr')){counter_session,1}]));
        mean_curr_F_early_high_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_FT_curr')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_early_high_rew_avail_rpe_FF_curr')){counter_session,1}]));

        mean_curr_TT_early_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_TT_curr')){counter_session,1}));
        mean_curr_FF_early_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_FF_curr')){counter_session,1}));
        mean_curr_FT_early_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_FT_curr')){counter_session,1}));
        mean_curr_TF_early_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_TF_curr')){counter_session,1}));
        mean_curr_T_early_low_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_TF_curr')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_TT_curr')){counter_session,1}]));
        mean_curr_F_early_low_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_FT_curr')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_early_low_rew_avail_rpe_FF_curr')){counter_session,1}]));

        % late harvests
        mean_prev_TT_late_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_TT_prev')){counter_session,1}));
        mean_prev_FF_late_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_FF_prev')){counter_session,1}));
        mean_prev_FT_late_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_FT_prev')){counter_session,1}));
        mean_prev_TF_late_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_TF_prev')){counter_session,1}));
        mean_prev_T_late_high_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_TF_prev')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_TT_prev')){counter_session,1}]));
        mean_prev_F_late_high_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_FT_prev')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_FF_prev')){counter_session,1}]));

        mean_prev_TT_late_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_TT_prev')){counter_session,1}));
        mean_prev_FF_late_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_FF_prev')){counter_session,1}));
        mean_prev_FT_late_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_FT_prev')){counter_session,1}));
        mean_prev_TF_late_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_TF_prev')){counter_session,1}));
        mean_prev_T_late_low_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_TT_prev')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_TT_prev')){counter_session,1}]));
        mean_prev_F_late_low_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_FT_prev')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_FF_prev')){counter_session,1}]));

        mean_curr_TT_late_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_TT_curr')){counter_session,1}));
        mean_curr_FF_late_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_FF_curr')){counter_session,1}));
        mean_curr_FT_late_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_FT_curr')){counter_session,1}));
        mean_curr_TF_late_high_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_TF_curr')){counter_session,1}));
        mean_curr_T_late_high_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_TF_curr')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_TT_curr')){counter_session,1}]));
        mean_curr_F_late_high_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_FT_curr')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_late_high_rew_avail_rpe_FF_curr')){counter_session,1}]));

        mean_curr_TT_late_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_TT_curr')){counter_session,1}));
        mean_curr_FF_late_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_FF_curr')){counter_session,1}));
        mean_curr_FT_late_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_FT_curr')){counter_session,1}));
        mean_curr_TF_late_low_rew_avail_(counter_session,:) = nanmean(nanmean(lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_TF_curr')){counter_session,1}));
        mean_curr_T_late_low_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_TF_curr')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_TT_curr')){counter_session,1}]));
        mean_curr_F_late_low_rew_avail_(counter_session,:) = nanmean(nanmean([lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_FT_curr')){counter_session,1}; lick.(strcat(variable_list(counter_variable), '_late_low_rew_avail_rpe_FF_curr')){counter_session,1}]));

    end
    num_event = 1;
    num_category = 5;
    for counter_event = 1 : num_event
        %%% all %%%
        s(counter_event) = subplot(3,2,[1]);
        hold on
        mean_ = [nanmean(mean_prev_TT_(:,counter_event))  nanmean(mean_curr_TT_(:,counter_event)); nanmean(mean_prev_FF_(:,counter_event))  nanmean(mean_curr_FF_(:,counter_event)); ...
            nanmean(mean_prev_FT_(:,counter_event))  nanmean(mean_curr_FT_(:,counter_event)); nanmean(mean_prev_TF_(:,counter_event))  nanmean(mean_curr_TF_(:,counter_event)); ...
            nanmean(mean_prev_T_(:,counter_event))  nanmean(mean_curr_T_(:,counter_event)); nanmean(mean_prev_F_(:,counter_event))  nanmean(mean_curr_F_(:,counter_event)) ];
        sem_ = [nanstd(mean_prev_TT_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_TT_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_FF_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_FF_(:,counter_event))/sqrt(count.num_session); ...
            nanstd(mean_prev_FT_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_FT_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_TF_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_TF_(:,counter_event))/sqrt(count.num_session); ...
            nanstd(mean_prev_T_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_T_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_F_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_F_(:,counter_event))/sqrt(count.num_session);];
        % bar(x_label,mean_, 'Facecolor', 'black')
        % add error bars
        [ngroups, nbars] = size(mean_);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, mean_(:,i), sem_(:,i), '.k', 'LineWidth', 2, 'MarkerSize', 30);
        end
        xlim([0.5 6.5])
        xticks([1 : 6])
        xticklabels(x_label)
        xline([1.5 2.5 3.5 4.5 5.5],'k', 'LineWidth',1)
        xlabel('RPE event')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        %      title(['all harvest'])
        title([' all harvest'])

        s2(counter_event) = subplot(3,2,[2]);
        hold on
        x_T_prev = cellstr(repmat(strcat(string(x_label(5)), 'prev'), length(mean_prev_T_), 1));
        x_T_curr = cellstr(repmat(strcat(string(x_label(5)), 'curr'), length(mean_curr_T_), 1));
        x_F_prev = cellstr(repmat(strcat(string(x_label(6)), 'prev'), length(mean_prev_T_), 1));
        x_F_curr = cellstr(repmat(strcat(string(x_label(6)), 'curr'), length(mean_curr_T_), 1));

        violinplot([mean_prev_T_(:,counter_event); mean_curr_T_(:,counter_event); mean_prev_F_(:,counter_event); mean_curr_F_(:,counter_event)], [x_T_prev;x_T_curr;x_F_prev;x_F_curr]);

        %%% early %%%
        % high rew
        s(counter_event + num_event*1) = subplot(3,2,3);
        hold on
        mean_ = [nanmean(mean_prev_TT_early_high_rew_avail_(:,counter_event))  nanmean(mean_curr_TT_early_high_rew_avail_(:,counter_event)); nanmean(mean_prev_FF_early_high_rew_avail_(:,counter_event))  nanmean(mean_curr_FF_early_high_rew_avail_(:,counter_event)); ...
            nanmean(mean_prev_FT_early_high_rew_avail_(:,counter_event))  nanmean(mean_curr_FT_early_high_rew_avail_(:,counter_event)); nanmean(mean_prev_TF_early_high_rew_avail_(:,counter_event))  nanmean(mean_curr_TF_early_high_rew_avail_(:,counter_event));  ...
            nanmean(mean_prev_T_early_high_rew_avail_(:,counter_event))  nanmean(mean_curr_T_early_high_rew_avail_(:,counter_event)); nanmean(mean_prev_F_early_high_rew_avail_(:,counter_event))  nanmean(mean_curr_F_early_high_rew_avail_(:,counter_event))];
        sem_ = [nanstd(mean_prev_TT_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_TT_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_FF_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_FF_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session); ...
            nanstd(mean_prev_FT_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_FT_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_TF_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_TF_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session);  ...
            nanstd(mean_prev_T_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_T_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_F_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_F_early_high_rew_avail_(:,counter_event))/sqrt(count.num_session)];
        %  bar(x_label,mean_, 'Facecolor', 'black')
        % add error bars
        [ngroups, nbars] = size(mean_);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, mean_(:,i), sem_(:,i), '.k', 'LineWidth', 2, 'MarkerSize', 30);
        end
        xlim([0.5 6.5])
        xticks([1 : 6])
        xticklabels(x_label)
        xline([1.5 2.5 3.5 4.5 5.5],'k', 'LineWidth',1)
        xlabel('RPE event')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        title(['early high reward'])

        % low rew
        s(counter_event + num_event*2) = subplot(3,2,4);
        hold on
        mean_ = [nanmean(mean_prev_TT_early_low_rew_avail_(:,counter_event))  nanmean(mean_curr_TT_early_low_rew_avail_(:,counter_event)); nanmean(mean_prev_FF_early_low_rew_avail_(:,counter_event))  nanmean(mean_curr_FF_early_low_rew_avail_(:,counter_event)); ...
            nanmean(mean_prev_FT_early_low_rew_avail_(:,counter_event))  nanmean(mean_curr_FT_early_low_rew_avail_(:,counter_event)); nanmean(mean_prev_TF_early_low_rew_avail_(:,counter_event))  nanmean(mean_curr_TF_early_low_rew_avail_(:,counter_event)) ; ...
            nanmean(mean_prev_T_early_low_rew_avail_(:,counter_event))  nanmean(mean_curr_T_early_low_rew_avail_(:,counter_event)); nanmean(mean_prev_F_early_low_rew_avail_(:,counter_event))  nanmean(mean_curr_F_early_low_rew_avail_(:,counter_event))];
        sem_ = [nanstd(mean_prev_TT_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_TT_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_FF_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_FF_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session); ...
            nanstd(mean_prev_FT_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_FT_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_TF_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_TF_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session) ; ...
            nanstd(mean_prev_T_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_T_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_F_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_F_early_low_rew_avail_(:,counter_event))/sqrt(count.num_session)];
        %     bar(x_label,mean_, 'Facecolor', 'black')
        % add error bars
        [ngroups, nbars] = size(mean_);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, mean_(:,i), sem_(:,i), '.k', 'LineWidth', 2, 'MarkerSize', 30);
        end
        xlim([0.5 6.5])
        xticks([1 : 6])
        xticklabels(x_label)
        xline([1.5 2.5 3.5 4.5 5.5],'k', 'LineWidth',1)
        xlabel('RPE event')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        title(['early low reward'])

        %%% late %%%
        % high rew
        s(counter_event + num_event*3) = subplot(3,2,5);
        hold on
        mean_ = [nanmean(mean_prev_TT_late_high_rew_avail_(:,counter_event))  nanmean(mean_curr_TT_late_high_rew_avail_(:,counter_event)); nanmean(mean_prev_FF_late_high_rew_avail_(:,counter_event))  nanmean(mean_curr_FF_late_high_rew_avail_(:,counter_event)); ...
            nanmean(mean_prev_FT_late_high_rew_avail_(:,counter_event))  nanmean(mean_curr_FT_late_high_rew_avail_(:,counter_event)); nanmean(mean_prev_TF_late_high_rew_avail_(:,counter_event))  nanmean(mean_curr_TF_late_high_rew_avail_(:,counter_event)) ; ...
            nanmean(mean_prev_T_late_high_rew_avail_(:,counter_event))  nanmean(mean_curr_T_late_high_rew_avail_(:,counter_event)); nanmean(mean_prev_F_late_high_rew_avail_(:,counter_event))  nanmean(mean_curr_F_late_high_rew_avail_(:,counter_event))];
        sem_ = [nanstd(mean_prev_TT_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_TT_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_FF_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_FF_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session); ...
            nanstd(mean_prev_FT_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_FT_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_TF_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_TF_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session) ; ...
            nanstd(mean_prev_T_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_T_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_F_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_F_late_high_rew_avail_(:,counter_event))/sqrt(count.num_session)];
        %      bar(x_label,mean_, 'Facecolor', 'black')
        % add error bars
        [ngroups, nbars] = size(mean_);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, mean_(:,i), sem_(:,i), '.k', 'LineWidth', 2, 'MarkerSize', 30);
        end
        xlim([0.5 6.5])
        xticks([1 : 6])
        xticklabels(x_label)
        xline([1.5 2.5 3.5 4.5 5.5],'k', 'LineWidth',1)
        xlabel('RPE event')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        title(['late high reward'])

        % low rew
        s(counter_event + num_event*4) = subplot(3,2,6);
        hold on
        mean_ = [nanmean(mean_prev_TT_late_low_rew_avail_(:,counter_event))  nanmean(mean_curr_TT_late_low_rew_avail_(:,counter_event)); nanmean(mean_prev_FF_late_low_rew_avail_(:,counter_event))  nanmean(mean_curr_FF_late_low_rew_avail_(:,counter_event)); ...
            nanmean(mean_prev_FT_late_low_rew_avail_(:,counter_event))  nanmean(mean_curr_FT_late_low_rew_avail_(:,counter_event)); nanmean(mean_prev_TF_late_low_rew_avail_(:,counter_event))  nanmean(mean_curr_TF_late_low_rew_avail_(:,counter_event)) ; ...
            nanmean(mean_prev_T_late_low_rew_avail_(:,counter_event))  nanmean(mean_curr_T_late_low_rew_avail_(:,counter_event)); nanmean(mean_prev_F_late_low_rew_avail_(:,counter_event))  nanmean(mean_curr_F_late_low_rew_avail_(:,counter_event))];
        sem_ = [nanstd(mean_prev_TT_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_TT_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_FF_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_FF_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session); ...
            nanstd(mean_prev_FT_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_FT_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_TF_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_TF_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session) ; ...
            nanstd(mean_prev_T_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_T_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session); nanstd(mean_prev_F_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session)  nanstd(mean_curr_F_late_low_rew_avail_(:,counter_event))/sqrt(count.num_session)];
        %     bar(x_label,mean_, 'Facecolor', 'black')
        % add error bars
        [ngroups, nbars] = size(mean_);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, mean_(:,i), sem_(:,i), '.k', 'LineWidth', 2, 'MarkerSize', 30);
        end
        xlim([0.5 6.5])
        xticks([1 : 6])
        xticklabels(x_label)
        xline([1.5 2.5 3.5 4.5 5.5],'k', 'LineWidth',1)
        xlabel('RPE event')
        ylabel(variable_list(counter_variable), 'Interpreter', 'none')
        title(['late low reward'])
    end
    linkaxes([s, s2], 'y')
    sgtitle([path.path_data_monkey_sorted(6:end) ' | sessions: ' num2str(count.num_session) ], 'Interpreter', 'none')
    ESN_Beautify_Plot
    eval(['fig' num2str(counter_variable) '.WindowState = ''maximized'';'])
    saveas(gcf, [path.out_path 'RPE' filesep path.path_data_monkey_sorted '_' char(variable_list(counter_variable))], 'pdf');
end
%% Figure - sac work
% variable_list = ["ISR","eye_vigor","eye_vm_max", "eye_dm_max", "reaction", "dist_err_eye_tgt", "det_var_err_eye_tgt"];
variable_list = ["ISR","eye_vigor","eye_vm_max", "eye_dm_max", "reaction", "dist_err_eye_tgt"];

for counter_variable = 1 : length(variable_list)
    clearvars mean_ err_ub err_lb
    x_err = 1:count.max_num_sac_work;

    eval(['fig' num2str(counter_variable) ' = figure;'])

    % average accross sessions
    s1 = subplot(2,3,1);
    hold on
    for counter_session = 1 : count.num_session
        mean_(counter_session,:) = nanmean(sac.(variable_list(counter_variable)){counter_session,1});
    end
    err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
    err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);

    shade(x_err, err_ub, 'k', x_err, err_lb, 'k', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','k', 'FillColor', 'k');
    plot(nanmean(mean_),'k', 'LineWidth', 1);
    xlabel('sac num. in work')
    ylabel(variable_list(counter_variable), 'Interpreter', 'none')
    xlim([0 8])
    xticks([1:7])
    title('all work')

    s2 = subplot(2,3,2);
    hold on
    for counter_session = 1 : count.num_session
        mean_(counter_session,:) = nanmean(sac.(strcat(variable_list(counter_variable), '_early')){counter_session,1});
    end
    err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
    err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);

    shade(x_err, err_ub, 'r', x_err, err_lb, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
    plot(nanmean(mean_),'r', 'LineWidth', 1);

    for counter_session = 1 : count.num_session
        mean_(counter_session,:) = nanmean(sac.(strcat(variable_list(counter_variable), '_late')){counter_session,1});
    end
    err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
    err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);

    shade(x_err, err_ub, 'b', x_err, err_lb, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
    plot(nanmean(mean_),'b', 'LineWidth', 1);
    xlabel('sac num. in work')
    ylabel(variable_list(counter_variable), 'Interpreter', 'none')
    xlim([0 8])
    xticks([1:7])
    title('early (r) vs late (b) work')

    s3 = subplot(2,3,3);
    % high vs low rew history
    hold on
    for counter_session = 1 : count.num_session
        mean_(counter_session,:) = nanmean(sac.(strcat(variable_list(counter_variable), '_high_rew_hist')){counter_session,1});
    end
    err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
    err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
    shade(x_err, err_ub, 'r', x_err, err_lb, 'r', 'LineWidth',  0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
    plot(nanmean(mean_),'r', 'LineWidth', 1);

    for counter_session = 1 : count.num_session
        mean_(counter_session,:) = nanmean(sac.(strcat(variable_list(counter_variable), '_low_rew_hist')){counter_session,1});
    end
    err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
    err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
    shade(x_err, err_ub, 'b', x_err, err_lb, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
    plot(nanmean(mean_),'b', 'LineWidth', 1);
    xlabel('sac num. in harvest')
    ylabel(variable_list(counter_variable), 'Interpreter', 'none')
    xlim([0 8])
    xticks([1:7])
    title('high rew history (r) vs low rew history (b)')

    s4 = subplot(2,3,4);
    % high vs low rew during early
    hold on
    for counter_session = 1 : count.num_session
        mean_(counter_session,:) = nanmean(sac.(strcat(variable_list(counter_variable), '_early_high_rew_hist')){counter_session,1});
    end
    err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
    err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
    shade(x_err, err_ub, 'r', x_err, err_lb, 'r', 'LineWidth',  0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
    plot(nanmean(mean_),'r', 'LineWidth', 1);

    for counter_session = 1 : count.num_session
        mean_(counter_session,:) = nanmean(sac.(strcat(variable_list(counter_variable), '_early_low_rew_hist')){counter_session,1});
    end
    err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
    err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
    shade(x_err, err_ub, 'b', x_err, err_lb, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
    plot(nanmean(mean_),'b', 'LineWidth', 1);
    xlabel('sac num. in harvest')
    ylabel(variable_list(counter_variable), 'Interpreter', 'none')
    xlim([0 8])
    xticks([1:7])
    title('early work: high rew history (r) vs low rew history (b)')

    % high vs low rew during late
    s5 = subplot(2,3,5);
    hold on
    for counter_session = 1 : count.num_session
        mean_(counter_session,:) = nanmean(sac.(strcat(variable_list(counter_variable), '_late_high_rew_hist')){counter_session,1});
    end
    err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
    err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
    shade(x_err, err_ub, 'r', x_err, err_lb, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
    plot(nanmean(mean_),'r', 'LineWidth', 1);

    for counter_session = 1 : count.num_session
        mean_(counter_session,:) = nanmean(sac.(strcat(variable_list(counter_variable), '_late_low_rew_hist')){counter_session,1});
    end
    err_ub = nanmean(mean_) + nanstd(mean_)/sqrt(count.num_session);
    err_lb = nanmean(mean_) - nanstd(mean_)/sqrt(count.num_session);
    shade(x_err, err_ub, 'b', x_err, err_lb, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
    plot(nanmean(mean_),'b', 'LineWidth', 1);
    xlabel('sac num. in harvest')
    ylabel(variable_list(counter_variable), 'Interpreter', 'none')
    xlim([0 8])
    xticks([1:7])
    title('late work: high rew history (r) vs low rew history (b)')

    linkaxes([s1 s2 s3 s4 s5], 'y')
    sgtitle([path.path_data_monkey_sorted(6:end) ' | sessions: ' num2str(count.num_session) ], 'Interpreter', 'none')
    ESN_Beautify_Plot
    eval(['fig' num2str(counter_variable) '.WindowState = ''maximized'';'])
    saveas(gcf, [path.out_path 'WORK' filesep  path.path_data_monkey_sorted '_' char(variable_list(counter_variable))], 'pdf');
end
end

