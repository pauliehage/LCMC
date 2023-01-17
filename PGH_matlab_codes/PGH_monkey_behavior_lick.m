%% Produces LICKS_ALL_DATA and EXPERIMENT_PARAMS
function [LICKS_ALL_DATA, EXPERIMENT_PARAMS] = PGH_monkey_behavior_lick(mat_file_address, flag_figure, params, funcs)
%% Get file_name and file_path

% if there is no inputs, then set pathnames to pwd
if nargin < 1
    [file_name,file_path] = uigetfile(pwd, 'Select DLC mat file');
    flag_figure = 1;
    params = [];
    funcs = [];
end
if nargin == 1
    [file_path,file_name,file_ext] = fileparts(mat_file_address);
    file_name = [file_name file_ext];
    flag_figure = 1;
    params = [];
    funcs = [];
end
if nargin >= 2
    [file_path,file_name,file_ext] = fileparts(mat_file_address);
    file_name = [file_name file_ext];

end
% add filesep ('/' or '\') to the end of file_path
if ~strcmp(file_path(end), filesep)
    file_path = [file_path filesep];
end

% parts = strsplit(file_path, ["/", "\"]);
% name = erase(parts(length(parts) - 2), '-');
% file_name = char(extractBetween(name, 3, 15));

EXPERIMENT_PARAMS.mat_FileName = file_name;
EXPERIMENT_PARAMS.mat_PathName = file_path;
EXPERIMENT_PARAMS.file_name = file_name;
EXPERIMENT_PARAMS.flag_figure = flag_figure;
EXPERIMENT_PARAMS.flag_figure_debug = 0;

%% Interp
DLC.FILE.interp = 1;
%% Load Data
fprintf('Loading: ')
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure mat_file_address params funcs;
filename = EXPERIMENT_PARAMS.mat_FileName;
path_name = EXPERIMENT_PARAMS.mat_PathName;
load([path_name filename], 'data');
DLC.data = data;
DLC.FILE.path_to_analyzed = path_name;

fprintf([EXPERIMENT_PARAMS.file_name ' --> Completed. \n'])

%% Find video sampling rate (FPS)
% check if the FPS file exists in analyzed figs folder
fprintf('Finding video FPS: ')
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
path_name = EXPERIMENT_PARAMS.mat_PathName;
if ~strcmp(path_name(end), filesep);path_name = [path_name filesep];end

datehour = EXPERIMENT_PARAMS.file_name(1:13);

path_to_raw_data = [path_name '..' filesep ...
    '..' filesep '..' filesep 'raw_data' filesep];
path_to_analyzed_figs_tongue = [path_name '..' filesep ...
    '..' filesep '..' filesep 'analyzed_figs' filesep 'behavior_data' filesep 'tongue' filesep];
if isempty(dir(path_to_analyzed_figs_tongue))
    mkdir(path_to_analyzed_figs_tongue);
end

dir_FPS = dir([path_name, '*_video.mat']);
if ~isempty(dir_FPS)
    load([path_name dir_FPS(1).name],'FPS', 'height', 'width', 'duration', 'num_frames')
else
    fprintf('FPS not found, computing now ... \n')
    [LED_FPS, FPS, height, width, duration, num_frames] = PGH_estimate_vid_fps(path_to_raw_data,flag_figure, 30, 150, 1);
    save([path_name datehour '_video.mat'],'FPS','LED_FPS', 'height', 'width', 'duration', 'num_frames');
    saveas(gcf,[path_to_analyzed_figs_tongue  datehour '_FPS'], 'pdf');

end

EXPERIMENT_PARAMS.FPS = FPS;
DLC.FILE.vid_height = height;
DLC.FILE.vid_width = width;
DLC.FILE.duration = duration;
DLC.FILE.num_frames = num_frames;
fprintf([num2str(EXPERIMENT_PARAMS.FPS) ' --> Completed. \n'])

%% Analyze DLC data
fprintf(['Analyzing: ' EXPERIMENT_PARAMS.mat_FileName '\n'])
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
flag_qa = 0;
if flag_qa == 1
    [DLC,EXPERIMENT_PARAMS] = quality_assurance(DLC,EXPERIMENT_PARAMS);
end
[DLC,EXPERIMENT_PARAMS] = analyze(DLC,EXPERIMENT_PARAMS, params, funcs);
[LICKS_ALL_DATA, EXPERIMENT_PARAMS] = build_LICKS_ALL_DATA(DLC, EXPERIMENT_PARAMS, params, funcs);
plot_lick_sorter(LICKS_ALL_DATA, EXPERIMENT_PARAMS, params, funcs);
end

%% MAIN FUNCTIONS %%
%% MAIN Function: Analyze
function [DLC,EXPERIMENT_PARAMS] = analyze(DLC,EXPERIMENT_PARAMS, params, funcs)
[DLC,EXPERIMENT_PARAMS] = extract_scale_shift(DLC,EXPERIMENT_PARAMS, params, funcs);
[DLC,EXPERIMENT_PARAMS] = detect_licks_and_bouts(DLC,EXPERIMENT_PARAMS, params, funcs);
[DLC,EXPERIMENT_PARAMS] = calculate_lick_kinematics(DLC,EXPERIMENT_PARAMS, params, funcs);
[DLC,EXPERIMENT_PARAMS] = geometrization(DLC,EXPERIMENT_PARAMS, params, funcs);
[DLC,EXPERIMENT_PARAMS] = sort_licks(DLC,EXPERIMENT_PARAMS, params, funcs);
[DLC,EXPERIMENT_PARAMS] = quantify_food(DLC,EXPERIMENT_PARAMS, params, funcs);
[DLC,EXPERIMENT_PARAMS] = detect_harvest(DLC,EXPERIMENT_PARAMS, params, funcs);
%PGH_gif_maker(DLC,562,'png');
PGH_plot_sample_trace(DLC);
%PGH_plot_lick_sorter_averages(DLC,'pdf','r')
end

%% MAIN Function: Quality assurance
function [DLC,EXPERIMENT_PARAMS] = quality_assurance(DLC,EXPERIMENT_PARAMS)
fprintf(['Performing QA: ' EXPERIMENT_PARAMS.file_name])
vid_obj = VideoReader(strcat(DLC.FILE.path_to_analyzed, EXPERIMENT_PARAMS.file_name(1:17), '.mp4'));
[DLC,EXPERIMENT_PARAMS] = qa_stationary(DLC,EXPERIMENT_PARAMS, vid_obj);
% Should be set to 1 if the lighting condition is bad
check_lighting_condition = 0;
[DLC,EXPERIMENT_PARAMS] = qa_food_and_tongue(DLC,EXPERIMENT_PARAMS, vid_obj, check_lighting_condition);
save_qa(DLC,EXPERIMENT_PARAMS);
fprintf(' --> Completed. \n')
end

%% MAIN Function: Build LICKS_ALL_DATA
function [LICKS_ALL_DATA, EXPERIMENT_PARAMS] = build_LICKS_ALL_DATA(DLC, EXPERIMENT_PARAMS, params, funcs)
fprintf(['Building LICK_DATA_ALL: ' EXPERIMENT_PARAMS.mat_FileName ' ... ' '\n'])
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
%% Build validity
validity = ones(1, DLC.IND.num_lick);
is_filter_outlier = 0;
if is_filter_outlier == 1
    validity(isoutlier(DLC.KINEMATIC.d_lick_max')) = 0;
    validity(isoutlier(DLC.KINEMATIC.v_lick_max')) = 0;
    validity(isoutlier(DLC.KINEMATIC.v_lick_min')) = 0;
    validity(isoutlier(DLC.TIME.time_lick_duration')) = 0;
end
LICKS_ALL_DATA.validity = validity;

%% Build tag and lick_tag_list
lick_tag_list = params.lick.tag_name_list;
lick_tag_bout_list = params.lick.tag_bout_name_list;
lick_tag_harvest_list = params.lick.tag_harvest_name_list;

is_groom= logical(DLC.CLASS.is_grooming_lick(LICKS_ALL_DATA.validity == 1));
is_inner_tube_success =  logical(DLC.CLASS.is_r_reward_inner_tube_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_reward_inner_tube_lick(LICKS_ALL_DATA.validity == 1));
is_inner_tube_fail = logical(DLC.CLASS.is_r_noreward_inner_tube_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_noreward_inner_tube_lick(LICKS_ALL_DATA.validity == 1));
is_outer_edge_success = logical(DLC.CLASS.is_r_reward_outer_edge_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_reward_outer_edge_lick(LICKS_ALL_DATA.validity == 1));
is_outer_edge_fail = logical(DLC.CLASS.is_r_noreward_outer_edge_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_noreward_outer_edge_lick(LICKS_ALL_DATA.validity == 1));
is_under_tube_success =  logical(DLC.CLASS.is_r_reward_under_tube_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_reward_under_tube_lick(LICKS_ALL_DATA.validity == 1));
is_under_tube_fail = logical(DLC.CLASS.is_r_noreward_under_tube_lick(LICKS_ALL_DATA.validity == 1) + DLC.CLASS.is_l_noreward_under_tube_lick(LICKS_ALL_DATA.validity == 1));
is_bout_start = ismember(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1) ,DLC.IND.ind_lick_onset_str_bout);
is_bout_end = ismember(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1) ,DLC.IND.ind_lick_onset_end_bout);
is_harvest_start = ismember(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1) ,DLC.IND.ind_lick_onset_str_harvest);
is_harvest_end = ismember(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1) ,DLC.IND.ind_lick_onset_end_harvest);

tag = zeros(1,length(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1)));
tag_bout = zeros(1,length(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1)));
tag_harvest = zeros(1,length(DLC.IND.ind_lick_onset(LICKS_ALL_DATA.validity == 1)));

tag(is_groom) = 1;
tag(is_inner_tube_success) = 2;
tag(is_inner_tube_fail) = 3;
tag(is_outer_edge_success) = 4;
tag(is_outer_edge_fail) = 5;
tag(is_under_tube_success) = 6;
tag(is_under_tube_fail) = 7;

LICKS_ALL_DATA.tag = tag;
EXPERIMENT_PARAMS.lick_tag_list = lick_tag_list;

tag_bout(is_bout_start) = 1;
tag_bout(is_bout_end) = 2;

LICKS_ALL_DATA.tag_bout = tag_bout;
EXPERIMENT_PARAMS.lick_tag_bout_list = lick_tag_bout_list;

tag_harvest(is_harvest_start) = 1;
tag_harvest(is_harvest_end) = 2;

LICKS_ALL_DATA.tag_harvest = tag_harvest;
EXPERIMENT_PARAMS.lick_tag_harvest_list = lick_tag_harvest_list;
%% Build time tags
LICKS_ALL_DATA.time_onset = DLC.TIME.time_lick_onset(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.time_vmax = DLC.TIME.time_v_lick_max_abs(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.time_dmax = DLC.TIME.time_d_lick_max_abs(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.time_vmin = DLC.TIME.time_v_lick_min_abs(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.time_offset = DLC.TIME.time_lick_offset(LICKS_ALL_DATA.validity == 1)';

%% Build kinematics
% max amp, vm+/-, ang
LICKS_ALL_DATA.tongue_dm_max = DLC.KINEMATIC.d_lick_max(LICKS_ALL_DATA.validity == 1)'; % mm
LICKS_ALL_DATA.tongue_vm_max = DLC.KINEMATIC.v_lick_max(LICKS_ALL_DATA.validity == 1)'; % mm/s
LICKS_ALL_DATA.tongue_vm_min = DLC.KINEMATIC.v_lick_min(LICKS_ALL_DATA.validity == 1)'; % mm/s
LICKS_ALL_DATA.tongue_ang_max = DLC.KINEMATIC.angle_lick_max(LICKS_ALL_DATA.validity == 1)';  % deg w.r.t face normal vector

% displacement, velocity, and angle of tongue during lick onset to
% offset
LICKS_ALL_DATA.tongue_dm = DLC.KINEMATIC.d_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_vm = DLC.KINEMATIC.v_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_ang = DLC.KINEMATIC.angle_lick(LICKS_ALL_DATA.validity == 1,:)';

% durations: licks, bout, harvest
LICKS_ALL_DATA.duration_lick = DLC.TIME.time_lick_duration(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.duration_bout = DLC.TIME.time_bout_duration';
LICKS_ALL_DATA.duration_harvest = DLC.TIME.time_harvest_duration';

% displacement, velocity, and angle stream data
LICKS_ALL_DATA.tongue_dm_stream = DLC.KINEMATIC.d_tip;
LICKS_ALL_DATA.tongue_vm_stream = DLC.KINEMATIC.v_tip;
LICKS_ALL_DATA.tongue_ang_stream = DLC.KINEMATIC.angle_midtip;

% time stream data
LICKS_ALL_DATA.time_1K_stream = DLC.TIME.time_1K';


%% Build positions (x, y)
% position: tongue
LICKS_ALL_DATA.tongue_tip_px = DLC.KINEMATIC.tip_tongue_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_tip_py = DLC.KINEMATIC.tip_tongue_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_r_px = DLC.KINEMATIC.r_tongue_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_r_py = DLC.KINEMATIC.r_tongue_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_l_px = DLC.KINEMATIC.l_tongue_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_l_py = DLC.KINEMATIC.l_tongue_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_mid_px = DLC.KINEMATIC.mid_tongue_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.tongue_mid_py = DLC.KINEMATIC.mid_tongue_y_lick(LICKS_ALL_DATA.validity == 1,:)';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.tongue_tip_px_onset = DLC.KINEMATIC.tip_tongue_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_py_onset = DLC.KINEMATIC.tip_tongue_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_px_onset = DLC.KINEMATIC.r_tongue_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_py_onset = DLC.KINEMATIC.r_tongue_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_px_onset = DLC.KINEMATIC.l_tongue_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_py_onset = DLC.KINEMATIC.l_tongue_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_px_onset = DLC.KINEMATIC.mid_tongue_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_py_onset = DLC.KINEMATIC.mid_tongue_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_px_vmax = DLC.KINEMATIC.tip_tongue_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_py_vmax = DLC.KINEMATIC.tip_tongue_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_px_vmax = DLC.KINEMATIC.r_tongue_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_py_vmax = DLC.KINEMATIC.r_tongue_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_px_vmax = DLC.KINEMATIC.l_tongue_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_py_vmax = DLC.KINEMATIC.l_tongue_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_px_vmax = DLC.KINEMATIC.mid_tongue_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_py_vmax = DLC.KINEMATIC.mid_tongue_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_px_dmax = DLC.KINEMATIC.tip_tongue_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_py_dmax = DLC.KINEMATIC.tip_tongue_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_px_dmax = DLC.KINEMATIC.r_tongue_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_py_dmax = DLC.KINEMATIC.r_tongue_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_px_dmax = DLC.KINEMATIC.l_tongue_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_py_dmax = DLC.KINEMATIC.l_tongue_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_px_dmax = DLC.KINEMATIC.mid_tongue_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_py_dmax = DLC.KINEMATIC.mid_tongue_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_px_vmin = DLC.KINEMATIC.tip_tongue_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_py_vmin = DLC.KINEMATIC.tip_tongue_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_px_vmin = DLC.KINEMATIC.r_tongue_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_py_vmin = DLC.KINEMATIC.r_tongue_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_px_vmin = DLC.KINEMATIC.l_tongue_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_py_vmin = DLC.KINEMATIC.l_tongue_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_px_vmin = DLC.KINEMATIC.mid_tongue_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_py_vmin = DLC.KINEMATIC.mid_tongue_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_px_offset = DLC.KINEMATIC.tip_tongue_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_tip_py_offset = DLC.KINEMATIC.tip_tongue_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_px_offset = DLC.KINEMATIC.r_tongue_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_r_py_offset = DLC.KINEMATIC.r_tongue_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_px_offset = DLC.KINEMATIC.l_tongue_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_l_py_offset = DLC.KINEMATIC.l_tongue_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_px_offset = DLC.KINEMATIC.mid_tongue_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.tongue_mid_py_offset = DLC.KINEMATIC.mid_tongue_y_offset(LICKS_ALL_DATA.validity == 1);

% position: nose
LICKS_ALL_DATA.nose_r_px = DLC.KINEMATIC.r_nose_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.nose_r_py = DLC.KINEMATIC.r_nose_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.nose_l_px = DLC.KINEMATIC.l_nose_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.nose_l_py = DLC.KINEMATIC.l_nose_y_lick(LICKS_ALL_DATA.validity == 1,:)';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.nose_r_px_onset = DLC.KINEMATIC.r_nose_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_py_onset = DLC.KINEMATIC.r_nose_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_px_onset = DLC.KINEMATIC.l_nose_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_py_onset = DLC.KINEMATIC.l_nose_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_px_vmax = DLC.KINEMATIC.r_nose_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_py_vmax = DLC.KINEMATIC.r_nose_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_px_vmax = DLC.KINEMATIC.l_nose_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_py_vmax = DLC.KINEMATIC.l_nose_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_px_dmax = DLC.KINEMATIC.r_nose_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_py_dmax = DLC.KINEMATIC.r_nose_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_px_dmax = DLC.KINEMATIC.l_nose_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_py_dmax = DLC.KINEMATIC.l_nose_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_px_vmin = DLC.KINEMATIC.r_nose_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_py_vmin = DLC.KINEMATIC.r_nose_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_px_vmin = DLC.KINEMATIC.l_nose_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_py_vmin = DLC.KINEMATIC.l_nose_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_px_offset = DLC.KINEMATIC.r_nose_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_r_py_offset = DLC.KINEMATIC.r_nose_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_px_offset = DLC.KINEMATIC.l_nose_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.nose_l_py_offset = DLC.KINEMATIC.l_nose_y_offset(LICKS_ALL_DATA.validity == 1);

% position: reward
LICKS_ALL_DATA.rew_r_px = DLC.KINEMATIC.r_food_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rew_r_py = DLC.KINEMATIC.r_food_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rew_l_px = DLC.KINEMATIC.l_food_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rew_l_py = DLC.KINEMATIC.l_food_y_lick(LICKS_ALL_DATA.validity == 1,:)';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.rew_r_px_onset = DLC.KINEMATIC.r_food_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_py_onset = DLC.KINEMATIC.r_food_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_px_onset = DLC.KINEMATIC.l_food_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_py_onset = DLC.KINEMATIC.l_food_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_px_vmax = DLC.KINEMATIC.r_food_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_py_vmax = DLC.KINEMATIC.r_food_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_px_vmax = DLC.KINEMATIC.l_food_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_py_vmax = DLC.KINEMATIC.l_food_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_px_dmax = DLC.KINEMATIC.r_food_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_py_dmax = DLC.KINEMATIC.r_food_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_px_dmax = DLC.KINEMATIC.l_food_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_py_dmax = DLC.KINEMATIC.l_food_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_px_vmin = DLC.KINEMATIC.r_food_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_py_vmin = DLC.KINEMATIC.r_food_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_px_vmin = DLC.KINEMATIC.l_food_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_py_vmin = DLC.KINEMATIC.l_food_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_px_offset = DLC.KINEMATIC.r_food_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_r_py_offset = DLC.KINEMATIC.r_food_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_px_offset = DLC.KINEMATIC.l_food_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rew_l_py_offset = DLC.KINEMATIC.l_food_y_offset(LICKS_ALL_DATA.validity == 1);

% position: tubes
LICKS_ALL_DATA.rtube_r_px = DLC.KINEMATIC.r_tube_r_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rtube_r_py = DLC.KINEMATIC.r_tube_r_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rtube_l_px = DLC.KINEMATIC.r_tube_l_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.rtube_l_py = DLC.KINEMATIC.r_tube_l_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.ltube_r_px = DLC.KINEMATIC.l_tube_r_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.ltube_r_py = DLC.KINEMATIC.l_tube_r_y_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.ltube_l_px = DLC.KINEMATIC.l_tube_l_x_lick(LICKS_ALL_DATA.validity == 1,:)';
LICKS_ALL_DATA.ltube_l_py = DLC.KINEMATIC.l_tube_l_y_lick(LICKS_ALL_DATA.validity == 1,:)';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.rtube_r_px_onset = DLC.KINEMATIC.r_tube_r_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_py_onset = DLC.KINEMATIC.r_tube_r_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_px_onset = DLC.KINEMATIC.r_tube_l_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_py_onset = DLC.KINEMATIC.r_tube_l_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_px_vmax = DLC.KINEMATIC.r_tube_r_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_py_vmax = DLC.KINEMATIC.r_tube_r_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_px_vmax = DLC.KINEMATIC.r_tube_l_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_py_vmax = DLC.KINEMATIC.r_tube_l_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_px_dmax = DLC.KINEMATIC.r_tube_r_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_py_dmax = DLC.KINEMATIC.r_tube_r_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_px_dmax = DLC.KINEMATIC.r_tube_l_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_py_dmax = DLC.KINEMATIC.r_tube_l_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_px_vmin = DLC.KINEMATIC.r_tube_r_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_py_vmin = DLC.KINEMATIC.r_tube_r_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_px_vmin = DLC.KINEMATIC.r_tube_l_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_py_vmin = DLC.KINEMATIC.r_tube_l_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_px_offset = DLC.KINEMATIC.r_tube_r_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_r_py_offset = DLC.KINEMATIC.r_tube_r_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_px_offset = DLC.KINEMATIC.r_tube_l_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.rtube_l_py_offset = DLC.KINEMATIC.r_tube_l_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_px_onset = DLC.KINEMATIC.l_tube_r_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_py_onset = DLC.KINEMATIC.l_tube_r_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_px_onset = DLC.KINEMATIC.l_tube_l_x_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_py_onset = DLC.KINEMATIC.l_tube_l_y_onset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_px_vmax = DLC.KINEMATIC.l_tube_r_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_py_vmax = DLC.KINEMATIC.l_tube_r_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_px_vmax = DLC.KINEMATIC.l_tube_l_x_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_py_vmax = DLC.KINEMATIC.l_tube_l_y_vmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_px_dmax = DLC.KINEMATIC.l_tube_r_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_py_dmax = DLC.KINEMATIC.l_tube_r_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_px_dmax = DLC.KINEMATIC.l_tube_l_x_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_py_dmax = DLC.KINEMATIC.l_tube_l_y_dmax(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_px_vmin = DLC.KINEMATIC.l_tube_r_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_py_vmin = DLC.KINEMATIC.l_tube_r_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_px_vmin = DLC.KINEMATIC.l_tube_l_x_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_py_vmin = DLC.KINEMATIC.l_tube_l_y_vmin(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_px_offset = DLC.KINEMATIC.l_tube_r_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_r_py_offset = DLC.KINEMATIC.l_tube_r_y_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_px_offset = DLC.KINEMATIC.l_tube_l_x_offset(LICKS_ALL_DATA.validity == 1);
LICKS_ALL_DATA.ltube_l_py_offset = DLC.KINEMATIC.l_tube_l_y_offset(LICKS_ALL_DATA.validity == 1);

%% Build reward-tube capacity
% rew capacity
LICKS_ALL_DATA.rew_capacity_r_lick_onset = DLC.FOOD.r_tube_food_lick_onset(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.rew_capacity_r_lick_offset = DLC.FOOD.r_tube_food_lick_offset(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.rew_capacity_l_lick_onset = DLC.FOOD.l_tube_food_lick_onset(LICKS_ALL_DATA.validity == 1)';
LICKS_ALL_DATA.rew_capacity_l_lick_offset = DLC.FOOD.l_tube_food_lick_offset(LICKS_ALL_DATA.validity == 1)';

LICKS_ALL_DATA.rew_capacity_r_bout_start = DLC.FOOD.r_tube_food_bout_start';
LICKS_ALL_DATA.rew_capacity_r_bout_end = DLC.FOOD.r_tube_food_bout_end';
LICKS_ALL_DATA.rew_capacity_l_bout_start = DLC.FOOD.l_tube_food_bout_start';
LICKS_ALL_DATA.rew_capacity_l_bout_end = DLC.FOOD.l_tube_food_bout_end';
end

%% MAIN Function: Plot lick sorter summary
function plot_lick_sorter(LICKS_ALL_DATA, EXPERIMENT_PARAMS, params, funcs)
fprintf(['Ploting: ' EXPERIMENT_PARAMS.mat_FileName '\n'])
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
%% Plot LICKS_ALL_DATA
if ~EXPERIMENT_PARAMS.flag_figure
    return;
end
params.cell_name = EXPERIMENT_PARAMS.file_name(1:13);
params.duration = EXPERIMENT_PARAMS.duration_video;
PGH_plot_lick_sorter(LICKS_ALL_DATA, params)
fprintf(' --> Completed. \n')

end

%% UTILITY FUNCTIONS %%
%% Function: Extract, scale, and shift DLC data
function [DLC, EXPERIMENT_PARAMS] = extract_scale_shift(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Extracting, scaling, and shifting DLC data ...')
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;

% specify scale for pixel to mm conversion
webcam_scale = 7/70;
blackfly_scale = 7/47;
scale = [webcam_scale blackfly_scale];
EXPERIMENT_PARAMS.scale = scale;
FPS = EXPERIMENT_PARAMS.FPS;

% load, interpolate (1K), scale, and rotate DLC data from table into variables
data_FPS = table2array(DLC.data);
time_vid = ((1/FPS) : (1/FPS) : size(data_FPS,1)/FPS)';

% interpolation from FPS to 1K
if DLC.FILE.interp == 1
    time_1K = (time_vid(1) : 0.001 : time_vid(end))';
    for counter_fields = 1 : size(data_FPS,2)
        data_1K(:,counter_fields) = interp1(time_vid,data_FPS(:,counter_fields),time_1K);
    end
else
    time_1K = time_vid; % for making gifs
    data_1K = data_FPS;
end

tip_tongue_x = scale(2)*data_1K(:,1);
tip_tongue_y = scale(2)*data_1K(:,2);
r_tongue_x = scale(2)*data_1K(:,3);
r_tongue_y = scale(2)*data_1K(:,4);
l_tongue_x = scale(2)*data_1K(:,5);
l_tongue_y = scale(2)*data_1K(:,6);
mid_tongue_x = scale(2)*data_1K(:,7);
mid_tongue_y = scale(2)*data_1K(:,8);
r_nose_x = scale(2)*data_1K(:,9);
r_nose_y = scale(2)*data_1K(:,10);
l_nose_x = scale(2)*data_1K(:,11);
l_nose_y = scale(2)*data_1K(:,12);
r_food_x = scale(2)*data_1K(:,13);
r_food_y = scale(2)*data_1K(:,14);
l_food_x = scale(2)*data_1K(:,15);
l_food_y = scale(2)*data_1K(:,16);
r_tube_r_x = scale(2)*data_1K(:,17);
r_tube_r_y = scale(2)*data_1K(:,18);
r_tube_l_x = scale(2)*data_1K(:,19);
r_tube_l_y = scale(2)*data_1K(:,20);
l_tube_r_x = scale(2)*data_1K(:,21);
l_tube_r_y = scale(2)*data_1K(:,22);
l_tube_l_x = scale(2)*data_1K(:,23);
l_tube_l_y = scale(2)*data_1K(:,24);

% find origin and shift DLC data
% clustering approach
% [~, cent_cluster_x] = kmeans(tip_tongue_x, 3);
% [~, cent_cluster_y] = kmeans(tip_tongue_y, 3);
% [x0, ind_min] = min(cent_cluster_x);
% if x0 < (mean(r_nose_x) + mean(l_nose_x))/2
%     [~, ind_max] = max(cent_cluster_x);
%     ind_mid = 1:3;
%     ind_mid(ind_max) = [];
%     ind_mid(ind_min) = [];
%     x0 = cent_cluster_x(ind_mid);
% end
%
% midpoint = (mean(r_nose_y) + mean(l_nose_y)) / 2;
% if abs(cent_cluster_y(1) - midpoint) < abs(cent_cluster_y(2) - midpoint) && abs(cent_cluster_y(1) - midpoint) < abs(cent_cluster_y(3) - midpoint)
%     y0 = cent_cluster_y(1);
% elseif abs(cent_cluster_y(2) - midpoint) < abs(cent_cluster_y(1) - midpoint) && abs(cent_cluster_y(2) - midpoint) < abs(cent_cluster_y(3) - midpoint)
%     y0 = cent_cluster_y(2);
% else
%     y0 = cent_cluster_y(3);
% end

% clustering approach
% nose marker approach
y0 = (mean(r_nose_y) + mean(l_nose_y)) / 2;
x0 = (mean(r_nose_x) + mean(l_nose_x)) / 2 + 5;

EXPERIMENT_PARAMS.duration_video = time_1K(end) - time_1K(1) ;
DLC.POINTS.tip_tongue_x = tip_tongue_x - x0;
DLC.POINTS.tip_tongue_y = tip_tongue_y - y0;
DLC.POINTS.r_tongue_x = r_tongue_x - x0;
DLC.POINTS.r_tongue_y = r_tongue_y - y0;
DLC.POINTS.l_tongue_x = l_tongue_x - x0;
DLC.POINTS.l_tongue_y = l_tongue_y - y0;
DLC.POINTS.mid_tongue_x = mid_tongue_x - x0;
DLC.POINTS.mid_tongue_y = mid_tongue_y - y0;
DLC.POINTS.r_nose_x = r_nose_x - x0;
DLC.POINTS.r_nose_y = r_nose_y - y0;
DLC.POINTS.l_nose_x = l_nose_x - x0;
DLC.POINTS.l_nose_y = l_nose_y - y0;
DLC.POINTS.l_food_x = l_food_x - x0;
DLC.POINTS.l_food_y = l_food_y - y0;
DLC.POINTS.l_tube_r_x = l_tube_r_x - x0;
DLC.POINTS.l_tube_r_y = l_tube_r_y - y0;
DLC.POINTS.l_tube_l_x = l_tube_l_x - x0;
DLC.POINTS.l_tube_l_y = l_tube_l_y - y0;
DLC.POINTS.r_food_x = r_food_x - x0;
DLC.POINTS.r_food_y = r_food_y - y0;
DLC.POINTS.r_tube_r_x = r_tube_r_x - x0;
DLC.POINTS.r_tube_r_y = r_tube_r_y - y0;
DLC.POINTS.r_tube_l_x = r_tube_l_x - x0;
DLC.POINTS.r_tube_l_y = r_tube_l_y - y0;
DLC.POINTS.x0 = x0;
DLC.POINTS.y0 = y0;
DLC.TIME.time_1K = time_1K;

fprintf(' --> Completed. \n')
end

%% Function: Detect licks and bouts
function [DLC, EXPERIMENT_PARAMS] = detect_licks_and_bouts(DLC,EXPERIMENT_PARAMS, params, funcs)
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
fprintf('Detecting licks and bouts ... ');
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
tip_tongue_x = DLC.POINTS.tip_tongue_x;
tip_tongue_y = DLC.POINTS.tip_tongue_y;

% set windows
if DLC.FILE.interp == 1
    smooth_window = 20;
    detection_window = 100;
    window_time = 500;
    bout_threshold = 1000;
else
    smooth_window = 2;
    detection_window = 10;
    window_time = 50;
    bout_threshold = 100;
end
prom_thresh = 0.85;

% compute and smooth displacement trace d_tip
d_tip = movmean(sqrt(tip_tongue_x.^2 + tip_tongue_y.^2),smooth_window)';

% compute negative signal and subtract within bout mins to take licks to 0
% while tongue in mouth
[peak_neg,ind_peak_neg,~,~]= findpeaks(-d_tip, 'MinPeakProminence', 1, 'MinPeakDistance',detection_window,'Annotate', 'extent');
is_del = isoutlier(peak_neg);
peak_neg(is_del) = [];
ind_peak_neg(is_del) = [];
d_tip_neg_interp = interp1(ind_peak_neg,peak_neg,1: length(d_tip));
d_tip = d_tip + d_tip_neg_interp;
d_tip(d_tip < 0) = 0;
d_tip(isnan(d_tip)) = 0;
if DLC.FILE.interp == 1
    v_tip = [0 diff(d_tip)]*1000;
else
    v_tip = [0 diff(d_tip)]*100;
end

% detect onset and offset of licks
[peak,ind_peak,~,prom]= findpeaks(d_tip, 'MinPeakProminence', 3, 'MinPeakDistance',detection_window,'Annotate', 'extent');
ind_lick_onset = nan(length(peak),1);
ind_lick_offset = nan(length(peak),1);
for counter_peak = 1 : length(peak)
    if ind_peak(counter_peak)-window_time < 1
        window_time = window_time + (ind_peak(counter_peak)-window_time) - 1;
    elseif ind_peak(counter_peak)+window_time > length(d_tip)
        window_time = window_time - (ind_peak(counter_peak) + window_time - length(d_tip));
    else
        if DLC.FILE.interp == 1
            window_time = 500;
        else
            window_time = 50;
        end
    end
    ind_lick_onset_ = find(d_tip(ind_peak(counter_peak)-window_time : ind_peak(counter_peak)) <= peak(counter_peak) - prom(counter_peak)*prom_thresh,1, 'last' ) + ind_peak(counter_peak)-window_time - 1;
    ind_lick_offset_ =  find(d_tip(ind_peak(counter_peak) : ind_peak(counter_peak)+window_time)  <= peak(counter_peak) - prom(counter_peak)*prom_thresh,1, 'first' ) + ind_peak(counter_peak);

    %recompute onset/offset using v_tip and values from d_tip
    ind_lick_onset_ = (find(v_tip(ind_lick_onset_-detection_window/2 : ind_lick_onset_+detection_window/2 ) <= 30,1, 'last' ) + ind_lick_onset_ - detection_window/2-1)-1;
    ind_lick_offset_ = (find(v_tip(ind_lick_offset_-detection_window/2 : ind_lick_offset_+detection_window/2 ) <= -30,1, 'last' ) + ind_lick_offset_ - detection_window/2-1) + 1;

    if isempty(ind_lick_onset_) || isempty(ind_lick_offset_)
        ind_lick_onset_ = nan;
        ind_lick_offset_ = nan;
    elseif d_tip(ind_lick_onset_)<prom(counter_peak) && d_tip(ind_lick_onset_) > 3
        ind_lick_onset_ = find(v_tip(ind_lick_onset_-detection_window : ind_lick_onset_) >= 30,1, 'first' ) + ind_lick_onset_ - detection_window-1;
    elseif d_tip(ind_lick_onset_)>prom(counter_peak) || d_tip(ind_lick_offset_)>prom(counter_peak)
        ind_lick_onset_ = nan;
        ind_lick_offset_ = nan;
    end

    if isempty(ind_lick_onset_) || isempty(ind_lick_offset_)
        ind_lick_onset_ = nan;
        ind_lick_offset_ = nan;
    end
    ind_lick_onset(counter_peak,1) = ind_lick_onset_;
    ind_lick_offset(counter_peak,1) =  ind_lick_offset_;
end
is_nan = isnan(ind_lick_onset) | isnan(ind_lick_offset) | (ind_lick_offset - ind_lick_onset)<=0;
ind_lick_onset(is_nan) = [];
ind_lick_offset(is_nan) = [];

is_outlier = isoutlier(d_tip(ind_lick_onset)) | isoutlier(d_tip(ind_lick_offset));
ind_lick_onset(is_outlier) = [];
ind_lick_offset(is_outlier) = [];

time_1K = DLC.TIME.time_1K;
time_lick_onset = time_1K(ind_lick_onset);
time_lick_offset =  time_1K(ind_lick_offset);
time_lick_duration = time_lick_offset - time_lick_onset;

validity = time_lick_duration <= 500;
ind_lick_onset(~validity) = [];
ind_lick_offset(~validity) = [];
time_lick_onset(~validity) = [];
time_lick_offset(~validity) = [];
time_lick_duration(~validity) = [];

num_lick = length(ind_lick_onset);

% 1s threshold for lick bouts
ind_lick_onset_str_bout_ = [1; find(diff(ind_lick_onset) >bout_threshold) + 1];
ind_lick_onset_str_bout = ind_lick_onset(ind_lick_onset_str_bout_);
ind_lick_onset_end_bout_ = [find(diff(ind_lick_onset) >bout_threshold); length(ind_lick_onset)];
ind_lick_onset_end_bout = ind_lick_onset(ind_lick_onset_end_bout_);

num_lick_bout = (ind_lick_onset_end_bout_ - ind_lick_onset_str_bout_) + 1;

validity = num_lick_bout>=3;
num_lick_bout(~validity) = [];
ind_lick_onset_str_bout(~validity) = [];
ind_lick_onset_end_bout(~validity) = [];

time_lick_onset_str_bout = time_1K(ind_lick_onset_str_bout);
time_lick_onset_end_bout = time_1K(ind_lick_onset_end_bout);
time_bout_duration = time_lick_onset_end_bout - time_lick_onset_str_bout;
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

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure
    hold on;
    plot(time_1K,d_tip,'k');
    plot(time_1K(ind_lick_onset_str_bout), d_tip(ind_lick_onset_str_bout), 'or', 'MarkerSize',10 );
    plot(time_1K(ind_lick_onset_end_bout), d_tip(ind_lick_onset_end_bout), 'ob', 'MarkerSize',10 );
    plot(time_1K(ind_lick_onset), d_tip(ind_lick_onset), '.r', 'MarkerSize',10 );
    plot(time_1K(ind_lick_offset), d_tip(ind_lick_onset), '.b', 'MarkerSize',10 );

    yyaxis right;
    hold on
    plot(time_1K,v_tip,'c');
    plot(time_1K(ind_lick_onset_str_bout), v_tip(ind_lick_onset_str_bout), 'or', 'MarkerSize',10 );
    plot(time_1K(ind_lick_onset_end_bout), v_tip(ind_lick_onset_end_bout), 'ob', 'MarkerSize',10 );
    plot(time_1K(ind_lick_onset), v_tip(ind_lick_onset), '.r', 'MarkerSize',10 );
    plot(time_1K(ind_lick_offset), v_tip(ind_lick_offset), '.b', 'MarkerSize',10 );

    xlabel('Time (s)');
    ylabel('Distance (mm)');
    title([EXPERIMENT_PARAMS.file_name ': ' num2str(num_lick) ' licks | ' num2str(num_bout) ' bouts'], 'interpreter', 'none')
    ESN_Beautify_Plot(gcf, [20 10])
end
fprintf(' --> Completed. \n')
end

%% Function: Calculate lick kinematics: d, v, a, angle, ILI(inter lick interval), ILR(instantaneous lick rate = 1/ILI)
function [DLC, EXPERIMENT_PARAMS] = calculate_lick_kinematics(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Calculating lick kinematics ...');
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
params.lick.length_trace = 500;
length_trace = params.lick.length_trace;
d_tip = DLC.KINEMATIC.d_tip;
time_lick_duration = DLC.TIME.time_lick_duration;

num_lick = DLC.IND.num_lick;
num_bout = DLC.IND.num_bout;
time_1K = DLC.TIME.time_1K;
ind_lick_onset = DLC.IND.ind_lick_onset;
ind_lick_offset = DLC.IND.ind_lick_offset;
ind_lick_onset_str_bout = DLC.IND.ind_lick_onset_str_bout;
ind_lick_onset_end_bout = DLC.IND.ind_lick_onset_end_bout;
tip_tongue_x = DLC.POINTS.tip_tongue_x;
tip_tongue_y = DLC.POINTS.tip_tongue_y;
r_tongue_x = DLC.POINTS.r_tongue_x;
r_tongue_y = DLC.POINTS.r_tongue_y;
l_tongue_x = DLC.POINTS.l_tongue_x;
l_tongue_y = DLC.POINTS.l_tongue_y;
mid_tongue_x = DLC.POINTS.mid_tongue_x;
mid_tongue_y = DLC.POINTS.mid_tongue_y;
l_tube_r_x = DLC.POINTS.l_tube_r_x;
l_tube_r_y = DLC.POINTS.l_tube_r_y;
l_tube_l_x = DLC.POINTS.l_tube_l_x;
l_tube_l_y = DLC.POINTS.l_tube_l_y;
l_food_x = DLC.POINTS.l_food_x;
l_food_y = DLC.POINTS.l_food_y;
r_tube_r_x = DLC.POINTS.r_tube_r_x;
r_tube_r_y = DLC.POINTS.r_tube_r_y;
r_tube_l_x = DLC.POINTS.r_tube_l_x;
r_tube_l_y = DLC.POINTS.r_tube_l_y;
r_food_x = DLC.POINTS.r_food_x;
r_food_y = DLC.POINTS.r_food_y;
r_nose_x = DLC.POINTS.r_nose_x;
r_nose_y = DLC.POINTS.r_nose_y;
l_nose_x = DLC.POINTS.l_nose_x;
l_nose_y = DLC.POINTS.l_nose_y;
midtip_y = tip_tongue_y - mid_tongue_y;
midtip_x = tip_tongue_x - mid_tongue_x;
time_lick_onset = DLC.TIME.time_lick_onset;

% build kinematicks _lick
for counter_lick =  1 : 1 : num_lick
    inds_ = (ind_lick_onset(counter_lick) ) : ind_lick_offset(counter_lick);
    if length(inds_) > 500
        inds_ = inds_(1:500);
        pad = [];
    elseif length(inds_) < 500
        % inds_ = padarray(inds_,[0 500-length(inds_)], inds_(end), 'post');
        pad = zeros(1,500-length(inds_));
    end

    % d_lick
    d_lick(counter_lick, :) = [d_tip(inds_) pad];
    [d_lick_max(counter_lick,1), ind_d_lick_max_local] = max(d_lick(counter_lick, :));
    ind_d_lick_max(counter_lick,1) = ind_lick_onset(counter_lick) + ind_d_lick_max_local - 1;

    % v_lick
    v_lick(counter_lick, :) = [0 (diff(d_lick(counter_lick,:))./ (time_1K(2) - time_1K(1)))];
    [v_lick_max(counter_lick,1), ind_v_lick_max_local] = max(v_lick(counter_lick,:));
    ind_v_lick_max(counter_lick,1) = ind_lick_onset(counter_lick) + ind_v_lick_max_local - 1;

    [v_lick_min(counter_lick,1), ind_v_lick_min_local] = min(v_lick(counter_lick,:));
    ind_v_lick_min(counter_lick,1) = ind_lick_onset(counter_lick) + ind_v_lick_min_local - 1;

    % angle_lick
    angle_lick(counter_lick, :) = [rad2deg(atan2(midtip_y(inds_),midtip_x(inds_)))' pad];

    % points
    ind_onset_(counter_lick) = ind_lick_onset(counter_lick);
    ind_vmax_(counter_lick) = ind_v_lick_max(counter_lick);
    ind_dmax_(counter_lick) = ind_d_lick_max(counter_lick);
    ind_vmin_(counter_lick) = ind_v_lick_min(counter_lick);
    ind_offset_(counter_lick) = ind_lick_offset(counter_lick);

    tip_tongue_x_lick(counter_lick, :) = [tip_tongue_x(inds_)' pad];
    tip_tongue_y_lick(counter_lick, :) = [tip_tongue_y(inds_)' pad];
    tip_tongue_x_onset(counter_lick) = tip_tongue_x(ind_onset_(counter_lick));
    tip_tongue_y_onset(counter_lick) = tip_tongue_y(ind_onset_(counter_lick));
    tip_tongue_x_vmax(counter_lick) = tip_tongue_x(ind_vmax_(counter_lick));
    tip_tongue_y_vmax(counter_lick) = tip_tongue_y(ind_vmax_(counter_lick));
    tip_tongue_x_dmax(counter_lick) = tip_tongue_x(ind_dmax_(counter_lick));
    tip_tongue_y_dmax(counter_lick) = tip_tongue_y(ind_dmax_(counter_lick));
    tip_tongue_x_vmin(counter_lick) = tip_tongue_x(ind_vmin_(counter_lick));
    tip_tongue_y_vmin(counter_lick) = tip_tongue_y(ind_vmin_(counter_lick));
    tip_tongue_x_offset(counter_lick) = tip_tongue_x(ind_offset_(counter_lick));
    tip_tongue_y_offset(counter_lick) = tip_tongue_y(ind_offset_(counter_lick));

    r_tongue_x_lick(counter_lick, :) = [r_tongue_x(inds_)' pad];
    r_tongue_y_lick(counter_lick, :) = [r_tongue_y(inds_)' pad];
    r_tongue_x_onset(counter_lick) = r_tongue_x(ind_onset_(counter_lick));
    r_tongue_y_onset(counter_lick) = r_tongue_y(ind_onset_(counter_lick));
    r_tongue_x_vmax(counter_lick) = r_tongue_x(ind_vmax_(counter_lick));
    r_tongue_y_vmax(counter_lick) = r_tongue_y(ind_vmax_(counter_lick));
    r_tongue_x_dmax(counter_lick) = r_tongue_x(ind_dmax_(counter_lick));
    r_tongue_y_dmax(counter_lick) = r_tongue_y(ind_dmax_(counter_lick));
    r_tongue_x_vmin(counter_lick) = r_tongue_x(ind_vmin_(counter_lick));
    r_tongue_y_vmin(counter_lick) = r_tongue_y(ind_vmin_(counter_lick));
    r_tongue_x_offset(counter_lick) = r_tongue_x(ind_offset_(counter_lick));
    r_tongue_y_offset(counter_lick) = r_tongue_y(ind_offset_(counter_lick));

    l_tongue_x_lick(counter_lick, :) = [l_tongue_x(inds_)' pad];
    l_tongue_y_lick(counter_lick, :) = [l_tongue_y(inds_)' pad];
    l_tongue_x_onset(counter_lick) = l_tongue_x(ind_onset_(counter_lick));
    l_tongue_y_onset(counter_lick) = l_tongue_y(ind_onset_(counter_lick));
    l_tongue_x_vmax(counter_lick) = l_tongue_x(ind_vmax_(counter_lick));
    l_tongue_y_vmax(counter_lick) = l_tongue_y(ind_vmax_(counter_lick));
    l_tongue_x_dmax(counter_lick) = l_tongue_x(ind_dmax_(counter_lick));
    l_tongue_y_dmax(counter_lick) = l_tongue_y(ind_dmax_(counter_lick));
    l_tongue_x_vmin(counter_lick) = l_tongue_x(ind_vmin_(counter_lick));
    l_tongue_y_vmin(counter_lick) = l_tongue_y(ind_vmin_(counter_lick));
    l_tongue_x_offset(counter_lick) = l_tongue_x(ind_offset_(counter_lick));
    l_tongue_y_offset(counter_lick) = l_tongue_y(ind_offset_(counter_lick));

    mid_tongue_x_lick(counter_lick, :) = [mid_tongue_x(inds_)' pad];
    mid_tongue_y_lick(counter_lick, :) = [mid_tongue_y(inds_)' pad];
    mid_tongue_x_onset(counter_lick) = mid_tongue_x(ind_onset_(counter_lick));
    mid_tongue_y_onset(counter_lick) = mid_tongue_y(ind_onset_(counter_lick));
    mid_tongue_x_vmax(counter_lick) = mid_tongue_x(ind_vmax_(counter_lick));
    mid_tongue_y_vmax(counter_lick) = mid_tongue_y(ind_vmax_(counter_lick));
    mid_tongue_x_dmax(counter_lick) = mid_tongue_x(ind_dmax_(counter_lick));
    mid_tongue_y_dmax(counter_lick) = mid_tongue_y(ind_dmax_(counter_lick));
    mid_tongue_x_vmin(counter_lick) = mid_tongue_x(ind_vmin_(counter_lick));
    mid_tongue_y_vmin(counter_lick) = mid_tongue_y(ind_vmin_(counter_lick));
    mid_tongue_x_offset(counter_lick) = mid_tongue_x(ind_offset_(counter_lick));
    mid_tongue_y_offset(counter_lick) = mid_tongue_y(ind_offset_(counter_lick));

    r_food_x_lick(counter_lick, :) = [r_food_x(inds_)' pad];
    r_food_y_lick(counter_lick, :) = [r_food_y(inds_)' pad];
    r_food_x_onset(counter_lick) = r_food_x(ind_onset_(counter_lick));
    r_food_y_onset(counter_lick) = r_food_y(ind_onset_(counter_lick));
    r_food_x_vmax(counter_lick) = r_food_x(ind_vmax_(counter_lick));
    r_food_y_vmax(counter_lick) = r_food_y(ind_vmax_(counter_lick));
    r_food_x_dmax(counter_lick) = r_food_x(ind_dmax_(counter_lick));
    r_food_y_dmax(counter_lick) = r_food_y(ind_dmax_(counter_lick));
    r_food_x_vmin(counter_lick) = r_food_x(ind_vmin_(counter_lick));
    r_food_y_vmin(counter_lick) = r_food_y(ind_vmin_(counter_lick));
    r_food_x_offset(counter_lick) = r_food_x(ind_offset_(counter_lick));
    r_food_y_offset(counter_lick) = r_food_y(ind_offset_(counter_lick));

    l_food_x_lick(counter_lick, :) = [l_food_x(inds_)' pad];
    l_food_y_lick(counter_lick, :) = [l_food_y(inds_)' pad];
    l_food_x_onset(counter_lick) = l_food_x(ind_onset_(counter_lick));
    l_food_y_onset(counter_lick) = l_food_y(ind_onset_(counter_lick));
    l_food_x_vmax(counter_lick) = l_food_x(ind_vmax_(counter_lick));
    l_food_y_vmax(counter_lick) = l_food_y(ind_vmax_(counter_lick));
    l_food_x_dmax(counter_lick) = l_food_x(ind_dmax_(counter_lick));
    l_food_y_dmax(counter_lick) = l_food_y(ind_dmax_(counter_lick));
    l_food_x_vmin(counter_lick) = l_food_x(ind_vmin_(counter_lick));
    l_food_y_vmin(counter_lick) = l_food_y(ind_vmin_(counter_lick));
    l_food_x_offset(counter_lick) = l_food_x(ind_offset_(counter_lick));
    l_food_y_offset(counter_lick) = l_food_y(ind_offset_(counter_lick));

    r_nose_x_lick(counter_lick, :) = [r_nose_x(inds_)' pad];
    r_nose_y_lick(counter_lick, :) = [r_nose_y(inds_)' pad];
    r_nose_x_onset(counter_lick) = r_nose_x(ind_onset_(counter_lick));
    r_nose_y_onset(counter_lick) = r_nose_y(ind_onset_(counter_lick));
    r_nose_x_vmax(counter_lick) = r_nose_x(ind_vmax_(counter_lick));
    r_nose_y_vmax(counter_lick) = r_nose_y(ind_vmax_(counter_lick));
    r_nose_x_dmax(counter_lick) = r_nose_x(ind_dmax_(counter_lick));
    r_nose_y_dmax(counter_lick) = r_nose_y(ind_dmax_(counter_lick));
    r_nose_x_vmin(counter_lick) = r_nose_x(ind_vmin_(counter_lick));
    r_nose_y_vmin(counter_lick) = r_nose_y(ind_vmin_(counter_lick));
    r_nose_x_offset(counter_lick) = r_nose_x(ind_offset_(counter_lick));
    r_nose_y_offset(counter_lick) = r_nose_y(ind_offset_(counter_lick));

    l_nose_x_lick(counter_lick, :) = [l_nose_x(inds_)' pad];
    l_nose_y_lick(counter_lick, :) = [l_nose_y(inds_)' pad];
    l_nose_x_onset(counter_lick) = l_nose_x(ind_onset_(counter_lick));
    l_nose_y_onset(counter_lick) = l_nose_y(ind_onset_(counter_lick));
    l_nose_x_vmax(counter_lick) = l_nose_x(ind_vmax_(counter_lick));
    l_nose_y_vmax(counter_lick) = l_nose_y(ind_vmax_(counter_lick));
    l_nose_x_dmax(counter_lick) = l_nose_x(ind_dmax_(counter_lick));
    l_nose_y_dmax(counter_lick) = l_nose_y(ind_dmax_(counter_lick));
    l_nose_x_vmin(counter_lick) = l_nose_x(ind_vmin_(counter_lick));
    l_nose_y_vmin(counter_lick) = l_nose_y(ind_vmin_(counter_lick));
    l_nose_x_offset(counter_lick) = l_nose_x(ind_offset_(counter_lick));
    l_nose_y_offset(counter_lick) = l_nose_y(ind_offset_(counter_lick));

    r_tube_r_x_lick(counter_lick, :) = [r_tube_r_x(inds_)' pad];
    r_tube_r_y_lick(counter_lick, :) = [r_tube_r_y(inds_)' pad];
    r_tube_r_x_onset(counter_lick) = r_tube_r_x(ind_onset_(counter_lick));
    r_tube_r_y_onset(counter_lick) = r_tube_r_y(ind_onset_(counter_lick));
    r_tube_r_x_vmax(counter_lick) = r_tube_r_x(ind_vmax_(counter_lick));
    r_tube_r_y_vmax(counter_lick) = r_tube_r_y(ind_vmax_(counter_lick));
    r_tube_r_x_dmax(counter_lick) = r_tube_r_x(ind_dmax_(counter_lick));
    r_tube_r_y_dmax(counter_lick) = r_tube_r_y(ind_dmax_(counter_lick));
    r_tube_r_x_vmin(counter_lick) = r_tube_r_x(ind_vmin_(counter_lick));
    r_tube_r_y_vmin(counter_lick) = r_tube_r_y(ind_vmin_(counter_lick));
    r_tube_r_x_offset(counter_lick) = r_tube_r_x(ind_offset_(counter_lick));
    r_tube_r_y_offset(counter_lick) = r_tube_r_y(ind_offset_(counter_lick));

    r_tube_l_x_lick(counter_lick, :) = [r_tube_l_x(inds_)' pad];
    r_tube_l_y_lick(counter_lick, :) = [r_tube_l_y(inds_)' pad];
    r_tube_l_x_onset(counter_lick) = r_tube_l_x(ind_onset_(counter_lick));
    r_tube_l_y_onset(counter_lick) = r_tube_l_y(ind_onset_(counter_lick));
    r_tube_l_x_vmax(counter_lick) = r_tube_l_x(ind_vmax_(counter_lick));
    r_tube_l_y_vmax(counter_lick) = r_tube_l_y(ind_vmax_(counter_lick));
    r_tube_l_x_dmax(counter_lick) = r_tube_l_x(ind_dmax_(counter_lick));
    r_tube_l_y_dmax(counter_lick) = r_tube_l_y(ind_dmax_(counter_lick));
    r_tube_l_x_vmin(counter_lick) = r_tube_l_x(ind_vmin_(counter_lick));
    r_tube_l_y_vmin(counter_lick) = r_tube_l_y(ind_vmin_(counter_lick));
    r_tube_l_x_offset(counter_lick) = r_tube_l_x(ind_offset_(counter_lick));
    r_tube_l_y_offset(counter_lick) = r_tube_l_y(ind_offset_(counter_lick));

    l_tube_r_x_lick(counter_lick, :) = [l_tube_r_x(inds_)' pad];
    l_tube_r_y_lick(counter_lick, :) = [l_tube_r_y(inds_)' pad];
    l_tube_r_x_onset(counter_lick) = l_tube_r_x(ind_onset_(counter_lick));
    l_tube_r_y_onset(counter_lick) = l_tube_r_y(ind_onset_(counter_lick));
    l_tube_r_x_vmax(counter_lick) = l_tube_r_x(ind_vmax_(counter_lick));
    l_tube_r_y_vmax(counter_lick) = l_tube_r_y(ind_vmax_(counter_lick));
    l_tube_r_x_dmax(counter_lick) = l_tube_r_x(ind_dmax_(counter_lick));
    l_tube_r_y_dmax(counter_lick) = l_tube_r_y(ind_dmax_(counter_lick));
    l_tube_r_x_vmin(counter_lick) = l_tube_r_x(ind_vmin_(counter_lick));
    l_tube_r_y_vmin(counter_lick) = l_tube_r_y(ind_vmin_(counter_lick));
    l_tube_r_x_offset(counter_lick) = l_tube_r_x(ind_offset_(counter_lick));
    l_tube_r_y_offset(counter_lick) = l_tube_r_y(ind_offset_(counter_lick));

    l_tube_l_x_lick(counter_lick, :) = [l_tube_l_x(inds_)' pad];
    l_tube_l_y_lick(counter_lick, :) = [l_tube_l_y(inds_)' pad];
    l_tube_l_x_onset(counter_lick) = l_tube_l_x(ind_onset_(counter_lick));
    l_tube_l_y_onset(counter_lick) = l_tube_l_y(ind_onset_(counter_lick));
    l_tube_l_x_vmax(counter_lick) = l_tube_l_x(ind_vmax_(counter_lick));
    l_tube_l_y_vmax(counter_lick) = l_tube_l_y(ind_vmax_(counter_lick));
    l_tube_l_x_dmax(counter_lick) = l_tube_l_x(ind_dmax_(counter_lick));
    l_tube_l_y_dmax(counter_lick) = l_tube_l_y(ind_dmax_(counter_lick));
    l_tube_l_x_vmin(counter_lick) = l_tube_l_x(ind_vmin_(counter_lick));
    l_tube_l_y_vmin(counter_lick) = l_tube_l_y(ind_vmin_(counter_lick));
    l_tube_l_x_offset(counter_lick) = l_tube_l_x(ind_offset_(counter_lick));
    l_tube_l_y_offset(counter_lick) = l_tube_l_y(ind_offset_(counter_lick));
end

% ILI bout
for counter_bout = 1 : 1 : num_bout
    ILI_bout(counter_bout,1) = mean(diff(time_lick_onset(find(ind_lick_onset >= ind_lick_onset_str_bout(counter_bout) & ind_lick_onset <= ind_lick_onset_end_bout(counter_bout)))));
end



v_tip = ([0 (diff(d_tip)./ (time_1K(2) - time_1K(1)))]);
angle_midtip = rad2deg(atan2(midtip_y, midtip_x))';
angle_lick_max = angle_midtip(ind_d_lick_max)';

ILR_bout = 1./ILI_bout;
ILI_bout((ILR_bout==Inf)) = nan;
ILR_bout((ILR_bout==Inf)) = nan;

ILI_lick = diff(time_lick_onset);
ILI_lick = [ILI_lick; nan];
ILR_lick = 1./ILI_lick;

time_d_lick_max_abs = time_1K(ind_d_lick_max);
time_d_lick_max_rel = (time_d_lick_max_abs - time_lick_onset);
time_v_lick_max_abs = time_1K(ind_v_lick_max);
time_v_lick_max_rel = (time_v_lick_max_abs - time_lick_onset);
time_v_lick_min_abs = time_1K(ind_v_lick_min);
time_v_lick_min_rel = (time_v_lick_min_abs - time_lick_onset);

DLC.KINEMATIC.d_tip = d_tip;
DLC.KINEMATIC.d_lick = d_lick;
DLC.KINEMATIC.d_lick_max = d_lick_max;
DLC.IND.ind_d_lick_max = ind_d_lick_max;
DLC.KINEMATIC.v_tip = v_tip;
DLC.KINEMATIC.angle_lick = angle_lick;
DLC.KINEMATIC.angle_midtip = angle_midtip;
DLC.KINEMATIC.angle_lick_max = angle_lick_max;
DLC.KINEMATIC.v_lick = v_lick;
DLC.KINEMATIC.v_lick_max = v_lick_max;
DLC.IND.ind_v_lick_max = ind_v_lick_max;
DLC.KINEMATIC.v_lick_min = v_lick_min;
DLC.IND.ind_v_lick_min = ind_v_lick_min;
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

DLC.KINEMATIC.tip_tongue_x_lick = tip_tongue_x_lick;
DLC.KINEMATIC.tip_tongue_y_lick = tip_tongue_y_lick ;
DLC.KINEMATIC.r_tongue_x_lick = r_tongue_x_lick ;
DLC.KINEMATIC.r_tongue_y_lick = r_tongue_y_lick ;
DLC.KINEMATIC.l_tongue_x_lick = l_tongue_x_lick ;
DLC.KINEMATIC.l_tongue_y_lick = l_tongue_y_lick ;
DLC.KINEMATIC.mid_tongue_x_lick = mid_tongue_x_lick ;
DLC.KINEMATIC.mid_tongue_y_lick = mid_tongue_y_lick ;
DLC.KINEMATIC.r_nose_x_lick = r_nose_x_lick ;
DLC.KINEMATIC.r_nose_y_lick = r_nose_y_lick ;
DLC.KINEMATIC.l_nose_x_lick = l_nose_x_lick ;
DLC.KINEMATIC.l_nose_y_lick = l_nose_y_lick ;
DLC.KINEMATIC.l_food_x_lick = l_food_x_lick ;
DLC.KINEMATIC.l_food_y_lick = l_food_y_lick ;
DLC.KINEMATIC.l_tube_r_x_lick = l_tube_r_x_lick ;
DLC.KINEMATIC.l_tube_r_y_lick = l_tube_r_y_lick ;
DLC.KINEMATIC.l_tube_l_x_lick = l_tube_l_x_lick ;
DLC.KINEMATIC.l_tube_l_y_lick = l_tube_l_y_lick ;
DLC.KINEMATIC.r_food_x_lick = r_food_x_lick ;
DLC.KINEMATIC.r_food_y_lick = r_food_y_lick ;
DLC.KINEMATIC.r_tube_r_x_lick = r_tube_r_x_lick ;
DLC.KINEMATIC.r_tube_r_y_lick = r_tube_r_y_lick ;
DLC.KINEMATIC.r_tube_l_x_lick = r_tube_l_x_lick ;
DLC.KINEMATIC.r_tube_l_y_lick = r_tube_l_y_lick ;

DLC.KINEMATIC.tip_tongue_x_onset = tip_tongue_x_onset;
DLC.KINEMATIC.tip_tongue_y_onset = tip_tongue_y_onset ;
DLC.KINEMATIC.r_tongue_x_onset = r_tongue_x_onset ;
DLC.KINEMATIC.r_tongue_y_onset = r_tongue_y_onset ;
DLC.KINEMATIC.l_tongue_x_onset = l_tongue_x_onset ;
DLC.KINEMATIC.l_tongue_y_onset = l_tongue_y_onset ;
DLC.KINEMATIC.mid_tongue_x_onset = mid_tongue_x_onset ;
DLC.KINEMATIC.mid_tongue_y_onset = mid_tongue_y_onset ;
DLC.KINEMATIC.r_nose_x_onset = r_nose_x_onset ;
DLC.KINEMATIC.r_nose_y_onset = r_nose_y_onset ;
DLC.KINEMATIC.l_nose_x_onset = l_nose_x_onset ;
DLC.KINEMATIC.l_nose_y_onset = l_nose_y_onset ;
DLC.KINEMATIC.l_food_x_onset = l_food_x_onset ;
DLC.KINEMATIC.l_food_y_onset = l_food_y_onset ;
DLC.KINEMATIC.l_tube_r_x_onset = l_tube_r_x_onset ;
DLC.KINEMATIC.l_tube_r_y_onset = l_tube_r_y_onset ;
DLC.KINEMATIC.l_tube_l_x_onset = l_tube_l_x_onset ;
DLC.KINEMATIC.l_tube_l_y_onset = l_tube_l_y_onset ;
DLC.KINEMATIC.r_food_x_onset = r_food_x_onset ;
DLC.KINEMATIC.r_food_y_onset = r_food_y_onset ;
DLC.KINEMATIC.r_tube_r_x_onset = r_tube_r_x_onset ;
DLC.KINEMATIC.r_tube_r_y_onset = r_tube_r_y_onset ;
DLC.KINEMATIC.r_tube_l_x_onset = r_tube_l_x_onset ;
DLC.KINEMATIC.r_tube_l_y_onset = r_tube_l_y_onset ;

DLC.KINEMATIC.tip_tongue_x_vmax = tip_tongue_x_vmax;
DLC.KINEMATIC.tip_tongue_y_vmax = tip_tongue_y_vmax ;
DLC.KINEMATIC.r_tongue_x_vmax = r_tongue_x_vmax ;
DLC.KINEMATIC.r_tongue_y_vmax = r_tongue_y_vmax ;
DLC.KINEMATIC.l_tongue_x_vmax = l_tongue_x_vmax ;
DLC.KINEMATIC.l_tongue_y_vmax = l_tongue_y_vmax ;
DLC.KINEMATIC.mid_tongue_x_vmax = mid_tongue_x_vmax ;
DLC.KINEMATIC.mid_tongue_y_vmax = mid_tongue_y_vmax ;
DLC.KINEMATIC.r_nose_x_vmax = r_nose_x_vmax ;
DLC.KINEMATIC.r_nose_y_vmax = r_nose_y_vmax ;
DLC.KINEMATIC.l_nose_x_vmax = l_nose_x_vmax ;
DLC.KINEMATIC.l_nose_y_vmax = l_nose_y_vmax ;
DLC.KINEMATIC.l_food_x_vmax = l_food_x_vmax ;
DLC.KINEMATIC.l_food_y_vmax = l_food_y_vmax ;
DLC.KINEMATIC.l_tube_r_x_vmax = l_tube_r_x_vmax ;
DLC.KINEMATIC.l_tube_r_y_vmax = l_tube_r_y_vmax ;
DLC.KINEMATIC.l_tube_l_x_vmax = l_tube_l_x_vmax ;
DLC.KINEMATIC.l_tube_l_y_vmax = l_tube_l_y_vmax ;
DLC.KINEMATIC.r_food_x_vmax = r_food_x_vmax ;
DLC.KINEMATIC.r_food_y_vmax = r_food_y_vmax ;
DLC.KINEMATIC.r_tube_r_x_vmax = r_tube_r_x_vmax ;
DLC.KINEMATIC.r_tube_r_y_vmax = r_tube_r_y_vmax ;
DLC.KINEMATIC.r_tube_l_x_vmax = r_tube_l_x_vmax ;
DLC.KINEMATIC.r_tube_l_y_vmax = r_tube_l_y_vmax ;

DLC.KINEMATIC.tip_tongue_x_dmax = tip_tongue_x_dmax;
DLC.KINEMATIC.tip_tongue_y_dmax = tip_tongue_y_dmax ;
DLC.KINEMATIC.r_tongue_x_dmax = r_tongue_x_dmax ;
DLC.KINEMATIC.r_tongue_y_dmax = r_tongue_y_dmax ;
DLC.KINEMATIC.l_tongue_x_dmax = l_tongue_x_dmax ;
DLC.KINEMATIC.l_tongue_y_dmax = l_tongue_y_dmax ;
DLC.KINEMATIC.mid_tongue_x_dmax = mid_tongue_x_dmax ;
DLC.KINEMATIC.mid_tongue_y_dmax = mid_tongue_y_dmax ;
DLC.KINEMATIC.r_nose_x_dmax = r_nose_x_dmax ;
DLC.KINEMATIC.r_nose_y_dmax = r_nose_y_dmax ;
DLC.KINEMATIC.l_nose_x_dmax = l_nose_x_dmax ;
DLC.KINEMATIC.l_nose_y_dmax = l_nose_y_dmax ;
DLC.KINEMATIC.l_food_x_dmax = l_food_x_dmax ;
DLC.KINEMATIC.l_food_y_dmax = l_food_y_dmax ;
DLC.KINEMATIC.l_tube_r_x_dmax = l_tube_r_x_dmax ;
DLC.KINEMATIC.l_tube_r_y_dmax = l_tube_r_y_dmax ;
DLC.KINEMATIC.l_tube_l_x_dmax = l_tube_l_x_dmax ;
DLC.KINEMATIC.l_tube_l_y_dmax = l_tube_l_y_dmax ;
DLC.KINEMATIC.r_food_x_dmax = r_food_x_dmax ;
DLC.KINEMATIC.r_food_y_dmax = r_food_y_dmax ;
DLC.KINEMATIC.r_tube_r_x_dmax = r_tube_r_x_dmax ;
DLC.KINEMATIC.r_tube_r_y_dmax = r_tube_r_y_dmax ;
DLC.KINEMATIC.r_tube_l_x_dmax = r_tube_l_x_dmax ;
DLC.KINEMATIC.r_tube_l_y_dmax = r_tube_l_y_dmax ;

DLC.KINEMATIC.tip_tongue_x_vmin = tip_tongue_x_vmin;
DLC.KINEMATIC.tip_tongue_y_vmin = tip_tongue_y_vmin ;
DLC.KINEMATIC.r_tongue_x_vmin = r_tongue_x_vmin ;
DLC.KINEMATIC.r_tongue_y_vmin = r_tongue_y_vmin ;
DLC.KINEMATIC.l_tongue_x_vmin = l_tongue_x_vmin ;
DLC.KINEMATIC.l_tongue_y_vmin = l_tongue_y_vmin ;
DLC.KINEMATIC.mid_tongue_x_vmin = mid_tongue_x_vmin ;
DLC.KINEMATIC.mid_tongue_y_vmin = mid_tongue_y_vmin ;
DLC.KINEMATIC.r_nose_x_vmin = r_nose_x_vmin ;
DLC.KINEMATIC.r_nose_y_vmin = r_nose_y_vmin ;
DLC.KINEMATIC.l_nose_x_vmin = l_nose_x_vmin ;
DLC.KINEMATIC.l_nose_y_vmin = l_nose_y_vmin ;
DLC.KINEMATIC.l_food_x_vmin = l_food_x_vmin ;
DLC.KINEMATIC.l_food_y_vmin = l_food_y_vmin ;
DLC.KINEMATIC.l_tube_r_x_vmin = l_tube_r_x_vmin ;
DLC.KINEMATIC.l_tube_r_y_vmin = l_tube_r_y_vmin ;
DLC.KINEMATIC.l_tube_l_x_vmin = l_tube_l_x_vmin ;
DLC.KINEMATIC.l_tube_l_y_vmin = l_tube_l_y_vmin ;
DLC.KINEMATIC.r_food_x_vmin = r_food_x_vmin ;
DLC.KINEMATIC.r_food_y_vmin = r_food_y_vmin ;
DLC.KINEMATIC.r_tube_r_x_vmin = r_tube_r_x_vmin ;
DLC.KINEMATIC.r_tube_r_y_vmin = r_tube_r_y_vmin ;
DLC.KINEMATIC.r_tube_l_x_vmin = r_tube_l_x_vmin ;
DLC.KINEMATIC.r_tube_l_y_vmin = r_tube_l_y_vmin ;

DLC.KINEMATIC.tip_tongue_x_offset = tip_tongue_x_offset;
DLC.KINEMATIC.tip_tongue_y_offset = tip_tongue_y_offset ;
DLC.KINEMATIC.r_tongue_x_offset = r_tongue_x_offset ;
DLC.KINEMATIC.r_tongue_y_offset = r_tongue_y_offset ;
DLC.KINEMATIC.l_tongue_x_offset = l_tongue_x_offset ;
DLC.KINEMATIC.l_tongue_y_offset = l_tongue_y_offset ;
DLC.KINEMATIC.mid_tongue_x_offset = mid_tongue_x_offset ;
DLC.KINEMATIC.mid_tongue_y_offset = mid_tongue_y_offset ;
DLC.KINEMATIC.r_nose_x_offset = r_nose_x_offset ;
DLC.KINEMATIC.r_nose_y_offset = r_nose_y_offset ;
DLC.KINEMATIC.l_nose_x_offset = l_nose_x_offset ;
DLC.KINEMATIC.l_nose_y_offset = l_nose_y_offset ;
DLC.KINEMATIC.l_food_x_offset = l_food_x_offset ;
DLC.KINEMATIC.l_food_y_offset = l_food_y_offset ;
DLC.KINEMATIC.l_tube_r_x_offset = l_tube_r_x_offset ;
DLC.KINEMATIC.l_tube_r_y_offset = l_tube_r_y_offset ;
DLC.KINEMATIC.l_tube_l_x_offset = l_tube_l_x_offset ;
DLC.KINEMATIC.l_tube_l_y_offset = l_tube_l_y_offset ;
DLC.KINEMATIC.r_food_x_offset = r_food_x_offset ;
DLC.KINEMATIC.r_food_y_offset = r_food_y_offset ;
DLC.KINEMATIC.r_tube_r_x_offset = r_tube_r_x_offset ;
DLC.KINEMATIC.r_tube_r_y_offset = r_tube_r_y_offset ;
DLC.KINEMATIC.r_tube_l_x_offset = r_tube_l_x_offset ;
DLC.KINEMATIC.r_tube_l_y_offset = r_tube_l_y_offset ;

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure;
    hold on;
    subplot(1,6,1)
    boxplot(d_lick_max, {'D max'});
    ylabel('mm');
    ylim([0 inf])
    subplot(1,6,2)
    boxplot(v_lick_max, {'V max'});
    ylabel('mm/s');
    ylim([0 inf])
    subplot(1,6,3)
    boxplot(v_lick_min, {'V min'});
    ylabel('mm/s');
    ylim([-inf 0])
    subplot(1,6,4)
    boxplot(time_lick_duration*1000, {'Duration'});
    ylabel('ms');
    ylim([0 inf])
    subplot(1,6,5)
    boxplot(ILI_bout*1000, {'ILI'});
    ylabel('ms');
    ylim([0 inf])
    subplot(1,6,6)
    boxplot(ILR_bout, {'ILR'});
    ylabel('hz');
    ylim([0 inf])
    sgtitle([EXPERIMENT_PARAMS.file_name], 'interpreter', 'none')
    ESN_Beautify_Plot(gcf, [20 10])
end

fprintf(' --> Completed. \n')
end

%% Function: Geometrization
function [DLC, EXPERIMENT_PARAMS] = geometrization(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Geometrizing frames ...');

tip_tongue_x = DLC.POINTS.tip_tongue_x;
tip_tongue_y = DLC.POINTS.tip_tongue_y;
r_tongue_x = DLC.POINTS.r_tongue_x;
r_tongue_y = DLC.POINTS.r_tongue_y;
l_tongue_x = DLC.POINTS.l_tongue_x;
l_tongue_y = DLC.POINTS.l_tongue_y;
mid_tongue_x = DLC.POINTS.mid_tongue_x;
mid_tongue_y = DLC.POINTS.mid_tongue_y;
l_tube_r_x = DLC.POINTS.l_tube_r_x;
l_tube_r_y = DLC.POINTS.l_tube_r_y;
l_tube_l_x = DLC.POINTS.l_tube_l_x;
l_tube_l_y = DLC.POINTS.l_tube_l_y;
l_food_x = DLC.POINTS.l_food_x;
l_food_y = DLC.POINTS.l_food_y;
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
        l_tongue_y(counter_frame),mid_tongue_y(counter_frame),r_tongue_y(counter_frame)], 'Simplify', false);

    geo_r_tube_empty = polyshape([r_food_x(counter_frame),r_tube_r_x(counter_frame),...
        r_tube_r_x(counter_frame),r_tube_l_x(counter_frame), r_tube_l_x(counter_frame)],[r_food_y(counter_frame),...
        r_food_y(counter_frame), r_tube_r_y(counter_frame),r_tube_l_y(counter_frame), r_food_y(counter_frame)], 'Simplify', false);

    geo_r_tube_full = polyshape([r_tube_r_x(counter_frame), r_tube_r_x(counter_frame),...
        r_tube_l_x(counter_frame), r_tube_l_x(counter_frame), r_food_x(counter_frame)],...
        [r_food_y(counter_frame), r_tube_r_y(counter_frame) + abs(max(r_food_y)-min(r_food_y)), ...
        r_tube_l_y(counter_frame) + abs(max(r_food_y)-min(r_food_y)), r_food_y(counter_frame), r_food_y(counter_frame)], 'Simplify', false);

    geo_inter_tongue_r_tube_empty = intersect(geo_tongue, geo_r_tube_empty);

    geo_inter_tongue_r_tube_full = intersect(geo_tongue, geo_r_tube_full);

    geo_all = [geo_tongue geo_r_tube_empty geo_r_tube_full ...
        geo_inter_tongue_r_tube_empty geo_inter_tongue_r_tube_full];

    geo_l_tube_empty = polyshape([l_food_x(counter_frame),l_tube_r_x(counter_frame),...
        l_tube_r_x(counter_frame),l_tube_l_x(counter_frame), l_tube_l_x(counter_frame)],[l_food_y(counter_frame),...
        l_food_y(counter_frame), l_tube_r_y(counter_frame),l_tube_l_y(counter_frame), l_food_y(counter_frame)], 'Simplify', false);

    geo_l_tube_full = polyshape([l_tube_r_x(counter_frame), l_tube_r_x(counter_frame),...
        l_tube_l_x(counter_frame), l_tube_l_x(counter_frame), l_food_x(counter_frame)],...
        [l_food_y(counter_frame), l_tube_r_y(counter_frame) - abs(max(l_food_y)-min(l_food_y)), ...
        l_tube_l_y(counter_frame) - abs(max(l_food_y)-min(l_food_y)), l_food_y(counter_frame), l_food_y(counter_frame)], 'Simplify', false);

    geo_inter_tongue_l_tube_empty = intersect(geo_tongue, geo_l_tube_empty);

    geo_inter_tongue_l_tube_full = intersect(geo_tongue, geo_l_tube_full);

    geo_all = [geo_tongue geo_r_tube_empty geo_r_tube_full geo_l_tube_empty geo_l_tube_full...
        geo_inter_tongue_r_tube_empty geo_inter_tongue_r_tube_full geo_inter_tongue_l_tube_empty...
        geo_inter_tongue_l_tube_full ];


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
    [cent_tongue_l_tube_empty_x(counter_lick, 1),cent_tongue_l_tube_empty_y(counter_lick, 1)]  = centroid(geo_inter_tongue_l_tube_empty);
    [cent_tongue_l_tube_full_x(counter_lick, 1),cent_tongue_l_tube_full_y(counter_lick, 1)] = centroid(geo_inter_tongue_l_tube_full);
    bool_tongue_l_tube_empty(counter_lick, 1)= bool_overlaps_all(1,4);
    bool_tongue_l_tube_full(counter_lick, 1) = bool_overlaps_all(1,5);
    area_l_tube_empty(counter_lick, 1) = area(geo_l_tube_empty);
    area_l_tube_full(counter_lick, 1) = area(geo_l_tube_full);
    area_inter_tongue_l_tube_empty(counter_lick, 1) = area(geo_inter_tongue_l_tube_empty);
    area_inter_tongue_l_tube_full(counter_lick, 1) = area(geo_inter_tongue_l_tube_full);

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

fprintf(' --> Completed. \n')
end

%% Function: Sort licks
function [DLC, EXPERIMENT_PARAMS] = sort_licks(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Sorting licks ...');

bool_tongue_r_tube_full = DLC.GEO.bool_tongue_r_tube_full;
bool_tongue_l_tube_full = DLC.GEO.bool_tongue_l_tube_full;
area_r_tube = DLC.GEO.area_r_tube_empty + DLC.GEO.area_r_tube_full;
area_l_tube = DLC.GEO.area_l_tube_empty + DLC.GEO.area_l_tube_full;
r_tube_food_perc = DLC.GEO.area_r_tube_full./area_r_tube;
l_tube_food_perc = DLC.GEO.area_l_tube_full./area_l_tube;

tip_tongue_x = DLC.POINTS.tip_tongue_x;
tip_tongue_y = DLC.POINTS.tip_tongue_y;
r_tongue_x = DLC.POINTS.r_tongue_x;
r_tongue_y = DLC.POINTS.r_tongue_y;
l_tongue_x = DLC.POINTS.l_tongue_x;
l_tongue_y = DLC.POINTS.l_tongue_y;

r_nose_x =  DLC.POINTS.r_nose_x;
r_nose_y = DLC.POINTS.r_nose_y;
l_nose_x =  DLC.POINTS.l_nose_x;
l_nose_y = DLC.POINTS.l_nose_y;

r_tube_r_x = DLC.POINTS.r_tube_r_x;
r_tube_r_y = DLC.POINTS.r_tube_r_y;
r_tube_l_x = DLC.POINTS.r_tube_l_x;
r_tube_l_y = DLC.POINTS.r_tube_l_y;
% r_food_y = DLC.POINTS.r_food_y;
% r_tube_food_perc = DLC.data.r_food_y./(DLC.FILE.vid_height - (DLC.data.r_tube_r_y + DLC.data.r_tube_l_y)/2);

l_tube_r_x = DLC.POINTS.l_tube_r_x;
l_tube_r_y = DLC.POINTS.l_tube_r_y;
l_tube_l_x = DLC.POINTS.l_tube_l_x;
l_tube_l_y = DLC.POINTS.l_tube_l_y;
% l_food_y = DLC.POINTS.l_food_y;
% l_tube_food_perc = DLC.data.l_food_y./((DLC.data.l_tube_r_y + DLC.data.l_tube_l_y)/2);

ind_d_lick_max = DLC.IND.ind_d_lick_max;
ind_lick_onset = DLC.IND.ind_lick_onset;
ind_lick_offset = DLC.IND.ind_lick_offset;
time_1K = DLC.TIME.time_1K';

% Classify inner, outer edge, under tube, and grooming
is_r_inner_tube_lick = false(size(ind_lick_onset));
is_l_inner_tube_lick = false(size(ind_lick_onset));
is_r_outer_edge_lick = false(size(ind_lick_onset));
is_l_outer_edge_lick = false(size(ind_lick_onset));
is_r_under_tube_lick = false(size(ind_lick_onset));
is_l_under_tube_lick = false(size(ind_lick_onset));
is_grooming_lick = false(size(ind_lick_onset));

% tongue markers relative to tube markers
tip_crosses_far_edge_r_tube = tip_tongue_x > r_tube_r_x + std(r_tube_r_x);
tip_crosses_far_edge_l_tube = tip_tongue_x > l_tube_r_x + std(l_tube_r_x);
tip_crosses_close_edge_r_tube = tip_tongue_x < r_tube_l_x + std(r_tube_l_x);
tip_crosses_close_edge_l_tube = tip_tongue_x < l_tube_l_x + std(l_tube_l_x);

tip_x_in_r_tube = tip_tongue_x < r_tube_r_x & tip_tongue_x > r_tube_l_x; 
r_x_in_r_tube = r_tongue_x < r_tube_r_x & r_tongue_x > r_tube_l_x;
l_x_in_r_tube = l_tongue_x < r_tube_r_x & l_tongue_x > r_tube_l_x; 
tip_x_in_l_tube = tip_tongue_x < l_tube_r_x & tip_tongue_x > l_tube_l_x; 
r_x_in_l_tube = r_tongue_x < l_tube_r_x & r_tongue_x > l_tube_l_x; 
l_x_in_l_tube = l_tongue_x < l_tube_r_x & l_tongue_x > l_tube_l_x;  
tip_y_in_r_tube = tip_tongue_y > (r_tube_r_y+r_tube_l_y)/2 - 1; 
r_y_in_r_tube = r_tongue_y > (r_tube_r_y+r_tube_l_y)/2 - 1; 
l_y_in_r_tube = l_tongue_y > (r_tube_r_y+r_tube_l_y)/2 - 1; 
tip_y_in_l_tube = tip_tongue_y < (l_tube_r_y+l_tube_l_y)/2 + 1; 
r_y_in_l_tube = r_tongue_y < (l_tube_r_y+l_tube_l_y)/2 + 1; 
l_y_in_l_tube = l_tongue_y < (l_tube_r_y+l_tube_l_y)/2 + 1; 

tip_in_r_tube = tip_x_in_r_tube & tip_y_in_r_tube; 
r_in_r_tube = r_x_in_r_tube & r_y_in_r_tube; 
l_in_r_tube = l_x_in_r_tube & l_y_in_r_tube; 
tip_in_l_tube = tip_x_in_l_tube & tip_y_in_l_tube;
r_in_l_tube = r_x_in_l_tube & r_y_in_l_tube; 
l_in_l_tube = l_x_in_l_tube & l_y_in_l_tube;

for i = 1:length(ind_lick_onset)
    done = false;
    ind = ind_lick_onset(i):ind_lick_offset(i);

%   INNER TUBE LICK DEFINITION
%   - there are frames where the tip, l, r tongue markers are all in 
%     the tube during the lick
%   - when tongue markers are in the y-coordinate range of the tube, 
%     they don't ever cross either of the edges of the tube

    % r inner tube lick
    tongue_in_tube = tip_in_r_tube(ind) & ...
        r_in_r_tube(ind) & ...
        l_in_r_tube(ind);
    tongue_crosses_edge = (tip_y_in_r_tube(ind) & ~tip_x_in_r_tube(ind)) | ...
        (r_y_in_r_tube(ind) & ~r_x_in_r_tube(ind)) | ...
        (l_y_in_r_tube(ind) & ~l_x_in_r_tube(ind));
    if isempty(find(tongue_crosses_edge, 1)) && ~isempty(find(tongue_in_tube, 1))
        is_r_inner_tube_lick(i) = true;
        done = true;
    end

    % l inner tube lick
    tongue_in_tube = tip_in_l_tube(ind) & ...
        r_in_l_tube(ind) & ...
        l_in_l_tube(ind);
    tongue_crosses_edge = (tip_y_in_l_tube(ind) & ~tip_x_in_l_tube(ind)) | ...
        (r_y_in_l_tube(ind) & ~r_x_in_l_tube(ind)) | ...
        (l_y_in_l_tube(ind) & ~l_x_in_l_tube(ind));
    if isempty(find(tongue_crosses_edge, 1)) && ~isempty(find(tongue_in_tube, 1))
        is_l_inner_tube_lick(i) = true;
        done = true;
    end

%   UNDER TUBE LICK DEFINITION
%   - not an inner tube
%   - the tip tongue marker crosses the close edge of the tube when it is
%     in the y_coordinate range of the tube
%   - otherwise, there are frames where the tip tongue marker is in the 
%     y_coordinate range of the tube but never crosses the far edge of the tube
    tongue_crosses_close_edge = tip_crosses_close_edge_r_tube(ind) & tip_y_in_r_tube(ind);
    tongue_crosses_far_edge = tip_crosses_far_edge_r_tube(ind) & tip_y_in_r_tube(ind);
    % r under tube lick
    if ~done && (~isempty(find(tongue_crosses_close_edge, 1))) || ...
            (~done && (isempty(find(tongue_crosses_far_edge, 1))) && ~isempty(find(tip_y_in_r_tube(ind), 1)))
        is_r_under_tube_lick(i) = true;
        done = true;
    end

    tongue_crosses_close_edge = tip_crosses_close_edge_l_tube(ind) & tip_y_in_l_tube(ind);
    tongue_crosses_far_edge = tip_crosses_far_edge_l_tube(ind) & tip_y_in_l_tube(ind);
    % l under tube lick
    if ~done && (~isempty(find(tongue_crosses_close_edge, 1))) || ...
            (~done && (isempty(find(tongue_crosses_far_edge, 1))) && ~isempty(find(tip_y_in_l_tube(ind), 1)))
         is_l_under_tube_lick(i) = true;
        done = true;
    end

%   OUTER EDGE LICK DEFINITION
%   - not an inner tube or under tube lick
%   - there are frames where the tip tongue marker has crossed the far 
%     edge when it is in the y-coordinate range of the tube

    % r outer edge lick
    tongue_crosses_far_edge = tip_crosses_far_edge_r_tube(ind) & tip_y_in_r_tube(ind);
    if ~done && ~isempty(find(tongue_crosses_far_edge, 1))
        is_r_outer_edge_lick(i) = true;
        done = true;
    end

     % l outer edge lick
    tongue_crosses_far_edge = tip_crosses_far_edge_l_tube(ind) & tip_y_in_l_tube(ind);
    if ~done && ~isempty(find(tongue_crosses_far_edge, 1))
        is_l_outer_edge_lick(i) = true;
        done = true;
    end

    % GROOMING LICK DEFINITION
    % not inner, outer edge, or under tube lick
    if ~done
        is_grooming_lick(i) = true;
    end
end

% grooming
ind_grooming_lick = find(is_grooming_lick);
% inner reward/noreward
is_r_reward_inner_tube_lick = bool_tongue_r_tube_full == 1 & is_r_inner_tube_lick;
is_r_noreward_inner_tube_lick = is_r_inner_tube_lick & ~is_r_reward_inner_tube_lick;
is_l_reward_inner_tube_lick = bool_tongue_l_tube_full == 1 & is_l_inner_tube_lick;
is_l_noreward_inner_tube_lick = is_l_inner_tube_lick & ~is_l_reward_inner_tube_lick;
ind_r_reward_inner_tube_lick = find(is_r_reward_inner_tube_lick);
ind_r_noreward_inner_tube_lick = find(is_r_noreward_inner_tube_lick);
ind_l_reward_inner_tube_lick = find(is_l_reward_inner_tube_lick);
ind_l_noreward_inner_tube_lick = find(is_l_noreward_inner_tube_lick);
% outer edge reward/noreward
is_r_reward_outer_edge_lick = bool_tongue_r_tube_full == 1 & r_tube_food_perc >= 0.90 & is_r_outer_edge_lick;
is_r_noreward_outer_edge_lick = is_r_outer_edge_lick & ~is_r_reward_outer_edge_lick;
is_l_reward_outer_edge_lick = bool_tongue_l_tube_full == 1 & l_tube_food_perc >= 0.90 & is_l_outer_edge_lick;
is_l_noreward_outer_edge_lick = is_l_outer_edge_lick & ~is_l_reward_outer_edge_lick;
ind_r_reward_outer_edge_lick = find(is_r_reward_outer_edge_lick);
ind_r_noreward_outer_edge_lick = find(is_r_noreward_outer_edge_lick);
ind_l_reward_outer_edge_lick = find(is_l_reward_outer_edge_lick);
ind_l_noreward_outer_edge_lick = find(is_l_noreward_outer_edge_lick);
% under tube reward/noreward
is_r_reward_under_tube_lick = bool_tongue_r_tube_full == 1 & r_tube_food_perc >= 0.90 & is_r_under_tube_lick;
is_r_noreward_under_tube_lick = is_r_under_tube_lick & ~is_r_reward_under_tube_lick;
is_l_reward_under_tube_lick = bool_tongue_l_tube_full == 1 & l_tube_food_perc >= 0.90 & is_l_under_tube_lick;
is_l_noreward_under_tube_lick = is_l_under_tube_lick & ~is_l_reward_under_tube_lick;
ind_r_reward_under_tube_lick = find(is_r_reward_under_tube_lick);
ind_r_noreward_under_tube_lick = find(is_r_noreward_under_tube_lick);
ind_l_reward_under_tube_lick = find(is_l_reward_under_tube_lick);
ind_l_noreward_under_tube_lick = find(is_l_noreward_under_tube_lick);

% inds for reward-driven licks
ind_lick_onset_r_reward = ind_lick_onset(sort([ind_r_reward_inner_tube_lick;ind_r_reward_outer_edge_lick;ind_r_reward_under_tube_lick; ind_r_noreward_inner_tube_lick;ind_r_noreward_outer_edge_lick;ind_r_noreward_under_tube_lick]));
time_lick_onset_r_reward = time_1K(ind_lick_onset_r_reward);
ind_lick_offset_r_reward = ind_lick_offset(sort([ind_r_reward_inner_tube_lick;ind_r_reward_outer_edge_lick;ind_r_reward_under_tube_lick; ind_r_noreward_inner_tube_lick;ind_r_noreward_outer_edge_lick;ind_r_noreward_under_tube_lick]));
time_lick_offset_r_reward = time_1K(ind_lick_offset_r_reward);
ind_lick_onset_l_reward = ind_lick_onset(sort([ind_l_reward_inner_tube_lick;ind_l_reward_outer_edge_lick;ind_l_reward_under_tube_lick; ind_l_noreward_inner_tube_lick;ind_l_noreward_outer_edge_lick;ind_l_noreward_under_tube_lick]));
time_lick_onset_l_reward = time_1K(ind_lick_onset_l_reward);
ind_lick_offset_l_reward = ind_lick_offset(sort([ind_l_reward_inner_tube_lick;ind_l_reward_outer_edge_lick;ind_l_reward_under_tube_lick; ind_l_noreward_inner_tube_lick;ind_l_noreward_outer_edge_lick;ind_l_noreward_under_tube_lick]));
time_lick_offset_l_reward = time_1K(ind_lick_offset_l_reward);

% add vars to DLC
DLC.CLASS.is_grooming_lick = is_grooming_lick;
DLC.CLASS.is_r_reward_inner_tube_lick = is_r_reward_inner_tube_lick;
DLC.CLASS.is_r_reward_outer_edge_lick = is_r_reward_outer_edge_lick;
DLC.CLASS.is_r_reward_under_tube_lick = is_r_reward_under_tube_lick;
DLC.CLASS.is_r_noreward_inner_tube_lick = is_r_noreward_inner_tube_lick;
DLC.CLASS.is_r_noreward_outer_edge_lick = is_r_noreward_outer_edge_lick;
DLC.CLASS.is_r_noreward_under_tube_lick = is_r_noreward_under_tube_lick;

DLC.CLASS.ind_lick_onset_r_reward = ind_lick_onset_r_reward;
DLC.TIME.time_lick_onset_r_reward = time_lick_onset_r_reward;
DLC.CLASS.ind_lick_offset_r_reward = ind_lick_offset_r_reward;
DLC.TIME.time_lick_offset_r_reward = time_lick_offset_r_reward;

DLC.CLASS.is_l_reward_inner_tube_lick = is_l_reward_inner_tube_lick;
DLC.CLASS.is_l_reward_outer_edge_lick = is_l_reward_outer_edge_lick;
DLC.CLASS.is_l_reward_under_tube_lick = is_l_reward_under_tube_lick;
DLC.CLASS.is_l_noreward_inner_tube_lick = is_l_noreward_inner_tube_lick;
DLC.CLASS.is_l_noreward_outer_edge_lick = is_l_noreward_outer_edge_lick;
DLC.CLASS.is_l_noreward_under_tube_lick = is_l_noreward_under_tube_lick;

DLC.CLASS.ind_lick_onset_l_reward = ind_lick_onset_l_reward;
DLC.TIME.time_lick_onset_l_reward = time_lick_onset_l_reward;
DLC.CLASS.ind_lick_offset_l_reward = ind_lick_offset_l_reward;
DLC.TIME.time_lick_offset_l_reward = time_lick_offset_l_reward;

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure;
    hold on;
    % plot(tip_tongue_x(ind_d_lick_max),tip_tongue_y(ind_d_lick_max), '.k');
    plot(tip_tongue_x(ind_d_lick_max(ind_grooming_lick)),tip_tongue_y(ind_d_lick_max(ind_grooming_lick)), 'ob');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_inner_tube_lick)), 'og');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_inner_tube_lick)), 'or');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_outer_edge_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_outer_edge_lick)), 'oc');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_outer_edge_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_outer_edge_lick)), 'om');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_under_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_under_tube_lick)), 'oy');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_under_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_under_tube_lick)), 'o', 'Color', [0.7 0.7 0.7]);
    plot(r_nose_x, r_nose_y,'ok');
    plot(l_nose_x, l_nose_y,'ok');
    plot(r_tube_r_x(ind_d_lick_max),r_tube_r_y(ind_d_lick_max),'sk');
    plot(r_tube_l_x(ind_d_lick_max),r_tube_l_y(ind_d_lick_max),'sk');
     plot(tip_tongue_x(ind_d_lick_max(ind_l_noreward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_noreward_inner_tube_lick)), 'og');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_reward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_reward_inner_tube_lick)), 'or');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_noreward_outer_edge_lick)),tip_tongue_y(ind_d_lick_max(ind_l_noreward_outer_edge_lick)), 'oc');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_reward_outer_edge_lick)),tip_tongue_y(ind_d_lick_max(ind_l_reward_outer_edge_lick)), 'om');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_noreward_under_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_noreward_under_tube_lick)), 'oy');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_reward_under_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_reward_under_tube_lick)), 'o', 'Color', [0.7 0.7 0.7]);
    plot(l_tube_r_x(ind_d_lick_max), l_tube_r_y(ind_d_lick_max),'sk');
    plot(l_tube_l_x(ind_d_lick_max), l_tube_l_y(ind_d_lick_max),'sk');
    xlabel('x position (mm)');
    ylabel('y position (mm)');
    set(gca, 'YDir','reverse')
    xlim([-10 15]);
    ylim([-25 25]);
    title([ EXPERIMENT_PARAMS.file_name ' | Lick sorting summary'], 'interpreter', 'none');
    ESN_Beautify_Plot(gcf, [20 10])
end

fprintf(' --> Completed. \n')
end

%% Function: Quantify food in tube
function [DLC, EXPERIMENT_PARAMS] = quantify_food(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Quantifying food in tube ...');
num_lick = DLC.IND.num_lick;
num_bout = DLC.IND.num_bout;
ind_lick_onset_str_bout = DLC.IND.ind_lick_onset_str_bout;
ind_lick_onset_end_bout = DLC.IND.ind_lick_onset_end_bout ;
ind_lick_onset = DLC.IND.ind_lick_onset;
ind_lick_offset = DLC.IND.ind_lick_offset;

r_tube_r_y = DLC.POINTS.r_tube_r_y;
r_food_x = DLC.POINTS.r_food_x;
r_food_y = DLC.POINTS.r_food_y;
r_food_x(r_food_y < 0) = nan;
r_food_y(r_food_y < 0) = nan;
r_food_max = max(r_food_y);
r_food_min = min(r_food_y);
r_tube_food = (r_food_y - r_food_max)./(r_tube_r_y - r_food_max);


l_tube_r_y = DLC.POINTS.l_tube_r_y;
l_food_x = DLC.POINTS.l_food_x;
l_food_y = DLC.POINTS.l_food_y;
l_food_x(l_food_y > 0) = nan;
l_food_y(l_food_y > 0) = nan;
l_food_max = max(l_food_y);
l_food_min = min(l_food_y);
l_tube_food = (l_food_y - l_food_min)./(l_tube_r_y - l_food_min);

for counter_lick =  1 : 1 : num_lick
    inds_ = ind_lick_onset(counter_lick):ind_lick_offset(counter_lick);
    ind_onset_(counter_lick) = DLC.IND.ind_lick_onset(counter_lick);
    ind_vmax_(counter_lick) = DLC.IND.ind_v_lick_max(counter_lick);
    ind_dmax_(counter_lick) = DLC.IND.ind_d_lick_max(counter_lick);
    ind_vmin_(counter_lick) = DLC.IND.ind_v_lick_min(counter_lick);
    ind_offset_(counter_lick) = DLC.IND.ind_lick_offset(counter_lick);
    for counter_inds_ = 1 : 1: length(inds_)
        r_tube_food_(counter_inds_) = r_tube_food(inds_(counter_inds_));
        r_tube_food_lick_onset(counter_lick,1) = r_tube_food_(1) ;
        r_tube_food_lick_offset(counter_lick,1) = r_tube_food_(end) ;
        l_tube_food_(counter_inds_) = l_tube_food(inds_(counter_inds_));
        l_tube_food_lick_onset(counter_lick,1) = l_tube_food_(1) ;
        l_tube_food_lick_offset(counter_lick,1) = l_tube_food_(end) ;

    end
end

for counter_bout =  1 : 1 : num_bout
    inds_ = ind_lick_onset_str_bout(counter_bout):ind_lick_onset_end_bout(counter_bout);
    for counter_inds_ = 1 : 1: length(inds_)
        r_tube_food_(counter_inds_) = r_tube_food(inds_(counter_inds_));
        r_tube_food_bout_start(counter_bout, 1) = r_tube_food_(1);
        r_tube_food_bout_end(counter_bout, 1) = r_tube_food_(end);
        l_tube_food_(counter_inds_) = l_tube_food(inds_(counter_inds_));
        l_tube_food_bout_start(counter_bout, 1) = l_tube_food_(1);
        l_tube_food_bout_end(counter_bout, 1) = l_tube_food_(end);
    end
end

DLC.FOOD.r_tube_food = r_tube_food;
DLC.FOOD.r_food_x = r_food_x;
DLC.FOOD.r_food_y = r_food_y;
DLC.FOOD.l_tube_food = l_tube_food;
DLC.FOOD.l_food_x = l_food_x;
DLC.FOOD.l_food_y = l_food_y;

DLC.FOOD.r_tube_food_lick_onset = r_tube_food_lick_onset;
DLC.FOOD.r_tube_food_lick_offset = r_tube_food_lick_offset;
DLC.FOOD.l_tube_food_lick_onset = l_tube_food_lick_onset;
DLC.FOOD.l_tube_food_lick_offset = l_tube_food_lick_offset;

DLC.FOOD.r_tube_food_bout_start = r_tube_food_bout_start;
DLC.FOOD.r_tube_food_bout_end = r_tube_food_bout_end;
DLC.FOOD.l_tube_food_bout_start = l_tube_food_bout_start;
DLC.FOOD.l_tube_food_bout_end = l_tube_food_bout_end;
fprintf(' --> Completed. \n')
end

%% Function: Detect harvest str, end, num, times, and duration
function [DLC, EXPERIMENT_PARAMS] = detect_harvest(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf('Detecting harvest ...');

d_tip = DLC.KINEMATIC.d_tip;
time_1K = DLC.TIME.time_1K';
ind_lick_onset = DLC.IND.ind_lick_onset;
is_grooming_lick = DLC.CLASS.is_grooming_lick;
ind_lick_onset_grooming = (ind_lick_onset(is_grooming_lick));
ind_lick_onset_reward = (ind_lick_onset(~is_grooming_lick));
ind_lick_onset_r_reward = DLC.CLASS.ind_lick_onset_r_reward;
time_lick_onset_r_reward = DLC.TIME.time_lick_onset_r_reward;
ind_lick_onset_l_reward = DLC.CLASS.ind_lick_onset_l_reward;
time_lick_onset_l_reward = DLC.TIME.time_lick_onset_l_reward;
r_tube_food = DLC.FOOD.r_tube_food;
l_tube_food = DLC.FOOD.l_tube_food;

ind_lick_onset_str_bout = DLC.IND.ind_lick_onset_str_bout;
ind_lick_onset_end_bout = DLC.IND.ind_lick_onset_end_bout;

% Check bouts that start or end with grooming lick
% ind_lick_onset_str_bout_grooming = ismember(ind_lick_onset_str_bout, ind_lick_onset_grooming);
% ind_lick_onset_end_bout_grooming = ismember(ind_lick_onset_end_bout, ind_lick_onset_grooming);

for counter_bout = 1 : DLC.IND.num_bout
    % harvest start
    if isempty(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset_grooming,1))
        ind_lick_onset_str_harvest(counter_bout,1) = ind_lick_onset_str_bout(counter_bout);
    elseif ~isempty(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset_grooming,1))
        shift_ind = 0;
        bool_shift = 1;
        while bool_shift == 1
            if ~isempty(find(ind_lick_onset(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset,1) + shift_ind) == ind_lick_onset_grooming, 1)) && ...
                    ind_lick_onset(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset,1) + shift_ind) ~= ind_lick_onset(end)
                shift_ind = shift_ind + 1;
            elseif isempty(find(ind_lick_onset(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset,1) + shift_ind) == ind_lick_onset_grooming, 1)) || ...
                    ind_lick_onset(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset,1) + shift_ind) == ind_lick_onset(end)
                bool_shift = 0;
            end
        end
        ind_lick_onset_str_harvest(counter_bout,1) = ind_lick_onset(find(ind_lick_onset_str_bout(counter_bout) == ind_lick_onset,1) + shift_ind);
    end
    % harvest end
    if isempty(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset_grooming,1))
        ind_lick_onset_end_harvest(counter_bout,1) = ind_lick_onset_end_bout(counter_bout);
    elseif ~isempty(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset_grooming,1))
        shift_ind = 0;
        bool_shift = 1;
        while bool_shift == 1
            if ~isempty(find(ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) == ind_lick_onset_grooming, 1)) && ...
                    ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) ~= ind_lick_onset_str_harvest(counter_bout) && ...
                    ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) ~= ind_lick_onset(1)
                shift_ind = shift_ind + 1;
            elseif isempty(find(ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) == ind_lick_onset_grooming, 1)) || ...
                    ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) == ind_lick_onset_str_harvest(counter_bout) || ...
                    ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind) == ind_lick_onset(1)
                bool_shift = 0;
            end
        end
        ind_lick_onset_end_harvest(counter_bout,1) = ind_lick_onset(find(ind_lick_onset_end_bout(counter_bout) == ind_lick_onset,1) - shift_ind);
    end
end

% Check harcest that start or end with grooming lick
% ind_lick_onset_str_harvest_grooming = ismember(ind_lick_onset_str_harvest, ind_lick_onset_grooming);
% ind_lick_onset_end_harvest_grooming = ismember(ind_lick_onset_end_harvest, ind_lick_onset_grooming);

% Count number of licks in harvest
for counter_harvest = 1 : length(ind_lick_onset_str_harvest)
    inds_ = ind_lick_onset_str_harvest(counter_harvest) : ind_lick_onset_end_harvest(counter_harvest);
    num_lick_harvest(counter_harvest,1) = length(find(inds_ == ind_lick_onset));
end

validity = num_lick_harvest>=3;
num_lick_harvest(~validity) = [];
ind_lick_onset_str_harvest(~validity) = [];
ind_lick_onset_end_harvest(~validity) = [];

time_lick_onset_str_harvest = time_1K(ind_lick_onset_str_harvest);
time_lick_onset_end_harvest = time_1K(ind_lick_onset_end_harvest);
time_harvest_duration = (time_1K(ind_lick_onset_end_harvest) - ...
    time_1K(ind_lick_onset_str_harvest))';

% Determine direction of bouts
num_r = [];
num_l = [];
num_g = [];

for i = 1:length(ind_lick_onset_str_bout)
    inds_lick_onset = ind_lick_onset_str_bout(i):ind_lick_onset_end_bout(i);
    num_r = [num_r; length(find(ismember(inds_lick_onset, ind_lick_onset_r_reward)))];
    num_l = [num_l; length(find(ismember(inds_lick_onset, ind_lick_onset_l_reward)))];
    num_g = [num_g; length(find(ismember(inds_lick_onset, ind_lick_onset_grooming)))];

end
is_bout_r = (num_r > num_l) & (num_r > num_g) & num_r > 2;
is_bout_l = (num_l > num_r) & (num_l > num_g) & num_l > 2;
is_bout_g = (num_g > num_r) & (num_g > num_l) & num_g > 2;

% Determine direction of harvest
num_r = [];
num_l = [];

for i = 1:length(ind_lick_onset_str_harvest)
    inds_lick_onset = ind_lick_onset_str_harvest(i):ind_lick_onset_end_harvest(i);
    num_r = [num_r; length(find(ismember(inds_lick_onset, ind_lick_onset_r_reward)))];
    num_l = [num_l; length(find(ismember(inds_lick_onset, ind_lick_onset_l_reward)))];
end
is_harvest_r = num_r > num_l;
is_harvest_l = num_l > num_r;


DLC.IND.ind_lick_onset_str_harvest = ind_lick_onset_str_harvest;
DLC.IND.ind_lick_onset_end_harvest = ind_lick_onset_end_harvest;
DLC.IND.num_lick_harvest = num_lick_harvest;
DLC.TIME.time_lick_onset_str_harvest = time_lick_onset_str_harvest;
DLC.TIME.time_lick_onset_end_harvest = time_lick_onset_end_harvest;
DLC.TIME.time_harvest_duration = time_harvest_duration;
DLC.CLASS.is_bout_r = is_bout_r;
DLC.CLASS.is_harvest_r = is_harvest_r;
DLC.CLASS.is_bout_l = is_bout_l;
DLC.CLASS.is_harvest_l = is_harvest_l;
DLC.CLASS.is_bout_g = is_bout_g;
if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure;
    hold on;
    plot(time_1K,d_tip,'-k');
    plot(time_lick_onset_str_harvest,20,'*g');
    plot(time_lick_onset_end_harvest,20,'*r');
    ylim([0 25])
    xlabel('Time (s)');
    ylabel('Displacement (mm)');
    if ~isempty(time_lick_onset_r_reward)
        plot(time_lick_onset_r_reward,19,'.r');
    end
    if ~isempty(time_lick_onset_l_reward)
        plot(time_lick_onset_l_reward,19,'.b');
    end
    yyaxis right
    plot(time_1K,r_tube_food,'.-r');
    plot(time_1K,l_tube_food,'.-b');
    ylabel('Reward capacity (0:Empty | 1:Full)')
    title([ EXPERIMENT_PARAMS.file_name ': Harvest | ' num2str(sum(DLC.CLASS.is_bout_l)) ' left bouts | ' num2str(sum(DLC.CLASS.is_bout_r)) ' right bouts'], 'interpreter', 'none');
    ESN_Beautify_Plot(gcf, [20 10])
end

fprintf(' --> Completed. \n')
end

%% Function: QA on stationary markers
function [DLC, EXPERIMENT_PARAMS] = qa_stationary(DLC,EXPERIMENT_PARAMS,vid_obj)
FPS = EXPERIMENT_PARAMS.FPS;
data_ = table2array(DLC.data);

% video
frame = read(vid_obj, 1);
vid_height = vid_obj.height;
vid_width = vid_obj.width;

% nose markers
fprintf("\nCorrecting nose markers...");
r_nose_x = data_(:,9);
r_nose_y = data_(:,10);
l_nose_x = data_(:,11);
l_nose_y = data_(:,12);
% Eliminate erroneous tracking from nose markers
[r_nose_x_qa, r_nose_y_qa, ~] = find_erroneous_markers(r_nose_x, r_nose_y, [0 0 vid_width vid_height], "r_nose", frame);
[l_nose_x_qa, l_nose_y_qa, ~] = find_erroneous_markers(l_nose_x, l_nose_y, [0 0 vid_width vid_height], "l_nose", frame);
r_nose_x_qa = fillmissing(r_nose_x_qa, 'linear');
r_nose_y_qa = fillmissing(r_nose_y_qa, 'linear');
l_nose_x_qa = fillmissing(l_nose_x_qa, 'linear');
l_nose_y_qa = fillmissing(l_nose_y_qa, 'linear');
DLC.data.r_nose_x = r_nose_x_qa;
DLC.data.r_nose_y = r_nose_y_qa;
DLC.data.l_nose_x = l_nose_x_qa;
DLC.data.l_nose_y = l_nose_y_qa;

% tube markers
fprintf("\nCorrecting tube markers...");
r_tube_r_x = data_(:,17);
r_tube_r_y = data_(:,18);
r_tube_l_x = data_(:,19);
r_tube_l_y = data_(:,20);
l_tube_r_x = data_(:,21);
l_tube_r_y = data_(:,22);
l_tube_l_x = data_(:,23);
l_tube_l_y = data_(:,24);
% check if the tube markers are at a reasonable position
l_l_length = median(l_tube_l_y);
l_r_length = median(l_tube_r_y);
if  abs(l_l_length - l_r_length) >= (l_l_length/4) || abs(l_l_length - l_r_length) >= (l_r_length/4)
    % manual correction
    [l_tube_l_x_qa, l_tube_l_y_qa, ~] = manual_correction(l_tube_l_x, l_tube_l_y, "l\_tube\_l", frame);
    [l_tube_r_x_qa, l_tube_r_y_qa, ~] = manual_correction(l_tube_r_x, l_tube_r_y, "l\_tube\_r", frame);
else
    [l_tube_r_x_qa, l_tube_r_y_qa, ~] = find_erroneous_markers(l_tube_r_x, l_tube_r_y, [0 0 vid_width vid_height/2], "l_tube_r", frame);
    [l_tube_l_x_qa, l_tube_l_y_qa, ~] = find_erroneous_markers(l_tube_l_x, l_tube_l_y, [0 0 vid_width vid_height/2], "l_tube_l", frame);
end
    
r_l_length = vid_height - median(r_tube_l_y);
r_r_length = vid_height - median(r_tube_r_y);
if abs(r_l_length - r_r_length) >= (r_l_length/4) || abs(r_l_length - r_r_length) >= (r_r_length/4)
    % manual correction
    [r_tube_l_x_qa, r_tube_l_y_qa, ~] = manual_correction(r_tube_l_x, r_tube_l_y, "r_tube_l", frame);
    [r_tube_r_x_qa, r_tube_r_y_qa, ~] = manual_correction(r_tube_r_x, r_tube_r_y, "r_tube_r", frame);
else
    [r_tube_r_x_qa, r_tube_r_y_qa, ~] = find_erroneous_markers(r_tube_r_x, r_tube_r_y, [0 vid_height/2 vid_width vid_height], "r_tube_r", frame);
    [r_tube_l_x_qa, r_tube_l_y_qa, ~] = find_erroneous_markers(r_tube_l_x, r_tube_l_y, [0 vid_height/2 vid_width vid_height], "r_tube_l", frame);
end
% Perform linear interpolation
l_tube_r_x_qa = fillmissing(l_tube_r_x_qa, 'linear');
l_tube_r_y_qa = fillmissing(l_tube_r_y_qa, 'linear');
l_tube_l_x_qa = fillmissing(l_tube_l_x_qa, 'linear');
l_tube_l_y_qa = fillmissing(l_tube_l_y_qa, 'linear');
r_tube_r_x_qa = fillmissing(r_tube_r_x_qa, 'linear');
r_tube_r_y_qa = fillmissing(r_tube_r_y_qa, 'linear');
r_tube_l_x_qa = fillmissing(r_tube_l_x_qa, 'linear');
r_tube_l_y_qa = fillmissing(r_tube_l_y_qa, 'linear');

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure
    hold on;
    set(gca,'Ydir','reverse')
    title("Stationary markers (Blue: DLC | Green: QA)")
    plot(l_tube_r_x, l_tube_r_y, 'b*')
    plot(l_tube_l_x, l_tube_l_y, 'b*')
    plot(r_tube_r_x, r_tube_r_y, 'b*')
    plot(r_tube_l_x, r_tube_l_y, 'b*')
    plot(r_nose_x, r_nose_y, 'b*')
    plot(l_nose_x, l_nose_y, 'b*')
    plot(l_tube_r_x_qa, l_tube_r_y_qa, 'g*')
    plot(l_tube_l_x_qa, l_tube_l_y_qa, 'g*')
    plot(r_tube_r_x_qa, r_tube_r_y_qa, 'g*')
    plot(r_tube_l_x_qa, r_tube_l_y_qa, 'g*')
    plot(r_nose_x_qa, r_nose_y_qa, 'g*')
    plot(l_nose_x_qa, l_nose_y_qa, 'g*')
end

DLC.data.l_tube_r_x = l_tube_r_x_qa;
DLC.data.l_tube_r_y = l_tube_r_y_qa;
DLC.data.l_tube_l_x = l_tube_l_x_qa;
DLC.data.l_tube_l_y = l_tube_l_y_qa;
DLC.data.r_tube_r_x = r_tube_r_x_qa;
DLC.data.r_tube_r_y = r_tube_r_y_qa;
DLC.data.r_tube_l_x = r_tube_l_x_qa;
DLC.data.r_tube_l_y = r_tube_l_y_qa;
fprintf(' --> Completed. \n')
    
% helper function to detect shift, outlier, and outside expected boundary
% rect = [x_lower, y_lower, x_upper, y_upper]      
function [x, y, ind_DLC_ERR] = find_erroneous_markers(x, y, rect, type, frame)
    mov_x = movmean(x, 1000);
    mov_y = movmean(y, 1000);
    
    if ~isempty(regexp(type, regexptranslate('wildcard', '*nose'))) || (max(mov_x) - min(mov_x) < 5 && max(mov_y) - min(mov_y) < 5)
        ind_DLC_ERR = [];
        % isoutlier(median) method
        x_err = isoutlier(x, "median");
        y_err = isoutlier(y, "median");
        ind_DLC_ERR = unique([find(x_err);find(y_err)]);

        % boundary method
        for i = 1:length(x)
            if x(i) < rect(1) || x(i) > rect(3) || y(i) < rect(2) || y(i) > rect(4)
                ind_DLC_ERR = [ind_DLC_ERR; i];
            end
        end
        ind_DLC_ERR = unique(ind_DLC_ERR);
        x(ind_DLC_ERR) = NaN;
        y(ind_DLC_ERR) = NaN;
    else
        % manual correction
        [x, y, ind_DLC_ERR] = manual_correction(x, y, type, frame);
    end

    % prepare data for interpolation
    if isnan(x(length(x)))
        x(length(x)) = x(find(~isnan(x), 1, "last"));
        y(length(y)) = y(find(~isnan(y), 1, "last"));
    end
    if isnan(x(1))
        x(1) = x(find(~isnan(x), 1, "first"));
        y(1) = y(find(~isnan(y), 1, "first"));
    end

end
function [x, y, ind_DLC_ERR] = manual_correction(x, y, type, frame)
    fig = figure;
    fig.WindowState = 'Maximize';
    imshow(frame)
    title(type, 'Interpreter', 'none')
    x_ref = [];
    y_ref = [];
    while isempty(x_ref)
        [x_ref, y_ref] = getpts;
        if length(x_ref) > 1
            x_ref = x_ref(end);
            y_ref = y_ref (end);
        end
    end
    ind_DLC_ERR = [];
    ind_DLC_ERR = [ind_DLC_ERR; find(abs(x - x_ref)>5)];
    ind_DLC_ERR = [ind_DLC_ERR; find(abs(y - y_ref)>5)];
    ind_DLC_ERR = unique(ind_DLC_ERR);
    if length(ind_DLC_ERR) == length(x)
        x = x_ref*ones(size(x));
        y = y_ref*ones(size(y));
    else
        x(ind_DLC_ERR) = NaN;
        y(ind_DLC_ERR) = NaN;
    end
    close(fig)
end
function [x, y] = kalman_filtering(x, y, FPS)
    %% define main variables for kalman filter
    dt = 1 / FPS;
    u = 0; % acceleration magnitude
    Q = [x(1); y(1); 0; 0]; % initialized state--it has four components: [positionX; positionY; velocityX; velocityY] of the hexbug
    Q_estimate = Q;  % estimate of initial location (what we are updating)
    accel_noise_mag = 2; % process noise: the variability in how fast the tracking is speeding up (stdv of acceleration: meters/sec^2)
    n_x = 10;  % measurement noise in the x direction
    n_y = 10;  % measurement noise in the y direction
    Ez = [n_x 0; 0 n_y];
    Ex = [dt^4/4 0 dt^3/2 0; ...
        0 dt^4/4 0 dt^3/2; ...
        dt^3/2 0 dt^2 0; ...
        0 dt^3/2 0 dt^2].*accel_noise_mag^2; % Ex convert the process noise (stdv) into covariance matrix
    P = Ex; % estimate of initial position variance (covariance matrix)
    %% Define update equations in 2D (Coefficent matrices): A physics based model for where we expect the marker to be [state transition (state + velocity)] + [input control (acceleration)]
    A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; % state update matrix
    B = [(dt^2/2); (dt^2/2); dt; dt]; % input control matrix
    C = [1 0 0 0; 0 1 0 0];  % measurement function C (applied to the state estimate Q to get the expected next new measurement)
    %% Initize result variables
    Q_loc_meas = []; % the path extracted by the tracking algo
    %% Initize estimation variables
    p_estimate = []; % position estimate
    v_estimate = []; % velocity estimate
    predic_state = [];
    predic_var = [];
    for t = 1:size(x)
        % load the given tracking
        Q_loc_meas(:,t) = [x(t); y(t)];
        %% do the kalman filter
        % Predict next state of the with the last state and predicted motion.
        Q_estimate = A * Q_estimate + B * u;
        predic_state = [predic_state; Q_estimate(1)] ;
        % Predict next covariance
        P = A * P * A' + Ex;
        predic_var = [predic_var; P] ;
        % Predicted measurement covariance
        % Kalman Gain
        K = P*C'*inv(C*P*C'+Ez);
        % Update the state estimate.
        if ~isnan(Q_loc_meas(:,t))
            Q_estimate = Q_estimate + K * (Q_loc_meas(:,t) - C * Q_estimate);
        end
        % update covariance estimation.
        P = (eye(4)-K*C)*P;
        %% Store data
        p_estimate = [p_estimate; [Q_estimate(1) Q_estimate(2)]];
        v_estimate = [v_estimate; [Q_estimate(3) Q_estimate(4)]];
    end
    x = p_estimate(:, 1);
    y = p_estimate(:, 2);
end
end

%% Function: QA on food and tongue
function [DLC, EXPERIMENT_PARAMS] = qa_food_and_tongue(DLC,EXPERIMENT_PARAMS,vid_obj,light_condition)
data_ = table2array(DLC.data);
path_to_raw = [EXPERIMENT_PARAMS.mat_PathName '..' filesep '..' filesep '..' filesep 'raw_data'];
% video
vid_height = vid_obj.height;
vid_width = vid_obj.width;

r_nose_x = DLC.data.r_nose_x;
r_nose_y = DLC.data.r_nose_y;
l_nose_x = DLC.data.l_nose_x;
l_nose_y = DLC.data.l_nose_y;
if light_condition == 1
    fprintf("Checking light condition in video ...");

    % 1) Lighting condition (for both tongue and food markers)
    avg_pi = [];
    mean_nose_x = round((mean(r_nose_x) + mean(l_nose_x))/2);
    mean_nose_y = round((mean(r_nose_y) + mean(l_nose_y))/2);

    for i = 1:length(r_nose_x)
        frame = read(vid_obj, i);
        avg_pi = [avg_pi; frame(mean_nose_y, mean_nose_x, :)];
    end
    % Find dark frames
    avg_pi = mean(avg_pi, 3);
    TF = isoutlier(avg_pi,'movmedian',length(avg_pi)/2);
    dark = [];
    for i = 1:length(TF)
        if TF(i)
            if avg_pi(i) < mean(avg_pi)
                dark = [dark; avg_pi(i)];
            else
                TF(i) = 0;
            end
        end
    end

    y0 = (mean(r_nose_y) + mean(l_nose_y ))/2;
    y_lower_lim = 10;
    y_upper_lim = vid_height - 10;
    x_upper_lim = vid_width;
    fprintf(' --> Completed. \n')
else
    TF = [];
    mean_nose_x = round((mean(r_nose_x) + mean(l_nose_x))/2);
    y0 = (mean(r_nose_y) + mean(l_nose_y ))/2;
    y_lower_lim = 10;
    y_upper_lim = vid_height - 10;
    x_upper_lim = vid_width;
end

% QA for food markers
fprintf("Correcting food markers...");
r_food_x = data_(:,13);
r_food_y = data_(:,14);
l_food_x = data_(:,15);
l_food_y = data_(:,16);
r_tube_r_x = DLC.data.r_tube_r_x;
r_tube_r_y = DLC.data.r_tube_r_y;
r_tube_l_x = DLC.data.r_tube_l_x;
r_tube_l_y = DLC.data.r_tube_l_y;
l_tube_r_x = DLC.data.l_tube_r_x;
l_tube_r_y = DLC.data.l_tube_r_y;
l_tube_l_x = DLC.data.l_tube_l_x;
l_tube_l_y = DLC.data.l_tube_l_y;

% 2) Boundary
TF_left = TF;
TF_right = TF;
mean_r_tube_l_x = mean(r_tube_l_x);
mean_r_tube_r_x = mean(r_tube_r_x);
mean_l_tube_l_x = mean(l_tube_l_x);
mean_l_tube_r_x = mean(l_tube_r_x);
for i = 1:length(l_food_y)
    if l_food_y(i) > y0 || l_food_x(i) < mean_l_tube_l_x ||  l_food_x(i) > mean_l_tube_r_x || l_food_y(i) <= y_lower_lim
        TF_left(i) = 1;
    end
    if r_food_y(i) < y0 || r_food_x(i) < mean_r_tube_l_x ||  r_food_x(i) > mean_r_tube_r_x || r_food_y(i) >= y_upper_lim
        TF_right(i) = 1;
    end
end
l_food_x_qa = l_food_x;
l_food_y_qa = l_food_y;
r_food_x_qa = r_food_x;
r_food_y_qa = r_food_y;
% l_food
l_food_incorrect_inds = find(TF_left);
l_food_x_qa(l_food_incorrect_inds) = NaN;
l_food_y_qa(l_food_incorrect_inds) = NaN;
if isnan(l_food_x_qa(1))
    l_food_x_qa(1) = l_food_x_qa(find(~TF_left, 1));
    l_food_y_qa(1) = l_food_y_qa(find(~TF_left, 1));
end
if isnan(l_food_x_qa(length(l_food_x_qa)))
    l_food_x_qa(length(l_food_x_qa)) = l_food_x_qa(find(~TF_left, 1, 'last'));
    l_food_y_qa(length(l_food_x_qa)) = l_food_y_qa(find(~TF_left, 1, 'last'));
end
% r_food
r_food_incorrect_inds = find(TF_right);
r_food_x_qa(r_food_incorrect_inds) = NaN;
r_food_y_qa(r_food_incorrect_inds) = NaN;
if isnan(r_food_x_qa(1))
    r_food_x_qa(1) = r_food_x_qa(find(~TF_right, 1));
    r_food_y_qa(1) = r_food_y_qa(find(~TF_right, 1));
end
if isnan(r_food_x_qa(length(r_food_x_qa)))
    r_food_x_qa(length(r_food_x_qa)) = r_food_x_qa(find(~TF_right, 1, 'last'));
    r_food_y_qa(length(r_food_x_qa)) = r_food_y_qa(find(~TF_right, 1, 'last'));
end
% Fixing all incorrection food markers with makima interpolation
l_food_x_qa = fillmissing(l_food_x_qa, 'makima');
l_food_y_qa = fillmissing(l_food_y_qa, 'makima');
r_food_x_qa = fillmissing(r_food_x_qa, 'makima');
r_food_y_qa = fillmissing(r_food_y_qa, 'makima');

% 3) Detect jumping markers
% l_food
max_diff = mean([median(l_tube_r_y), median(l_tube_l_y)]) / 4;
dy = diff(l_food_y_qa);
jump_inds = zeros(size(l_food_y_qa));
rise = NaN;
for i = 1:length(dy)
    if isnan(rise)
        % find jumping start
        if abs(dy(i)) > max_diff
            if dy(i) > 0
                rise = true;
            else
                rise = false;
            end
                jump_inds(i+1) = 1;
        end
    else
        % find jumping end
        if abs(dy(i)) > max_diff
            if (rise && dy(i) < 0) || (~rise && dy(i) > 0)
                jump_inds(i+1) = -1;
                rise = NaN;
            end
        end
    end
end
jumping = false;
for i = 1:length(jump_inds)
    window_ = 300;
    if i+window_ > length(jump_inds)
        window_ = window_ - (length(jump_inds) - i + window_);
    end
    if jump_inds(i) == 1 && ~isempty(find(jump_inds(i:i+window_)==-1, 1))
            jumping = true;
    end
    if jump_inds(i) == -1
        jumping = false;
    end
    if jumping
        l_food_y_qa(i) = NaN;
        l_food_x_qa(i) = NaN;
        l_food_incorrect_inds = [l_food_incorrect_inds i];
    end
end

% r_food
max_diff = mean([median(r_tube_r_y), median(r_tube_l_y)]) / 4;
dy = diff(r_food_y_qa);
jump_inds = zeros(size(r_food_y_qa));
rise = NaN;
for i = 1:length(dy)
    if isnan(rise)
        % find jumping start
        if abs(dy(i)) > max_diff
            if dy(i) > 0
                rise = true;
            else
                rise = false;
            end
                jump_inds(i+1) = 1;
        end
    else
        % find jumping end
        if abs(dy(i)) > max_diff
            if (rise && dy(i) < 0) || (~rise && dy(i) > 0)
                jump_inds(i+1) = -1;
                rise = NaN;
            end
        end
    end
end
jumping = false;
for i = 1:length(jump_inds)
    window_ = 300;
    if i+window_ > length(jump_inds)
        window_ = window_ - (length(jump_inds) - i + window_);
    end
    if jump_inds(i) == 1 && ~isempty(find(jump_inds(i:i+window_)==-1, 1))
            jumping = true;
    end
    if jump_inds(i) == -1
        jumping = false;
    end
    if jumping
        r_food_y_qa(i) = NaN;
        r_food_x_qa(i) = NaN;
        r_food_incorrect_inds = [r_food_incorrect_inds i];
    end
end

l_food_x_qa = fillmissing(l_food_x_qa, 'makima');
l_food_y_qa = fillmissing(l_food_y_qa, 'makima');
r_food_x_qa = fillmissing(r_food_x_qa, 'makima');
r_food_y_qa = fillmissing(r_food_y_qa, 'makima');

l_food_incorrect_inds = unique(sort(l_food_incorrect_inds));
r_food_incorrect_inds = unique(sort(r_food_incorrect_inds));
if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure
    hold on
    set(gca,'Ydir','reverse')
    title("Food markers (Blue/Red: QA-left/right | *: DLC)")
    plot(l_food_incorrect_inds, l_food_y(l_food_incorrect_inds), 'b*')
    plot(r_food_incorrect_inds, r_food_y(r_food_incorrect_inds), 'r*')
    plot(l_food_y_qa, 'b-')
    plot(r_food_y_qa, 'r-')
end
% Save variables
DLC.data.l_food_x = l_food_x_qa;
DLC.data.l_food_y = l_food_y_qa;
DLC.data.r_food_x = r_food_x_qa;
DLC.data.r_food_y = r_food_y_qa;
fprintf(' --> Completed. \n')

% QA for tongue markers
fprintf("Correcting tongue markers...");
tip_tongue_x = data_(:,1);
tip_tongue_y = data_(:,2);
r_tongue_x = data_(:,3);
r_tongue_y = data_(:,4);
l_tongue_x = data_(:,5);
l_tongue_y = data_(:,6);
mid_tongue_x = data_(:,7);
mid_tongue_y = data_(:,8);

TF_tip = TF;
TF_l = TF;
TF_r = TF;
TF_mid = TF;
% 1) Boundary
for i = 2:length(tip_tongue_y)
    if tip_tongue_y(i) <= y_lower_lim || tip_tongue_y(i) >= y_upper_lim || abs(angle(i) - angle(i-1)) > 80
        TF_tip(i) = 1;
    end
    if l_tongue_y(i) <= y_lower_lim || l_tongue_y(i) >= y_upper_lim
        TF_l(i) = 1;
    end
    if r_tongue_y(i) <= y_lower_lim || r_tongue_y(i) >= y_upper_lim
        TF_r(i) = 1;
    end
    if mid_tongue_y(i) <= y_lower_lim || mid_tongue_x(i) >= x_upper_lim || mid_tongue_y(i) >= y_upper_lim
        TF_mid(i) = 1;
    end
end

% 2) Abrupt change in distance between the 4 markers
% 3) Abnormal shape (sum of angles not equal to 360)
tl_lengths = zeros(size(tip_tongue_x));
lm_lengths = zeros(size(tip_tongue_x));
mr_lengths = zeros(size(tip_tongue_x));
rt_lengths = zeros(size(tip_tongue_x));
tlm_angles = zeros(size(tip_tongue_x));
lmr_angles = zeros(size(tip_tongue_x));
mrt_angles = zeros(size(tip_tongue_x));
rtl_angles = zeros(size(tip_tongue_x));
for i = 1:length(tip_tongue_x)
    tip = [tip_tongue_x(i), tip_tongue_y(i)];
    mid = [mid_tongue_x(i), mid_tongue_y(i)];
    l = [l_tongue_x(i), l_tongue_y(i)];
    r = [r_tongue_x(i), r_tongue_y(i)];
    
    % populate length
    tl_lengths(i) = norm(tip-l);
    lm_lengths(i) = norm(l-mid);
    mr_lengths(i) = norm(mid-r);
    rt_lengths(i) = norm(r-tip);
    % populate angles
    tlm_angles(i) = find_angle(tip, l, mid);
    lmr_angles(i) = find_angle(l, mid, r);
    mrt_angles(i) = find_angle(mid, r, tip);
    rtl_angles(i) = find_angle(r, tip, l);
end

inds = [];
[~, ind] = findpeaks(diff(tl_lengths), "MinPeakProminence", 30);
inds = [inds; ind];
[~, ind] = findpeaks(diff(lm_lengths), "MinPeakProminence", 30);
inds = [inds; ind];
[~, ind] = findpeaks(diff(mr_lengths), "MinPeakProminence", 30);
inds = [inds; ind];
[~, ind] = findpeaks(diff(rt_lengths), "MinPeakProminence", 30);
inds = [inds; ind];
angles = rtl_angles + mrt_angles + lmr_angles + tlm_angles;
inds = [inds; find(angles < 359 | angles > 361)];

% Fixing all incorrect tongue markers using makima interpolation
tip_tongue_x_qa = tip_tongue_x;
tip_tongue_y_qa = tip_tongue_y;
l_tongue_x_qa = l_tongue_x;
l_tongue_y_qa = l_tongue_y;
mid_tongue_x_qa = mid_tongue_x;
mid_tongue_y_qa = mid_tongue_y;
r_tongue_x_qa = r_tongue_x;
r_tongue_y_qa = r_tongue_y;

tip_incorrect_inds = unique(sort([find(TF_tip)'; inds]));
l_incorrect_inds = unique(sort([find(TF_l)'; inds]));
r_incorrect_inds = unique(sort([find(TF_r)'; inds]));
mid_incorrect_inds = unique(sort([find(TF_mid)'; inds]));
tip_tongue_x_qa(tip_incorrect_inds) = NaN;
tip_tongue_y_qa(tip_incorrect_inds) = NaN;
l_tongue_x_qa(l_incorrect_inds) = NaN;
l_tongue_y_qa(l_incorrect_inds) = NaN;
r_tongue_x_qa(r_incorrect_inds) = NaN;
r_tongue_y_qa(r_incorrect_inds) = NaN;
mid_tongue_x_qa(mid_incorrect_inds) = NaN;
mid_tongue_y_qa(mid_incorrect_inds) = NaN;
tip_tongue_x_qa = fillmissing(tip_tongue_x_qa, 'makima');
tip_tongue_y_qa = fillmissing(tip_tongue_y_qa, 'makima');
l_tongue_x_qa = fillmissing(l_tongue_x_qa, 'makima');
l_tongue_y_qa = fillmissing(l_tongue_y_qa, 'makima');
r_tongue_x_qa = fillmissing(r_tongue_x_qa, 'makima');
r_tongue_y_qa = fillmissing(r_tongue_y_qa, 'makima');
mid_tongue_x_qa = fillmissing(mid_tongue_x_qa, 'makima');
mid_tongue_y_qa = fillmissing(mid_tongue_y_qa, 'makima');

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    tip_dist = sqrt((tip_tongue_x).^2 + (tip_tongue_y).^2);
    l_dist = sqrt((l_tongue_x).^2 + (l_tongue_y).^2);
    r_dist = sqrt((r_tongue_x).^2 + (r_tongue_y).^2);
    mid_dist = sqrt((mid_tongue_x).^2 + (mid_tongue_y).^2);

    tip_dist_qa = sqrt((tip_tongue_x_qa).^2 + (tip_tongue_y_qa).^2);
    l_dist_qa = sqrt((l_tongue_x_qa).^2 + (l_tongue_y_qa).^2);
    r_dist_qa = sqrt((r_tongue_x_qa).^2 + (r_tongue_y_qa).^2);
    mid_dist_qa = sqrt((mid_tongue_x_qa).^2 + (mid_tongue_y_qa).^2);

    figure
    tip_inds = find(abs(tip_dist(tip_incorrect_inds)-tip_dist_qa(tip_incorrect_inds)) > 5);
    l_inds = find(abs(l_dist(l_incorrect_inds)-l_dist_qa(l_incorrect_inds)) > 5);
    r_inds = find(abs(r_dist(r_incorrect_inds)-r_dist_qa(r_incorrect_inds)) > 5);
    mid_inds = find(abs(mid_dist(mid_incorrect_inds)-mid_dist_qa(mid_incorrect_inds)) > 5);
    tip_incorrect_inds = tip_incorrect_inds(tip_inds);
    l_incorrect_inds = l_incorrect_inds(l_inds);
    r_incorrect_inds = r_incorrect_inds(r_inds);
    mid_incorrect_inds = mid_incorrect_inds(mid_inds);

    subplot(4, 2, 1)
    plot(tip_tongue_x_qa, '-')
    hold on
    plot(tip_incorrect_inds, tip_tongue_x(tip_incorrect_inds), '*')
    title("Tip tongue x QA (*: DLC)")

    subplot(4, 2, 2)
    plot(tip_tongue_y_qa, '-')
    hold on
    plot(tip_incorrect_inds, tip_tongue_y(tip_incorrect_inds), '*')
    title("Tip tongue y QA (*: DLC)")

    subplot(4, 2, 3)
    plot(l_tongue_x_qa, '-')
    hold on
    plot(l_incorrect_inds, l_tongue_x(l_incorrect_inds), '*')
    title("Left tongue x QA (*: DLC)")

    subplot(4, 2, 4)
    plot(l_tongue_y_qa, '-')
    hold on
    plot(l_incorrect_inds, l_tongue_y(l_incorrect_inds), '*')
    title("Left tongue y QA (*: DLC)")

    subplot(4, 2, 5)
    plot(r_tongue_x_qa, '-')
    hold on
    plot(r_incorrect_inds, l_tongue_x(r_incorrect_inds), '*')
    title("Right tongue x QA (*: DLC)")

    subplot(4, 2, 6)
    plot(r_tongue_y_qa, '-')
    hold on
    plot(r_incorrect_inds, r_tongue_y(r_incorrect_inds), '*')
    title("Right tongue y QA (*: DLC)")

    subplot(4, 2, 7)
    plot(mid_tongue_x_qa, '-')
    hold on
    plot(mid_incorrect_inds, mid_tongue_x(mid_incorrect_inds), '*')
    title("Mid tongue x QA (*: DLC)")

    subplot(4, 2, 8)
    plot(mid_tongue_y_qa, '-')
    hold on
    plot(mid_incorrect_inds, mid_tongue_y(mid_incorrect_inds), '*')
    title("Mid tongue y QA (*: DLC)")

%     inds = sort(unique([tip_incorrect_inds; l_incorrect_inds; r_incorrect_inds; mid_incorrect_inds]));
%     for i = 1:length(inds)
%         imwrite(read(vid_obj, inds(i)), strcat(num2str(inds(i)), '.png'))
%     end
end
DLC.data.tip_tongue_x = tip_tongue_x_qa;
DLC.data.tip_tongue_y = tip_tongue_y_qa;
DLC.data.r_tongue_x = r_tongue_x_qa;
DLC.data.r_tongue_y = r_tongue_y_qa;
DLC.data.l_tongue_x = l_tongue_x_qa;
DLC.data.l_tongue_y = l_tongue_y_qa;
DLC.data.mid_tongue_x = mid_tongue_x_qa;
DLC.data.mid_tongue_y = mid_tongue_y_qa;
fprintf(' --> Completed. \n')

function angle = find_angle(a, b, c)
x1 = b(1);
y1 = b(2);
x2 = a(1);
y2 = a(2);
x3 = c(1);
y3 = c(2);
angle = atan2(abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)), ...
            (x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)) * 180/pi;
end
end

%% Function: Update .mat and csv after QA
function save_qa(DLC, EXPERIMENT_PARAMS)
% update csv file
update_csv = 0;
if update_csv == 1
    csv_file = dir(strcat(DLC.FILE.path_to_analyzed, '*.csv'));
    [csv_table, t, ~] = xlsread(strcat(DLC.FILE.path_to_analyzed, csv_file.name));
    csv_table(:,2) = DLC.data.tip_tongue_x;
    csv_table(:,3) = DLC.data.tip_tongue_y;
    csv_table(:,5) = DLC.data.r_tongue_x;
    csv_table(:,6) = DLC.data.r_tongue_y;
    csv_table(:,8) = DLC.data.l_tongue_x;
    csv_table(:,9) = DLC.data.l_tongue_y;
    csv_table(:,11) = DLC.data.mid_tongue_x;
    csv_table(:,12) = DLC.data.mid_tongue_y;
    csv_table(:,14) = DLC.data.r_nose_x;
    csv_table(:,15) = DLC.data.r_nose_y;
    csv_table(:,17) = DLC.data.l_nose_x;
    csv_table(:,18) = DLC.data.l_nose_y;
    csv_table(:,20) = DLC.data.r_food_x;
    csv_table(:,21) = DLC.data.r_food_y;
    csv_table(:,23) = DLC.data.l_food_x;
    csv_table(:,24) = DLC.data.l_food_y;
    csv_table(:,26) = DLC.data.r_tube_r_x;
    csv_table(:,27) = DLC.data.r_tube_r_y;
    csv_table(:,29) = DLC.data.r_tube_l_x;
    csv_table(:,30) = DLC.data.r_tube_l_y;
    csv_table(:,32) = DLC.data.l_tube_r_x;
    csv_table(:,33) = DLC.data.l_tube_r_y;
    csv_table(:,35) = DLC.data.l_tube_l_x;
    csv_table(:,36) = DLC.data.l_tube_l_y;
    writecell(t, strcat(DLC.FILE.path_to_analyzed, csv_file.name));
    writematrix(csv_table,strcat(DLC.FILE.path_to_analyzed, csv_file.name), 'WriteMode', 'append');
end
% update .mat file
update_mat = 1;
if update_mat == 1
    data = DLC.data;
    save(strcat(DLC.FILE.path_to_analyzed, EXPERIMENT_PARAMS.mat_FileName), 'data')
end
end


%% Function: PGH_plot_lick_sorter_averages
function PGH_plot_lick_sorter_averages(DLC,format,direction)
is_IT = DLC.CLASS.(['is_' direction '_reward_inner_tube_lick']) | DLC.CLASS.(['is_' direction '_noreward_inner_tube_lick']);
is_OE = DLC.CLASS.(['is_' direction '_reward_outer_edge_lick']) | DLC.CLASS.(['is_' direction '_noreward_outer_edge_lick']);
is_UT = DLC.CLASS.(['is_' direction '_reward_under_tube_lick']) | DLC.CLASS.(['is_' direction '_noreward_under_tube_lick']);

figure
subplot(3,4,1)
hold on
line([nanmean(DLC.KINEMATIC.r_tube_r_x_dmax(:,is_IT)) nanmean(DLC.KINEMATIC.r_tube_r_x_dmax(:,is_IT))] ...
    ,[nanmean(DLC.KINEMATIC.r_tube_r_y_dmax(:,is_IT)) 20], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.r_tube_l_x_dmax(:,is_IT)) nanmean(DLC.KINEMATIC.r_tube_l_x_dmax(:,is_IT))] ...
    ,[nanmean(DLC.KINEMATIC.r_tube_l_y_dmax(:,is_IT)) 20], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.r_tube_r_x_dmax(:,is_IT)) nanmean(DLC.KINEMATIC.r_tube_l_x_dmax(:,is_IT))] ...
    , [nanmean(DLC.KINEMATIC.r_tube_r_y_dmax(:,is_IT)) nanmean(DLC.KINEMATIC.r_tube_l_y_dmax(:,is_IT))], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.l_tube_r_x_dmax(:,is_IT)) nanmean(DLC.KINEMATIC.l_tube_r_x_dmax(:,is_IT))] ...
    ,[-20 nanmean(DLC.KINEMATIC.l_tube_r_y_dmax(:,is_IT))], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.l_tube_l_x_dmax(:,is_IT)) nanmean(DLC.KINEMATIC.l_tube_l_x_dmax(:,is_IT))] ...
    ,[-20 nanmean(DLC.KINEMATIC.l_tube_l_y_dmax(:,is_IT))], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.l_tube_r_x_dmax(:,is_IT)) nanmean(DLC.KINEMATIC.l_tube_l_x_dmax(:,is_IT))] ...
    , [nanmean(DLC.KINEMATIC.l_tube_r_y_dmax(:,is_IT)) nanmean(DLC.KINEMATIC.l_tube_l_y_dmax(:,is_IT))], 'Color', 'black');

plot(DLC.KINEMATIC.tip_tongue_x_dmax(is_IT), DLC.KINEMATIC.tip_tongue_y_dmax(is_IT), '.b','MarkerSize',20)
plot(DLC.KINEMATIC.r_food_x_dmax(is_IT), DLC.KINEMATIC.r_food_y_dmax(is_IT), '.y','MarkerSize',20)
plot(DLC.KINEMATIC.l_food_x_dmax(is_IT), DLC.KINEMATIC.l_food_y_dmax(is_IT), '.y','MarkerSize',20)
plot(DLC.KINEMATIC.l_nose_x_dmax(is_IT), DLC.KINEMATIC.l_nose_y_dmax(is_IT), '.k','MarkerSize',20)
plot(DLC.KINEMATIC.r_nose_x_dmax(is_IT), DLC.KINEMATIC.r_nose_y_dmax(is_IT), '.k','MarkerSize',20)
plot(DLC.KINEMATIC.r_tube_r_x_dmax(is_IT), DLC.KINEMATIC.r_tube_r_y_dmax(is_IT), '.g','MarkerSize',20)
plot(DLC.KINEMATIC.r_tube_l_x_dmax(is_IT), DLC.KINEMATIC.r_tube_l_y_dmax(is_IT), '.g','MarkerSize',20)
plot(DLC.KINEMATIC.l_tube_r_x_dmax(is_IT), DLC.KINEMATIC.l_tube_r_y_dmax(is_IT), '.g','MarkerSize',20)
plot(DLC.KINEMATIC.l_tube_l_x_dmax(is_IT), DLC.KINEMATIC.l_tube_l_y_dmax(is_IT), '.g','MarkerSize',20)

% plot(DLC.KINEMATIC.tip_tongue_x_lick(is_IT,:), DLC.KINEMATIC.tip_tongue_y_lick(is_IT,:), 'k','LineWidth', 0.1)
% plot(nanmean(DLC.KINEMATIC.tip_tongue_x_lick(is_IT,:)), nanmean(DLC.KINEMATIC.tip_tongue_y_lick(is_IT,:)), 'k','LineWidth', 2)

mean_x_ = nanmean(DLC.KINEMATIC.tip_tongue_x_lick(is_IT,:));
sem_x_ = nanstd(DLC.KINEMATIC.tip_tongue_x_lick(is_IT,:))/sqrt(sum(is_IT));
mean_y_ = nanmean(DLC.KINEMATIC.tip_tongue_y_lick(is_IT,:));
sem_y_ = nanstd(DLC.KINEMATIC.tip_tongue_y_lick(is_IT,:))/sqrt(sum(is_IT));
% shade(mean_x_ + sem_x_, mean_y_ + sem_y_, 'r',  mean_x_ - sem_x_, mean_y_ - sem_y_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_x_+sem_x_, mean_y_+sem_y_, 'r','LineWidth', 1)
plot(mean_x_-sem_x_, mean_y_-sem_y_, 'r','LineWidth', 1)
plot(mean_x_, mean_y_, 'k','LineWidth', 2)

title(['Inner tube lick (mean+/-sem), n = ' num2str(sum(is_IT))], 'interpret', 'none');
xlim([-10, 20]);
ylim([-20, 20]);
xlabel('x (mm)')
xlabel('y (mm)')
set(gca, 'YDir','reverse')

subplot(3,4,2)
hold on
mean_ = nanmean(DLC.KINEMATIC.d_lick(is_IT,:));
sem_ = nanstd(DLC.KINEMATIC.d_lick(is_IT,:))/sqrt(sum(is_IT));
shade(1:length(mean_), mean_ + sem_, 'r', 1:length(mean_), mean_ - sem_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_, 'k', 'LineWidth', 1)

title(['Displacement'], 'interpret', 'none');
xlim([-50, 550]);
ylim([0, inf]);
xlabel('time (ms)')
ylabel('displacement (mm)')

subplot(3,4,3)
hold on
mean_ = nanmean(DLC.KINEMATIC.v_lick(is_IT,:));
sem_ = nanstd(DLC.KINEMATIC.v_lick(is_IT,:))/sqrt(sum(is_IT));
shade(1:length(mean_), mean_ + sem_, 'r', 1:length(mean_), mean_ - sem_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_, 'k', 'LineWidth', 1)

title(['Velocity'], 'interpret', 'none');
xlim([-50, 550]);
ylim([-inf, inf]);
xlabel('time (ms)')
ylabel('velocity (mm/s)')

subplot(3,4,4)
hold on
mean_ = nanmean(DLC.KINEMATIC.angle_lick(is_IT,:));
sem_ = nanstd(DLC.KINEMATIC.angle_lick(is_IT,:))/sqrt(sum(is_IT));
shade(1:length(mean_), mean_ + sem_, 'r', 1:length(mean_), mean_ - sem_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_, 'k', 'LineWidth', 1)

title(['Angle'], 'interpret', 'none');
xlim([-50, 550]);
ylim([-inf, inf]);
xlabel('time (ms)')
ylabel('angle (deg)')

subplot(3,4,5)
hold on
line([nanmean(DLC.KINEMATIC.r_tube_r_x_dmax(:,is_OE)) nanmean(DLC.KINEMATIC.r_tube_r_x_dmax(:,is_OE))] ...
    ,[nanmean(DLC.KINEMATIC.r_tube_r_y_dmax(:,is_OE)) 20], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.r_tube_l_x_dmax(:,is_OE)) nanmean(DLC.KINEMATIC.r_tube_l_x_dmax(:,is_OE))] ...
    ,[nanmean(DLC.KINEMATIC.r_tube_l_y_dmax(:,is_OE)) 20], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.r_tube_r_x_dmax(:,is_OE)) nanmean(DLC.KINEMATIC.r_tube_l_x_dmax(:,is_OE))] ...
    , [nanmean(DLC.KINEMATIC.r_tube_r_y_dmax(:,is_OE)) nanmean(DLC.KINEMATIC.r_tube_l_y_dmax(:,is_OE))], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.l_tube_r_x_dmax(:,is_OE)) nanmean(DLC.KINEMATIC.l_tube_r_x_dmax(:,is_OE))] ...
    ,[-20 nanmean(DLC.KINEMATIC.l_tube_r_y_dmax(:,is_OE))], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.l_tube_l_x_dmax(:,is_OE)) nanmean(DLC.KINEMATIC.l_tube_l_x_dmax(:,is_OE))] ...
    ,[-20 nanmean(DLC.KINEMATIC.l_tube_l_y_dmax(:,is_OE))], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.l_tube_r_x_dmax(:,is_OE)) nanmean(DLC.KINEMATIC.l_tube_l_x_dmax(:,is_OE))] ...
    , [nanmean(DLC.KINEMATIC.l_tube_r_y_dmax(:,is_OE)) nanmean(DLC.KINEMATIC.l_tube_l_y_dmax(:,is_OE))], 'Color', 'black');

plot(DLC.KINEMATIC.tip_tongue_x_dmax(is_OE), DLC.KINEMATIC.tip_tongue_y_dmax(is_OE), '.b','MarkerSize',20)
plot(DLC.KINEMATIC.r_food_x_dmax(is_OE), DLC.KINEMATIC.r_food_y_dmax(is_OE), '.y','MarkerSize',20)
plot(DLC.KINEMATIC.l_food_x_dmax(is_OE), DLC.KINEMATIC.l_food_y_dmax(is_OE), '.y','MarkerSize',20)
plot(DLC.KINEMATIC.l_nose_x_dmax(is_OE), DLC.KINEMATIC.l_nose_y_dmax(is_OE), '.k','MarkerSize',20)
plot(DLC.KINEMATIC.r_nose_x_dmax(is_OE), DLC.KINEMATIC.r_nose_y_dmax(is_OE), '.k','MarkerSize',20)
plot(DLC.KINEMATIC.r_tube_r_x_dmax(is_OE), DLC.KINEMATIC.r_tube_r_y_dmax(is_OE), '.g','MarkerSize',20)
plot(DLC.KINEMATIC.r_tube_l_x_dmax(is_OE), DLC.KINEMATIC.r_tube_l_y_dmax(is_OE), '.g','MarkerSize',20)
plot(DLC.KINEMATIC.l_tube_r_x_dmax(is_OE), DLC.KINEMATIC.l_tube_r_y_dmax(is_OE), '.g','MarkerSize',20)
plot(DLC.KINEMATIC.l_tube_l_x_dmax(is_OE), DLC.KINEMATIC.l_tube_l_y_dmax(is_OE), '.g','MarkerSize',20)

mean_x_ = nanmean(DLC.KINEMATIC.tip_tongue_x_lick(is_OE,:));
sem_x_ = nanstd(DLC.KINEMATIC.tip_tongue_x_lick(is_OE,:))/sqrt(sum(is_OE));
mean_y_ = nanmean(DLC.KINEMATIC.tip_tongue_y_lick(is_OE,:));
sem_y_ = nanstd(DLC.KINEMATIC.tip_tongue_y_lick(is_OE,:))/sqrt(sum(is_OE));
% shade(mean_x_ + sem_x_, mean_y_ + sem_y_, 'r',  mean_x_ - sem_x_, mean_y_ - sem_y_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_x_+sem_x_, mean_y_+sem_y_, 'r','LineWidth', 1)
plot(mean_x_-sem_x_, mean_y_-sem_y_, 'r','LineWidth', 1)
plot(mean_x_, mean_y_, 'k','LineWidth', 2)

title(['Outer edge lick (mean+/-sem), n = ' num2str(sum(is_OE))], 'interpret', 'none');
xlim([-10, 20]);
ylim([-20, 20]);
xlabel('x (mm)')
xlabel('y (mm)')
set(gca, 'YDir','reverse')

subplot(3,4,6)
hold on
mean_ = nanmean(DLC.KINEMATIC.d_lick(is_OE,:));
sem_ = nanstd(DLC.KINEMATIC.d_lick(is_OE,:))/sqrt(sum(is_OE));
shade(1:length(mean_), mean_ + sem_, 'r', 1:length(mean_), mean_ - sem_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_, 'k', 'LineWidth', 1)

title(['Displacement'], 'interpret', 'none');
xlim([-50, 550]);
ylim([0, inf]);
xlabel('time (ms)')
ylabel('displacement (mm)')

subplot(3,4,7)
hold on
mean_ = nanmean(DLC.KINEMATIC.v_lick(is_OE,:));
sem_ = nanstd(DLC.KINEMATIC.v_lick(is_OE,:))/sqrt(sum(is_OE));
shade(1:length(mean_), mean_ + sem_, 'r', 1:length(mean_), mean_ - sem_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_, 'k', 'LineWidth', 1)

title(['Velocity'], 'interpret', 'none');
xlim([-50, 550]);
ylim([-inf, inf]);
xlabel('time (ms)')
ylabel('velocity (mm/s)')

subplot(3,4,8)
hold on
mean_ = nanmean(DLC.KINEMATIC.angle_lick(is_OE,:));
sem_ = nanstd(DLC.KINEMATIC.angle_lick(is_OE,:))/sqrt(sum(is_OE));
shade(1:length(mean_), mean_ + sem_, 'r', 1:length(mean_), mean_ - sem_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_, 'k', 'LineWidth', 1)

title(['Angle'], 'interpret', 'none');
xlim([-50, 550]);
ylim([-inf, inf]);
xlabel('time (ms)')
ylabel('angle (deg)')

subplot(3,4,9)
hold on
line([nanmean(DLC.KINEMATIC.r_tube_r_x_dmax(:,is_UT)) nanmean(DLC.KINEMATIC.r_tube_r_x_dmax(:,is_UT))] ...
    ,[nanmean(DLC.KINEMATIC.r_tube_r_y_dmax(:,is_UT)) 20], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.r_tube_l_x_dmax(:,is_UT)) nanmean(DLC.KINEMATIC.r_tube_l_x_dmax(:,is_UT))] ...
    ,[nanmean(DLC.KINEMATIC.r_tube_l_y_dmax(:,is_UT)) 20], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.r_tube_r_x_dmax(:,is_UT)) nanmean(DLC.KINEMATIC.r_tube_l_x_dmax(:,is_UT))] ...
    , [nanmean(DLC.KINEMATIC.r_tube_r_y_dmax(:,is_UT)) nanmean(DLC.KINEMATIC.r_tube_l_y_dmax(:,is_UT))], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.l_tube_r_x_dmax(:,is_UT)) nanmean(DLC.KINEMATIC.l_tube_r_x_dmax(:,is_UT))] ...
    ,[-20 nanmean(DLC.KINEMATIC.l_tube_r_y_dmax(:,is_UT))], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.l_tube_l_x_dmax(:,is_UT)) nanmean(DLC.KINEMATIC.l_tube_l_x_dmax(:,is_UT))] ...
    ,[-20 nanmean(DLC.KINEMATIC.l_tube_l_y_dmax(:,is_UT))], 'Color', 'black');
line([nanmean(DLC.KINEMATIC.l_tube_r_x_dmax(:,is_UT)) nanmean(DLC.KINEMATIC.l_tube_l_x_dmax(:,is_UT))] ...
    , [nanmean(DLC.KINEMATIC.l_tube_r_y_dmax(:,is_UT)) nanmean(DLC.KINEMATIC.l_tube_l_y_dmax(:,is_UT))], 'Color', 'black');

plot(DLC.KINEMATIC.tip_tongue_x_dmax(is_UT), DLC.KINEMATIC.tip_tongue_y_dmax(is_UT), '.b','MarkerSize',20)
plot(DLC.KINEMATIC.r_food_x_dmax(is_UT), DLC.KINEMATIC.r_food_y_dmax(is_UT), '.y','MarkerSize',20)
plot(DLC.KINEMATIC.l_food_x_dmax(is_UT), DLC.KINEMATIC.l_food_y_dmax(is_UT), '.y','MarkerSize',20)
plot(DLC.KINEMATIC.l_nose_x_dmax(is_UT), DLC.KINEMATIC.l_nose_y_dmax(is_UT), '.k','MarkerSize',20)
plot(DLC.KINEMATIC.r_nose_x_dmax(is_UT), DLC.KINEMATIC.r_nose_y_dmax(is_UT), '.k','MarkerSize',20)
plot(DLC.KINEMATIC.r_tube_r_x_dmax(is_UT), DLC.KINEMATIC.r_tube_r_y_dmax(is_UT), '.g','MarkerSize',20)
plot(DLC.KINEMATIC.r_tube_l_x_dmax(is_UT), DLC.KINEMATIC.r_tube_l_y_dmax(is_UT), '.g','MarkerSize',20)
plot(DLC.KINEMATIC.l_tube_r_x_dmax(is_UT), DLC.KINEMATIC.l_tube_r_y_dmax(is_UT), '.g','MarkerSize',20)
plot(DLC.KINEMATIC.l_tube_l_x_dmax(is_UT), DLC.KINEMATIC.l_tube_l_y_dmax(is_UT), '.g','MarkerSize',20)

mean_x_ = nanmean(DLC.KINEMATIC.tip_tongue_x_lick(is_UT,:));
sem_x_ = nanstd(DLC.KINEMATIC.tip_tongue_x_lick(is_UT,:))/sqrt(sum(is_UT));
mean_y_ = nanmean(DLC.KINEMATIC.tip_tongue_y_lick(is_UT,:));
sem_y_ = nanstd(DLC.KINEMATIC.tip_tongue_y_lick(is_UT,:))/sqrt(sum(is_UT));
% shade(mean_x_ + sem_x_, mean_y_ + sem_y_, 'r',  mean_x_ - sem_x_, mean_y_ - sem_y_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_x_+sem_x_, mean_y_+sem_y_, 'r','LineWidth', 1)
plot(mean_x_-sem_x_, mean_y_-sem_y_, 'r','LineWidth', 1)
plot(mean_x_, mean_y_, 'k','LineWidth', 2)

title(['Under tube lick (mean+/-sem), n = ' num2str(sum(is_UT))], 'interpret', 'none');
xlim([-10, 20]);
ylim([-20, 20]);
xlabel('x (mm)')
xlabel('y (mm)')
set(gca, 'YDir','reverse')

subplot(3,4,10)
hold on
mean_ = nanmean(DLC.KINEMATIC.d_lick(is_UT,:));
sem_ = nanstd(DLC.KINEMATIC.d_lick(is_UT,:))/sqrt(sum(is_UT));
shade(1:length(mean_), mean_ + sem_, 'r', 1:length(mean_), mean_ - sem_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_, 'k', 'LineWidth', 1)

title(['Displacement'], 'interpret', 'none');
xlim([-50, 550]);
ylim([0, inf]);
xlabel('time (ms)')
ylabel('displacement (mm)')

subplot(3,4,11)
hold on
mean_ = nanmean(DLC.KINEMATIC.v_lick(is_UT,:));
sem_ = nanstd(DLC.KINEMATIC.v_lick(is_UT,:))/sqrt(sum(is_UT));
shade(1:length(mean_), mean_ + sem_, 'r', 1:length(mean_), mean_ - sem_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_, 'k', 'LineWidth', 1)

title(['Velocity'], 'interpret', 'none');
xlim([-50, 550]);
ylim([-inf, inf]);
xlabel('time (ms)')
ylabel('velocity (mm/s)')

subplot(3,4,12)
hold on
mean_ = nanmean(DLC.KINEMATIC.angle_lick(is_UT,:));
sem_ = nanstd(DLC.KINEMATIC.angle_lick(is_UT,:))/sqrt(sum(is_UT));
shade(1:length(mean_), mean_ + sem_, 'r', 1:length(mean_), mean_ - sem_, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
plot(mean_, 'k', 'LineWidth', 1)

title(['Angle'], 'interpret', 'none');
xlim([-50, 550]);
ylim([-inf, inf]);
xlabel('time (ms)')
ylabel('angle (deg)')

% get name
file_name_dir = dir([DLC.FILE.path_to_analyzed '*_ANALYZED.mat']);
file_name_mat = file_name_dir.name(1:13);
sgtitle([file_name_mat '_' direction], 'Interpreter', 'none')
ESN_Beautify_Plot(gcf,[18 16])
saveas(gcf, [file_name_mat '_lick_sorter_' direction], format);
end
%% Function: PGH_plot_sample_trace
function PGH_plot_sample_trace(DLC)
clearvars -except DLC
% load alignment file for sample diff
path_to_analyzed_tongue = DLC.FILE.path_to_analyzed;
align_file = dir([path_to_analyzed_tongue '*aligned.mat']);
file_name = align_file.name(1:13);

load([path_to_analyzed_tongue align_file.name], 'align_PD');
sample_diff = align_PD.sample_diff;
exp_start_time = align_PD.BEHAVE_PD_xcorr_time_1K(1,1);

% load SACS_ALL_DATA for trial data
path_to_analyzed_eye = [DLC.FILE.path_to_analyzed '..' filesep 'eye' filesep];
sac_file = dir([path_to_analyzed_eye '*RECAL.mat']);
load([path_to_analyzed_eye sac_file.name], 'TRIALS_DATA', 'SACS_ALL_DATA');
time_end_trial = TRIALS_DATA.time_end' - exp_start_time + sample_diff/1000;
tag_sac = SACS_ALL_DATA.tag';
time_onset = SACS_ALL_DATA.time_onset' - exp_start_time + sample_diff/1000; time_onset(tag_sac >= 11) = [];
eye_vm_ = SACS_ALL_DATA.eye_r_vm'; eye_vm_(tag_sac >= 11,:) = [];
eye_vm = smooth(reshape(eye_vm_',[],1),5); 

% build eye time
window = (-49:1:100)*0.001;
for counter_sac = 1 : length(time_onset)
    time_eye_vm_(counter_sac,:) = time_onset(counter_sac) + window;
end
time_eye_vm = reshape(time_eye_vm_',[],1);

% filter
eye_vm(eye_vm < 50) = 0;

% define time interval 
x_lim = [1683 1703];

% plot
figure;
hold on;
plot(DLC.TIME.time_1K,DLC.KINEMATIC.d_tip,'-k');
plot(DLC.TIME.time_d_lick_max_abs,25,'.k', 'MarkerSize', 20);
plot(time_end_trial,25,'.r', 'MarkerSize', 20)
ylim([0 26])
xlim(x_lim)
xlabel('Time (s)');
ylabel('Lick displacement (mm)');
set(gca, 'Ycolor','k')
yyaxis right
plot(time_eye_vm,eye_vm, '-b')
ylim([0 1000])
ylabel('Saccade velocity (deg/s)')
%plot(DLC.TIME.time_1K,ESN_smooth(DLC.FOOD.r_tube_food),'-r');
%plot(DLC.TIME.time_1K,DLC.TIME.l_tube_food,'.-b');
%ylim([0 1.5])
%ylabel('Reward capacity (0:Empty | 1:Full)')
set(gca, 'Ycolor','k')
title([file_name '_' num2str(x_lim(1)) '-' num2str(x_lim(2))], 'interpreter', 'none');
ESN_Beautify_Plot(gcf, [20 10])
end