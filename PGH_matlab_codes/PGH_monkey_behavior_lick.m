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

%% Load Data
fprintf('Loading: ')
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure mat_file_address params funcs;
filename = EXPERIMENT_PARAMS.mat_FileName;
path_name = EXPERIMENT_PARAMS.mat_PathName;
load([path_name filename], 'data');
DLC.data = data;
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
    fprintf('FPS not found, computing now ...')
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
flag_qa = 1;
if flag_qa == 1
    [DLC,EXPERIMENT_PARAMS] = quality_assurance(DLC,EXPERIMENT_PARAMS, params, funcs);
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
end

%% MAIN Function: Quality assurance
function [DLC,EXPERIMENT_PARAMS] = quality_assurance(DLC,EXPERIMENT_PARAMS, params, funcs)
fprintf(['Performing QA: ' EXPERIMENT_PARAMS.file_name])
[DLC,EXPERIMENT_PARAMS] = qa_stationary(DLC,EXPERIMENT_PARAMS, params, funcs);
[DLC,EXPERIMENT_PARAMS] = qa_food_and_tongue(DLC,EXPERIMENT_PARAMS, params, funcs);
fprintf(' --> Completed. \n')

end

%% MAIN Function: Build LICKS_ALL_DATA
function [LICKS_ALL_DATA, EXPERIMENT_PARAMS] = build_LICKS_ALL_DATA(DLC, EXPERIMENT_PARAMS, params, funcs)
fprintf(['Building LICK_DATA_ALL: ' EXPERIMENT_PARAMS.mat_FileName ' ... ' '\n'])
clearvars -except DLC EXPERIMENT_PARAMS LICKS_ALL_DATA flag_figure params funcs;
%% Build tag and lick_tag_list

lick_tag_list = params.lick.tag_name_list;

is_groom= logical(DLC.CLASS.is_grooming_lick);
is_inner_tube_success =  logical(DLC.CLASS.is_r_reward_inner_tube_lick + DLC.CLASS.is_l_reward_inner_tube_lick);
is_inner_tube_fail = logical(DLC.CLASS.is_r_noreward_inner_tube_lick + DLC.CLASS.is_l_noreward_inner_tube_lick);
is_outer_tube_success = logical(DLC.CLASS.is_r_reward_outer_tube_lick + DLC.CLASS.is_l_reward_outer_tube_lick);
is_outer_tube_fail = logical(DLC.CLASS.is_r_noreward_outer_tube_lick + DLC.CLASS.is_l_noreward_outer_tube_lick);
is_bout_start = ismember(DLC.IND.ind_lick_onset ,DLC.IND.ind_lick_onset_str_bout);
is_bout_end = ismember(DLC.IND.ind_lick_onset ,DLC.IND.ind_lick_onset_end_bout);
is_harvest_start = ismember(DLC.IND.ind_lick_onset ,DLC.IND.ind_lick_onset_str_harvest);
is_harvest_end = ismember(DLC.IND.ind_lick_onset ,DLC.IND.ind_lick_onset_end_harvest);


tag(is_groom) = 1;
tag(is_inner_tube_success) = 2;
tag(is_inner_tube_fail) = 3;
tag(is_outer_tube_success) = 4;
tag(is_outer_tube_fail) = 5;
tag(is_bout_start) = 6;
tag(is_bout_end) = 7;
tag(is_harvest_start) = 8;
tag(is_harvest_end) = 9;

LICKS_ALL_DATA.tag = tag;
EXPERIMENT_PARAMS.lick_tag_list = lick_tag_list;

is_bout_start_r = ismember(DLC.IND.ind_lick_onset ,DLC.IND.ind_lick_onset_str_bout(DLC.CLASS.is_bout_r));
is_bout_start_l = ismember(DLC.IND.ind_lick_onset ,DLC.IND.ind_lick_onset_str_bout(DLC.CLASS.is_bout_l));
is_bout_start_g = ismember(DLC.IND.ind_lick_onset ,DLC.IND.ind_lick_onset_str_bout(DLC.CLASS.is_bout_g));
is_bout_end_r = ismember(DLC.IND.ind_lick_onset ,DLC.IND.ind_lick_onset_end_bout(DLC.CLASS.is_bout_r));
is_bout_end_l = ismember(DLC.IND.ind_lick_onset ,DLC.IND.ind_lick_onset_end_bout(DLC.CLASS.is_bout_l));
is_bout_end_g = ismember(DLC.IND.ind_lick_onset ,DLC.IND.ind_lick_onset_end_bout(DLC.CLASS.is_bout_g));

tag_bout(is_bout_start_r) = 1;
tag_bout(is_bout_start_l) = 2;
tag_bout(is_bout_start_g) = 3;
tag_bout(is_bout_end_r) = 4;
tag_bout(is_bout_end_l) = 5;
tag_bout(is_bout_end_g) = 6;

% LICKS_ALL_DATA.is_bout_str = is_bout_start';
% LICKS_ALL_DATA.is_bout_end = is_bout_end';
% LICKS_ALL_DATA.is_harvest_start = is_harvest_start';
% LICKS_ALL_DATA.is_harvest_end = is_harvest_end';
LICKS_ALL_DATA.tag_bout = tag_bout;


%% Build time tags
LICKS_ALL_DATA.time_onset = DLC.TIME.time_lick_onset';
LICKS_ALL_DATA.time_vmax = DLC.TIME.time_v_lick_max_abs';
LICKS_ALL_DATA.time_dmax = DLC.TIME.time_d_lick_max_abs';
LICKS_ALL_DATA.time_vmin = DLC.TIME.time_v_lick_min_abs';
LICKS_ALL_DATA.time_offset = DLC.TIME.time_lick_offset';

%% Build kinematics
% max amp, vm+/-, ang
LICKS_ALL_DATA.tongue_dm_max = DLC.KINEMATIC.d_lick_max'; % mm
LICKS_ALL_DATA.tongue_vm_max = DLC.KINEMATIC.v_lick_max'; % mm/s
LICKS_ALL_DATA.tongue_vm_min = DLC.KINEMATIC.v_lick_min'; % mm/s
LICKS_ALL_DATA.tongue_ang_max = DLC.KINEMATIC.angle_lick_max';  % deg w.r.t face normal vector

% displacement, velocity, and angle of tongue during lick onset to
% offset
LICKS_ALL_DATA.tongue_dm = DLC.KINEMATIC.d_lick';
LICKS_ALL_DATA.tongue_vm = DLC.KINEMATIC.v_lick';
LICKS_ALL_DATA.tongue_ang = DLC.KINEMATIC.angle_lick';

% durations: licks, bout, harvest
LICKS_ALL_DATA.duration_lick = DLC.TIME.time_lick_duration';
LICKS_ALL_DATA.duration_bout = DLC.TIME.time_bout_duration';
LICKS_ALL_DATA.duration_harvest = DLC.TIME.time_harvest_duration';

% displacement, velocity, and angle stream data
LICKS_ALL_DATA.tongue_dm_stream = DLC.KINEMATIC.d_tip';
LICKS_ALL_DATA.tongue_vm_stream = DLC.KINEMATIC.v_tip';
LICKS_ALL_DATA.tongue_ang_stream = DLC.KINEMATIC.angle_midtip';

% time stream data
LICKS_ALL_DATA.time_1K_stream = DLC.TIME.time_1K';

%% Build validity
validity = ones(1, DLC.IND.num_lick);
validity(LICKS_ALL_DATA.duration_lick > 0.5) = 0;
validity(LICKS_ALL_DATA.tongue_dm_max > 25) = 0;
LICKS_ALL_DATA.validity = validity;
%% Build positions (x, y)
% position: tongue
LICKS_ALL_DATA.tongue_tip_px = DLC.KINEMATIC.tip_tongue_x_lick';
LICKS_ALL_DATA.tongue_tip_py = DLC.KINEMATIC.tip_tongue_y_lick';
LICKS_ALL_DATA.tongue_r_px = DLC.KINEMATIC.r_tongue_x_lick';
LICKS_ALL_DATA.tongue_r_py = DLC.KINEMATIC.r_tongue_y_lick';
LICKS_ALL_DATA.tongue_l_px = DLC.KINEMATIC.l_tongue_x_lick';
LICKS_ALL_DATA.tongue_l_py = DLC.KINEMATIC.l_tongue_y_lick';
LICKS_ALL_DATA.tongue_mid_px = DLC.KINEMATIC.mid_tongue_x_lick';
LICKS_ALL_DATA.tongue_mid_py = DLC.KINEMATIC.mid_tongue_y_lick';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.tongue_tip_px_onset = DLC.KINEMATIC.tip_tongue_x_onset;
LICKS_ALL_DATA.tongue_tip_py_onset = DLC.KINEMATIC.tip_tongue_y_onset;
LICKS_ALL_DATA.tongue_r_px_onset = DLC.KINEMATIC.r_tongue_x_onset;
LICKS_ALL_DATA.tongue_r_py_onset = DLC.KINEMATIC.r_tongue_y_onset;
LICKS_ALL_DATA.tongue_l_px_onset = DLC.KINEMATIC.l_tongue_x_onset;
LICKS_ALL_DATA.tongue_l_py_onset = DLC.KINEMATIC.l_tongue_y_onset;
LICKS_ALL_DATA.tongue_mid_px_onset = DLC.KINEMATIC.mid_tongue_x_onset;
LICKS_ALL_DATA.tongue_mid_py_onset = DLC.KINEMATIC.mid_tongue_y_onset;
LICKS_ALL_DATA.tongue_tip_px_vmax = DLC.KINEMATIC.tip_tongue_x_vmax;
LICKS_ALL_DATA.tongue_tip_py_vmax = DLC.KINEMATIC.tip_tongue_y_vmax;
LICKS_ALL_DATA.tongue_r_px_vmax = DLC.KINEMATIC.r_tongue_x_vmax;
LICKS_ALL_DATA.tongue_r_py_vmax = DLC.KINEMATIC.r_tongue_y_vmax;
LICKS_ALL_DATA.tongue_l_px_vmax = DLC.KINEMATIC.l_tongue_x_vmax;
LICKS_ALL_DATA.tongue_l_py_vmax = DLC.KINEMATIC.l_tongue_y_vmax;
LICKS_ALL_DATA.tongue_mid_px_vmax = DLC.KINEMATIC.mid_tongue_x_vmax;
LICKS_ALL_DATA.tongue_mid_py_vmax = DLC.KINEMATIC.mid_tongue_y_vmax;
LICKS_ALL_DATA.tongue_tip_px_dmax = DLC.KINEMATIC.tip_tongue_x_dmax;
LICKS_ALL_DATA.tongue_tip_py_dmax = DLC.KINEMATIC.tip_tongue_y_dmax;
LICKS_ALL_DATA.tongue_r_px_dmax = DLC.KINEMATIC.r_tongue_x_dmax;
LICKS_ALL_DATA.tongue_r_py_dmax = DLC.KINEMATIC.r_tongue_y_dmax;
LICKS_ALL_DATA.tongue_l_px_dmax = DLC.KINEMATIC.l_tongue_x_dmax;
LICKS_ALL_DATA.tongue_l_py_dmax = DLC.KINEMATIC.l_tongue_y_dmax;
LICKS_ALL_DATA.tongue_mid_px_dmax = DLC.KINEMATIC.mid_tongue_x_dmax;
LICKS_ALL_DATA.tongue_mid_py_dmax = DLC.KINEMATIC.mid_tongue_y_dmax;
LICKS_ALL_DATA.tongue_tip_px_vmin = DLC.KINEMATIC.tip_tongue_x_vmin;
LICKS_ALL_DATA.tongue_tip_py_vmin = DLC.KINEMATIC.tip_tongue_y_vmin;
LICKS_ALL_DATA.tongue_r_px_vmin = DLC.KINEMATIC.r_tongue_x_vmin;
LICKS_ALL_DATA.tongue_r_py_vmin = DLC.KINEMATIC.r_tongue_y_vmin;
LICKS_ALL_DATA.tongue_l_px_vmin = DLC.KINEMATIC.l_tongue_x_vmin;
LICKS_ALL_DATA.tongue_l_py_vmin = DLC.KINEMATIC.l_tongue_y_vmin;
LICKS_ALL_DATA.tongue_mid_px_vmin = DLC.KINEMATIC.mid_tongue_x_vmin;
LICKS_ALL_DATA.tongue_mid_py_vmin = DLC.KINEMATIC.mid_tongue_y_vmin;
LICKS_ALL_DATA.tongue_tip_px_offset = DLC.KINEMATIC.tip_tongue_x_offset;
LICKS_ALL_DATA.tongue_tip_py_offset = DLC.KINEMATIC.tip_tongue_y_offset;
LICKS_ALL_DATA.tongue_r_px_offset = DLC.KINEMATIC.r_tongue_x_offset;
LICKS_ALL_DATA.tongue_r_py_offset = DLC.KINEMATIC.r_tongue_y_offset;
LICKS_ALL_DATA.tongue_l_px_offset = DLC.KINEMATIC.l_tongue_x_offset;
LICKS_ALL_DATA.tongue_l_py_offset = DLC.KINEMATIC.l_tongue_y_offset;
LICKS_ALL_DATA.tongue_mid_px_offset = DLC.KINEMATIC.mid_tongue_x_offset;
LICKS_ALL_DATA.tongue_mid_py_offset = DLC.KINEMATIC.mid_tongue_y_offset;

% position: nose
LICKS_ALL_DATA.nose_r_px = DLC.KINEMATIC.r_nose_x_lick';
LICKS_ALL_DATA.nose_r_py = DLC.KINEMATIC.r_nose_y_lick';
LICKS_ALL_DATA.nose_l_px = DLC.KINEMATIC.l_nose_x_lick';
LICKS_ALL_DATA.nose_l_py = DLC.KINEMATIC.l_nose_y_lick';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.nose_r_px_onset = DLC.KINEMATIC.r_nose_x_onset;
LICKS_ALL_DATA.nose_r_py_onset = DLC.KINEMATIC.r_nose_y_onset;
LICKS_ALL_DATA.nose_l_px_onset = DLC.KINEMATIC.l_nose_x_onset;
LICKS_ALL_DATA.nose_l_py_onset = DLC.KINEMATIC.l_nose_y_onset;
LICKS_ALL_DATA.nose_r_px_vmax = DLC.KINEMATIC.r_nose_x_vmax;
LICKS_ALL_DATA.nose_r_py_vmax = DLC.KINEMATIC.r_nose_y_vmax;
LICKS_ALL_DATA.nose_l_px_vmax = DLC.KINEMATIC.l_nose_x_vmax;
LICKS_ALL_DATA.nose_l_py_vmax = DLC.KINEMATIC.l_nose_y_vmax;
LICKS_ALL_DATA.nose_r_px_dmax = DLC.KINEMATIC.r_nose_x_dmax;
LICKS_ALL_DATA.nose_r_py_dmax = DLC.KINEMATIC.r_nose_y_dmax;
LICKS_ALL_DATA.nose_l_px_dmax = DLC.KINEMATIC.l_nose_x_dmax;
LICKS_ALL_DATA.nose_l_py_dmax = DLC.KINEMATIC.l_nose_y_dmax;
LICKS_ALL_DATA.nose_r_px_vmin = DLC.KINEMATIC.r_nose_x_vmin;
LICKS_ALL_DATA.nose_r_py_vmin = DLC.KINEMATIC.r_nose_y_vmin;
LICKS_ALL_DATA.nose_l_px_vmin = DLC.KINEMATIC.l_nose_x_vmin;
LICKS_ALL_DATA.nose_l_py_vmin = DLC.KINEMATIC.l_nose_y_vmin;
LICKS_ALL_DATA.nose_r_px_offset = DLC.KINEMATIC.r_nose_x_offset;
LICKS_ALL_DATA.nose_r_py_offset = DLC.KINEMATIC.r_nose_y_offset;
LICKS_ALL_DATA.nose_l_px_offset = DLC.KINEMATIC.l_nose_x_offset;
LICKS_ALL_DATA.nose_l_py_offset = DLC.KINEMATIC.l_nose_y_offset;

% position: reward
LICKS_ALL_DATA.rew_r_px = DLC.KINEMATIC.r_food_x_lick';
LICKS_ALL_DATA.rew_r_py = DLC.KINEMATIC.r_food_y_lick';
LICKS_ALL_DATA.rew_l_px = DLC.KINEMATIC.l_food_x_lick';
LICKS_ALL_DATA.rew_l_py = DLC.KINEMATIC.l_food_y_lick';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.rew_r_px_onset = DLC.KINEMATIC.r_food_x_onset;
LICKS_ALL_DATA.rew_r_py_onset = DLC.KINEMATIC.r_food_y_onset;
LICKS_ALL_DATA.rew_l_px_onset = DLC.KINEMATIC.l_food_x_onset;
LICKS_ALL_DATA.rew_l_py_onset = DLC.KINEMATIC.l_food_y_onset;
LICKS_ALL_DATA.rew_r_px_vmax = DLC.KINEMATIC.r_food_x_vmax;
LICKS_ALL_DATA.rew_r_py_vmax = DLC.KINEMATIC.r_food_y_vmax;
LICKS_ALL_DATA.rew_l_px_vmax = DLC.KINEMATIC.l_food_x_vmax;
LICKS_ALL_DATA.rew_l_py_vmax = DLC.KINEMATIC.l_food_y_vmax;
LICKS_ALL_DATA.rew_r_px_dmax = DLC.KINEMATIC.r_food_x_dmax;
LICKS_ALL_DATA.rew_r_py_dmax = DLC.KINEMATIC.r_food_y_dmax;
LICKS_ALL_DATA.rew_l_px_dmax = DLC.KINEMATIC.l_food_x_dmax;
LICKS_ALL_DATA.rew_l_py_dmax = DLC.KINEMATIC.l_food_y_dmax;
LICKS_ALL_DATA.rew_r_px_vmin = DLC.KINEMATIC.r_food_x_vmin;
LICKS_ALL_DATA.rew_r_py_vmin = DLC.KINEMATIC.r_food_y_vmin;
LICKS_ALL_DATA.rew_l_px_vmin = DLC.KINEMATIC.l_food_x_vmin;
LICKS_ALL_DATA.rew_l_py_vmin = DLC.KINEMATIC.l_food_y_vmin;
LICKS_ALL_DATA.rew_r_px_offset = DLC.KINEMATIC.r_food_x_offset;
LICKS_ALL_DATA.rew_r_py_offset = DLC.KINEMATIC.r_food_y_offset;
LICKS_ALL_DATA.rew_l_px_offset = DLC.KINEMATIC.l_food_x_offset;
LICKS_ALL_DATA.rew_l_py_offset = DLC.KINEMATIC.l_food_y_offset;

% position: tubes
LICKS_ALL_DATA.rtube_r_px = DLC.KINEMATIC.r_tube_r_x_lick';
LICKS_ALL_DATA.rtube_r_py = DLC.KINEMATIC.r_tube_r_y_lick';
LICKS_ALL_DATA.rtube_l_px = DLC.KINEMATIC.r_tube_l_x_lick';
LICKS_ALL_DATA.rtube_l_py = DLC.KINEMATIC.r_tube_l_y_lick';
LICKS_ALL_DATA.ltube_r_px = DLC.KINEMATIC.l_tube_r_x_lick';
LICKS_ALL_DATA.ltube_r_py = DLC.KINEMATIC.l_tube_r_y_lick';
LICKS_ALL_DATA.ltube_l_px = DLC.KINEMATIC.l_tube_l_x_lick';
LICKS_ALL_DATA.ltube_l_py = DLC.KINEMATIC.l_tube_l_y_lick';
% position at specific kinematic event (onset, vmax, dmax, vmin, offset)
LICKS_ALL_DATA.rtube_r_px_onset = DLC.KINEMATIC.r_tube_r_x_onset;
LICKS_ALL_DATA.rtube_r_py_onset = DLC.KINEMATIC.r_tube_r_y_onset;
LICKS_ALL_DATA.rtube_l_px_onset = DLC.KINEMATIC.r_tube_l_x_onset;
LICKS_ALL_DATA.rtube_l_py_onset = DLC.KINEMATIC.r_tube_l_y_onset;
LICKS_ALL_DATA.rtube_r_px_vmax = DLC.KINEMATIC.r_tube_r_x_vmax;
LICKS_ALL_DATA.rtube_r_py_vmax = DLC.KINEMATIC.r_tube_r_y_vmax;
LICKS_ALL_DATA.rtube_l_px_vmax = DLC.KINEMATIC.r_tube_l_x_vmax;
LICKS_ALL_DATA.rtube_l_py_vmax = DLC.KINEMATIC.r_tube_l_y_vmax;
LICKS_ALL_DATA.rtube_r_px_dmax = DLC.KINEMATIC.r_tube_r_x_dmax;
LICKS_ALL_DATA.rtube_r_py_dmax = DLC.KINEMATIC.r_tube_r_y_dmax;
LICKS_ALL_DATA.rtube_l_px_dmax = DLC.KINEMATIC.r_tube_l_x_dmax;
LICKS_ALL_DATA.rtube_l_py_dmax = DLC.KINEMATIC.r_tube_l_y_dmax;
LICKS_ALL_DATA.rtube_r_px_vmin = DLC.KINEMATIC.r_tube_r_x_vmin;
LICKS_ALL_DATA.rtube_r_py_vmin = DLC.KINEMATIC.r_tube_r_y_vmin;
LICKS_ALL_DATA.rtube_l_px_vmin = DLC.KINEMATIC.r_tube_l_x_vmin;
LICKS_ALL_DATA.rtube_l_py_vmin = DLC.KINEMATIC.r_tube_l_y_vmin;
LICKS_ALL_DATA.rtube_r_px_offset = DLC.KINEMATIC.r_tube_r_x_offset;
LICKS_ALL_DATA.rtube_r_py_offset = DLC.KINEMATIC.r_tube_r_y_offset;
LICKS_ALL_DATA.rtube_l_px_offset = DLC.KINEMATIC.r_tube_l_x_offset;
LICKS_ALL_DATA.rtube_l_py_offset = DLC.KINEMATIC.r_tube_l_y_offset;
LICKS_ALL_DATA.ltube_r_px_onset = DLC.KINEMATIC.l_tube_r_x_onset;
LICKS_ALL_DATA.ltube_r_py_onset = DLC.KINEMATIC.l_tube_r_y_onset;
LICKS_ALL_DATA.ltube_l_px_onset = DLC.KINEMATIC.l_tube_l_x_onset;
LICKS_ALL_DATA.ltube_l_py_onset = DLC.KINEMATIC.l_tube_l_y_onset;
LICKS_ALL_DATA.ltube_r_px_vmax = DLC.KINEMATIC.l_tube_r_x_vmax;
LICKS_ALL_DATA.ltube_r_py_vmax = DLC.KINEMATIC.l_tube_r_y_vmax;
LICKS_ALL_DATA.ltube_l_px_vmax = DLC.KINEMATIC.l_tube_l_x_vmax;
LICKS_ALL_DATA.ltube_l_py_vmax = DLC.KINEMATIC.l_tube_l_y_vmax;
LICKS_ALL_DATA.ltube_r_px_dmax = DLC.KINEMATIC.l_tube_r_x_dmax;
LICKS_ALL_DATA.ltube_r_py_dmax = DLC.KINEMATIC.l_tube_r_y_dmax;
LICKS_ALL_DATA.ltube_l_px_dmax = DLC.KINEMATIC.l_tube_l_x_dmax;
LICKS_ALL_DATA.ltube_l_py_dmax = DLC.KINEMATIC.l_tube_l_y_dmax;
LICKS_ALL_DATA.ltube_r_px_vmin = DLC.KINEMATIC.l_tube_r_x_vmin;
LICKS_ALL_DATA.ltube_r_py_vmin = DLC.KINEMATIC.l_tube_r_y_vmin;
LICKS_ALL_DATA.ltube_l_px_vmin = DLC.KINEMATIC.l_tube_l_x_vmin;
LICKS_ALL_DATA.ltube_l_py_vmin = DLC.KINEMATIC.l_tube_l_y_vmin;
LICKS_ALL_DATA.ltube_r_px_offset = DLC.KINEMATIC.l_tube_r_x_offset;
LICKS_ALL_DATA.ltube_r_py_offset = DLC.KINEMATIC.l_tube_r_y_offset;
LICKS_ALL_DATA.ltube_l_px_offset = DLC.KINEMATIC.l_tube_l_x_offset;
LICKS_ALL_DATA.ltube_l_py_offset = DLC.KINEMATIC.l_tube_l_y_offset;

%% Build reward-tube capacity
% rew capacity
LICKS_ALL_DATA.rew_capacity_r = DLC.FOOD.r_tube_food_lick';
LICKS_ALL_DATA.rew_capacity_l = DLC.FOOD.l_tube_food_lick';




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
time_vid = (1/FPS) : (1/FPS) : size(data_FPS,1)/FPS';

% interpolation from FPS to 1K
time_1K = time_vid(1) : 0.001 : time_vid(end);
for counter_fields = 1 : size(data_FPS,2)
    data_1K(:,counter_fields) = interp1(time_vid,data_FPS(:,counter_fields),time_1K);
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

% find origin, shift, and rotate DLC data
[~, cent_cluster_x] = kmeans(tip_tongue_x, 3);
[~, cent_cluster_y] = kmeans(tip_tongue_y, 3);
[x0, ind_min] = min(cent_cluster_x);
if x0 < (mean(r_nose_x) + mean(l_nose_x))/2
    [~, ind_max] = max(cent_cluster_x);
    ind_mid = 1:3;
    ind_mid(ind_max) = [];
    ind_mid(ind_min) = [];
    x0 = cent_cluster_x(ind_mid);
end

midpoint = (mean(r_nose_y) + mean(l_nose_y)) / 2;
if abs(cent_cluster_y(1) - midpoint) < abs(cent_cluster_y(2) - midpoint) && abs(cent_cluster_y(1) - midpoint) < abs(cent_cluster_y(3) - midpoint)
    y0 = cent_cluster_y(1);
elseif abs(cent_cluster_y(2) - midpoint) < abs(cent_cluster_y(1) - midpoint) && abs(cent_cluster_y(2) - midpoint) < abs(cent_cluster_y(3) - midpoint)
    y0 = cent_cluster_y(2);
else
    y0 = cent_cluster_y(3);
end

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

%         mid_tongue_x = DLC.POINTS.mid_tongue_x;
%         mid_tongue_y = DLC.POINTS.mid_tongue_y;
%         d_mid = sqrt(mid_tongue_x.^2 + mid_tongue_y.^2);
%         d_mid_tip = sqrt((tip_tongue_x-mid_tongue_x).^2 + (tip_tongue_y-mid_tongue_y).^2);
%
%         d_tip_1 = d_mid+d_mid_tip;
d_tip = sqrt(tip_tongue_x.^2 + tip_tongue_y.^2);

time_1K = DLC.TIME.time_1K;

[pks,locs,~, ~]= findpeaks(d_tip, 'MinPeakHeight', 1.5, 'MinPeakProminence', 1);
onset = zeros(size(locs));
offset = zeros(size(locs));
for i = 1:length(locs)
    % find onset ind
    if i == 1
        min_amp = min(d_tip(1:locs(i)));
        threshold = (pks(i) - min_amp) * 0.1 + min_amp;
        ind = find(d_tip(1:locs(i)) < threshold) - 1;
        onset(i) = ind(length(ind)) - 2;
    else
        min_amp = min(d_tip(locs(i-1):locs(i)));
        threshold = (pks(i) - min_amp) * 0.1 + min_amp;
        ind = find(d_tip(locs(i-1):locs(i)) < threshold) - 1;
        onset(i) = ind(length(ind)) + locs(i - 1) - 2;
    end

    % find offset ind
    if i == length(locs)
        min_amp = min(d_tip(locs(i):length(d_tip)));
        threshold = (pks(i) - min_amp) * 0.1 + min_amp;
        ind = find(d_tip(locs(i):length(d_tip)) < threshold, 1) - 1;
        offset(i) = ind + locs(i) + 2;
    else
        min_amp = min(d_tip(locs(i):locs(i+1)));
        threshold = (pks(i) - min_amp) * 0.1 + min_amp;
        ind = find(d_tip(locs(i):locs(i+1)) < threshold, 1) - 1;
        offset(i) = ind + locs(i) + 2;
    end
end
onset = onset - 3;
offset = offset + 1;
del = find(onset <= 0 | offset >= length(d_tip));
onset(del) = [];
offset(del) = [];
locs(del)= [];
pks(del)= [];
onset_medians = movmedian(d_tip(onset), 300);

% find indices to delete (multiple peaks in a lick)
del = [];
for i = 2:length(locs)
    trough = min(d_tip(locs(i - 1):(locs(i)))) - onset_medians(i);
    if trough > 0
        if trough > onset_medians(i) + 1 || ((trough > ((pks(i-1) - onset_medians(i) + 1) * 0.2) && trough > (pks(i) - onset_medians(i) + 1) * 0.2))
            if ~isempty(del)
                del_size = size(del);
                if del(del_size(1), 2) == (i-1)
                    del(del_size(1), 2) = i;
                else
                    del = [del; [i-1 i]];
                end
            else
                del = [del; [i-1 i]];
            end
        end
    end
end

del_size = size(del);
for i = 1:del_size(1)
    if del(i, 1) == 1
        starting = onset(1);
    else
        starting_found = 0;
        starting = del(i, 1) - 1;
        while ~starting_found
            if isnan(offset(starting))
                starting = starting - 1;
            else
                starting_found = 1;
            end
        end
        starting = offset(starting);
    end

    if del(i, 2) == length(offset)
        ending = offset(length(offset));
    else
        ending_found = 0;
        ending = del(i, 2) + 1;
        while ~ending_found
            if isnan(ending)
                ending = ending + 1;
            else
                ending_found = 1;
            end
        end
        ending = onset(ending);
    end
    locs(del(i, 1):del(i, 2))= NaN;
    pks(del(i, 1):del(i, 2))= NaN;
    [val, pos] = max(d_tip(starting:ending));
    pks(del(i, 1)) = val;
    locs(del(i, 1)) = pos + starting;
    onset(del(i, 1)) = starting;
    onset(del(i, 1) + 1:del(i, 2)) = NaN;
    offset(del(i, 1)) = ending;
    offset(del(i, 1) + 1:del(i, 2)) = NaN;
end

pks = pks(~isnan(pks));
locs = locs(~isnan(locs));
onset = onset(~isnan(onset));
offset = offset(~isnan(offset));
onset_medians = movmedian(d_tip(onset), 300);
offset_medians = movmedian(d_tip(offset), 300);

% Compute onsets and offsets based on filtered locs
ind_lick_onset = zeros(size(locs));
ind_lick_offset = zeros(size(locs));
for i = 1:length(locs)
    % find onset ind
    if i == 1
        min_amp = min(d_tip(1:locs(i)));
        if onset_medians(i) > min_amp
            min_amp = onset_medians(i);
        end
        threshold = (pks(i) - min_amp) * 0.1 + min_amp;
        ind = find(d_tip(1:locs(i)) < threshold) -1 ;
        ind_lick_onset(i) = ind(length(ind)) ;
    else
        min_amp = min(d_tip(locs(i-1):locs(i)));
        if onset_medians(i) > min_amp
            min_amp = onset_medians(i);
        end
        threshold = (pks(i) - min_amp) * 0.1 + min_amp;
        ind = find(d_tip(locs(i-1):locs(i)) < threshold) -1 ;
        ind_lick_onset(i) = ind(length(ind)) + locs(i - 1) ;
    end

    % ensure onset occurs after the previous offset
    if i > 1
        while ind_lick_onset(i) < ind_lick_offset(i-1)
            if d_tip(ind_lick_offset(i-1) - 1) - d_tip(ind_lick_offset(i-1)) < d_tip(ind_lick_onset(i) + 1) - d_tip(ind_lick_onset(i))
                ind_lick_offset(i-1) = ind_lick_offset(i-1) - 1;
            else
                ind_lick_onset(i) = ind_lick_onset(i) + 1;
            end
        end
    end

    % find offset ind
    if i == length(locs)
        min_amp = min(d_tip(locs(i):length(d_tip)));
        if offset_medians(i) > min_amp
            min_amp = offset_medians(i);
        end
        threshold = (pks(i) - min_amp) * 0.1 + min_amp;
        ind = find(d_tip(locs(i):length(d_tip)) < threshold, 1) - 1;
        if ind + locs(i) + 1 > length(d_tip)
            ind_lick_offset(i) = length(d_tip);
        else
            ind_lick_offset(i) = ind + locs(i) ;
        end
    else
        min_amp = min(d_tip(locs(i):locs(i+1)));
        if offset_medians(i) > min_amp
            min_amp = offset_medians(i);
        end
        threshold = (pks(i) - min_amp) * 0.1 + min_amp;
        ind = find(d_tip(locs(i):locs(i+1)) < threshold, 1) - 1;
        ind_lick_offset(i) = ind + locs(i) ;
    end
end

% detect and remove d_tip bias caused by noise during selection of origin
bias = inf;
while bias > 0.1
    for counter_ind = 1:length(ind_lick_onset) - 1
        mean_bias_window(counter_ind) = mean( d_tip(ind_lick_offset(counter_ind):ind_lick_onset(counter_ind+1)) );
    end
    bias = mean(mean_bias_window);
    d_tip = d_tip - bias;
    d_tip(d_tip < 0) = 0;
end

% recalculate ind_lick_onset and ind_lick_offset using new d_tip
% with bias removed
for counter_lick = 1 : length(ind_lick_onset)
    inds_ = ind_lick_onset(counter_lick) : ind_lick_offset(counter_lick);
    if( ~isempty(inds_(find(d_tip(inds_) > 0, 1, 'first'))-1 ))
        ind_lick_onset_new(counter_lick,1) = inds_(find(d_tip(inds_) > 0, 1, 'first'))-1;
        ind_lick_offset_new(counter_lick,1) = inds_(find(d_tip(inds_) > 0, 1, 'last'))+1;
    else
        continue;
    end
end

ind_lick_onset_0 = find(ind_lick_onset_new == 0);
ind_lick_offset_0 = find(ind_lick_offset_new == 0);

if ~isempty(ind_lick_onset_0)
    ind_lick_onset_new(ind_lick_onset_0) = [];
    ind_lick_offset_new(ind_lick_onset_0) = [];
end

ind_lick_onset = ind_lick_onset_new;
ind_lick_offset = ind_lick_offset_new;

num_lick = length(ind_lick_onset);

time_lick_onset = time_1K(ind_lick_onset);
time_lick_offset =  time_1K(ind_lick_offset);
time_lick_duration = time_lick_offset - time_lick_onset;

ind_lick_onset_str_bout_ = [1; find(diff(ind_lick_onset) >1000) + 1];
ind_lick_onset_str_bout = ind_lick_onset(ind_lick_onset_str_bout_);
ind_lick_onset_end_bout_ = [find(diff(ind_lick_onset) >1000); length(ind_lick_onset)];
ind_lick_onset_end_bout = ind_lick_onset(ind_lick_onset_end_bout_);

time_lick_onset_str_bout = time_1K(ind_lick_onset_str_bout);
time_lick_onset_end_bout = time_1K(ind_lick_onset_end_bout);
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
DLC.TIME.time_lick_onset_str_bout = time_lick_onset_str_bout';
DLC.TIME.time_lick_onset_end_bout = time_lick_onset_end_bout';
DLC.TIME.time_lick_onset = time_lick_onset';
DLC.TIME.time_lick_offset = time_lick_offset';
DLC.TIME.time_lick_duration = time_lick_duration';
DLC.TIME.time_bout_duration = time_bout_duration';

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure
    hold on;
    plot(time_1K,d_tip,'k');
    plot(time_lick_onset_str_bout,21,'*g');
    plot(time_lick_onset_end_bout,21,'*r');
    plot(time_lick_onset, 20, '.g');
    plot(time_lick_offset, 20, '.r');

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

d_tip = DLC.KINEMATIC.d_tip;
time_lick_duration = DLC.TIME.time_lick_duration;

num_lick = DLC.IND.num_lick;
num_bout = DLC.IND.num_bout;
time_1K = DLC.TIME.time_1K';
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
% d_lick
for counter_lick =  1 : 1 : num_lick
    inds_ = (ind_lick_onset(counter_lick) ) : ind_lick_offset(counter_lick);
    for counter_inds_ = 1 : 1: length(inds_)
        d_lick_(counter_inds_) = d_tip(inds_(counter_inds_));
        d_lick(counter_lick, counter_inds_) = d_lick_(counter_inds_);
    end
    [d_lick_max_local, ind_d_lick_max_local] = max(d_tip(inds_));
    ind_d_lick_max(counter_lick,1) = ind_lick_onset(counter_lick) + ind_d_lick_max_local - 1;
    d_lick_max(counter_lick,1) =  d_lick_max_local;
end

% v_lick
for counter_lick =  1 : 1 : num_lick
    inds_ = ind_lick_onset(counter_lick):ind_lick_offset(counter_lick);
    v_lick_ = [0 (diff(d_tip(inds_))./ (time_1K(2) - time_1K(1)))'];
    v_lick(counter_lick, 1:length(inds_)) = v_lick_;
    [v_lick_max_local, ind_v_lick_max_local] = max(v_lick(counter_lick,:));
    ind_v_lick_max(counter_lick,1) = ind_lick_onset(counter_lick) + ind_v_lick_max_local - 1;
    v_lick_max(counter_lick,1) =  v_lick_max_local;

    [v_lick_min_local, ind_v_lick_min_local] = min(v_lick(counter_lick,:));
    ind_v_lick_min(counter_lick,1) = ind_lick_onset(counter_lick) + ind_v_lick_min_local - 1;
    v_lick_min(counter_lick,1) =  v_lick_min_local;
end

% a_lick
for counter_lick =  1 : 1 : num_lick
    inds_ = ind_lick_onset(counter_lick):ind_lick_offset(counter_lick) ;
    a_lick_ = [0 0 (diff(d_tip(inds_),2)./ (time_1K(2) - time_1K(1)))'];
    a_lick(counter_lick, 1:length(inds_)) = a_lick_;
    [a_lick_max_local, ind_a_lick_max_local] = max(a_lick(counter_lick,:));
    ind_a_lick_max(counter_lick,1) = ind_lick_onset(counter_lick) + ind_a_lick_max_local - 1;
    a_lick_max(counter_lick,1) =  a_lick_max_local;
end

% angle_lick
for counter_lick = 1 : 1 : num_lick
    inds_ = ind_lick_onset(counter_lick):ind_lick_offset(counter_lick);
    for counter_inds_ = 1 : 1: length(inds_)
        angle_lick_(counter_inds_) = rad2deg(atan2(midtip_y(inds_(counter_inds_)),...
            midtip_x(inds_(counter_inds_))));
        angle_lick(counter_lick, counter_inds_) = angle_lick_(counter_inds_);
    end
end

% ILI bout
for counter_bout = 1 : 1 : num_bout
    ILI_bout(counter_bout,1) = mean(diff(time_lick_onset(find(ind_lick_onset >= ind_lick_onset_str_bout(counter_bout) & ind_lick_onset <= ind_lick_onset_end_bout(counter_bout)))));
end

% points_reduced within lick onset and offset only, and at specific
% kinematic events
for counter_lick =  1 : 1 : num_lick
    inds_ = ind_lick_onset(counter_lick):ind_lick_offset(counter_lick);
    ind_onset_(counter_lick) = ind_lick_onset(counter_lick);
    ind_vmax_(counter_lick) = ind_v_lick_max(counter_lick);
    ind_dmax_(counter_lick) = ind_d_lick_max(counter_lick);
    ind_vmin_(counter_lick) = ind_v_lick_min(counter_lick);
    ind_offset_(counter_lick) = ind_lick_offset(counter_lick);
    if length(inds_) > 450
        inds_ = inds_(1:450);
        inds_ = padarray(inds_,[0 50], inds_(1), 'pre');
    elseif length(inds_) <= 450
        inds_ = padarray(inds_,[0 50], inds_(1), 'pre');
        inds_ = padarray(inds_,[0 500-length(inds_)], inds_(end), 'post');
    end
    for counter_inds_ = 1 : 1: length(inds_)
        tip_tongue_x_(counter_inds_) = tip_tongue_x(inds_(counter_inds_));
        tip_tongue_x_lick(counter_lick, counter_inds_) = tip_tongue_x_(counter_inds_);
        tip_tongue_y_(counter_inds_) = tip_tongue_y(inds_(counter_inds_));
        tip_tongue_y_lick(counter_lick, counter_inds_) = tip_tongue_y_(counter_inds_);
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

        r_tongue_x_(counter_inds_) = r_tongue_x(inds_(counter_inds_));
        r_tongue_x_lick(counter_lick, counter_inds_) = r_tongue_x_(counter_inds_);
        r_tongue_y_(counter_inds_) = r_tongue_y(inds_(counter_inds_));
        r_tongue_y_lick(counter_lick, counter_inds_) = r_tongue_y_(counter_inds_);
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

        l_tongue_x_(counter_inds_) = l_tongue_x(inds_(counter_inds_));
        l_tongue_x_lick(counter_lick, counter_inds_) = l_tongue_x_(counter_inds_);
        l_tongue_y_(counter_inds_) = l_tongue_y(inds_(counter_inds_));
        l_tongue_y_lick(counter_lick, counter_inds_) = l_tongue_y_(counter_inds_);
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

        mid_tongue_x_(counter_inds_) = mid_tongue_x(inds_(counter_inds_));
        mid_tongue_x_lick(counter_lick, counter_inds_) = mid_tongue_x_(counter_inds_);
        mid_tongue_y_(counter_inds_) = mid_tongue_y(inds_(counter_inds_));
        mid_tongue_y_lick(counter_lick, counter_inds_) = mid_tongue_y_(counter_inds_);
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

        r_food_x_(counter_inds_) = r_food_x(inds_(counter_inds_));
        r_food_x_lick(counter_lick, counter_inds_) = r_food_x_(counter_inds_);
        r_food_y_(counter_inds_) = r_food_y(inds_(counter_inds_));
        r_food_y_lick(counter_lick, counter_inds_) = r_food_y_(counter_inds_);
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

        l_food_x_(counter_inds_) = l_food_x(inds_(counter_inds_));
        l_food_x_lick(counter_lick, counter_inds_) = l_food_x_(counter_inds_);
        l_food_y_(counter_inds_) = l_food_y(inds_(counter_inds_));
        l_food_y_lick(counter_lick, counter_inds_) = l_food_y_(counter_inds_);
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

        r_nose_x_(counter_inds_) = r_nose_x(inds_(counter_inds_));
        r_nose_x_lick(counter_lick, counter_inds_) = r_nose_x_(counter_inds_);
        r_nose_y_(counter_inds_) = r_nose_y(inds_(counter_inds_));
        r_nose_y_lick(counter_lick, counter_inds_) = r_nose_y_(counter_inds_);
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

        l_nose_x_(counter_inds_) = l_nose_x(inds_(counter_inds_));
        l_nose_x_lick(counter_lick, counter_inds_) = l_nose_x_(counter_inds_);
        l_nose_y_(counter_inds_) = l_nose_y(inds_(counter_inds_));
        l_nose_y_lick(counter_lick, counter_inds_) = l_nose_y_(counter_inds_);
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

        r_tube_r_x_(counter_inds_) = r_tube_r_x(inds_(counter_inds_));
        r_tube_r_x_lick(counter_lick, counter_inds_) = r_tube_r_x_(counter_inds_);
        r_tube_r_y_(counter_inds_) = r_tube_r_y(inds_(counter_inds_));
        r_tube_r_y_lick(counter_lick, counter_inds_) = r_tube_r_y_(counter_inds_);
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

        r_tube_l_x_(counter_inds_) = r_tube_l_x(inds_(counter_inds_));
        r_tube_l_x_lick(counter_lick, counter_inds_) = r_tube_l_x_(counter_inds_);
        r_tube_l_y_(counter_inds_) = r_tube_l_y(inds_(counter_inds_));
        r_tube_l_y_lick(counter_lick, counter_inds_) = r_tube_l_y_(counter_inds_);
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

        l_tube_r_x_(counter_inds_) = l_tube_r_x(inds_(counter_inds_));
        l_tube_r_x_lick(counter_lick, counter_inds_) = l_tube_r_x_(counter_inds_);
        l_tube_r_y_(counter_inds_) = l_tube_r_y(inds_(counter_inds_));
        l_tube_r_y_lick(counter_lick, counter_inds_) = l_tube_r_y_(counter_inds_);
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

        l_tube_l_x_(counter_inds_) = l_tube_l_x(inds_(counter_inds_));
        l_tube_l_x_lick(counter_lick, counter_inds_) = l_tube_l_x_(counter_inds_);
        l_tube_l_y_(counter_inds_) = l_tube_l_y(inds_(counter_inds_));
        l_tube_l_y_lick(counter_lick, counter_inds_) = l_tube_l_y_(counter_inds_);
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
end

v_tip = ([0 (diff(d_tip)./ (time_1K(2) - time_1K(1)))'])';
a_tip = ([0 0 (diff(d_tip, 2)./(time_1K(2) - time_1K(1)))'])';
angle_midtip = rad2deg(atan2(midtip_y, midtip_x));
angle_lick_max = angle_midtip(ind_d_lick_max);

ILR_bout = 1./ILI_bout;
ILI_bout(find(ILR_bout==Inf)) = nan;
ILR_bout(find(ILR_bout==Inf)) = nan;

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
r_food_y = DLC.POINTS.r_food_y;

l_tube_r_x = DLC.POINTS.l_tube_r_x;
l_tube_r_y = DLC.POINTS.l_tube_r_y;
l_tube_l_x = DLC.POINTS.l_tube_l_x;
l_tube_l_y = DLC.POINTS.l_tube_l_y;
l_food_y = DLC.POINTS.l_food_y;

ind_d_lick_max = DLC.IND.ind_d_lick_max;
ind_lick_onset = DLC.IND.ind_lick_onset;
ind_lick_offset = DLC.IND.ind_lick_offset;
time_1K = DLC.TIME.time_1K';

% Classify outer, inner, and grooming
is_r_reward_outer_tube_lick = false(size(ind_lick_onset));
is_r_noreward_outer_tube_lick = false(size(ind_lick_onset));
is_l_reward_outer_tube_lick = false(size(ind_lick_onset));
is_l_noreward_outer_tube_lick = false(size(ind_lick_onset));
is_r_inner_tube_lick = false(size(ind_lick_onset));
is_l_inner_tube_lick = false(size(ind_lick_onset));
is_grooming_lick = false(size(ind_lick_onset));

% tongue in tube
tip_towards_r_tube = tip_tongue_x < r_tube_r_x & tip_tongue_x > r_tube_l_x;
r_towards_r_tube = r_tongue_x < r_tube_r_x & r_tongue_x > r_tube_l_x;
l_towards_r_tube = l_tongue_x < r_tube_r_x & l_tongue_x > r_tube_l_x;
tip_towards_l_tube = tip_tongue_x < l_tube_r_x & tip_tongue_x > l_tube_l_x;
r_towards_l_tube = r_tongue_x < l_tube_r_x & r_tongue_x > l_tube_l_x;
l_towards_l_tube = l_tongue_x < l_tube_r_x & l_tongue_x > l_tube_l_x;
tip_crossed_r_tube = tip_tongue_y > (r_tube_r_y+r_tube_l_y)/2 - 1;
r_crossed_r_tube = r_tongue_y > (r_tube_r_y+r_tube_l_y)/2 - 1;
l_crossed_r_tube = l_tongue_y > (r_tube_r_y+r_tube_l_y)/2 - 1;
tip_crossed_l_tube = tip_tongue_y < (l_tube_r_y+l_tube_l_y)/2 + 1;
r_crossed_l_tube = r_tongue_y < (l_tube_r_y+l_tube_l_y)/2 + 1;
l_crossed_l_tube = l_tongue_y < (l_tube_r_y+l_tube_l_y)/2 + 1;
tip_in_r_tube = tip_towards_r_tube & tip_crossed_r_tube;
r_in_r_tube = r_towards_r_tube & r_crossed_r_tube;
l_in_r_tube = l_towards_r_tube & l_crossed_r_tube;
tip_in_l_tube = tip_towards_l_tube & tip_crossed_l_tube;
r_in_l_tube = r_towards_l_tube & r_crossed_l_tube;
l_in_l_tube = l_towards_l_tube & l_crossed_l_tube;

for i = 1:length(ind_lick_onset)
    done = false;
    ind = ind_lick_onset(i):ind_lick_offset(i);

    r_outer = tip_crossed_r_tube(ind) & r_crossed_r_tube(ind) & l_crossed_r_tube(ind)...
        & ((r_in_r_tube(ind) & ~l_in_r_tube(ind)) | (l_in_r_tube(ind) & ~r_in_r_tube(ind)));
    r_inner = tip_in_r_tube(ind) & r_in_r_tube(ind) & l_in_r_tube(ind);
    r_food_overflow = r_food_y(ind) < (r_tube_r_y(ind) + r_tube_l_y(ind))/2;
    if ~isempty(find(r_outer, 1))
        % r_outer
        if bool_tongue_r_tube_full(i) && ~isempty(find(r_outer & r_food_overflow, 1))
            is_r_reward_outer_tube_lick(i) = true;
        else
            is_r_noreward_outer_tube_lick(i) = true;
        end
        done = true;
    elseif ~isempty(find(r_inner, 1))
        % r_inner
        is_r_inner_tube_lick(i) = true;
        done = true;
    end

    if ~done
        l_outer = tip_crossed_l_tube(ind) & r_crossed_l_tube(ind) & l_crossed_l_tube(ind)...
            & ((r_in_l_tube(ind) & ~l_in_l_tube(ind)) | (l_in_l_tube(ind) & ~r_in_l_tube(ind)));
        l_inner = tip_in_l_tube(ind) & r_in_l_tube(ind) & l_in_l_tube(ind);
        l_food_overflow = l_food_y(ind) > (l_tube_r_y(ind) + l_tube_l_y(ind))/2;
        if ~isempty(find(l_outer, 1))
            % l_outer
            if bool_tongue_l_tube_full(i) && ~isempty(find(l_outer & l_food_overflow, 1))
                is_l_reward_outer_tube_lick(i) = true;
            else
                is_l_noreward_outer_tube_lick(i) = true;
            end
            done = true;
        elseif ~isempty(find(l_inner, 1))
            % l_inner
            is_l_inner_tube_lick(i) = true;
            done = true;
        end
    end

    if ~done
        if ~isempty(find(tip_crossed_r_tube(ind), 1))
            % r_outer
            if bool_tongue_r_tube_full(i) && ~isempty(find(tip_crossed_r_tube(ind) & r_food_overflow, 1))
                is_r_reward_outer_tube_lick(i) = true;
            else
                is_r_noreward_outer_tube_lick(i) = true;
            end
        elseif ~isempty(find(tip_crossed_l_tube(ind), 1))
            % l_outer
            if bool_tongue_l_tube_full(i) && ~isempty(find(tip_crossed_l_tube(ind) & l_food_overflow, 1))
                is_l_reward_outer_tube_lick(i) = true;
            else
                is_l_noreward_outer_tube_lick(i) = true;
            end
        else
            is_grooming_lick(i) = true;
        end
    end
end

is_r_reward_inner_tube_lick = bool_tongue_r_tube_full == 1 & is_r_inner_tube_lick;
is_r_noreward_inner_tube_lick = is_r_inner_tube_lick & ~is_r_reward_inner_tube_lick;
is_l_reward_inner_tube_lick = bool_tongue_l_tube_full == 1 & is_l_inner_tube_lick;
is_l_noreward_inner_tube_lick = is_l_inner_tube_lick & ~is_l_reward_inner_tube_lick;
ind_grooming_lick = find(is_grooming_lick);
ind_r_reward_inner_tube_lick = find(is_r_reward_inner_tube_lick);
ind_r_noreward_inner_tube_lick = find(is_r_noreward_inner_tube_lick);
ind_r_reward_outer_tube_lick = find(is_r_reward_outer_tube_lick);
ind_r_noreward_outer_tube_lick = find(is_r_noreward_outer_tube_lick);
ind_l_reward_inner_tube_lick = find(is_l_reward_inner_tube_lick);
ind_l_noreward_inner_tube_lick = find(is_l_noreward_inner_tube_lick);
ind_l_reward_outer_tube_lick = find(is_l_reward_outer_tube_lick);
ind_l_noreward_outer_tube_lick = find(is_l_noreward_outer_tube_lick);

ind_lick_onset_r_reward = ind_lick_onset(sort([ind_r_reward_inner_tube_lick;ind_r_reward_outer_tube_lick;ind_r_noreward_inner_tube_lick;ind_r_noreward_outer_tube_lick]));
time_lick_onset_r_reward = time_1K(ind_lick_onset_r_reward);
ind_lick_offset_r_reward = ind_lick_offset(sort([ind_r_reward_inner_tube_lick;ind_r_reward_outer_tube_lick;ind_r_noreward_inner_tube_lick;ind_r_noreward_outer_tube_lick]));
time_lick_offset_r_reward = time_1K(ind_lick_offset_r_reward);
ind_lick_onset_l_reward = ind_lick_onset(sort([ind_l_reward_inner_tube_lick;ind_l_reward_outer_tube_lick;ind_l_noreward_inner_tube_lick;ind_l_noreward_outer_tube_lick]));
time_lick_onset_l_reward = time_1K(ind_lick_onset_l_reward);
ind_lick_offset_l_reward = ind_lick_offset(sort([ind_l_reward_inner_tube_lick;ind_l_reward_outer_tube_lick;ind_l_noreward_inner_tube_lick;ind_l_noreward_outer_tube_lick]));
time_lick_offset_l_reward = time_1K(ind_lick_offset_l_reward);

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

if EXPERIMENT_PARAMS.flag_figure_debug == 1
    figure;
    hold on;
    % plot(tip_tongue_x(ind_d_lick_max),tip_tongue_y(ind_d_lick_max), '.k');
    plot(tip_tongue_x(ind_d_lick_max(ind_grooming_lick)),tip_tongue_y(ind_d_lick_max(ind_grooming_lick)), 'ob');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_inner_tube_lick)), 'og');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_inner_tube_lick)), 'or');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_noreward_outer_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_noreward_outer_tube_lick)), 'oc');
    plot(tip_tongue_x(ind_d_lick_max(ind_r_reward_outer_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_r_reward_outer_tube_lick)), 'om');
    plot(r_nose_x, r_nose_y,'ok');
    plot(l_nose_x, l_nose_y,'ok');
    plot(r_tube_r_x(ind_d_lick_max),r_tube_r_y(ind_d_lick_max),'sk');
    plot(r_tube_l_x(ind_d_lick_max),r_tube_l_y(ind_d_lick_max),'sk');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_noreward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_noreward_inner_tube_lick)), 'og');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_reward_inner_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_reward_inner_tube_lick)), 'or');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_noreward_outer_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_noreward_outer_tube_lick)), 'oc');
    plot(tip_tongue_x(ind_d_lick_max(ind_l_reward_outer_tube_lick)),tip_tongue_y(ind_d_lick_max(ind_l_reward_outer_tube_lick)), 'om');
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
        r_tube_food_lick(counter_lick, counter_inds_) = r_tube_food_(counter_inds_);
        l_tube_food_(counter_inds_) = l_tube_food(inds_(counter_inds_));
        l_tube_food_lick(counter_lick, counter_inds_) = l_tube_food_(counter_inds_);
    end
end

DLC.FOOD.r_tube_food = r_tube_food;
DLC.FOOD.r_food_x = r_food_x;
DLC.FOOD.r_food_y = r_food_y;
DLC.FOOD.l_tube_food = l_tube_food;
DLC.FOOD.l_food_x = l_food_x;
DLC.FOOD.l_food_y = l_food_y;

DLC.FOOD.r_tube_food_lick = r_tube_food_lick;
DLC.FOOD.l_tube_food_lick = l_tube_food_lick;

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
ind_lick_onset_r_reward = DLC.CLASS.ind_lick_onset_r_reward;
time_lick_onset_r_reward = DLC.TIME.time_lick_onset_r_reward;
r_tube_food = DLC.FOOD.r_tube_food;
ind_lick_onset_l_reward = DLC.CLASS.ind_lick_onset_l_reward;
time_lick_onset_l_reward = DLC.TIME.time_lick_onset_l_reward;
l_tube_food = DLC.FOOD.l_tube_food;

ind_lick_onset_str_bout = DLC.IND.ind_lick_onset_str_bout;
ind_lick_onset_end_bout = DLC.IND.ind_lick_onset_end_bout;
ind_lick_onset = ind_lick_onset(~is_grooming_lick);
ind_lick_onset_str_harvest_ = [1; find(diff(ind_lick_onset) >1000) + 1];
ind_lick_onset_str_harvest = ind_lick_onset(ind_lick_onset_str_harvest_);
ind_lick_onset_end_harvest_ = [find(diff(ind_lick_onset) >1000); length(ind_lick_onset)];
ind_lick_onset_end_harvest = ind_lick_onset(ind_lick_onset_end_harvest_);

num_lick_harvest = (ind_lick_onset_end_harvest_ - ind_lick_onset_str_harvest_) + 1;

time_lick_onset_str_harvest = time_1K(ind_lick_onset_str_harvest);
time_lick_onset_end_harvest = time_1K(ind_lick_onset_end_harvest);
time_harvest_duration = time_1K(ind_lick_onset_end_harvest) - ...
    time_1K(ind_lick_onset_str_harvest);

inds_del_harvest = num_lick_harvest<3;
num_lick_harvest(inds_del_harvest) = [];
time_harvest_duration(inds_del_harvest) = [];
ind_lick_onset_str_harvest(inds_del_harvest) = [];
ind_lick_onset_end_harvest(inds_del_harvest) = [];
time_lick_onset_str_harvest(inds_del_harvest) = [];
time_lick_onset_end_harvest(inds_del_harvest) = [];

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
        plot(time_lick_onset_r_reward,19,'.c');
    end
    if ~isempty(time_lick_onset_l_reward)
        plot(time_lick_onset_l_reward,19,'.b');
    end
    yyaxis right
    plot(time_1K,r_tube_food,'.-c');
    plot(time_1K,l_tube_food,'.-b');
    ylabel('Reward capacity (0:Empty | 1:Full)')
    title([ EXPERIMENT_PARAMS.file_name ': Harvest | ' num2str(sum(DLC.CLASS.is_bout_l)) ' left bouts | ' num2str(sum(DLC.CLASS.is_bout_r)) ' right bouts'], 'interpreter', 'none');
    ESN_Beautify_Plot(gcf, [20 10])
end

fprintf(' --> Completed. \n')
end

%% Function: QA on stationary markers
function [DLC, EXPERIMENT_PARAMS] = qa_stationary(DLC,EXPERIMENT_PARAMS, params, funcs)
FPS = EXPERIMENT_PARAMS.FPS;
data_ = table2array(DLC.data);

vid_height = DLC.FILE.vid_height;
% kalman filter for staionary markers (nose and tube markers)
fprintf("\nCorrecting nose markers...");
r_nose_x = data_(:,9);
r_nose_y = data_(:,10);
l_nose_x = data_(:,11);
l_nose_y = data_(:,12);
% Eliminate erroneous tracking from nose markers
[r_nose_x_qa, r_nose_y_qa, ~] = find_erroneous_markers(r_nose_x, r_nose_y, "r", vid_height/2);
[l_nose_x_qa, l_nose_y_qa, ~] = find_erroneous_markers(l_nose_x, l_nose_y, "l",vid_height/2);
% Perform kalman filtering
%         [r_nose_x_kalman, r_nose_y_kalman] = kalman_filtering(r_nose_x_qa, r_nose_y_qa, FPS);
%         [l_nose_x_kalman, l_nose_y_kalman] = kalman_filtering(l_nose_x_qa, l_nose_y_qa, FPS);
% Perform linear interpolation
r_nose_x_kalman = fillmissing(r_nose_x_qa, 'linear');
r_nose_y_kalman = fillmissing(r_nose_y_qa, 'linear');
l_nose_x_kalman = fillmissing(l_nose_x_qa, 'linear');
l_nose_y_kalman = fillmissing(l_nose_y_qa, 'linear');
% Save variables to DLC
%         DLC.ERR.MEASURED.r_nose = [r_nose_x(DLC.ERR.IND.r_nose) r_nose_y(DLC.ERR.IND.r_nose)];
%         DLC.ERR.MEASURED.l_nose = [l_nose_x(DLC.ERR.IND.l_nose) l_nose_y(DLC.ERR.IND.l_nose)];
%         DLC.ERR.POST_QA.r_nose = [r_nose_x_kalman(DLC.ERR.IND.r_nose) r_nose_y_kalman(DLC.ERR.IND.r_nose)];
%         DLC.ERR.POST_QA.l_nose = [l_nose_x_kalman(DLC.ERR.IND.l_nose) l_nose_y_kalman(DLC.ERR.IND.l_nose)];


data_(:, 9) = r_nose_x_kalman;
data_(:, 10) = r_nose_y_kalman;
data_(:, 11) = l_nose_x_kalman;
data_(:, 12) = l_nose_y_kalman;
DLC.data.r_nose_x = r_nose_x_kalman;
DLC.data.r_nose_y = r_nose_y_kalman;
DLC.data.l_nose_x = l_nose_x_kalman;
DLC.data.l_nose_y = l_nose_y_kalman;
fprintf(' --> Completed. \n')
% kalman filter for tube markers
fprintf("Correcting tube markers...");
r_tube_r_x = data_(:,17);
r_tube_r_y = data_(:,18);
r_tube_l_x = data_(:,19);
r_tube_l_y = data_(:,20);
l_tube_r_x = data_(:,21);
l_tube_r_y = data_(:,22);
l_tube_l_x = data_(:,23);
l_tube_l_y = data_(:,24);
results = verify_food_tube_markers_median(median(l_tube_l_y), median(l_tube_r_y), median(r_tube_l_y), median(r_tube_r_y), vid_height);
if ~isempty(results)
    % find the marker that is incorrect, (whichever further from
    % the nose marker on that side is wrong.
    % then set it to expected position
    if results == 'l'
        if abs(l_tube_l_y - l_nose_y) > abs(l_tube_r_y - l_nose_y)
            % l_tube_l_y is wrong
            l_tube_l_y = l_tube_r_y;
            l_tube_l_x = l_tube_r_x - abs(r_tube_l_x - r_tube_r_x);
            %DLC.ERR.IND.l_tube_l = 1:length(l_tube_l_y);
        else
            % l_tube_r_y is wrong
            l_tube_r_y = l_tube_l_x + abs(r_tube_l_x - r_tube_r_x);
            l_tube_r_x = l_tube_l_y;
            %DLC.ERR.IND.l_tube_r = 1:length(l_tube_r_y);
        end
    else
        if abs(r_tube_l_y - r_nose_y) > abs(r_tube_r_y - r_nose_y)
            % r_tube_l_y is wrong
            r_tube_l_y = r_tube_r_y;
            r_tube_l_x = r_tube_r_x - abs(l_tube_l_x - l_tube_r_x);
        else
            % r_tube_r_y is wrong
            r_tube_r_y = r_tube_l_y;
            r_tube_r_x = r_tube_l_x - abs(l_tube_l_x - l_tube_r_x);
        end
    end
end
% Eliminate erroneous tracking from food tube markers
[l_tube_r_x_qa, l_tube_r_y_qa, ~] = find_erroneous_markers(l_tube_r_x, l_tube_r_y, "l", median(l_nose_y));
[l_tube_l_x_qa, l_tube_l_y_qa, ~] = find_erroneous_markers(l_tube_l_x, l_tube_l_y, "l", median(l_nose_y));
[r_tube_r_x_qa, r_tube_r_y_qa, ~] = find_erroneous_markers(r_tube_r_x, r_tube_r_y, "r", median(r_nose_y));
[r_tube_l_x_qa, r_tube_l_y_qa, ~] = find_erroneous_markers(r_tube_l_x, r_tube_l_y, "r", median(r_nose_y));
% Perform kalman filtering
%             [l_tube_r_x_kalman, l_tube_r_y_kalman] = kalman_filtering(l_tube_r_x_qa, l_tube_r_y_qa, FPS);
%             [l_tube_l_x_kalman, l_tube_l_y_kalman] = kalman_filtering(l_tube_l_x_qa, l_tube_l_y_qa, FPS);
%             [r_tube_r_x_kalman, r_tube_r_y_kalman] = kalman_filtering(r_tube_r_x_qa, r_tube_r_y_qa, FPS);
%             [r_tube_l_x_kalman, r_tube_l_y_kalman] = kalman_filtering(r_tube_l_x_qa, r_tube_l_y_qa, FPS);
% Perform linear interpolation
l_tube_r_x_kalman = fillmissing(l_tube_r_x_qa, 'linear');
l_tube_r_y_kalman = fillmissing(l_tube_r_y_qa, 'linear');
l_tube_l_x_kalman = fillmissing(l_tube_l_x_qa, 'linear');
l_tube_l_y_kalman = fillmissing(l_tube_l_y_qa, 'linear');
r_tube_r_x_kalman = fillmissing(r_tube_r_x_qa, 'linear');
r_tube_r_y_kalman = fillmissing(r_tube_r_y_qa, 'linear');
r_tube_l_x_kalman = fillmissing(r_tube_l_x_qa, 'linear');
r_tube_l_y_kalman = fillmissing(r_tube_l_y_qa, 'linear');
% Save variables
%             DLC.ERR.MEASURED.r_tube_r = [r_tube_r_x(DLC.ERR.IND.r_tube_r) r_tube_r_y(DLC.ERR.IND.r_tube_r)];
%             DLC.ERR.MEASURED.r_tube_l = [r_tube_l_x(DLC.ERR.IND.r_tube_l) r_tube_l_y(DLC.ERR.IND.r_tube_l)];
%             DLC.ERR.MEASURED.l_tube_r = [l_tube_r_x(DLC.ERR.IND.l_tube_r) l_tube_r_y(DLC.ERR.IND.l_tube_r)];
%             DLC.ERR.MEASURED.l_tube_l = [l_tube_l_x(DLC.ERR.IND.l_tube_l) l_tube_l_y(DLC.ERR.IND.l_tube_l)];
%             DLC.ERR.POST_QA.r_tube_r = [r_tube_r_x_kalman(DLC.ERR.IND.r_tube_r) r_tube_r_y_kalman(DLC.ERR.IND.r_tube_r)];
%             DLC.ERR.POST_QA.r_tube_l = [r_tube_l_x_kalman(DLC.ERR.IND.r_tube_l) r_tube_l_y_kalman(DLC.ERR.IND.r_tube_l)];
%             DLC.ERR.POST_QA.l_tube_r = [l_tube_r_x_kalman(DLC.ERR.IND.l_tube_r) l_tube_r_y_kalman(DLC.ERR.IND.l_tube_r)];
%             DLC.ERR.POST_QA.l_tube_l = [l_tube_l_x_kalman(DLC.ERR.IND.l_tube_l) l_tube_l_y_kalman(DLC.ERR.IND.l_tube_l)];
data_(:,17) = r_tube_r_x_kalman;
data_(:,18) = r_tube_r_y_kalman;
data_(:,19) = r_tube_l_x_kalman;
data_(:,20) = r_tube_l_y_kalman;
data_(:,21) = l_tube_r_x_kalman;
data_(:,22) = l_tube_r_y_kalman;
data_(:,23) = l_tube_l_x_kalman;
data_(:,24) = l_tube_l_y_kalman;
DLC.data.l_tube_r_x = l_tube_r_x_kalman;
DLC.data.l_tube_r_y = l_tube_r_y_kalman;
DLC.data.l_tube_l_x = l_tube_l_x_kalman;
DLC.data.l_tube_l_y = l_tube_l_y_kalman;
DLC.data.r_tube_r_x = r_tube_r_x_kalman;
DLC.data.r_tube_r_y = r_tube_r_y_kalman;
DLC.data.r_tube_l_x = r_tube_l_x_kalman;
DLC.data.r_tube_l_y = r_tube_l_y_kalman;
fprintf(' --> Completed. \n')
    function results = verify_food_tube_markers_median(l_tube_l_y, l_tube_r_y, r_tube_l_y, r_tube_r_y, vid_height)
        results = [];
        l_l_length = l_tube_l_y;
        l_r_length = l_tube_r_y;
        r_l_length = vid_height - r_tube_l_y;
        r_r_length = vid_height - r_tube_r_y;
        if  abs(l_l_length - l_r_length) >= (l_l_length/4) || abs(l_l_length - l_r_length) >= (l_r_length/4)
            results = 'l';
        elseif abs(r_l_length - r_r_length) >= (r_l_length/4) || abs(r_l_length - r_r_length) >= (r_r_length/4)
            results = 'r';
        end
    end
    function [x, y, ind_DLC_ERR] = find_erroneous_markers(x, y, type, ref)
        % define height and width of the boundaries beyond which the
        % stationary point tracking is considered erroneous

        x_ref = movmean(x, 300);
        y_ref = movmean(y, 300);

        % determine if the median of tube markers are reliable 
        klist=1:3;
        myfunc = @(X,K)(kmeans(X, K));
        eva = evalclusters(y, myfunc, 'CalinskiHarabasz', 'klist', klist);
        if eva.OptimalK ~= 1
            warning('off')
            [idx, c] = kmeans(y, eva.OptimalK);      
            if sum(abs(diff(sort(c))) < 5) ~= (eva.OptimalK - 1)
                dist = c - ref;
                if type== "l"
                    dif = min(abs(dist(find(dist < 0))));
                    y_ref = ones(size(y))*(-dif + ref);                    
                elseif type == "r"
                    y_ref = ones(size(y))*min(abs(dist(find(dist > 0)))) + ref;
                end
                cluster_ind = find(c == y_ref(1));
                x_ref = ones(size(y))*median(x(find(idx==cluster_ind)));
            end
        end

        ind_DLC_ERR = [];    
        for i = 1:length(x)
            rect = [x_ref(i) - 5 y_ref(i) - 5 10 10]; % [x_lower, y_lower, x_upper, y_upper]
            if x(i) < rect(1) || x(i) > rect(1) + rect(3) || y(i) < rect(2) || y(i) > rect(2) + rect(4)
                ind_DLC_ERR = [ind_DLC_ERR; i];
            end
        end
        x(ind_DLC_ERR) = NaN;
        y(ind_DLC_ERR) = NaN;

        % prepare data for interpolation
        if isnan(x(length(x)))
            x(length(x)) = x_ref(length(x));
            y(length(y)) = y_ref(length(y));
        end
        if isnan(x(1))
            x(1) = x_ref(1);
            y(1) = y_ref(1);
        end
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
function [DLC, EXPERIMENT_PARAMS] = qa_food_and_tongue(DLC,EXPERIMENT_PARAMS, params, funcs)
data_ = table2array(DLC.data);
path_to_raw = [EXPERIMENT_PARAMS.mat_PathName '..' filesep '..' filesep '..' filesep 'raw_data'];


qa_light_condition = 0;
if qa_light_condition == 1
    fprintf("Checking light condition in video ...");

    r_nose_x = data_(:,9);
    r_nose_y = data_(:,10);
    l_nose_x = data_(:,11);
    l_nose_y = data_(:,12);

    % 1) Lighting condition (for both tongue and food markers)
    % Read video
    %     dir_info = dir(path_to_raw);
    %     for i = 1:length(dir_info)
    %         if regexp(dir_info(i).name, regexptranslate('wildcard', '*.mp4'))
    %             file_name = dir_info(i).name;
    %         end
    %     end
    %     vid = VideoReader(fullfile(path_to_raw, file_name));
    % Define which pixel to look at
    avg_pi = [];
    mean_nose_x = round((mean(r_nose_x) + mean(l_nose_x))/2);
    mean_nose_y = round((mean(r_nose_y) + mean(l_nose_y))/2);

    for i = 1:length(r_nose_x)
        frame = read(vid, i);
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
    y_upper_lim = DLC.FILE.vid_height - 10;
    x_upper_lim = DLC.FILE.vid_width;
    fprintf(' --> Completed. \n')
else
    %     dir_info = dir(path_to_raw);
    %     for i = 1:length(dir_info)
    %         if regexp(dir_info(i).name, regexptranslate('wildcard', '*.mp4'))
    %             file_name = dir_info(i).name;
    %         end
    %     end
    %     vid = VideoReader(fullfile(path_to_raw, file_name));


    TF = [];
    r_nose_x = data_(:,9);
    r_nose_y = data_(:,10);
    l_nose_x = data_(:,11);
    l_nose_y = data_(:,12);
    mean_nose_x = round((mean(r_nose_x) + mean(l_nose_x))/2);
    mean_nose_y = round((mean(r_nose_y) + mean(l_nose_y))/2);
    y0 = (mean(r_nose_y) + mean(l_nose_y ))/2;
    y_lower_lim = 10;
    y_upper_lim = DLC.FILE.vid_height - 10;
    x_upper_lim = DLC.FILE.vid_width;
end

% QA for food markers
fprintf("Correcting food markers...");
r_food_x = data_(:,13);
r_food_y = data_(:,14);
l_food_x = data_(:,15);
l_food_y = data_(:,16);
r_tube_r_x = data_(:,17);
r_tube_l_x = data_(:,19);
l_tube_r_x = data_(:,21);
l_tube_l_x = data_(:,23);

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
% 3) Likelihood from DLC
% l_food
l_food_incorrect_inds = find(TF_left);
l_food_x(l_food_incorrect_inds) = NaN;
l_food_y(l_food_incorrect_inds) = NaN;
if isnan(l_food_x(1))
    l_food_x(1) = l_food_x(find(~TF_left, 1));
    l_food_y(1) = l_food_y(find(~TF_left, 1));
end
if isnan(l_food_x(length(l_food_x)))
    l_food_x(length(l_food_x)) = l_food_x(find(~TF_left, 1, 'last'));
    l_food_y(length(l_food_x)) = l_food_y(find(~TF_left, 1, 'last'));
end
% r_food
r_food_incorrect_inds = find(TF_right);
r_food_x(r_food_incorrect_inds) = NaN;
r_food_y(r_food_incorrect_inds) = NaN;
if isnan(r_food_x(1))
    r_food_x(1) = r_food_x(find(~TF_right, 1));
    r_food_y(1) = r_food_y(find(~TF_right, 1));
end
if isnan(r_food_x(length(r_food_x)))
    r_food_x(length(r_food_x)) = r_food_x(find(~TF_right, 1, 'last'));
    r_food_y(length(r_food_x)) = r_food_y(find(~TF_right, 1, 'last'));
end
% Fixing all incorrection food markers with makima interpolation
l_food_new_x = fillmissing(l_food_x, 'makima');
l_food_new_y = fillmissing(l_food_y, 'makima');
r_food_new_x = fillmissing(r_food_x, 'makima');
r_food_new_y = fillmissing(r_food_y, 'makima');
% Save variables
%         DLC.ERR.IND.l_food = l_food_incorrect_inds;
%         DLC.ERR.IND.r_food = r_food_incorrect_inds;
%         DLC.ERR.MEASURED.l_food = [l_food_old_x(DLC.ERR.IND.l_food) l_food_old_y(DLC.ERR.IND.l_food)];
%         DLC.ERR.MEASURED.r_food = [r_food_old_x(DLC.ERR.IND.r_food) r_food_old_y(DLC.ERR.IND.r_food)];
%         DLC.ERR.POST_QA.l_food = [l_food_new_x(DLC.ERR.IND.l_food) l_food_new_y(DLC.ERR.IND.l_food)];
%         DLC.ERR.POST_QA.r_food = [r_food_new_x(DLC.ERR.IND.r_food) r_food_new_y(DLC.ERR.IND.r_food)];
data_(:,13) = l_food_new_x;
data_(:,14) = l_food_new_y;
data_(:,15) = r_food_new_x;
data_(:,16) = r_food_new_y;
DLC.data.l_food_x = l_food_new_x;
DLC.data.l_food_y = l_food_new_y;
DLC.data.r_food_x = r_food_new_x;
DLC.data.r_food_y = r_food_new_y;

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
% Find incorrect tongue markers
for i = 2:length(tip_tongue_y)
    if tip_tongue_y(i) <= y_lower_lim || tip_tongue_x(i) < mean_nose_x || tip_tongue_y(i) >= y_upper_lim || abs(angle(i) - angle(i-1)) > 80
        TF_tip(i) = 1;
    end
    if l_tongue_y(i) <= y_lower_lim || l_tongue_x(i) < mean_nose_x || l_tongue_y(i) >= y_upper_lim
        TF_l(i) = 1;
    end

    if r_tongue_y(i) <= y_lower_lim || r_tongue_x(i) < mean_nose_x || r_tongue_y(i) >= y_upper_lim
        TF_r(i) = 1;
    end

    if mid_tongue_y(i) <= y_lower_lim || mid_tongue_x(i) < mean_nose_x || mid_tongue_x(i) >= x_upper_lim || mid_tongue_y(i) >= y_upper_lim
        TF_mid = 1;
    end
end

% 3) Likelihood
tip_incorrect_inds = find(TF_tip);
l_incorrect_inds = find(TF_l);
r_incorrect_inds = find(TF_r);
mid_incorrect_inds = find(TF_mid);

tip_tongue_x(tip_incorrect_inds) = NaN;
tip_tongue_y(tip_incorrect_inds) = NaN;
l_tongue_x(l_incorrect_inds) = NaN;
l_tongue_y(l_incorrect_inds) = NaN;
r_tongue_x(r_incorrect_inds) = NaN;
r_tongue_y(r_incorrect_inds) = NaN;
mid_tongue_x(mid_incorrect_inds) = NaN;
mid_tongue_y(mid_incorrect_inds) = NaN;
% Fixing all incorrect tongue markers using makima interpolation
tip_tongue_new_x = fillmissing(tip_tongue_x, 'makima');
tip_tongue_new_y = fillmissing(tip_tongue_y, 'makima');
l_tongue_new_x = fillmissing(l_tongue_x, 'makima');
l_tongue_new_y = fillmissing(l_tongue_y, 'makima');
r_tongue_new_x = fillmissing(r_tongue_x, 'makima');
r_tongue_new_y = fillmissing(r_tongue_y, 'makima');
mid_tongue_new_x = fillmissing(mid_tongue_x, 'makima');
mid_tongue_new_y = fillmissing(mid_tongue_y, 'makima');
% Save variables
%         DLC.ERR.IND.tip_tongue = tip_incorrect_inds;
%         DLC.ERR.IND.l_tongue = l_incorrect_inds;
%         DLC.ERR.IND.r_tongue = r_incorrect_inds;
%         DLC.ERR.IND.mid_tongue = mid_incorrect_inds;
%         DLC.ERR.MEASURED.tip_tongue = [tip_tongue_old_x(DLC.ERR.IND.tip_tongue) tip_tongue_old_y(DLC.ERR.IND.tip_tongue)];
%         DLC.ERR.MEASURED.l_tongue = [l_tongue_old_x(DLC.ERR.IND.l_tongue) l_tongue_old_y(DLC.ERR.IND.l_tongue)];
%         DLC.ERR.MEASURED.r_tongue = [r_tongue_old_x(DLC.ERR.IND.r_tongue) r_tongue_old_y(DLC.ERR.IND.r_tongue)];
%         DLC.ERR.MEASURED.mid_tongue = [mid_tongue_old_x(DLC.ERR.IND.mid_tongue) mid_tongue_old_y(DLC.ERR.IND.mid_tongue)];
%         DLC.ERR.POST_QA.tip_tongue = [tip_tongue_new_x(DLC.ERR.IND.tip_tongue) tip_tongue_new_y(DLC.ERR.IND.tip_tongue)];
%         DLC.ERR.POST_QA.l_tongue = [l_tongue_new_x(DLC.ERR.IND.l_tongue) l_tongue_new_y(DLC.ERR.IND.l_tongue)];
%         DLC.ERR.POST_QA.r_tongue = [r_tongue_new_x(DLC.ERR.IND.r_tongue) r_tongue_new_y(DLC.ERR.IND.r_tongue)];
%         DLC.ERR.POST_QA.mid_tongue = [mid_tongue_new_x(DLC.ERR.IND.mid_tongue) mid_tongue_new_y(DLC.ERR.IND.mid_tongue)];
data_(:,1) = tip_tongue_new_x;
data_(:,2) = tip_tongue_new_y;
data_(:,3) = r_tongue_new_x;
data_(:,4) = r_tongue_new_y;
data_(:,5) = l_tongue_new_x;
data_(:,6) = l_tongue_new_y;
data_(:,7) = mid_tongue_new_x;
data_(:,8) = mid_tongue_new_y;
DLC.data.tip_tongue_x = tip_tongue_new_x;
DLC.data.tip_tongue_y = tip_tongue_new_y;
DLC.data.r_tongue_x = r_tongue_new_x;
DLC.data.r_tongue_y = r_tongue_new_y;
DLC.data.l_tongue_x = l_tongue_new_x;
DLC.data.l_tongue_y = l_tongue_new_y;
DLC.data.mid_tongue_x = mid_tongue_new_x;
DLC.data.mid_tongue_y = mid_tongue_new_y;
fprintf(' --> Completed. \n')
end
