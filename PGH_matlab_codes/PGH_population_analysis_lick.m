%% MASTER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function PGH_population_analysis
function PGH_population_analysis_lick
% clc; clear; close all;
tic
path_data_monkey_sorted = {'data_125d','data_59d'};

%% params
params.lick.length_trace    = 500;
params.lick.inds_span       = ((-(params.lick.length_trace/2)+1) : 1 : (params.lick.length_trace/2))';
params.lick.ang_step        = 45;
params.lick.ang_edges       = -90 - params.lick.ang_step /2 : params.lick.ang_step  : 90 + params.lick.ang_step /2 ;
params.lick.ang_values      = -90 : params.lick.ang_step : 90;
params.lick.amp_step    = 5;
params.lick.amp_edges    = 0 : params.lick.amp_step : 25;
params.lick.vel_step    = 150;
params.lick.vel_edges    = 0 : params.lick.vel_step : 750;
params.lick.dur_step = 100;
params.lick.dur_edges = 0 : params.lick.dur_step : 500;
params.lick.tags_CS_ang_avg  = [1:9];
params.lick.tag_name_list = {  ...
    'groom', ... % tag 1
    'inner_tube_success', ... % tag 2
    'inner_tube_fail', ... % tag 3
    'outer_tube_success' ..., % tag 4
    'outer_tube_fail', ... % tag 5
    'bout_start', ... % tag 6
    'bout_end', ... % tag 7
    'harvest_start',  ...% tag 8
    'harvest_end', ... % tag 9
    };
params.lick.tag_bout_name_list = {  ...
    'bout_start_r', ... % tag 1
    'bout_start_l', ... % tag 2
    'bout_start_g', ... % tag 3
    'bout_end_r' ..., % tag 4
    'bout_end_l', ... % tag 5
    'bout_end_g', ... % tag 6
    };



%% Build functions
%extract_population_data(path_data_monkey_sorted);
%build_population_data
%% Plot Functions
plot_session_summary(params)
toc
end

%% UTIL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function extract_population_data
function extract_population_data(path_data_monkey_sorted)
num_subj = numel(path_data_monkey_sorted);
out_path = [pwd  filesep 'LICK' filesep 'session_data' filesep];
mkdir(out_path);

for counter_subject = 1:num_subj

    current_monkey_path_ = path_data_monkey_sorted{counter_subject};
    fprintf(['### ' 'Analyzing subject ', current_monkey_path_, ' no. ',...
        num2str(counter_subject), ' / ' num2str(num_subj) ' ### \n']);

    if ~strcmp(current_monkey_path_(end), filesep)
        current_monkey_path = [current_monkey_path_ filesep];
    end

    session_list = dir([current_monkey_path filesep '20*' filesep '20*']);
    session_list([session_list.isdir] == 0) = [];
    sess_list = {session_list.name};

    num_sess = length(sess_list);
    % Loop over sessions
    for counter_sess = 1 : 1 : num_sess
        current_sess = sess_list{counter_sess};
        fprintf(['  ### ' 'sess ', current_sess, ' no. ',...
            num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

        sess_path = [current_monkey_path, current_sess(1:7), filesep, current_sess, filesep];

        sess_name = regexprep(current_sess(3:end),'-','');
        sess_meta_data = readtable([sess_path sess_name '.xls']);

        rec_list = sess_meta_data.folder_name(logical...
            (sess_meta_data.ephys .* sess_meta_data.eye));

        n_recs = length(rec_list);

        for counter_rec = 1:n_recs
            path_to_rec_tongue = [sess_path rec_list{counter_rec} filesep 'analyzed_data' filesep 'behavior_data' filesep 'tongue' filesep];
            current_rec_data = dir([path_to_rec_tongue '*_ANALYZED.mat']);
            current_rec_data_name = current_rec_data.name;

            data_recording(counter_rec) = load([path_to_rec_tongue current_rec_data_name], 'LICKS_ALL_DATA','EXPERIMENT_PARAMS');
            data_recording(counter_rec).LICKS_ALL_DATA = rmfield(data_recording(counter_rec).LICKS_ALL_DATA, ...
                {'time_1K_stream', 'tongue_dm_stream', 'tongue_vm_stream','tongue_ang_stream', 'tongue_dm', 'tongue_vm', 'tongue_ang'});

            data_recording(counter_rec).LICKS_ALL_DATA.rec_duration = data_recording(counter_rec).EXPERIMENT_PARAMS.duration_video;

        end
        [LICKS_ALL_DATA] = PGH_combine_recs(data_recording);
        id = [current_sess '_combine_' num2str(n_recs)];
        LICKS_ALL_DATA.id = id;
        LICKS_ALL_DATA.subject = current_monkey_path_;
        LICKS_ALL_DATA.rec_list = rec_list;
        LICKS_ALL_DATA.rec_num = n_recs;

        % Save data
        save([out_path id '.mat'], 'LICKS_ALL_DATA','-v7.3')

    end
    fprintf('### ALL DONE. ###\n')
end

end


%% function build_population_data
function build_population_data
%% Load data
path_LICK_session_data = 'Z:\video_10TB\Paul\LICK\session_data';
out_path = [pwd  filesep 'LICK' filesep];

if ~strcmp(path_LICK_session_data(end), filesep);path_LICK_session_data = [path_LICK_session_data filesep];end
session_list = dir([path_LICK_session_data '*mat']);
num_session = numel(session_list);

%% Loop over sessions
for counter_session = 1 : num_session
    fprintf(['### ' 'Analyzing session no. ' num2str(counter_session), ' / ' num2str(num_session) ' ### \n']);
    data_session = load([session_list(counter_session).folder filesep session_list(counter_session).name]);
    population_lick_data(counter_session) = data_session.LICKS_ALL_DATA;
end
%% Save data
fprintf(['Saving .mat file' ' ...'])
save([out_path 'population_lick_data' '.mat'], 'population_lick_data', '-v7.3');

fprintf('### ALL DONE. ###\n')

end
%% function combine_recs
function [LICKS_ALL_DATA] = PGH_combine_recs(data_recording)

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


end

%% PLOT FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_session_summary(params)
end
