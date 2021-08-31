%% MASTER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function PH_population_coding_sac_sorter
function PH_population_coding_sac_sorter
%% Steps to prepare the combined eye & neuro data
clc; clear; close all;
tic
% (1) 
re_run_ESN_sac_sorter(); % this func will re-analyze the _ANALYZED files and add the SACS_ALL_DATA with tags to them
% (2) re_run_add_ephys_sac_sorter(); % this func will fuse the eye & neuro data and save the results with _sac name
% (3) combine_sac_files(); % combine the _sac files and form one file per cell. save the data in the ALL_PCELL_### folder.
% (3) combine_pair_files(); % This is for pair of simultaneous pCells 
% (4) CS_on_analysis(); % load cell_data files (_combine_) and add the CS_on_analysis  to them.
% (5) build_population_data(); % load cell_data files (_combine_) and form population data.
% (5) build_synchrony_data(); % form joint probability of cell_1 and cell_2.
% (6) build_neural_properties(); % load cell_data files (_combine_) and form properties data.
% (7) sac_modulation_index(); % load cell_data files (_combine_) and add the sac_modulation to them.
toc

%% Plot functions
clc; clear;
tic
% (1) plot_neural_properties(1); % load population_neural_properties and plot neural_properties
% (2) plot_CS_on_properties(2); % load population_neural_properties and plot CS_on properties
% (3) plot_population_data_iteratively(3);
% (3) plot_population_data_single_condition();
% (4.1) 
plot_synch_ratio();
% (4.2) plot_synch_MI
% (5) plot_modulation_z_score(4);
% (6.1) plot_single_session_modulation('210507_174229_03_combine_1');

% (6.2)
% plot the figs for the sessions requested by reza
% reza_list_8_16(1, 1:16) = {'210407_172124_02_combine_2', '210412_153217_03_2_combine_2', '210412_153217_03_combine_2', '210412_164800_03_combine_1', '210414_124944_02_combine_5', '210414_124944_03_combine_5', '210414_124944_04_2_combine_5', '210415_151027_03_combine_3', '210510_154635_03_2_combine_4', '210510_162010_04_combine_3', '210510_164823_03_combine_2', '210512_143906_02_2_combine_5', '210512_143906_02_combine_5', '210512_143906_03_combine_5', '210512_143906_07_combine_5', '210615_165101_01_2_combine_2'};
% for counter_file = 1:length(reza_list_8_16)
%  plot_single_session_modulation(char(reza_list_8_16(counter_file)));
% end

% (6.3)
%plot the figs for all sessions in cell_data
% dir_list = dir('*.mat');
% full_list = {dir_list.name}';
% for counter_file = 1:length(full_list)
% file_name = char(full_list(counter_file));
% file_name = file_name(1:end-4);
% plot_single_session_modulation(file_name);
% end


toc
end

%% UTIL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function build_pCell_ids()
function pCell_ids = build_pCell_ids(flag_pair_list)
JSP_global_variables(flag_pair_list);
pCell_list = PH_build_pCell_list(flag_pair_list);
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
pCell_ids = cell(num_pCells, 1);
for counter_pCell = 1 : num_pCells
    file_name_cell = pCell_list{counter_pCell, 1};
    if file_name_cell(18) == 's'
        id_          = file_name_cell(1:16);
    elseif isstrprop(file_name_cell(18),'digit')
        id_          = file_name_cell(1:18);
    else
        error('Build plot_data_compress: cell id does not follow the standards')
    end
    num_recording = nansum(pCell_list_isstr(counter_pCell, :));
    cell_file_name = [id_ '_' 'combine' '_' num2str(num_recording)];
    pCell_ids{counter_pCell, 1} = cell_file_name;
end
end
%% function idx_mirza_ramon()
function [idx_mirza, idx_ramon] = idx_mirza_ramon(flag_pair_list)
pCell_ids = build_pCell_ids(flag_pair_list);
num_pCells = size(pCell_ids, 1);
idx_mirza = false(num_pCells, 1);
idx_ramon = false(num_pCells, 1);
if flag_pair_list
    if num_pCells ~= 80
        error('mirza_ramon_idx: number of P-cells is not 80. Please modify the code.')
    end
    idx_mirza(1:54,1)  = true;
    idx_ramon(55:80,1) = true;
else
    if num_pCells ~= 143
        error('mirza_ramon_idx: number of P-cells is not 143. Please modify the code.')
    end
    idx_mirza(1:90,1)   = true;
    idx_ramon(91:143,1) = true;
end

end
%% function expand_index_event_data(event_data)
function event_data_expand = expand_index_event_data(event_data, dim)
if nargin < 2
    dim = 1;
end

JSP_global_variables();
global expand_index
if isempty(expand_index)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
event_data = logical(event_data);
event_data_expand = logical(event_data);
for counter = (-expand_index) : 1 : (expand_index)
    event_data_expand = event_data_expand | circshift(event_data, counter, dim);
end

end
%% function mutual_info(event_1, event_2)
function MI = mutual_info(event_1, event_2)
%%
JSP_global_variables();
global length_trace
if size(event_1, 1) == length_trace
    event_1 = event_1';
elseif size(event_1, 2) == length_trace
    % DO nothing
else
    fprintf('mutual_info: event_1 incorrect size.\n')
end
if size(event_2, 1) == length_trace
    event_2 = event_2';
elseif size(event_2, 2) == length_trace
    % DO nothing
else
    fprintf('mutual_info: event_1 incorrect size.\n')
end
MI = zeros(1, length_trace);
event_1 = logical(event_1);
event_2 = logical(event_2);
for state_1 = [true false]
    for state_2 = [true false]
        P_1 = (event_1 == state_1);
        P_2 = (event_2 == state_2);
        P_12 = P_1 & P_2;
        P_1 = mean(P_1);
        P_2 = mean(P_2);
        P_12 = mean(P_12);
        MI_ = ( P_12 .* log2(P_12 ./ P_1 ./ P_2) );
        MI_(isnan(MI_)) = 0;
        MI = MI + MI_;
    end
end

end
%% function estimate_pair_distance()
function pCell_contact_distance = estimate_pair_distance()
%% get pair_list
flag_pair_list = true;   % This should be true. DO NOT change it to false
JSP_global_variables(flag_pair_list);
pCell_ids = build_pCell_ids(flag_pair_list);
num_pCells = size(pCell_ids, 1);
path_data_monkey_sorted = uigetdir;
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

%% extract ch_num and electrode_type
pCell_ch_num = nan(num_pCells, 1);
pCell_electrode_ch_num = nan(num_pCells, 1);
pCell_ch_level = nan(num_pCells, 1);
% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    %% build address
    file_name_cell = pCell_ids{counter_pCell,1}(1:13);
    year_ = file_name_cell(1:2);
    month_ = file_name_cell(3:4);
    day_ = file_name_cell(5:6);
    hour_ = file_name_cell(8:9);
    minute_ = file_name_cell(10:11);
    second_ = file_name_cell(12:13);
    subFolder_month = ['20' year_ '-' month_ filesep];
    subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
    subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
    subFolder_data = ['analyzed_data' filesep];
    %% load meta_data
    file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_data];
    file_name = [file_name_cell '_meta_data.json'];
    fname = [file_path file_name];
    fid = fopen(fname);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    meta_data_ = jsondecode(str);
    num_of_channels = meta_data_.num_of_channels;
    contact_ch_num = str2double(pCell_ids{counter_pCell,1}(15:16));
    ch_level_ = 0;
    if num_of_channels==4
        if contact_ch_num==2
            ch_level_ = 1;
        elseif (contact_ch_num==1)||(contact_ch_num==3)||(contact_ch_num==4)
            ch_level_ = 2;
        else
            error('estimate_pair_distance: wrong ch_num for Tetrode.');
        end
    elseif num_of_channels==7
        if contact_ch_num==4
            ch_level_ = 1;
        elseif (contact_ch_num==3)||(contact_ch_num==5)||(contact_ch_num==6)
            ch_level_ = 2;
        elseif (contact_ch_num==1)||(contact_ch_num==2)||(contact_ch_num==7)
            ch_level_ = 3;
        else
            error('estimate_pair_distance: wrong ch_num for Heptode.')
        end
    elseif num_of_channels==64
        ch_level_ = 1;
    else
        error('estimate_pair_distance: unknown electrode type.')
    end
    %% store results
    pCell_electrode_ch_num(counter_pCell, 1) = num_of_channels;
    pCell_ch_num(counter_pCell, 1) = contact_ch_num;
    pCell_ch_level(counter_pCell, 1) = ch_level_;
end

%% estimate contact_distance
pCell_contact_distance = nan(num_pCells, 1);
for counter_pair = 1 : num_pCells/2
%     electrode_ch_num_1 = pCell_electrode_ch_num(2*counter_pair-1, 1);
%     electrode_ch_num_2 = pCell_electrode_ch_num(2*counter_pair,   1);
%     ch_num_1 = pCell_ch_num(2*counter_pair-1, 1);
%     ch_num_2 = pCell_ch_num(2*counter_pair,   1);
    ch_level_1_ = pCell_ch_level(2*counter_pair-1, 1);
    ch_level_2_ = pCell_ch_level(2*counter_pair,   1);
    ch_level_1 = min([ch_level_1_, ch_level_2_]);
    ch_level_2 = max([ch_level_1_, ch_level_2_]);
    contact_distance_ = 0;
    if ch_level_1 == ch_level_2
        contact_distance_ = 50;
    elseif (ch_level_1==1) && (ch_level_2==2)
        contact_distance_ = 100;
    elseif (ch_level_1==1) && (ch_level_2==3)
        contact_distance_ = 200;
    elseif (ch_level_1==2) && (ch_level_2==3)
        contact_distance_ = 110;
    else
        error('estimate_pair_distance: undefined condition.')
    end
    pCell_contact_distance(2*counter_pair-1, 1) = contact_distance_;
    pCell_contact_distance(2*counter_pair,   1) = contact_distance_;
end

end

%% RE-RUN RAW FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function re_run_ESN_sac_sorter()
function re_run_ESN_sac_sorter()
clc; close all;
flag_pair_list = false; % This should be false. DO NOT change it to true
pCell_list = PH_build_pCell_list(flag_pair_list);
path_data_monkey_sorted = uigetdir;

if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    num_recording = nansum(pCell_list_isstr(counter_pCell, :));
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name_cell = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name_cell(1:2);
        month_ = file_name_cell(3:4);
        day_ = file_name_cell(5:6);
        hour_ = file_name_cell(8:9);
        minute_ = file_name_cell(10:11);
        second_ = file_name_cell(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_data' filesep];
        %% load BEHAVE data
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        file_name = [file_name_cell(1:13) '_ANALYZED.mat'];
        if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
        BEHAVE = load([file_path file_name], 'EXPERIMENT_PARAMS', 'TRIALS_DATA');
        %% Build SACS_ALL_DATA using ESN_Sac_Sorter
        flag_session_figure = true;
        [SACS_ALL_DATA, TRIALS_DATA, EXPERIMENT_PARAMS] = JSP_Sac_Sorter(BEHAVE.TRIALS_DATA, BEHAVE.EXPERIMENT_PARAMS, flag_session_figure);

        %% Save _ANALYZED.mat Data to disk
        save([file_path file_name_cell(1:13) '_ANALYZED.mat'], ...
            'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_ALL_DATA', '-append', '-v7.3');

        %% Save _REDUCED.mat Data to disk
        rmfields_list = {'eye_l_vm_filt', 'eye_l_vy_filt', 'eye_l_vx_filt', 'eye_l_py_filt', 'eye_l_px_filt', ...
            'eye_r_vm_filt', 'eye_r_vy_filt', 'eye_r_vx_filt', 'eye_r_py_filt', 'eye_r_px_filt', ...
            'time', 'time_1K', 'target_visible', 'reward', 'tgt_py', 'tgt_px', 'time_tgt', ...
            'eye_l_vm', 'eye_r_vm', 'eye_l_vy', 'eye_l_vx', 'eye_r_vy', 'eye_r_vx', ...
            'eye_l_py', 'eye_l_px', 'eye_r_py', 'eye_r_px', 'time_eyelink', 'inds_invalid', 'inds_trial'};
        TRIALS_DATA = rmfield(TRIALS_DATA,rmfields_list);
        save([file_path file_name_cell(1:13) '_REDUCED.mat'], ...
            'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_ALL_DATA', '-append', '-v7.3');

        %% Save Fig
        hFig_ = gcf;
        file_name_plot_ = file_name_cell(1:13);
        file_path_plot_ = [file_path '..' filesep 'analyzed_figs' filesep];
        saveas(hFig_,[file_path_plot_ file_name_plot_ '_sac_sorter'], 'pdf');
        saveas(hFig_,[file_path_plot_ file_name_plot_ '_sac_sorter'], 'png');
        close(hFig_);
    end
end
fprintf('### ALL DONE. ###\n')
end
%% function re_run_add_ephys_sac_sorter()
function re_run_add_ephys_sac_sorter()
clc; close all;
flag_pair_list = false; % This should be false. DO NOT change it to true
pCell_list = PH_build_pCell_list(flag_pair_list);
path_data_monkey_sorted = uigetdir;

if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    num_recording = nansum(pCell_list_isstr(counter_pCell, :));
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name_cell = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name_cell(1:2);
        month_ = file_name_cell(3:4);
        day_ = file_name_cell(5:6);
        hour_ = file_name_cell(8:9);
        minute_ = file_name_cell(10:11);
        second_ = file_name_cell(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_data' filesep];
        %% RE-RUN add_ephys_sac_sorter
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
        file_name = file_name_cell;
        [SACS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS] = add_ephys_sac_sorter(file_path, file_name);
        %% Save _sac Data to disk
        save([file_path file_name_cell(1:end-4) '_sac' file_name_cell(end-3:end) '.mat'], ...
            'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS', '-v7.3');
        
    end
end
fprintf('### ALL DONE. ###\n')
end
%% function add_ephys_sac_sorter(file_path, file_name)
function [SACS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS] = add_ephys_sac_sorter(file_path, file_name)
%% load EPHYS sorted DATA
if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
file_name = [file_name '.psort'];
fprintf(['Loading ', file_name, ' ... ']);
DATA_PSORT = Psort_read_psort([file_path file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% load EPHYS EVENT DATA
file_name = [EPHYS.CH_sorted_file_name(1:13) '_EVE1_aligned.mat'];
EPHYS.CH_EVE = load([file_path file_name]);
if isfield(EPHYS.CH_EVE, 'EPHYS_time_15K')
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_15K(:);
else
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_30K(:);
end
EPHYS.CH_EVE.EPHYS_time_1K  = reshape(EPHYS.CH_EVE.EPHYS_time_1K ,[], 1);
EPHYS.CH_EVE.BEHAVE_time_1K = reshape(EPHYS.CH_EVE.BEHAVE_time_1K,[], 1);

%% load BEHAVE DATA
file_name = [EPHYS.CH_sorted_file_name(1:13) '_ANALYZED.mat'];
BEHAVE = load([file_path file_name]);
BEHAVE.EXPERIMENT_PARAMS.EPHYS_file_name = EPHYS.CH_sorted_file_name;
BEHAVE.EXPERIMENT_PARAMS.EPHYS_file_path = EPHYS.CH_sorted_file_path;

%% build EPHYS.CH_sorted from DATA_PSORT
JSP_global_variables();
global waveform_inds_span
if isempty(waveform_inds_span)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
ch_data = double(DATA_PSORT.topLevel_data.ch_data);
ch_time = double(DATA_PSORT.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT.topLevel_data.ss_index)));
CS_index = find(logical(double(DATA_PSORT.topLevel_data.cs_index)));
SS_time = ch_time(SS_index);
CS_time = ch_time(CS_index);

SS_inds = repmat(waveform_inds_span(:)', length(SS_index), 1) + repmat(SS_index(:), 1, length(waveform_inds_span));
SS_inds(SS_inds < 1) = 1;
SS_inds(SS_inds > length(ch_data)) = length(ch_data);
CS_inds = repmat(waveform_inds_span(:)', length(CS_index), 1) + repmat(CS_index(:), 1, length(waveform_inds_span));
CS_inds(CS_inds < 1) = 1;
CS_inds(CS_inds > length(ch_data)) = length(ch_data);
SS_waveform = ch_data(SS_inds);
CS_waveform = ch_data(CS_inds);

SS_waveform  = reshape(SS_waveform,[], length(waveform_inds_span));
CS_waveform  = reshape(CS_waveform,[], length(waveform_inds_span));

EPHYS.CH_sorted.SS_data.SS_ind = SS_index;
EPHYS.CH_sorted.CS_data.CS_ind = CS_index;
EPHYS.CH_sorted.SS_data.SS_time = SS_time;
EPHYS.CH_sorted.CS_data.CS_time = CS_time;
EPHYS.CH_sorted.SS_data.SS_waveform = SS_waveform;
EPHYS.CH_sorted.CS_data.CS_waveform = CS_waveform;
EPHYS.CH_sorted.duration = ch_time(end) - ch_time(1);

%% build SSxSS AUTO PROBABILITY
clearvars -except EPHYS BEHAVE
SS_time   = EPHYS.CH_sorted.SS_data.SS_time;
CS_time   = EPHYS.CH_sorted.CS_data.CS_time;
Corr_data = ESN_correlogram(SS_time, CS_time);
EPHYS.CH_sorted.Corr_data = Corr_data;

%% SS & CS train_aligned
clearvars -except EPHYS BEHAVE
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K;
length_time_ = length(EPHYS_time_1K);
CS_time = EPHYS.CH_sorted.CS_data.CS_time;
if isempty(CS_time)
    CS_time = EPHYS_time_1K(1);
end
CS_time(end+1) = max([EPHYS_time_1K(end), CS_time(end)])+1;
SS_time = EPHYS.CH_sorted.SS_data.SS_time;
if isempty(SS_time)
    SS_time = EPHYS_time_1K(1);
end
SS_time(end+1) = max([EPHYS_time_1K(end), SS_time(end)])+1;
EPHYS_CS_train_1K = false(size(EPHYS_time_1K));
EPHYS_SS_train_1K = false(size(EPHYS_time_1K));
counter_CS = find(CS_time >= EPHYS_time_1K(1), 1, 'first');
counter_SS = find(SS_time >= EPHYS_time_1K(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = EPHYS_time_1K(counter_time_point);
    if time_ponit_>=CS_time(counter_CS)
        EPHYS_CS_train_1K(counter_time_point) = true;
        counter_CS = counter_CS + 1;
    end
    if time_ponit_>=SS_time(counter_SS)
        EPHYS_SS_train_1K(counter_time_point) = true;
        counter_SS = counter_SS + 1;
    end
end
EPHYS.CH_EVE.EPHYS_CS_train_1K = EPHYS_CS_train_1K;
EPHYS.CH_EVE.EPHYS_SS_train_1K = EPHYS_SS_train_1K;

%% extract BEHAVE_ind from BEHAVE_EB_xcorr_time_1K
clearvars -except EPHYS BEHAVE
JSP_global_variables();
global data_type_eye_list data_type_BEHAVE_list data_type_neuro_list data_type_EPHYS_list event_type_list inds_span length_trace
if isempty(data_type_eye_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
fprintf(['Analyzing ', EPHYS.CH_sorted_file_name, ' ... ']);
REFRENCE_TIME = EPHYS.CH_EVE.BEHAVE_time_1K;
length_time_  = length(REFRENCE_TIME);
num_sacs      = length(BEHAVE.SACS_ALL_DATA.tag);

% Init variables
% use nan for eye data and logical false for neuro data
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_data_type_eye = 1 : length(data_type_eye_list)
        data_type_eye_name = data_type_eye_list{counter_data_type_eye};
        BEHAVE.SACS_ALL_DATA.([data_type_eye_name '_' event_type_name]) = nan(length_trace, num_sacs);
    end
    for counter_data_type_neuro = 1 : length(data_type_neuro_list)
        data_type_neuro_name = data_type_neuro_list{counter_data_type_neuro};
        BEHAVE.SACS_ALL_DATA.([data_type_neuro_name '_' event_type_name]) = false(length_trace, num_sacs);
    end
end

% Build stream dataset. 
% concatenate the data using the cell2mat and then interpolate the potential missing data
BEHAVE.stream = struct;
time_1K_cell2mat = cell2mat(BEHAVE.TRIALS_DATA.('time_1K')(:))';
BEHAVE.stream.('time_1K') = REFRENCE_TIME;
for counter_variable = 1 : 1 : length(data_type_BEHAVE_list)
    variable_name = data_type_BEHAVE_list{counter_variable};
    variable_cell2mat = cell2mat(BEHAVE.TRIALS_DATA.(variable_name)(:))';
    BEHAVE.stream.(variable_name) = interp1(time_1K_cell2mat, variable_cell2mat, REFRENCE_TIME, 'nearest', 'extrap');
end

% Loop over saccades and add the eye traces to the SACS_ALL_DATA
for counter_sac = 1 : num_sacs
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        event_time = BEHAVE.SACS_ALL_DATA.(['time' '_' event_type_name])(1,counter_sac);
        if isnan(event_time)
            % if the event_time is nan, then skip the event.
            continue;
        end
        event_ind = find(REFRENCE_TIME >= event_time, 1, 'first');
        if isempty(event_ind)
            % if the event_ind is empty, then skip the event.
            BEHAVE.SACS_ALL_DATA.validity(1, counter_sac) = false;
            continue;
        end
        event_inds = repmat( event_ind, 1, length(inds_span)) + repmat(inds_span(:)', length(event_ind), 1);
        event_inds( event_inds < 1 ) = 1;
        event_inds( event_inds > length_time_ ) = length_time_;
        for counter_data_type_eye = 1 : length(data_type_eye_list)
            data_type_eye_name = data_type_eye_list{counter_data_type_eye};
            data_type_BEHAVE_name = data_type_BEHAVE_list{counter_data_type_eye};
            BEHAVE.SACS_ALL_DATA.([data_type_eye_name '_' event_type_name])(:,counter_sac) = ...
                reshape(BEHAVE.stream.(data_type_BEHAVE_name)(event_inds), length_trace, 1);
        end
    end
end

% Loop over saccades and add neuro_SS & neuro_CS to the SACS_ALL_DATA
REFRENCE_TIME = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_time_1K(:); % the initial ind will be drawn from BEHAVE_EB_xcorr_time_1K
length_time_  = length(EPHYS.CH_EVE.EPHYS_time_1K); % the converted ind will be applied to EPHYS_time_1K
for counter_sac = 1 : num_sacs
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        event_time = BEHAVE.SACS_ALL_DATA.(['time' '_' event_type_name])(1,counter_sac);
        if isnan(event_time)
            % if the event_time is nan, then skip the event.
            continue;
        end
        event_ind = find(REFRENCE_TIME >= event_time, 1, 'first');
        if isempty(event_ind)
            % if the event_ind is empty, then skip the event.
            BEHAVE.SACS_ALL_DATA.validity(1, counter_sac) = false;
            continue;
        end
        % set the event_ind_converted based on align_states
        event_ind_converted = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(event_ind);
        % re-write the event_ind_converted to be based on align_photodiode if
        % the event is a cue presentation
        if strcmp(event_type_name,'visual')
            tag_ = BEHAVE.SACS_ALL_DATA.tag(1, counter_sac);
            if (tag_ == 1) || (tag_ == 2) || (tag_ == 3) || (tag_ == 6) || (tag_ == 7)
                % 'prim_success' tag 1 % 'prim_attempt' tag 2 % 'prim_fail' tag 3 % 'back_center_success' tag 6 % 'back_center_prim' tag 7
                event_ind_converted = EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_1K(event_ind);
            end
        end
        if isempty(event_ind_converted)
            % if the event_ind_converted is empty, then skip the event.
            BEHAVE.SACS_ALL_DATA.validity(1, counter_sac) = false;
            continue;
        end
        event_inds_converted = repmat( event_ind_converted, 1, length(inds_span)) + repmat(inds_span(:)', length(event_ind_converted), 1);
        event_inds_converted( event_inds_converted < 1 ) = 1;
        event_inds_converted( event_inds_converted > length_time_ ) = length_time_;
        for counter_data_type_eye = 1 : length(data_type_eye_list)
            data_type_neuro_name = data_type_neuro_list{counter_data_type_eye};
            data_type_EPHYS_name = data_type_EPHYS_list{counter_data_type_eye};
            BEHAVE.SACS_ALL_DATA.([data_type_neuro_name '_' event_type_name])(:,counter_sac) = ...
                logical(reshape(EPHYS.CH_EVE.(data_type_EPHYS_name)(event_inds_converted), length_trace, 1));
        end
    end
end

% Add neuro_CS_count_visual & neuro_CS_count_auditory to BEHAVE.SACS_ALL_DATA
% count the number of the the CS after the visual/auditory event.
% from the event till 200ms or sac onset, whichever happen first
BEHAVE.SACS_ALL_DATA.neuro_CS_count_visual   = zeros(size(BEHAVE.SACS_ALL_DATA.tag));
BEHAVE.SACS_ALL_DATA.neuro_CS_count_auditory = zeros(size(BEHAVE.SACS_ALL_DATA.tag));
for counter_sac = 1 : num_sacs
    ind_start     = (length_trace/2);
    time_onset    = BEHAVE.SACS_ALL_DATA.time_onset(   1, counter_sac);
    time_visual   = BEHAVE.SACS_ALL_DATA.time_visual(  1, counter_sac);
    time_auditory = BEHAVE.SACS_ALL_DATA.time_auditory(1, counter_sac);
    reaction_visual   = round((time_onset - time_visual  )*1000.0);
    reaction_auditory = round((time_onset - time_auditory)*1000.0);
    
    if (reaction_visual > 0) && (reaction_visual < 200)
        ind_end_visual = ind_start + reaction_visual;
    elseif (reaction_visual >= 200)
        ind_end_visual = ind_start + 200;
    elseif (reaction_visual <= 0)
        ind_end_visual = ind_start;
    else
        % this condition covers the nan values
        ind_end_visual = ind_start;
    end
    
    if (reaction_auditory > 0) && (reaction_auditory < 200)
        ind_end_auditory = ind_start + reaction_auditory;
    elseif (reaction_auditory >= 200)
        ind_end_auditory = ind_start + 200;
    elseif (reaction_auditory <= 0)
        ind_end_auditory = ind_start;
    else
        % this condition covers the nan values
        ind_end_auditory = ind_start;
    end
    
    BEHAVE.SACS_ALL_DATA.neuro_CS_count_visual(  1, counter_sac) = nansum(BEHAVE.SACS_ALL_DATA.neuro_CS_visual(  ind_start:ind_end_visual,   counter_sac));
    BEHAVE.SACS_ALL_DATA.neuro_CS_count_auditory(1, counter_sac) = nansum(BEHAVE.SACS_ALL_DATA.neuro_CS_auditory(ind_start:ind_end_auditory, counter_sac));
end

fprintf(' --> Completed. \n')

%% Add Neural_Properties
clearvars -except EPHYS BEHAVE
EPHYS.Neural_Properties = struct();
if length(EPHYS.CH_sorted.SS_data.SS_time)>1
    EPHYS.Neural_Properties.SS_num = length(EPHYS.CH_sorted.SS_data.SS_time);
    EPHYS.Neural_Properties.SS_duration = EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.SS_firing_rate = EPHYS.Neural_Properties.SS_num / EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.SS_time = EPHYS.CH_sorted.SS_data.SS_time;
    EPHYS.Neural_Properties.SS_waveform = nanmean(EPHYS.CH_sorted.SS_data.SS_waveform);
else
    EPHYS.Neural_Properties.SS_num = 0;
    EPHYS.Neural_Properties.SS_duration = EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.SS_firing_rate = 0;
    EPHYS.Neural_Properties.SS_time = [];
    EPHYS.Neural_Properties.SS_waveform = zeros(1, size(EPHYS.CH_sorted.CS_data.CS_waveform, 2));
end

if length(EPHYS.CH_sorted.CS_data.CS_time)>1
    EPHYS.Neural_Properties.CS_num = length(EPHYS.CH_sorted.CS_data.CS_time);
    EPHYS.Neural_Properties.CS_firing_rate = EPHYS.Neural_Properties.CS_num / EPHYS.CH_sorted.duration;
    EPHYS.Neural_Properties.CS_time = EPHYS.CH_sorted.CS_data.CS_time;
    EPHYS.Neural_Properties.CS_waveform = nanmean(EPHYS.CH_sorted.CS_data.CS_waveform);
else
    EPHYS.Neural_Properties.CS_num = 0;
    EPHYS.Neural_Properties.CS_firing_rate = 0;
    EPHYS.Neural_Properties.CS_time = [];
    EPHYS.Neural_Properties.CS_waveform = zeros(1, size(EPHYS.CH_sorted.SS_data.SS_waveform, 2));
end
EPHYS.Neural_Properties.Corr_data_CS_inds_span     = nanmean(EPHYS.CH_sorted.Corr_data.CS_inds_span);
EPHYS.Neural_Properties.Corr_data_CS_bin_size_time = nanmean(EPHYS.CH_sorted.Corr_data.CS_bin_size_time);
EPHYS.Neural_Properties.Corr_data_SS_inds_span     = nanmean(EPHYS.CH_sorted.Corr_data.SS_inds_span);
EPHYS.Neural_Properties.Corr_data_SS_bin_size_time = nanmean(EPHYS.CH_sorted.Corr_data.SS_bin_size_time);
EPHYS.Neural_Properties.Corr_data_SS_SSxSS_AUTO    = nanmean(EPHYS.CH_sorted.Corr_data.SS_SSxSS_AUTO);
EPHYS.Neural_Properties.Corr_data_CS_CSxSS_AUTO    = nanmean(EPHYS.CH_sorted.Corr_data.CS_CSxSS_AUTO);

%% outputs
SACS_ALL_DATA = BEHAVE.SACS_ALL_DATA;
Neural_Properties = EPHYS.Neural_Properties;
EXPERIMENT_PARAMS = BEHAVE.EXPERIMENT_PARAMS;

end
%% function combine_sac_files()
function combine_sac_files()
clc; close all;
flag_pair_list = false; % This should be false. DO NOT change it to true
pCell_list = PH_build_pCell_list(flag_pair_list);
path_data_monkey_sorted = uigetdir;
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
ALL_PCELL_name = ['ALL_PCELL_' num2str(num_pCells)];
if ~exist([path_data_monkey_sorted ALL_PCELL_name filesep 'cell_data'], 'dir')
    mkdir([path_data_monkey_sorted ALL_PCELL_name filesep 'cell_data']);
end
if ~exist([path_data_monkey_sorted ALL_PCELL_name filesep 'cell_figs_behave'], 'dir')
    mkdir([path_data_monkey_sorted ALL_PCELL_name filesep 'cell_figs_behave']);
end
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    num_recording = nansum(pCell_list_isstr(counter_pCell, :));
    clearvars data_recordings
    %% Loop over recordings
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name_cell = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name_cell(1:2);
        month_ = file_name_cell(3:4);
        day_ = file_name_cell(5:6);
        hour_ = file_name_cell(8:9);
        minute_ = file_name_cell(10:11);
        second_ = file_name_cell(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_data' filesep];
        %% load recording data
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
        data_recording = load([file_path file_name_cell(1:end-4) '_sac' file_name_cell(end-3:end) '.mat'], ...
            'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS');
        if file_name_cell(18) == 's'
            data_recording.id          = file_name_cell(1:16);
        elseif isstrprop(file_name_cell(18),'digit')
            data_recording.id          = file_name_cell(1:18);
        else
            error('Build plot_data_compress: cell id does not follow the standards')
        end
        data_recordings(counter_recording) = data_recording;
    end
    
    %% Init data_cell
    data_cell = struct;
    
    %% id
    data_cell.id = cell(num_recording, 1);
    for counter_recording = 1 : 1 : num_recording
        data_cell.id{counter_recording, 1} = data_recordings(counter_recording).id;
    end
    
    %% EXPERIMENT_PARAMS
    field_names_EXPERIMENT_PARAMS = fieldnames(data_recordings(1).EXPERIMENT_PARAMS);
    for counter_field = 1 : 1 : length(field_names_EXPERIMENT_PARAMS)
        field_name_EXPERIMENT_PARAMS = field_names_EXPERIMENT_PARAMS{counter_field};
        data_cell.EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS) = cell(num_recording, 1);
        for counter_recording = 1 : 1 : num_recording
            data_cell.EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS){counter_recording, 1} = data_recordings(counter_recording).EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS);
        end
    end
    
    %% Neural_Properties
    % Init variables
    data_cell.Neural_Properties.SS_num = 0;
    data_cell.Neural_Properties.SS_duration = 0;
    data_cell.Neural_Properties.SS_firing_rate = 0;
    data_cell.Neural_Properties.CS_num = 0;
    data_cell.Neural_Properties.CS_firing_rate = 0;
    data_cell.Neural_Properties.SS_time = [];
    data_cell.Neural_Properties.CS_time = [];
    variable_names = {'SS_waveform', 'CS_waveform', ...
        'Corr_data_CS_inds_span', 'Corr_data_CS_bin_size_time', 'Corr_data_SS_inds_span', 'Corr_data_SS_bin_size_time', ...
        'Corr_data_SS_SSxSS_AUTO','Corr_data_CS_CSxSS_AUTO'};
    for counter_variable = 1 : length(variable_names)
        variable_name = variable_names{counter_variable};
        data_cell.Neural_Properties.(variable_name) = zeros(size(data_recordings(1).Neural_Properties.(variable_name)));
    end
    % Loop over recordings
    for counter_recording = 1 : 1 : num_recording
        data_cell.Neural_Properties.SS_num = data_cell.Neural_Properties.SS_num + ...
            data_recordings(counter_recording).Neural_Properties.SS_num;
        data_cell.Neural_Properties.SS_duration = data_cell.Neural_Properties.SS_duration + ...
            data_recordings(counter_recording).Neural_Properties.SS_duration;
        data_cell.Neural_Properties.CS_num = data_cell.Neural_Properties.CS_num + ...
            data_recordings(counter_recording).Neural_Properties.CS_num;
        
        data_cell.Neural_Properties.SS_time = vertcat(data_cell.Neural_Properties.SS_time ,...
            data_recordings(counter_recording).Neural_Properties.SS_time);
        data_cell.Neural_Properties.CS_time = vertcat(data_cell.Neural_Properties.CS_time, ...
            data_recordings(counter_recording).Neural_Properties.CS_time);
        % compute weighted average for these variables
        for counter_variable = 1 : length(variable_names)
            variable_name = variable_names{counter_variable};
            if contains(variable_name, 'Corr_data_SS') || contains(variable_name, 'SS_waveform')
                num_spike = data_recordings(counter_recording).Neural_Properties.SS_num;
            elseif contains(variable_name, 'Corr_data_CS') || contains(variable_name, 'CS_waveform')
                num_spike = data_recordings(counter_recording).Neural_Properties.CS_num;
            end
            data_cell.Neural_Properties.(variable_name) = data_cell.Neural_Properties.(variable_name) + ...
                ( data_recordings(counter_recording).Neural_Properties.(variable_name) .* num_spike);
        end
    end
    data_cell.Neural_Properties.SS_firing_rate = ...
            data_cell.Neural_Properties.SS_num ./ data_cell.Neural_Properties.SS_duration;
    data_cell.Neural_Properties.CS_firing_rate = ...
            data_cell.Neural_Properties.CS_num ./ data_cell.Neural_Properties.SS_duration;
    % devide by number of events to compute the weigted average
    for counter_variable = 1 : length(variable_names)
        variable_name = variable_names{counter_variable};
        if contains(variable_name, 'Corr_data_SS') || contains(variable_name, 'SS_waveform')
            num_spike = data_cell.Neural_Properties.SS_num;
        elseif contains(variable_name, 'Corr_data_CS') || contains(variable_name, 'CS_waveform')
            num_spike = data_cell.Neural_Properties.CS_num;
        end
        data_cell.Neural_Properties.(variable_name) = data_cell.Neural_Properties.(variable_name) ./ num_spike;
    end
    
    %% SACS_ALL_DATA
    field_names_SACS_ALL_DATA = fieldnames(data_recordings(1).SACS_ALL_DATA);
    for counter_field = 1 : 1 : length(field_names_SACS_ALL_DATA)
        field_name_SACS_ALL_DATA = field_names_SACS_ALL_DATA{counter_field};
        data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = [];
        for counter_recording = 1 : 1 : num_recording
            data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = horzcat(...
                data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA), ...
                data_recordings(counter_recording).SACS_ALL_DATA.(field_name_SACS_ALL_DATA));
        end
    end
    
    %% Save data_cell
    id = data_cell.id;
    EXPERIMENT_PARAMS = data_cell.EXPERIMENT_PARAMS;
    Neural_Properties = data_cell.Neural_Properties;
    SACS_ALL_DATA = data_cell.SACS_ALL_DATA;
    cell_name = [data_cell.id{1} '_' 'combine' '_' num2str(length(data_cell.id))];
    fprintf([cell_name ': Saving .mat file ...'])
    save([path_data_monkey_sorted ALL_PCELL_name filesep 'cell_data' filesep cell_name '.mat'],...
        'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS', 'id','-v7.3')
    fprintf(' --> Completed. \n')
    
    %% plot_sac_sorter
    params.cell_name    = cell_name;
    params.duration     = Neural_Properties.SS_duration;
    params.sac_tag_list = EXPERIMENT_PARAMS.sac_tag_list{1};
    params.num_trials   = nansum(cell2mat(EXPERIMENT_PARAMS.num_trials));
    fprintf([cell_name ': Saving .png plot ...'])
    plot_sac_sorter(SACS_ALL_DATA, params)
    hFig_ = gcf;
%     saveas(hFig_,[path_data_monkey_sorted ALL_PCELL_name filesep 'cell_figs_behave' filesep cell_name], 'pdf');
    saveas(hFig_,[path_data_monkey_sorted ALL_PCELL_name filesep 'cell_figs_behave' filesep cell_name], 'png');
    close(hFig_)
    fprintf(' --> Completed. \n')
    
end
fprintf('### ALL DONE. ###\n')
end
%% function combine_pair_files()
function combine_pair_files()
clc; close all;
flag_pair_list = true; % This should be true. DO NOT change it to false
pCell_list = PH_build_pCell_list(flag_pair_list);
path_data_monkey_sorted = uigetdir;
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);
ALL_PCELL_name = ['ALL_PCELL_' num2str(num_pCells)];
if ~exist([path_data_monkey_sorted ALL_PCELL_name filesep 'pair_data'], 'dir')
    mkdir([path_data_monkey_sorted ALL_PCELL_name filesep 'pair_data']);
end
if ~exist([path_data_monkey_sorted ALL_PCELL_name filesep 'pair_figs_behave'], 'dir')
    mkdir([path_data_monkey_sorted ALL_PCELL_name filesep 'pair_figs_behave']);
end
%% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    num_recording = nansum(pCell_list_isstr(counter_pCell, :));
    clearvars data_recordings
    %% Loop over recordings
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name_cell = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        year_ = file_name_cell(1:2);
        month_ = file_name_cell(3:4);
        day_ = file_name_cell(5:6);
        hour_ = file_name_cell(8:9);
        minute_ = file_name_cell(10:11);
        second_ = file_name_cell(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_figs = ['analyzed_data' filesep];
        %% load recording data
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_figs];
        if ~strcmp(file_path(end), filesep);file_path = [file_path filesep];end
        data_recording = load([file_path file_name_cell(1:end-4) '_sac' file_name_cell(end-3:end) '.mat'], ...
            'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS');
        if file_name_cell(18) == 's'
            data_recording.id          = file_name_cell(1:16);
        elseif file_name_cell(18) == '2'
            data_recording.id          = file_name_cell(1:18);
        else
            error('Build plot_data_compress: cell id does not follow the standards')
        end
        data_recordings(counter_recording) = data_recording;
    end
    
    %% Init data_cell
    data_cell = struct;
    
    %% id
    data_cell.id = cell(num_recording, 1);
    for counter_recording = 1 : 1 : num_recording
        data_cell.id{counter_recording, 1} = data_recordings(counter_recording).id;
    end
    
    %% EXPERIMENT_PARAMS
    field_names_EXPERIMENT_PARAMS = fieldnames(data_recordings(1).EXPERIMENT_PARAMS);
    for counter_field = 1 : 1 : length(field_names_EXPERIMENT_PARAMS)
        field_name_EXPERIMENT_PARAMS = field_names_EXPERIMENT_PARAMS{counter_field};
        data_cell.EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS) = cell(num_recording, 1);
        for counter_recording = 1 : 1 : num_recording
            data_cell.EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS){counter_recording, 1} = data_recordings(counter_recording).EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS);
        end
    end
    
    %% Neural_Properties
    % Init variables
    data_cell.Neural_Properties.SS_num = 0;
    data_cell.Neural_Properties.SS_duration = 0;
    data_cell.Neural_Properties.SS_firing_rate = 0;
    data_cell.Neural_Properties.CS_num = 0;
    data_cell.Neural_Properties.CS_firing_rate = 0;
    data_cell.Neural_Properties.SS_time = [];
    data_cell.Neural_Properties.CS_time = [];
    variable_names = {'SS_waveform', 'CS_waveform', ...
        'Corr_data_CS_inds_span', 'Corr_data_CS_bin_size_time', 'Corr_data_SS_inds_span', 'Corr_data_SS_bin_size_time', ...
        'Corr_data_SS_SSxSS_AUTO','Corr_data_CS_CSxSS_AUTO'};
    for counter_variable = 1 : length(variable_names)
        variable_name = variable_names{counter_variable};
        data_cell.Neural_Properties.(variable_name) = zeros(size(data_recordings(1).Neural_Properties.(variable_name)));
    end
    % Loop over recordings
    for counter_recording = 1 : 1 : num_recording
        data_cell.Neural_Properties.SS_num = data_cell.Neural_Properties.SS_num + ...
            data_recordings(counter_recording).Neural_Properties.SS_num;
        data_cell.Neural_Properties.SS_duration = data_cell.Neural_Properties.SS_duration + ...
            data_recordings(counter_recording).Neural_Properties.SS_duration;
        data_cell.Neural_Properties.CS_num = data_cell.Neural_Properties.CS_num + ...
            data_recordings(counter_recording).Neural_Properties.CS_num;
        
        data_cell.Neural_Properties.SS_time = vertcat(data_cell.Neural_Properties.SS_time ,...
            data_recordings(counter_recording).Neural_Properties.SS_time);
        data_cell.Neural_Properties.CS_time = vertcat(data_cell.Neural_Properties.CS_time, ...
            data_recordings(counter_recording).Neural_Properties.CS_time);
        % compute weighted average for these variables
        for counter_variable = 1 : length(variable_names)
            variable_name = variable_names{counter_variable};
            if contains(variable_name, 'Corr_data_SS') || contains(variable_name, 'SS_waveform')
                num_spike = data_recordings(counter_recording).Neural_Properties.SS_num;
            elseif contains(variable_name, 'Corr_data_CS') || contains(variable_name, 'CS_waveform')
                num_spike = data_recordings(counter_recording).Neural_Properties.CS_num;
            end
            data_cell.Neural_Properties.(variable_name) = data_cell.Neural_Properties.(variable_name) + ...
                ( data_recordings(counter_recording).Neural_Properties.(variable_name) .* num_spike);
        end
    end
    data_cell.Neural_Properties.SS_firing_rate = ...
            data_cell.Neural_Properties.SS_num ./ data_cell.Neural_Properties.SS_duration;
    data_cell.Neural_Properties.CS_firing_rate = ...
            data_cell.Neural_Properties.CS_num ./ data_cell.Neural_Properties.SS_duration;
    % devide by number of events to compute the weigted average
    for counter_variable = 1 : length(variable_names)
        variable_name = variable_names{counter_variable};
        if contains(variable_name, 'Corr_data_SS') || contains(variable_name, 'SS_waveform')
            num_spike = data_cell.Neural_Properties.SS_num;
        elseif contains(variable_name, 'Corr_data_CS') || contains(variable_name, 'CS_waveform')
            num_spike = data_cell.Neural_Properties.CS_num;
        end
        data_cell.Neural_Properties.(variable_name) = data_cell.Neural_Properties.(variable_name) ./ num_spike;
    end
    
    %% SACS_ALL_DATA
    field_names_SACS_ALL_DATA = fieldnames(data_recordings(1).SACS_ALL_DATA);
    for counter_field = 1 : 1 : length(field_names_SACS_ALL_DATA)
        field_name_SACS_ALL_DATA = field_names_SACS_ALL_DATA{counter_field};
        data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = [];
        for counter_recording = 1 : 1 : num_recording
            data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = horzcat(...
                data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA), ...
                data_recordings(counter_recording).SACS_ALL_DATA.(field_name_SACS_ALL_DATA));
        end
    end
    
    %% Save data_cell
    id = data_cell.id;
    EXPERIMENT_PARAMS = data_cell.EXPERIMENT_PARAMS;
    Neural_Properties = data_cell.Neural_Properties;
    SACS_ALL_DATA = data_cell.SACS_ALL_DATA;
    cell_name = [data_cell.id{1} '_' 'combine' '_' num2str(length(data_cell.id))];
    fprintf([cell_name ': Saving .mat file ...'])
    save([path_data_monkey_sorted ALL_PCELL_name filesep 'pair_data' filesep cell_name '.mat'],...
        'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS', 'id','-v7.3')
    fprintf(' --> Completed. \n')
    
    %% plot_sac_sorter
    params.cell_name    = cell_name;
    params.duration     = Neural_Properties.SS_duration;
    params.sac_tag_list = EXPERIMENT_PARAMS.sac_tag_list{1};
    params.num_trials   = nansum(cell2mat(EXPERIMENT_PARAMS.num_trials));
    fprintf([cell_name ': Saving .png plot ...'])
    plot_sac_sorter(SACS_ALL_DATA, params)
    hFig_ = gcf;
%     saveas(hFig_,[path_data_monkey_sorted ALL_PCELL_name filesep 'pair_figs_behave' filesep cell_name], 'pdf');
    saveas(hFig_,[path_data_monkey_sorted ALL_PCELL_name filesep 'pair_figs_behave' filesep cell_name], 'png');
    close(hFig_)
    fprintf(' --> Completed. \n')
    
end
fprintf('### ALL DONE. ###\n')
end

%% RE-RUN PCELL FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function CS_on_analysis()
function CS_on_analysis()
clc; close all;
flag_pair_list = false;
JSP_global_variables(flag_pair_list);
global ang_edges ang_values range_cell_with_4dir_behave
if isempty(ang_edges)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
pCell_ids = build_pCell_ids(flag_pair_list);
num_pCells = size(pCell_ids, 1);
clearvars CS_on_population
%% Loop over pCells
for counter_pCell = 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    %% load SACS_ALL_DATA
    cell_file_name = pCell_ids{counter_pCell, 1};
    load([path_cell_data cell_file_name], 'SACS_ALL_DATA');
    
    %% compute CS-on
    visual_px_offset = SACS_ALL_DATA.visual_px_offset;
    visual_py_offset = SACS_ALL_DATA.visual_py_offset;
    eye_r_px_onset   = SACS_ALL_DATA.eye_r_px_onset;
    eye_r_py_onset   = SACS_ALL_DATA.eye_r_py_onset;
    delta_x = visual_px_offset - eye_r_px_onset;
    delta_y = visual_py_offset - eye_r_py_onset;
    visual_ang = wrapTo360(atan2d(delta_y, delta_x));
%     visual_ang = wrapTo360(SACS_ALL_DATA.eye_r_ang);
    neuro_CS_count = nansum((SACS_ALL_DATA.neuro_SS_visual(250:325,:)));
%     neuro_CS_count = SACS_ALL_DATA.neuro_CS_count_visual;

    visual_ang_bin = discretize(visual_ang, ang_edges);
    
    last_bin_id = length(ang_edges) - 1;
    visual_ang_bin(visual_ang_bin == last_bin_id) = 1; % wrap the circle around
    % 1: 0deg % 2: 45deg % 3: 90deg % 4: 135deg % 5: 180deg % 6: 225deg % 7: 270deg % 8: 315deg
    num_ang_bin = length(ang_edges) - 2;
    tag_bin = unique(SACS_ALL_DATA.tag);
    num_tag_bin = length(tag_bin);
    CS_count  = zeros(num_tag_bin, num_ang_bin);
    sac_count = zeros(num_tag_bin, num_ang_bin);
    CS_prob   = zeros(num_tag_bin, num_ang_bin);
    for counter_tag = 1 : num_tag_bin
        for counter_ang = 1 : num_ang_bin
            idx_tag = (SACS_ALL_DATA.tag == tag_bin(counter_tag));
            idx_ang = (visual_ang_bin == counter_ang);
            flag_last_cue = SACS_ALL_DATA.flag_last_cue;
            idx_ = idx_tag&idx_ang&flag_last_cue;
            CS_count(counter_tag, counter_ang) = nansum(neuro_CS_count(1, idx_));
            sac_count(counter_tag, counter_ang) = nansum(idx_);
            CS_prob(counter_tag, counter_ang) = CS_count(counter_tag, counter_ang) ./ sac_count(counter_tag, counter_ang);
        end
    end
    
    r = nansum(CS_prob.* repmat(exp(1i*deg2rad(ang_values)), num_tag_bin, 1) , 2); % compute weighted sum of cos and sin of angles
    CS_ang =  wrapTo360(rad2deg(angle(r))); % Computes the mean direction for circular data.
    CS_prob_sum = nansum(CS_prob,2); % sum of weights
    CS_rho = abs(r) ./ CS_prob_sum; % Computes mean resultant vector length for circular data.
    
%     CS_count_avg  = CS_count( 1, :) + CS_count( 4, :) + CS_count( 8, :); % CS_count( 8, :); % 
%     sac_count_avg = sac_count(1, :) + sac_count(4, :) + sac_count(8, :); % sac_count(8, :); % 
    CS_count_avg  = CS_count( 1, :); % 
    sac_count_avg = sac_count(1, :); % 
    
    
    % 'prim_success' tag 1 % 'corr_success' tag 4 % 'back_center_irrelev' tag 8
    CS_prob_avg = CS_count_avg ./ sac_count_avg;
    
    r_avg = nansum(CS_prob_avg.* exp(1i*deg2rad(ang_values)) , 2); % compute weighted sum of cos and sin of angles
    CS_ang_avg =  wrapTo360(rad2deg(angle(r_avg)));
%      CS_ang_avg =  135; % hard code CS on

    CS_rho_avg = abs(r_avg) ./ nansum(CS_prob_avg,2);
    
    if (counter_pCell >= range_cell_with_4dir_behave(1)) && (counter_pCell <= range_cell_with_4dir_behave(2))
        idx_CS_on_dir = discretize( ESN_Round(CS_ang_avg, 90.0, 'round'), ang_edges);
    else
        idx_CS_on_dir = discretize(CS_ang_avg, ang_edges);
    end
    
    if idx_CS_on_dir == last_bin_id; idx_CS_on_dir = 1; end
    
    idx_ = idx_CS_on_dir - 1; % make it 0-index format
    if (idx_ == 8); idx_ = 0; end
    idx_CS_tuned = mod((idx_ : 1 : idx_+7), 8) + 1;
    CS_prob_avg_tuned = CS_prob_avg(idx_CS_tuned);
    
    %% Von Mises
    if CS_rho_avg < 0.53
        vonMises_kappa = 2*CS_rho_avg + CS_rho_avg^3 + 5*CS_rho_avg^5/6;
    elseif CS_rho_avg>=0.53 && CS_rho_avg<0.85
        vonMises_kappa = -.4 + 1.39*CS_rho_avg + 0.43/(1-CS_rho_avg);
    else
        vonMises_kappa = 1/(CS_rho_avg^3 - 4*CS_rho_avg^2 + 3*CS_rho_avg);
    end
    
    % evaluate pdf
    vonMises_C = 1/(2*pi*besseli(0,vonMises_kappa));
    vonMises_pdf = vonMises_C * exp(vonMises_kappa*cosd(ang_values-CS_ang_avg));
    vonMises_var = 1 - (besseli(1,vonMises_kappa) / besseli(0,vonMises_kappa));
    vonMises_std = wrapTo360(rad2deg(sqrt(vonMises_var)));
    
    %% Build CS_on_data
    CS_on_data.sac_count = sac_count;
    CS_on_data.CS_count  = CS_count;
    CS_on_data.CS_prob = CS_prob;
    CS_on_data.CS_ang  = CS_ang;
    CS_on_data.CS_rho  = CS_rho;
    CS_on_data.CS_prob_avg = CS_prob_avg;
    CS_on_data.CS_prob_avg_tuned = CS_prob_avg_tuned;
    CS_on_data.CS_ang_avg  = CS_ang_avg;
    CS_on_data.CS_rho_avg  = CS_rho_avg;
    CS_on_data.vonMises_kappa  = vonMises_kappa;
    CS_on_data.vonMises_var  = vonMises_var;
    CS_on_data.vonMises_std  = vonMises_std;
    CS_on_data.vonMises_pdf  = vonMises_pdf;
    CS_on_data.idx_CS_on_dir  = idx_CS_on_dir;
    CS_on_data.idx_CS_tuned   = idx_CS_tuned;
    CS_on_data.visual_ang_bin       = visual_ang_bin;
    
    %% Append CS_on_data to cell_data
    save([path_cell_data cell_file_name], 'CS_on_data', '-append');
    CS_on_population(counter_pCell) = CS_on_data;
end
fprintf('### ALL DONE. ###\n')
save([path_cell_data '..' filesep 'CS_on_population.mat'], 'CS_on_population');
end
%% function sac_modulation_index()
function sac_modulation_index()
%% Description
% For statistical testing, we Z-scored the spike counts. 
% This is done by first subtracting the observed spike count by expected spike count then dividing the resulting value by expected standard deviation in spike count. 
% Assuming the spike count in each bin is described by a binomial distribution, if the probability of synchronized spikes by chance is p=N1*N2/(T/dT)^2, the expected standard deviation in spike count is sqrt(p(1-p)*T/dT).
% When assessing statistical significance of the synchrony using the average spike count within ±1ms, we further adjusted this standard deviation used for the Z-score by a factor of sqrt(dT/2). 
% Pairs with an average spike count within ±1 ms from origin above 3 in a Z-scored spike count (corresponding to an alpha-value of 0.00135) were considered to exhibit a significant degree of synchrony.

%%
clc; close all;
flag_pair_list = true;
JSP_global_variables(flag_pair_list);
global length_trace
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
pCell_ids = build_pCell_ids(flag_pair_list);
num_pCells = size(pCell_ids, 1);
counter_not_modulated = 1;
%% Loop over pCells
pCell_z_score = nan(num_pCells, 1);
for counter_pCell = 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    %% load SACS_ALL_DATA
    cell_file_name = pCell_ids{counter_pCell, 1};
    load([path_cell_data cell_file_name], 'SACS_ALL_DATA', 'Neural_Properties');
    %% Compute SS_baseline, SS_baseline_stdv
    SS_time = Neural_Properties.SS_time;
    SS_ISI = diff(SS_time);
    SS_ISI(abs(SS_ISI)>5)=0;
    SS_baseline = (length(SS_ISI)+1) ./ sum(SS_ISI) * 0.001; % 1ms probability
    
    %% Compute modulation_z_score
    neuro_SS_onset = SACS_ALL_DATA.neuro_SS_onset;
    idx_tag = false(size(SACS_ALL_DATA.tag));
    idx_tag = idx_tag | (SACS_ALL_DATA.tag == 1); % 'prim_success' tag 1
    idx_tag = idx_tag | (SACS_ALL_DATA.tag == 4); % 'corr_success' tag 4
    idx_tag = idx_tag | (SACS_ALL_DATA.tag == 6); % 'back_center_success' tag 6
    idx_tag = idx_tag | (SACS_ALL_DATA.tag == 7); % 'back_center_prim' tag 7
    idx_tag = idx_tag | (SACS_ALL_DATA.tag == 8); % 'back_center_irrelev' tag 8
    eye_r_amp_m = SACS_ALL_DATA.eye_r_amp_m;
    median_amp = median(eye_r_amp_m(idx_tag));
    idx_amp = (eye_r_amp_m > median_amp);
    idx_ = idx_tag & idx_amp;
    
    if length_trace ~= 500
        error('sac_modulation_index: length_trace is not 500. Please modify the code.')
    end
    
    SS_train = neuro_SS_onset(:, idx_);
    
    num_sacs = sum(idx_);
    num_perm = 2000;
    modulation_range_perm = nan(1, num_perm);
    modulation_trace_perm = nan(length_trace, num_perm);
    for counter_perm = 1 : num_perm
        idx_perm_ = randi(num_sacs, 1, num_sacs);
        modulation_trace_perm_ = smooth(nanmean(SS_train(:,idx_perm_),2), 3, 'sgolay', 2);
        modulation_range_perm_ = max(modulation_trace_perm_(201:350,1)) - min(modulation_trace_perm_(201:350,1));
        modulation_trace_perm(:, counter_perm) = modulation_trace_perm_;
        modulation_range_perm(:, counter_perm) = modulation_range_perm_;
    end
    
    modulation_trace_mean = ESN_smooth(nanmean(SS_train,2));
    modulation_trace_stdv = nanstd(modulation_trace_perm,0,2);
    modulation_range_mean = max(modulation_trace_mean(201:350,1)) - min(modulation_trace_mean(201:350,1));
    modulation_range_stdv = nanstd(modulation_range_perm,0,2);
    modulation_z_score = modulation_range_mean ./ modulation_range_stdv;
    modulation_baseline = nanmean(nanmean(SS_train(26:125,:),2));
    
    modulation_thresold = 3.0;
    is_modulated = (modulation_z_score > modulation_thresold);
    
    if ~is_modulated
        fprintf([num2str(counter_not_modulated) '. ' 'z-score: ' num2str(modulation_z_score) '\n']);
        counter_not_modulated = counter_not_modulated + 1;
    end
    
    %% Build sac_modulation
    sac_modulation.SS_baseline           = SS_baseline;
    sac_modulation.modulation_trace_mean = modulation_trace_mean;
    sac_modulation.modulation_trace_stdv = modulation_trace_stdv;
    sac_modulation.modulation_range_mean = modulation_range_mean;
    sac_modulation.modulation_range_stdv = modulation_range_stdv;
    sac_modulation.modulation_baseline   = modulation_baseline;
    sac_modulation.modulation_z_score    = modulation_z_score;
    sac_modulation.is_modulated          = is_modulated;
    pCell_z_score(counter_pCell, 1) = modulation_z_score;
    %% Append sac_modulation to cell_data
    save([path_cell_data cell_file_name], 'sac_modulation', '-append');
    
end
fprintf('### ALL DONE. ###\n')
figure();
histogram(pCell_z_score, 0:1:(max(pCell_z_score)+1));
end

%% BUILD FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function build_neural_properties()
function build_neural_properties()
clc; close all;
flag_pair_list = false;
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
pCell_ids = build_pCell_ids(flag_pair_list);
num_pCells = size(pCell_ids, 1);

%% Init variables
cell_file_name = pCell_ids{1, 1};
load([path_cell_data cell_file_name], 'Neural_Properties', 'CS_on_data');
population_neural_properties = struct;
field_names_Neural_Properties = fieldnames(Neural_Properties);
for counter_field_Neural_Properties = 1 : length(field_names_Neural_Properties)
    field_name_Neural_Properties = field_names_Neural_Properties{counter_field_Neural_Properties};
    if strcmp(field_name_Neural_Properties, 'SS_time') || strcmp(field_name_Neural_Properties, 'CS_time')
        % skip the SS_time and CS_time since their size are varibale
        continue;
    end
    population_neural_properties.(field_name_Neural_Properties) = nan(num_pCells, size(Neural_Properties.(field_name_Neural_Properties), 2) );
end
field_names_CS_on = fieldnames(CS_on_data);
for counter_field_CS_on = 1 : length(field_names_CS_on)
    field_name_CS_on = field_names_CS_on{counter_field_CS_on};
    if strcmp(field_name_CS_on, 'visual_ang_bin') || strcmp(field_name_CS_on, 'visual_ang_bin_tuned')
        % skip the visual_ang_bin and visual_ang_bin_tuned since their size are varibale
        continue;
    end
    if size(CS_on_data.(field_name_CS_on), 1) > 1
        population_neural_properties.(field_name_CS_on) = nan(num_pCells, size(CS_on_data.(field_name_CS_on), 2), size(CS_on_data.(field_name_CS_on), 1));
    elseif size(CS_on_data.(field_name_CS_on), 1) == 1
        population_neural_properties.(field_name_CS_on) = nan(num_pCells, size(CS_on_data.(field_name_CS_on), 2));
    end
end
%% Loop over pCells
for counter_pCell = 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    %% load SACS_ALL_DATA
    cell_file_name = pCell_ids{counter_pCell, 1};
    load([path_cell_data cell_file_name], 'Neural_Properties', 'CS_on_data');
    %%
    for counter_field_Neural_Properties = 1 : length(field_names_Neural_Properties)
        field_name_Neural_Properties = field_names_Neural_Properties{counter_field_Neural_Properties};
        if strcmp(field_name_Neural_Properties, 'SS_time') || strcmp(field_name_Neural_Properties, 'CS_time')
            % skip the SS_time and CS_time since their size are varibale
            continue;
        end
        population_neural_properties.(field_name_Neural_Properties)(counter_pCell, :) = ...
            Neural_Properties.(field_name_Neural_Properties);
    end
    for counter_field_CS_on = 1 : length(field_names_CS_on)
        field_name_CS_on = field_names_CS_on{counter_field_CS_on};
        if strcmp(field_name_CS_on, 'visual_ang_bin') || strcmp(field_name_CS_on, 'visual_ang_bin_tuned')
            % skip the visual_ang_bin and visual_ang_bin_tuned since their size are varibale
            continue;
        end
        
%         account for missing tag
        if size(population_neural_properties.(field_name_CS_on)(counter_pCell, :, :), 3) ~= ...
                size(CS_on_data.(field_name_CS_on),1)
          CS_on_data.(field_name_CS_on) = [CS_on_data.(field_name_CS_on) ;zeros(1, size(CS_on_data.(field_name_CS_on),2))];          
        end  
    
        if size(CS_on_data.(field_name_CS_on), 1) > 1
            population_neural_properties.(field_name_CS_on)(counter_pCell, :, :) = ...
                CS_on_data.(field_name_CS_on)';
        elseif size(CS_on_data.(field_name_CS_on), 1) == 1
            population_neural_properties.(field_name_CS_on)(counter_pCell, :) = ...
                CS_on_data.(field_name_CS_on);
        end
    end
end

%% Add CS_suupersion_time
CS_xprob_pCells = population_neural_properties.Corr_data_CS_CSxSS_AUTO;
CS_xprob_suppression = CS_xprob_pCells ./ nanmean(CS_xprob_pCells(:,20:50), 2);
[~,idx] = max(CS_xprob_suppression(:,56:end)>0.63, [], 2);
population_neural_properties.CS_suppression_time = idx+5;

%% Save data
fprintf(['Saving .mat files' ' ...'])
save([path_cell_data '..' filesep 'population_neural_properties' '.mat'], 'population_neural_properties', '-v7.3');
fprintf(' --> Completed. \n')

end
%% function build_population_data()
function build_population_data()
clc; close all;
flag_pair_list = false;    % This should be false. DO NOT change it to true
JSP_global_variables(flag_pair_list);
global event_type_list amp_edges vel_edges length_trace %ang_values ang_edges
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
pCell_ids = build_pCell_ids(flag_pair_list);
num_pCells = size(pCell_ids, 1);
cell_file_name = pCell_ids{1, 1};
load([path_cell_data cell_file_name], 'SACS_ALL_DATA', 'CS_on_data');
tag_bin = unique(SACS_ALL_DATA.tag);
num_tag_bin = length(tag_bin);
% num_ang_bin = length(ang_values);
num_ang_bin = length(CS_on_data.idx_CS_tuned);
num_amp_bin = length(amp_edges) - 1;
num_vel_bin = length(vel_edges) - 1;
flag_build_absol = true;
flag_build_tuned = false;
if ~xor(flag_build_absol, flag_build_tuned)
    fprintf('><ERROR><: build_population_data, Please select either flag_build_absol or flag_build_tuned. Not both. The code will run much faster this way.\n')
    return;
end

%% Init variables
% STRUCT (SS / CS / VM) -> 10x1 tag struct (amp / vel) -> variable name (onset / vmax / offset / visual) -> 6x8 cell (ampXang / velXang) -> 138x500 double (pCellXtrace)
fprintf(['Initializing population_data' ' ...'])
for counter_tag = 1 : num_tag_bin
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        if flag_build_absol
        SS_population_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        CS_population_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        VM_population_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        VT_population_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        num_sac_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        
        SS_population_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        CS_population_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        VM_population_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        VT_population_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        num_sac_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        end
        if flag_build_tuned
        SS_population_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        CS_population_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        VM_population_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        VT_population_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        num_sac_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        
        SS_population_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        CS_population_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        VM_population_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        VT_population_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        num_sac_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        end
        for counter_ang = 1 : num_ang_bin
            for counter_amp = 1 : num_amp_bin
                if flag_build_absol
                SS_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                CS_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                VM_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                VT_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                num_sac_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = zeros(num_pCells, 1);
                end
                if flag_build_tuned
                SS_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                CS_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                VM_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                VT_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                num_sac_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = zeros(num_pCells, 1);
                end
            end
            for counter_vel = 1 : num_vel_bin
                if flag_build_absol
                SS_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                CS_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                VM_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                VT_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                num_sac_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = zeros(num_pCells, 1);
                end
                if flag_build_tuned
                SS_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                CS_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                VM_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                VT_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                num_sac_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = zeros(num_pCells, 1);
                end
            end
        end
    end
end
fprintf(' --> Completed. \n')

%% Loop over pCells
for counter_pCell = 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    %% load SACS_ALL_DATA
    cell_file_name = pCell_ids{counter_pCell, 1};
    load([path_cell_data cell_file_name], 'SACS_ALL_DATA', 'CS_on_data');
    
    %% Compute data
    SACS_amp_bin = discretize(SACS_ALL_DATA.eye_r_amp_m,  amp_edges);
    SACS_vel_bin = discretize(SACS_ALL_DATA.eye_r_vm_max, vel_edges);
    eye_amp_x = repmat(SACS_ALL_DATA.eye_r_amp_x, length_trace, 1);
    eye_amp_y = repmat(SACS_ALL_DATA.eye_r_amp_y, length_trace, 1);
    eye_amp_m = repmat(SACS_ALL_DATA.eye_r_amp_m, length_trace, 1);
    
%     visual_px_offset = SACS_ALL_DATA.visual_px_offset;
%     visual_py_offset = SACS_ALL_DATA.visual_py_offset;
%     eye_r_px_onset   = SACS_ALL_DATA.eye_r_px_onset;
%     eye_r_py_onset   = SACS_ALL_DATA.eye_r_py_onset;
%     delta_x = visual_px_offset - eye_r_px_onset;
%     delta_y = visual_py_offset - eye_r_py_onset;
%     visual_ang = wrapTo360(atan2d(delta_y, delta_x));
%     visual_ang_bin = discretize(visual_ang, ang_edges);
%     last_bin_id = length(ang_edges) - 1;
%     visual_ang_bin(visual_ang_bin == last_bin_id) = 1; % wrap the circle around
    
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        eye_vx = SACS_ALL_DATA.(['eye_vx' '_' event_type_name]);
        eye_vy = SACS_ALL_DATA.(['eye_vy' '_' event_type_name]);
        eye_vm = sqrt(eye_vx.^2 + eye_vy.^2 );
        eye_vt = ( (eye_vx.*eye_amp_x) + (eye_vy.*eye_amp_y) ) ./ eye_amp_m;
        SACS_ALL_DATA.(['eye_vm' '_' event_type_name]) = eye_vm;
        SACS_ALL_DATA.(['eye_vt' '_' event_type_name]) = eye_vt;
        for counter_tag = 1 : num_tag_bin
            idx_tag = (SACS_ALL_DATA.tag == tag_bin(counter_tag));
            for counter_ang = 1 : num_ang_bin
%                  idx_ang_absol = (visual_ang_bin == counter_ang);
                 idx_ang_absol = (CS_on_data.visual_ang_bin == counter_ang);
                 idx_ang_tuned = (CS_on_data.visual_ang_bin == CS_on_data.idx_CS_tuned(counter_ang));
                for counter_amp = 1 : num_amp_bin
                    idx_amp = (SACS_amp_bin == counter_amp);
                    
                    if flag_build_absol
                    idx_absol = idx_tag & idx_amp & idx_ang_absol;
                    SS_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    CS_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    VM_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vm' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    VT_population_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vt' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    num_sac_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        nansum(idx_absol);
                    end
                    if flag_build_tuned
                    idx_tuned = idx_tag & idx_amp & idx_ang_tuned;
                    SS_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    CS_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    VM_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vm' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    VT_population_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vt' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    num_sac_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell, :) = ...
                        nansum(idx_tuned);
                    end
                end
                for counter_vel = 1 : num_vel_bin
                    idx_vel = (SACS_vel_bin == counter_vel);
                    if flag_build_absol
                    idx_absol = idx_tag & idx_vel & idx_ang_absol;
                    SS_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    CS_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    VM_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vm' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    VT_population_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vt' '_' event_type_name])(:,idx_absol), 2), 1, length_trace);
                    num_sac_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        nansum(idx_absol);
                    end
                    if flag_build_tuned
                    idx_tuned = idx_tag & idx_vel & idx_ang_tuned;
                    SS_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    CS_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    VM_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vm' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    VT_population_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        reshape(nanmean( SACS_ALL_DATA.(['eye_vt' '_' event_type_name])(:,idx_tuned), 2), 1, length_trace);
                    num_sac_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell, :) = ...
                        nansum(idx_tuned);
                    end
                end
            end
        end
    end

end
fprintf('### ALL DONE. ###\n')

%% Save data
fprintf(['Saving .mat files' ' ...'])
if flag_build_absol
save([path_cell_data '..' filesep 'SS_population_absol' '.mat'], 'SS_population_absol', '-v7.3');
save([path_cell_data '..' filesep 'CS_population_absol' '.mat'], 'CS_population_absol', '-v7.3');
save([path_cell_data '..' filesep 'VM_population_absol' '.mat'], 'VM_population_absol', '-v7.3');
save([path_cell_data '..' filesep 'VT_population_absol' '.mat'], 'VT_population_absol', '-v7.3');
save([path_cell_data '..' filesep 'num_sac_absol' '.mat'], 'num_sac_absol', '-v7.3');
end
if flag_build_tuned
save([path_cell_data '..' filesep 'SS_population_tuned' '.mat'], 'SS_population_tuned', '-v7.3');
save([path_cell_data '..' filesep 'CS_population_tuned' '.mat'], 'CS_population_tuned', '-v7.3');
save([path_cell_data '..' filesep 'VM_population_tuned' '.mat'], 'VM_population_tuned', '-v7.3');
save([path_cell_data '..' filesep 'VT_population_tuned' '.mat'], 'VT_population_tuned', '-v7.3');
save([path_cell_data '..' filesep 'num_sac_tuned' '.mat'], 'num_sac_tuned', '-v7.3');
end
fprintf(' --> Completed. \n')

end
%% function build_synchrony_data()
function build_synchrony_data()
clc; close all;
flag_pair_list = true;   % This should be true. DO NOT change it to false
JSP_global_variables(flag_pair_list);
global event_type_list amp_edges vel_edges length_trace ang_values range_cell_with_4dir_behave ang_edges expand_index
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
pCell_ids = build_pCell_ids(flag_pair_list);
num_pCells = size(pCell_ids, 1);
cell_file_name = pCell_ids{1, 1};
load([path_cell_data cell_file_name], 'SACS_ALL_DATA', 'CS_on_data');
tag_bin = unique(SACS_ALL_DATA.tag);
num_tag_bin = length(tag_bin);
num_ang_bin = length(CS_on_data.idx_CS_tuned);
num_amp_bin = length(amp_edges) - 1;
num_vel_bin = length(vel_edges) - 1;
flag_build_absol = true;
flag_build_tuned = false;
if ~xor(flag_build_absol, flag_build_tuned)
    fprintf('><ERROR><: build_synchrony_data, Please select either flag_build_absol or flag_build_tuned. Not both. The code will run much faster this way.\n')
    return;
end

%% Init variables
% STRUCT (SS / CS) -> 10x1 tag struct (amp / vel) -> variable name (onset / vmax / offset / visual) -> 6x8 cell (ampXang / velXang) -> 138x500 double (pCellXtrace)
fprintf(['Initializing synchrony_data' ' ...'])
for counter_tag = 1 : num_tag_bin
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        if flag_build_absol
        SS_synch_joint_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        SS_synch_cross_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        CS_synch_joint_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        CS_synch_cross_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        SS_synch_margn_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        CS_synch_margn_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        num_synch_absol.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        
        SS_synch_joint_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        SS_synch_cross_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        CS_synch_joint_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        CS_synch_cross_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        SS_synch_margn_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        CS_synch_margn_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        num_synch_absol.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        end
        if flag_build_tuned
        SS_synch_joint_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        SS_synch_cross_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        CS_synch_joint_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        CS_synch_cross_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        SS_synch_margn_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        CS_synch_margn_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        num_synch_tuned.amp(counter_tag).(event_type_name) = cell(num_amp_bin, num_ang_bin);
        
        SS_synch_joint_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        SS_synch_cross_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        CS_synch_joint_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        CS_synch_cross_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        SS_synch_margn_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        CS_synch_margn_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        num_synch_tuned.vel(counter_tag).(event_type_name) = cell(num_vel_bin, num_ang_bin);
        end
        for counter_ang = 1 : num_ang_bin
            for counter_amp = 1 : num_amp_bin
                if flag_build_absol
                SS_synch_joint_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                SS_synch_cross_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_joint_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_cross_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                SS_synch_margn_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_margn_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                num_synch_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = zeros(num_pCells, 1);
                end
                if flag_build_tuned
                SS_synch_joint_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                SS_synch_cross_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_joint_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_cross_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                SS_synch_margn_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_margn_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = nan(num_pCells, length_trace);
                num_synch_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang} = zeros(num_pCells, 1);
                end
            end
            for counter_vel = 1 : num_vel_bin
                if flag_build_absol
                SS_synch_joint_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                SS_synch_cross_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_joint_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_cross_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                SS_synch_margn_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_margn_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                num_synch_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = zeros(num_pCells, 1);
                end
                if flag_build_tuned
                SS_synch_joint_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                SS_synch_cross_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_joint_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_cross_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                SS_synch_margn_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                CS_synch_margn_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = nan(num_pCells, length_trace);
                num_synch_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang} = zeros(num_pCells, 1);
                end
            end
        end
    end
end
fprintf(' --> Completed. \n')

%% Loop over pCells
for counter_pCell = 1 : 2 : (num_pCells-1)
    fprintf(['### ' 'Analyzing pCell no. ', num2str((counter_pCell+1)/2), ' / ' num2str(num_pCells/2) ' ###' '\n']);
    %% load SACS_ALL_DATA
    cell_file_name_1 = pCell_ids{counter_pCell  , 1};
    cell_file_name_2 = pCell_ids{counter_pCell+1, 1};
    cell_1 = load([path_cell_data cell_file_name_1], 'SACS_ALL_DATA', 'CS_on_data');
    cell_2 = load([path_cell_data cell_file_name_2], 'SACS_ALL_DATA', 'CS_on_data');
    
    %% Re-calculate CS-on
    CS_count_avg_1  = cell_1.CS_on_data.CS_count( 1, :) + cell_1.CS_on_data.CS_count( 4, :) + cell_1.CS_on_data.CS_count( 8, :);
    CS_count_avg_2  = cell_2.CS_on_data.CS_count( 1, :) + cell_2.CS_on_data.CS_count( 4, :) + cell_2.CS_on_data.CS_count( 8, :);
    sac_count_avg_1 = cell_1.CS_on_data.sac_count(1, :) + cell_1.CS_on_data.sac_count(4, :) + cell_1.CS_on_data.sac_count(8, :);
    sac_count_avg_2 = cell_2.CS_on_data.sac_count(1, :) + cell_2.CS_on_data.sac_count(4, :) + cell_2.CS_on_data.sac_count(8, :);
    CS_count_pair    = CS_count_avg_1  + CS_count_avg_2;
    sac_count_pair   = sac_count_avg_1 + sac_count_avg_2;
    % 'prim_success' tag 1 % 'corr_success' tag 4
    CS_prob_pair = CS_count_pair ./ sac_count_pair;
    
    r_pair = nansum(CS_prob_pair.* exp(1i*deg2rad(ang_values)) , 2); % compute weighted sum of cos and sin of angles
    CS_ang_pair =  wrapTo360(rad2deg(angle(r_pair)));
    
    CS_rho_pair = abs(r_pair) ./ nansum(CS_prob_pair,2);
    
    if (counter_pCell >= range_cell_with_4dir_behave(1)) && (counter_pCell <= range_cell_with_4dir_behave(2))
        idx_CS_on_pair = discretize( ESN_Round(CS_ang_pair, 90.0, 'round'), ang_edges);
    else
        idx_CS_on_pair = discretize(CS_ang_pair, ang_edges);
    end
    
    last_bin_id = length(ang_edges) - 1;
    if idx_CS_on_pair == last_bin_id; idx_CS_on_pair = 1; end
    
    idx_ = idx_CS_on_pair - 1; % make it 0-index format
    if (idx_ == 8); idx_ = 0; end
    idx_CS_pair_tuned = mod((idx_ : 1 : idx_+7), 8) + 1;
    CS_prob_pair_tuned = CS_prob_pair(idx_CS_pair_tuned);
    
    cell_1.CS_on_data.CS_prob_pair = CS_prob_pair;
    cell_1.CS_on_data.CS_prob_pair_tuned = CS_prob_pair_tuned;
    cell_1.CS_on_data.CS_ang_pair  = CS_ang_pair;
    cell_1.CS_on_data.CS_rho_pair  = CS_rho_pair;
    cell_1.CS_on_data.idx_CS_on_pair  = idx_CS_on_pair;
    cell_1.CS_on_data.idx_CS_pair_tuned   = idx_CS_pair_tuned;
    
    cell_2.CS_on_data.CS_prob_pair = CS_prob_pair;
    cell_2.CS_on_data.CS_prob_pair_tuned = CS_prob_pair_tuned;
    cell_2.CS_on_data.CS_ang_pair  = CS_ang_pair;
    cell_2.CS_on_data.CS_rho_pair  = CS_rho_pair;
    cell_2.CS_on_data.idx_CS_on_pair  = idx_CS_on_pair;
    cell_2.CS_on_data.idx_CS_pair_tuned   = idx_CS_pair_tuned;
    
    %% Save CS_on_pair results
    CS_on_data = cell_1.CS_on_data;
    save([path_cell_data cell_file_name_1], 'CS_on_data', '-append');
    CS_on_data = cell_2.CS_on_data;
    save([path_cell_data cell_file_name_2], 'CS_on_data', '-append');
    
    %% Compute data
    SACS_amp_bin_1 = discretize(cell_1.SACS_ALL_DATA.eye_r_amp_m,  amp_edges);
    SACS_vel_bin_1 = discretize(cell_1.SACS_ALL_DATA.eye_r_vm_max, vel_edges);
    SACS_amp_bin_2 = discretize(cell_2.SACS_ALL_DATA.eye_r_amp_m,  amp_edges);
    SACS_vel_bin_2 = discretize(cell_2.SACS_ALL_DATA.eye_r_vm_max, vel_edges);
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        for counter_tag = 1 : num_tag_bin
            idx_tag_1 = (cell_1.SACS_ALL_DATA.tag == tag_bin(counter_tag));
            idx_tag_2 = (cell_2.SACS_ALL_DATA.tag == tag_bin(counter_tag));
            for counter_ang = 1 : num_ang_bin
                %%
                idx_ang_absol_1 = (cell_1.CS_on_data.visual_ang_bin == counter_ang);
                idx_ang_absol_2 = (cell_2.CS_on_data.visual_ang_bin == counter_ang);
                idx_ang_tuned_1 = (cell_1.CS_on_data.visual_ang_bin == cell_1.CS_on_data.idx_CS_pair_tuned(counter_ang));
                idx_ang_tuned_2 = (cell_2.CS_on_data.visual_ang_bin == cell_2.CS_on_data.idx_CS_pair_tuned(counter_ang));
                for counter_amp = 1 : num_amp_bin
                    idx_amp_1 = (SACS_amp_bin_1 == counter_amp);
                    idx_amp_2 = (SACS_amp_bin_2 == counter_amp);               
                    if flag_build_absol
                    idx_absol_1 = idx_tag_1 & idx_amp_1 & idx_ang_absol_1;
                    idx_absol_2 = idx_tag_2 & idx_amp_2 & idx_ang_absol_2;
                    event_SS_1 = cell_1.SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_absol_1);
                    event_CS_1 = cell_1.SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_absol_1);
                    num_synch_1 = nansum(idx_absol_1);
                    event_SS_2 = cell_2.SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_absol_2);
                    event_CS_2 = cell_2.SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_absol_2);
                    num_synch_2 = nansum(idx_absol_2);
%                     event_SS = reshape(nanmean( logical(event_SS_1) & logical(event_SS_2), 2), 1, length_trace);
%                     event_CS = reshape(nanmean( logical(event_CS_1) & logical(event_CS_2), 2), 1, length_trace);                    
%                     SS_synchrony_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = event_SS;
%                     CS_synchrony_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = event_CS;
%                     num_synch_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :)    = num_synch_1;
%                     SS_synchrony_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = event_SS;
%                     CS_synchrony_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = event_CS;
%                     num_synch_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :)    = num_synch_2;
                    event_SS_1 = expand_index_event_data(event_SS_1, 1);
                    event_SS_2 = expand_index_event_data(event_SS_2, 1);
                    event_CS_1 = expand_index_event_data(event_CS_1, 1);
                    event_CS_2 = expand_index_event_data(event_CS_2, 1);
                    event_SS_p1p2 = reshape(nanmean( ( logical(event_SS_1)) & ( logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_n1n2 = reshape(nanmean( (~logical(event_SS_1)) & (~logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_p1n2 = reshape(nanmean( ( logical(event_SS_1)) & (~logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_n1p2 = reshape(nanmean( (~logical(event_SS_1)) & ( logical(event_SS_2)), 2), 1, length_trace);
                    event_CS_p1p2 = reshape(nanmean( ( logical(event_CS_1)) & ( logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_n1n2 = reshape(nanmean( (~logical(event_CS_1)) & (~logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_p1n2 = reshape(nanmean( ( logical(event_CS_1)) & (~logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_n1p2 = reshape(nanmean( (~logical(event_CS_1)) & ( logical(event_CS_2)), 2), 1, length_trace);
                    %%
                    SS_synch_joint_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = event_SS_p1p2;
                    SS_synch_cross_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = event_SS_p1n2;
                    CS_synch_joint_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = event_CS_p1p2;
                    CS_synch_cross_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = event_CS_p1n2;
                    num_synch_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :)    = num_synch_1;
                    SS_synch_joint_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = event_SS_n1n2;
                    SS_synch_cross_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = event_SS_n1p2;
                    CS_synch_joint_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = event_CS_n1n2;
                    CS_synch_cross_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = event_CS_n1p2;
                    num_synch_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :)    = num_synch_2;
                    SS_synch_margn_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = ...
                        reshape(nanmean( event_SS_1, 2), 1, length_trace);
                    CS_synch_margn_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = ...
                        reshape(nanmean( event_CS_1, 2), 1, length_trace);
                    SS_synch_margn_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = ...
                        reshape(nanmean( event_SS_2, 2), 1, length_trace);
                    CS_synch_margn_absol.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = ...
                        reshape(nanmean( event_CS_2, 2), 1, length_trace);                   
                    end
                    if flag_build_tuned
                    idx_tuned_1 = idx_tag_1 & idx_amp_1 & idx_ang_tuned_1;
                    idx_tuned_2 = idx_tag_2 & idx_amp_2 & idx_ang_tuned_2;
                    event_SS_1 = cell_1.SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_tuned_1);
                    event_CS_1 = cell_1.SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_tuned_1);
                    num_synch_1 = nansum(idx_tuned_1);
                    event_SS_2 = cell_2.SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_tuned_2);
                    event_CS_2 = cell_2.SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_tuned_2);
                    num_synch_2 = nansum(idx_tuned_2);
                    event_SS_1 = expand_index_event_data(event_SS_1, 1);  % dim=1, expand along column. event_ is a 500xn matrix
                    event_SS_2 = expand_index_event_data(event_SS_2, 1);
                    event_CS_1 = expand_index_event_data(event_CS_1, 1);
                    event_CS_2 = expand_index_event_data(event_CS_2, 1);
                    event_SS_p1p2 = reshape(nanmean( ( logical(event_SS_1)) & ( logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_n1n2 = reshape(nanmean( (~logical(event_SS_1)) & (~logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_p1n2 = reshape(nanmean( ( logical(event_SS_1)) & (~logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_n1p2 = reshape(nanmean( (~logical(event_SS_1)) & ( logical(event_SS_2)), 2), 1, length_trace);
                    event_CS_p1p2 = reshape(nanmean( ( logical(event_CS_1)) & ( logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_n1n2 = reshape(nanmean( (~logical(event_CS_1)) & (~logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_p1n2 = reshape(nanmean( ( logical(event_CS_1)) & (~logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_n1p2 = reshape(nanmean( (~logical(event_CS_1)) & ( logical(event_CS_2)), 2), 1, length_trace);
                    SS_synch_joint_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = event_SS_p1p2;
                    SS_synch_cross_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = event_SS_p1n2;
                    CS_synch_joint_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = event_CS_p1p2;
                    CS_synch_cross_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = event_CS_p1n2;
                    num_synch_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :)    = num_synch_1;
                    SS_synch_joint_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = event_SS_n1n2;
                    SS_synch_cross_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = event_SS_n1p2;
                    CS_synch_joint_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = event_CS_n1n2;
                    CS_synch_cross_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = event_CS_n1p2;
                    num_synch_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :)    = num_synch_2;
                    SS_synch_margn_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = ...
                        reshape(nanmean( event_SS_1, 2), 1, length_trace);
                    CS_synch_margn_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell,   :) = ...
                        reshape(nanmean( event_CS_1, 2), 1, length_trace);
                    SS_synch_margn_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = ...
                        reshape(nanmean( event_SS_2, 2), 1, length_trace);
                    CS_synch_margn_tuned.amp(counter_tag).(event_type_name){counter_amp, counter_ang}(counter_pCell+1, :) = ...
                        reshape(nanmean( event_CS_2, 2), 1, length_trace);
                    end
                end
                for counter_vel = 1 : num_vel_bin
                    idx_vel_1 = (SACS_vel_bin_1 == counter_vel);
                    idx_vel_2 = (SACS_vel_bin_2 == counter_vel);
                    if flag_build_absol
                    idx_absol_1 = idx_tag_1 & idx_vel_1 & idx_ang_absol_1;
                    idx_absol_2 = idx_tag_2 & idx_vel_2 & idx_ang_absol_2;
                    event_SS_1 = cell_1.SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_absol_1);
                    event_CS_1 = cell_1.SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_absol_1);
                    num_synch_1 = nansum(idx_absol_1);
                    event_SS_2 = cell_2.SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_absol_2);
                    event_CS_2 = cell_2.SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_absol_2);
                    num_synch_2 = nansum(idx_absol_2);
                    event_SS_1 = expand_index_event_data(event_SS_1, 1);
                    event_SS_2 = expand_index_event_data(event_SS_2, 1);
                    event_CS_1 = expand_index_event_data(event_CS_1, 1);
                    event_CS_2 = expand_index_event_data(event_CS_2, 1);
                    event_SS_p1p2 = reshape(nanmean( ( logical(event_SS_1)) & ( logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_n1n2 = reshape(nanmean( (~logical(event_SS_1)) & (~logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_p1n2 = reshape(nanmean( ( logical(event_SS_1)) & (~logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_n1p2 = reshape(nanmean( (~logical(event_SS_1)) & ( logical(event_SS_2)), 2), 1, length_trace);
                    event_CS_p1p2 = reshape(nanmean( ( logical(event_CS_1)) & ( logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_n1n2 = reshape(nanmean( (~logical(event_CS_1)) & (~logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_p1n2 = reshape(nanmean( ( logical(event_CS_1)) & (~logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_n1p2 = reshape(nanmean( (~logical(event_CS_1)) & ( logical(event_CS_2)), 2), 1, length_trace);
                    SS_synch_joint_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = event_SS_p1p2;
                    SS_synch_cross_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = event_SS_p1n2;
                    CS_synch_joint_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = event_CS_p1p2;
                    CS_synch_cross_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = event_CS_p1n2;
                    num_synch_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :)    = num_synch_1;
                    SS_synch_joint_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = event_SS_n1n2;
                    SS_synch_cross_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = event_SS_n1p2;
                    CS_synch_joint_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = event_CS_n1n2;
                    CS_synch_cross_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = event_CS_n1p2;
                    num_synch_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :)    = num_synch_2;
                    SS_synch_margn_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = ...
                        reshape(nanmean( event_SS_1, 2), 1, length_trace);
                    CS_synch_margn_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = ...
                        reshape(nanmean( event_CS_1, 2), 1, length_trace);
                    SS_synch_margn_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = ...
                        reshape(nanmean( event_SS_2, 2), 1, length_trace);
                    CS_synch_margn_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = ...
                        reshape(nanmean( event_CS_2, 2), 1, length_trace);
                    
%                     SS_synchrony_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = event_SS;
%                     CS_synchrony_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = event_CS;
%                     num_synch_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :)    = num_synch_1;
%                     SS_synchrony_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = event_SS;
%                     CS_synchrony_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = event_CS;
%                     num_synch_absol.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :)    = num_synch_2;
                    end
                    if flag_build_tuned
                    idx_tuned_1 = idx_tag_1 & idx_vel_1 & idx_ang_tuned_1;
                    idx_tuned_2 = idx_tag_2 & idx_vel_2 & idx_ang_tuned_2;
                    event_SS_1 = cell_1.SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_tuned_1);
                    event_CS_1 = cell_1.SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_tuned_1);
                    num_synch_1 = nansum(idx_tuned_1);
                    event_SS_2 = cell_2.SACS_ALL_DATA.(['neuro_SS' '_' event_type_name])(:,idx_tuned_2);
                    event_CS_2 = cell_2.SACS_ALL_DATA.(['neuro_CS' '_' event_type_name])(:,idx_tuned_2);
                    num_synch_2 = nansum(idx_tuned_2);
                    event_SS_1 = expand_index_event_data(event_SS_1, 1);
                    event_SS_2 = expand_index_event_data(event_SS_2, 1);
                    event_CS_1 = expand_index_event_data(event_CS_1, 1);
                    event_CS_2 = expand_index_event_data(event_CS_2, 1);
                    event_SS_p1p2 = reshape(nanmean( ( logical(event_SS_1)) & ( logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_n1n2 = reshape(nanmean( (~logical(event_SS_1)) & (~logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_p1n2 = reshape(nanmean( ( logical(event_SS_1)) & (~logical(event_SS_2)), 2), 1, length_trace);
                    event_SS_n1p2 = reshape(nanmean( (~logical(event_SS_1)) & ( logical(event_SS_2)), 2), 1, length_trace);
                    event_CS_p1p2 = reshape(nanmean( ( logical(event_CS_1)) & ( logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_n1n2 = reshape(nanmean( (~logical(event_CS_1)) & (~logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_p1n2 = reshape(nanmean( ( logical(event_CS_1)) & (~logical(event_CS_2)), 2), 1, length_trace);
                    event_CS_n1p2 = reshape(nanmean( (~logical(event_CS_1)) & ( logical(event_CS_2)), 2), 1, length_trace);
                    SS_synch_joint_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = event_SS_p1p2;
                    SS_synch_cross_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = event_SS_p1n2;
                    CS_synch_joint_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = event_CS_p1p2;
                    CS_synch_cross_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = event_CS_p1n2;
                    num_synch_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :)    = num_synch_1;
                    SS_synch_joint_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = event_SS_n1n2;
                    SS_synch_cross_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = event_SS_n1p2;
                    CS_synch_joint_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = event_CS_n1n2;
                    CS_synch_cross_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = event_CS_n1p2;
                    num_synch_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :)    = num_synch_2;
                    SS_synch_margn_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = ...
                        reshape(nanmean( event_SS_1, 2), 1, length_trace);
                    CS_synch_margn_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell,   :) = ...
                        reshape(nanmean( event_CS_1, 2), 1, length_trace);
                    SS_synch_margn_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = ...
                        reshape(nanmean( event_SS_2, 2), 1, length_trace);
                    CS_synch_margn_tuned.vel(counter_tag).(event_type_name){counter_vel, counter_ang}(counter_pCell+1, :) = ...
                        reshape(nanmean( event_CS_2, 2), 1, length_trace);
                    end
                end
            end
        end
    end

end
fprintf('### ALL DONE. ###\n')

%% Save data
fprintf(['Saving .mat files' ' ...'])
if flag_build_absol
win_len=[num2str((expand_index*2) + 1) 'ms'];
save([path_cell_data '..' filesep 'SS_synch_joint_' win_len '_absol' '.mat'], 'SS_synch_joint_absol', '-v7.3');
save([path_cell_data '..' filesep 'SS_synch_cross_' win_len '_absol' '.mat'], 'SS_synch_cross_absol', '-v7.3');
save([path_cell_data '..' filesep 'CS_synch_joint_' win_len '_absol' '.mat'], 'CS_synch_joint_absol', '-v7.3');
save([path_cell_data '..' filesep 'CS_synch_cross_' win_len '_absol' '.mat'], 'CS_synch_cross_absol', '-v7.3');
save([path_cell_data '..' filesep 'SS_synch_margn_' win_len '_absol' '.mat'], 'SS_synch_margn_absol', '-v7.3');
save([path_cell_data '..' filesep 'CS_synch_margn_' win_len '_absol' '.mat'], 'CS_synch_margn_absol', '-v7.3');
save([path_cell_data '..' filesep 'num_synch_' win_len '_absol' '.mat'], 'num_synch_absol', '-v7.3');
end
if flag_build_tuned
win_len=[num2str((expand_index*2) + 1) 'ms'];
save([path_cell_data '..' filesep 'SS_synch_joint_' win_len '_tuned' '.mat'], 'SS_synch_joint_tuned', '-v7.3');
save([path_cell_data '..' filesep 'SS_synch_cross_' win_len '_tuned' '.mat'], 'SS_synch_cross_tuned', '-v7.3');
save([path_cell_data '..' filesep 'CS_synch_joint_' win_len '_tuned' '.mat'], 'CS_synch_joint_tuned', '-v7.3');
save([path_cell_data '..' filesep 'CS_synch_cross_' win_len '_tuned' '.mat'], 'CS_synch_cross_tuned', '-v7.3');
save([path_cell_data '..' filesep 'SS_synch_margn_' win_len '_tuned' '.mat'], 'SS_synch_margn_tuned', '-v7.3');
save([path_cell_data '..' filesep 'CS_synch_margn_' win_len '_tuned' '.mat'], 'CS_synch_margn_tuned', '-v7.3');
save([path_cell_data '..' filesep 'num_synch_' win_len '_tuned' '.mat'], 'num_synch_tuned', '-v7.3');
end
fprintf(' --> Completed. \n')

end

%% AVG FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function population_data_avg_over_levels
function [population_avg_levels, num_sac_avg_levels] = population_data_avg_over_levels(population_data, num_sac_data, variable, dim, levels)
%% Handle inputs
% variable = 'amp' or 'vel'
% dim = 1 -> average over different varibale amp/vel levels
% dim = 2 -> average over different ang levels
if nargin < 5
    levels = [];
end
if nargin < 4
    dim = 1;
end
if nargin < 3
    variable = 'amp';
end
if nargin < 2
    population_avg_levels = [];
    num_sac_avg_levels = [];
    return;
end

%% Calc necessary variables
JSP_global_variables();
global event_type_list inds_span
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
length_trace = length(inds_span);
num_pCells  = size(population_data.(variable)(1).onset{1, 1}, 1);
num_tag_bin = length(population_data.(variable));
num_ang_bin = size(population_data.(variable)(1).onset, 2);
num_var_bin = size(population_data.(variable)(1).onset, 1);
vec_tag_bin = 1:1:num_tag_bin;
vec_ang_bin = 1:1:num_ang_bin;
vec_var_bin = 1:1:num_var_bin;

%% Init population_avg_levels
if dim == 1
    if isempty(levels)
        levels = vec_var_bin;
    end
    num_var_bin     = length(levels);
    num_ang_bin_avg = num_ang_bin;
    num_var_bin_avg = 1;
    vec_var_bin     = levels;
    vec_ang_bin_avg = vec_ang_bin;
    vec_var_bin_avg = 1;
elseif dim == 2
    if isempty(levels)
        levels = vec_ang_bin;
    end
    num_ang_bin     = length(levels);
    num_ang_bin_avg = 1;
    num_var_bin_avg = num_var_bin;
    vec_ang_bin     = levels;
    vec_ang_bin_avg = 1;
    vec_var_bin_avg = vec_var_bin;
else
    if isempty(levels)
        levels = vec_ang_bin;
    end
    dim = 1;
    num_var_bin     = length(levels);
    num_ang_bin_avg = num_ang_bin;
    num_var_bin_avg = 1;
    vec_var_bin     = levels;
    vec_ang_bin_avg = vec_ang_bin;
    vec_var_bin_avg = 1;
end
population_avg_levels = struct;
num_sac_avg_levels    = struct;
for counter_tag = vec_tag_bin
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        population_avg_levels.(variable)(counter_tag).(event_type_name) = cell(num_var_bin_avg, num_ang_bin_avg);
        num_sac_avg_levels.(   variable)(counter_tag).(event_type_name) = cell(num_var_bin_avg, num_ang_bin_avg);
        for counter_ang = vec_ang_bin_avg
            for counter_var = vec_var_bin_avg
                population_avg_levels.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = nan(num_pCells, length_trace);
                num_sac_avg_levels.(   variable)(counter_tag).(event_type_name){counter_var, counter_ang} = zeros(num_pCells, 1);
            end
        end
    end
end

%% Compute population_avg_levels
fprintf(['population_data_avg_over_levels' ' ...'])
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_tag = vec_tag_bin
        if dim == 1
            for counter_ang = vec_ang_bin
                event_data_avg_ = zeros(num_pCells, length_trace);
                num_sac_avg_    = zeros(num_pCells, length_trace);
                for counter_var = vec_var_bin
                    event_data_ = ...
                        population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                    event_data_(isnan(event_data_)) = 0;
                    num_sac_ = ...
                           num_sac_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                    num_sac_ = repmat(num_sac_, 1, length_trace);
                    event_data_avg_ = event_data_avg_ + (event_data_ .* num_sac_);
                    num_sac_avg_    = num_sac_avg_    + num_sac_;
                end
                population_avg_levels.(variable)(counter_tag).(event_type_name){1, counter_ang} = event_data_avg_ ./ num_sac_avg_;
                num_sac_avg_levels.(variable)(counter_tag).(event_type_name){1, counter_ang} = num_sac_avg_(:,1);
            end
        elseif dim == 2
            for counter_var = vec_var_bin
                event_data_avg_ = zeros(num_pCells, length_trace);
                num_sac_avg_    = zeros(num_pCells, length_trace);
                for counter_ang = vec_ang_bin
                    event_data_ = ...
                        population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                    event_data_(isnan(event_data_)) = 0;
                    num_sac_ = ...
                           num_sac_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                    num_sac_ = repmat(num_sac_, 1, length_trace);
                    event_data_avg_ = event_data_avg_ + (event_data_ .* num_sac_);
                    num_sac_avg_    = num_sac_avg_    + num_sac_;
                end
                population_avg_levels.(variable)(counter_tag).(event_type_name){counter_var, 1} = event_data_avg_ ./ num_sac_avg_;
                num_sac_avg_levels.(   variable)(counter_tag).(event_type_name){counter_var, 1} = num_sac_avg_(:,1);
            end
        end
    end
end
fprintf(' --> Completed. \n')
end
%% function population_data_combine_levels
function [population_combined_levels, num_sac_combined_levels] = population_data_combine_levels(population_data, num_sac_data, variable, dim, idx_levels)
% The outputs of the function has the same structure and size as the inputs.
% The code will loop over idx_levels and then store the resulted avg data of those levels in the idx_levels(1).
% The data for other idx_levels, other than idx_levels(1), will be set to NaN.
%% Handle inputs
% variable = 'amp' or 'vel'
% dim = 1 -> average over different varibale amp/vel levels
% dim = 2 -> average over different ang levels
% idx_levels should be an array of integers
if nargin < 5
    population_combined_levels = population_data;
    num_sac_combined_levels = num_sac_data;
    return;
end
if nargin < 4
    dim = 1;
end
if nargin < 3
    variable = 'amp';
end
if nargin < 2
    population_combined_levels = [];
    num_sac_combined_levels = [];
    return;
end

%% Calc necessary variables
JSP_global_variables();
global event_type_list inds_span
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
length_trace = length(inds_span);
num_pCells  = size(population_data.(variable)(1).onset{1, 1}, 1);
num_tag_bin = length(population_data.(variable));
num_ang_bin = size(population_data.(variable)(1).onset, 2);
num_var_bin = size(population_data.(variable)(1).onset, 1);

%% Init population_avg_levels
population_combined_levels = population_data;
num_sac_combined_levels = num_sac_data;

%% Compute population_avg_levels
fprintf(['population_data_combine_levels' ' ...'])
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_tag = 1 : num_tag_bin
        if dim == 1
            for counter_ang = 1 : num_ang_bin
                event_data_avg_ = zeros(num_pCells, length_trace);
                num_sac_avg_    = zeros(num_pCells, 1);
                for counter_var = idx_levels
                    event_data_ = ...
                        population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                    event_data_(isnan(event_data_)) = 0;
                    num_sac_ = ...
                           num_sac_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                    event_data_avg_ = event_data_avg_ + (event_data_ .* repmat(num_sac_, 1, length_trace) );
                    num_sac_avg_    = num_sac_avg_    + num_sac_;
                    population_combined_levels.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = nan(size(event_data_));
                    num_sac_combined_levels.(   variable)(counter_tag).(event_type_name){counter_var, counter_ang} = zeros(size(num_sac_));
                end
                population_combined_levels.(variable)(counter_tag).(event_type_name){idx_levels(1), counter_ang} = event_data_avg_ ./ num_sac_avg_;
                num_sac_combined_levels.(   variable)(counter_tag).(event_type_name){idx_levels(1), counter_ang} = num_sac_avg_(:,1);
            end
        elseif dim == 2
            for counter_var = 1 : num_var_bin
                event_data_avg_ = zeros(num_pCells, length_trace);
                num_sac_avg_    = zeros(num_pCells, 1);
                for counter_ang = idx_levels
                    event_data_ = ...
                        population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                    event_data_(isnan(event_data_)) = 0;
                    num_sac_ = ...
                           num_sac_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                    event_data_avg_ = event_data_avg_ + (event_data_ .* repmat(num_sac_, 1, length_trace) );
                    num_sac_avg_    = num_sac_avg_    + num_sac_;
                    population_combined_levels.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = nan(size(event_data_));
                    num_sac_combined_levels.(   variable)(counter_tag).(event_type_name){counter_var, counter_ang} = zeros(size(num_sac_));
                end
                population_combined_levels.(variable)(counter_tag).(event_type_name){counter_var, idx_levels(1)} = event_data_avg_ ./ num_sac_avg_;
                num_sac_combined_levels.(   variable)(counter_tag).(event_type_name){counter_var, idx_levels(1)} = num_sac_avg_(:,1);
            end
        end
    end
end
fprintf(' --> Completed. \n')
end
%% function population_data_combine_tags
function [population_combined_tags, num_sac_combined_tags] = population_data_combine_tags(population_data, num_sac_data, variable, idx_tags)
% The outputs of the function has the same structure and size as the inputs.
% The code will loop over idx_tags and then store the resulted avg data of those tags in the idx_tags(1).
% The data for other tags, other than idx_tags(1), will be set to NaN.
%% Handle inputs
% variable = 'amp' or 'vel'
% idx_tags should be an array of integers
if nargin < 4
    population_combined_tags = population_data;
    num_sac_combined_tags    = num_sac_data;
    return;
end
if nargin < 3
    variable = 'amp';
end
if nargin < 2
    population_combined_tags = [];
    num_sac_combined_tags = [];
    return;
end

%% Calc necessary variables
JSP_global_variables();
global event_type_list inds_span
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
length_trace = length(inds_span);
num_pCells  = size(population_data.(variable)(1).onset{1, 1}, 1);
num_tag_bin = length(population_data.(variable));
num_ang_bin = size(population_data.(variable)(1).onset, 2);
num_var_bin = size(population_data.(variable)(1).onset, 1);

%% Init population_avg_levels
population_combined_tags = population_data;
num_sac_combined_tags    = num_sac_data;

%% Compute population_avg_levels
fprintf(['population_data_combine_tags' ' ...'])
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_ang = 1 : num_ang_bin
        for counter_var = 1 : num_var_bin
            event_data_avg_ = zeros(num_pCells, length_trace);
            num_sac_avg_    = zeros(num_pCells, 1);
            for counter_tag = idx_tags
                event_data_ = ...
                    population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                event_data_(isnan(event_data_)) = 0;
                num_sac_ = ...
                    num_sac_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                event_data_avg_ = event_data_avg_ + (event_data_ .* repmat(num_sac_, 1, length_trace) );
                num_sac_avg_    = num_sac_avg_    + num_sac_;
                population_combined_tags.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = nan(size(event_data_));
                num_sac_combined_tags.(   variable)(counter_tag).(event_type_name){counter_var, counter_ang} = zeros(size(num_sac_));
            end
            population_combined_tags.(variable)(idx_tags(1)).(event_type_name){counter_var, counter_ang} = event_data_avg_ ./ num_sac_avg_;
            num_sac_combined_tags.(   variable)(idx_tags(1)).(event_type_name){counter_var, counter_ang} = num_sac_avg_(:,1);
        end
    end
end
fprintf(' --> Completed. \n')
end
%% function population_data_avg_over_pCells
function [population_avg_pCells, population_sem_pCells] = population_data_avg_over_pCells(population_data, num_sac_data, variable)
%% Handle inputs
if nargin < 3
    variable = 'amp';
end
if nargin < 2
    population_avg_pCells = [];
    population_sem_pCells = [];
    return;
end

%% Calc necessary variables
JSP_global_variables();
global event_type_list inds_span
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
length_trace = length(inds_span);
num_pCells  = size(population_data.(variable)(1).onset{1, 1}, 1);
num_tag_bin = length(population_data.(variable));
num_ang_bin = size(population_data.(variable)(1).onset, 2);
num_var_bin = size(population_data.(variable)(1).onset, 1);

%% Init population_avg_pCells population_sem_pCells
population_avg_pCells = struct;
population_sem_pCells = struct;
for counter_tag = 1 : num_tag_bin
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        population_avg_pCells.(variable)(counter_tag).(event_type_name) = cell(num_var_bin, num_ang_bin);
        population_sem_pCells.(variable)(counter_tag).(event_type_name) = cell(num_var_bin, num_ang_bin);
        for counter_ang = 1 : num_ang_bin
            for counter_var = 1 : num_var_bin
                population_avg_pCells.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = nan(1, length_trace);
                population_sem_pCells.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = nan(1, length_trace);
            end
        end
    end
end

%% Compute population_avg_pCells
fprintf(['population_data_avg_over_pCells' ' ...'])
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_tag = 1 : num_tag_bin
        for counter_ang = 1 : num_ang_bin
            for counter_var = 1 : num_var_bin
                event_data_ = ...
                    population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                event_data_(isnan(event_data_)) = 0;
                num_sac_ = ...
                       num_sac_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                num_sac_norm1_ = num_sac_ ./ nansum(num_sac_, 1);
                avg_pCells_ = nansum(event_data_ .* repmat(num_sac_norm1_, 1, length_trace), 1);
                sem_pCells_ = sqrt( var(event_data_,num_sac_norm1_, 1, 'omitnan') ) ./ sqrt(sum(num_sac_>0));
                population_avg_pCells.(variable)(counter_tag).(event_type_name){counter_var, counter_ang}(1,:) = avg_pCells_;
                population_sem_pCells.(variable)(counter_tag).(event_type_name){counter_var, counter_ang}(1,:) = sem_pCells_;
            end
        end
    end
end
fprintf(' --> Completed. \n')

end
%% function population_data_perm_over_pCells
function [population_avg_pCells, population_sem_pCells] = population_data_perm_over_pCells(population_data, num_sac_data, variable)
%% Handle inputs
if nargin < 3
    variable = 'amp';
end
if nargin < 2
    population_avg_pCells = [];
    population_sem_pCells = [];
    return;
end

%% Calc necessary variables
JSP_global_variables();
global event_type_list inds_span
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
length_trace = length(inds_span);
num_pCells  = size(population_data.(variable)(1).onset{1, 1}, 1);
num_tag_bin = length(population_data.(variable));
num_ang_bin = size(population_data.(variable)(1).onset, 2);
num_var_bin = size(population_data.(variable)(1).onset, 1);
num_iterations = 1000;

%% Init population_std_sacs
population_avg_pCells = struct;
population_sem_pCells = struct;
for counter_tag = 1 : num_tag_bin
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        population_avg_pCells.(variable)(counter_tag).(event_type_name) = cell(num_var_bin, num_ang_bin);
        population_sem_pCells.(variable)(counter_tag).(event_type_name) = cell(num_var_bin, num_ang_bin);
        for counter_ang = 1 : num_ang_bin
            for counter_var = 1 : num_var_bin
                population_avg_pCells.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = nan(1, length_trace);
                population_sem_pCells.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = nan(1, length_trace);
            end
        end
    end
end

%% Compute population_avg_pCells
fprintf(['population_data_perm_over_pCells' ' ...'])
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_tag = 1 : num_tag_bin
        for counter_ang = 1 : num_ang_bin
            for counter_var = 1 : num_var_bin
                event_data_ = ...
                    population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                event_data_(isnan(event_data_)) = 0;
                num_sac_ = ...
                       num_sac_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                num_sac_matrix_ = repmat(num_sac_, 1, length_trace);
                   
                % SEM based on bootstrapping on trials
                %{
                % DISCONTINUED, DO NOT USE
                num_sac_total_ = nansum(num_sac_);
                idx_edges_ = [1; num_sac_];
                idx_edges_ = cumsum(idx_edges_);
                count_data_ = event_data_ .* num_sac_matrix_;
                event_data_perm_ = nan(num_iterations, length_trace);
                for counter_iteration = 1 : num_iterations
                    idx_iteration_ = randi(num_sac_total_, 1, num_sac_total_);
                    [N_idx_edges_,~] = histcounts(idx_iteration_,idx_edges_);
                    weight_cell_ = N_idx_edges_' ./ num_sac_;
                    weight_cell_matrix_ = repmat(weight_cell_, 1, length_trace);
                    event_data_iteration_ = count_data_ .* weight_cell_matrix_;
                    event_data_iteration_ = nansum(event_data_iteration_) ./ num_sac_total_;
                    event_data_perm_(counter_iteration, :) = event_data_iteration_;
                end
                %}
                
                % SEM based on bootstrapping on pCells
                event_data_perm_ = nan(num_iterations, length_trace);
                for counter_iteration = 1 : num_iterations
                    idx_iteration_ = randi(num_pCells, 1, num_pCells);
                    event_data_iteration_ = event_data_(idx_iteration_, :);
                    num_sac_iteration_    = num_sac_matrix_(   idx_iteration_, :);
                    avg_data_iteration_ = nansum(event_data_iteration_ .* num_sac_iteration_) ./ nansum(num_sac_iteration_);
                    event_data_perm_(counter_iteration, :) = avg_data_iteration_;
                end
                event_data_perm_avg = nanmean(event_data_perm_);
                event_data_perm_sem = nanstd(event_data_perm_) ./ sqrt(sum(num_sac_>0));
                population_avg_pCells.(variable)(counter_tag).(event_type_name){counter_var, counter_ang}(1,:) = event_data_perm_avg;
                population_sem_pCells.(variable)(counter_tag).(event_type_name){counter_var, counter_ang}(1,:) = event_data_perm_sem;
            end
        end
    end
end
fprintf(' --> Completed. \n')

end
%% function population_data_idx_pCells
function [population_idx_pCells, num_sac_idx_pCells] = population_data_idx_pCells(population_data, num_sac_data, variable, idx_pCells)
% idx_pCells: is a boolean array. 1 for including a pCell, and 0 for exluding a pCell
%% Handle inputs
if nargin < 3
    variable = 'amp';
end
if nargin < 2
    population_idx_pCells = [];
    num_sac_idx_pCells = [];
    return;
end

%% Calc necessary variables
JSP_global_variables();
global event_type_list inds_span
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
length_trace = length(inds_span);
num_pCells  = size(population_data.(variable)(1).onset{1, 1}, 1);
num_idx_pCells = sum(idx_pCells);
num_tag_bin = length(population_data.(variable));
num_ang_bin = size(population_data.(variable)(1).onset, 2);
num_var_bin = size(population_data.(variable)(1).onset, 1);

%% Init population_idx_pCells
population_idx_pCells = struct;
num_sac_idx_pCells    = struct;
for counter_tag = 1 : num_tag_bin
    for counter_event_type = 1 : length(event_type_list)
        event_type_name = event_type_list{counter_event_type};
        population_idx_pCells.(variable)(counter_tag).(event_type_name) = cell(num_var_bin, num_ang_bin);
        num_sac_idx_pCells.(variable)(counter_tag).(event_type_name) = cell(num_var_bin, num_ang_bin);
        for counter_ang = 1 : num_ang_bin
            for counter_var = 1 : num_var_bin
                population_idx_pCells.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = nan(num_idx_pCells, length_trace);
                num_sac_idx_pCells.(   variable)(counter_tag).(event_type_name){counter_var, counter_ang} = zeros(num_idx_pCells, 1);
            end
        end
    end
end

%% Compute population_avg_pCells
fprintf(['population_data_idx_pCells' ' ...'])
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_tag = 1 : num_tag_bin
        for counter_ang = 1 : num_ang_bin
            for counter_var = 1 : num_var_bin
                event_data_ = ...
                    population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                num_sac_ = ...
                       num_sac_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                event_data_idx_pCells_ = event_data_(idx_pCells,:);
                num_sac_idx_pCells_ = num_sac_(idx_pCells,:);
                population_idx_pCells.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = event_data_idx_pCells_;
                num_sac_idx_pCells.(   variable)(counter_tag).(event_type_name){counter_var, counter_ang} = num_sac_idx_pCells_;
            end
        end
    end
end
fprintf(' --> Completed. \n')
end
%% function population_data_smooth_pCells
function population_smooth = population_data_smooth_pCells(population_data, num_sac_data, variable)
%% Handle inputs
if nargin < 3
    variable = 'amp';
end
if nargin < 2
    population_smooth = [];
    return;
end

%% Calc necessary variables
JSP_global_variables();
global event_type_list inds_span
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
length_trace = length(inds_span);
num_pCells  = size(population_data.(variable)(1).onset{1, 1}, 1);
num_tag_bin = length(population_data.(variable));
num_ang_bin = size(population_data.(variable)(1).onset, 2);
num_var_bin = size(population_data.(variable)(1).onset, 1);

%% Init population_avg_pCells population_sem_pCells
population_smooth = population_data;

%% Compute population_avg_pCells
fprintf(['population_data_smooth_pCells' ' ...'])
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_tag = 1 : num_tag_bin
        for counter_ang = 1 : num_ang_bin
            for counter_var = 1 : num_var_bin
                event_data_ = ...
                    population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
                event_data_(isnan(event_data_)) = 0;
                event_data_smooth = ESN_smooth(event_data_, 2);
                population_smooth.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = event_data_smooth;
            end
        end
    end
end
fprintf(' --> Completed. \n')
end
%% function population_data_subtract_baseline
function population_data_baseline = population_data_subtract_baseline(population_data, firing_rate, variable)
%% Handle inputs

% variable = 'amp' or 'vel'
% dim = 1 -> average over different varibale amp/vel levels
% dim = 2 -> average over different ang levels
if nargin < 3
    variable = 'amp';
end
if nargin < 2
    population_data_baseline = [];
    return;
end

%% Calc necessary variables
JSP_global_variables();
global event_type_list inds_span
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
length_trace = length(inds_span);
num_pCells  = size(population_data.(variable)(1).onset{1, 1}, 1);
num_tag_bin = length(population_data.(variable));
num_ang_bin = size(population_data.(variable)(1).onset, 2);
num_var_bin = size(population_data.(variable)(1).onset, 1);

%% Compute population_data_baseline
fprintf(['population_data_subtract_baseline' ' ...'])
population_data_baseline = struct;
population_data_baseline.(variable) = population_data.(variable);
firing_rate = repmat(firing_rate(:), 1, length_trace) ./ 1000.0; % to convert Hz back to probability we should devide by 1000.
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_tag = 1 : num_tag_bin
        for counter_ang = 1 : num_ang_bin
            for counter_var = 1 : num_var_bin
                event_data_ = ...
                    population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
%                 event_data_baseline_ = event_data_ - firing_rate;
                  event_data_baseline_ = event_data_ - 0;

                population_data_baseline.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = event_data_baseline_;
            end
        end
    end
end
fprintf(' --> Completed. \n')

end
%% function synchrony_data_ratio
function synch_ratio = synchrony_data_ratio(synch_joint, synch_margn, variable)
%% Handle inputs
if nargin < 3
    variable = 'amp';
end
if nargin < 2
    synch_ratio = [];
    return;
end

%% Calc necessary variables
JSP_global_variables();
global event_type_list inds_span
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
length_trace = length(inds_span);
num_pCells  = size(synch_joint.(variable)(1).onset{1, 1}, 1);
num_tag_bin = length(synch_joint.(variable));
num_ang_bin = size(synch_joint.(variable)(1).onset, 2);
num_var_bin = size(synch_joint.(variable)(1).onset, 1);

%% Init population_avg_pCells population_sem_pCells
synch_ratio = synch_joint;

%% Compute population_avg_pCells
fprintf(['synchrony_data_ratio' ' ...'])
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_tag = 1 : num_tag_bin
        for counter_ang = 1 : num_ang_bin
            for counter_var = 1 : num_var_bin
                event_data_joint_ = ...
                    synch_joint.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
%                 event_data_joint_(isnan(event_data_joint_)) = 0;
                event_data_margn_ = ...
                    synch_margn.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
%                 event_data_margn_(isnan(event_data_margn_)) = 0;
                event_p1p2 = event_data_joint_(1:2:(num_pCells-1), :);
                event_p1   = event_data_margn_(1:2:(num_pCells-1), :);
                event_p2   = event_data_margn_(2:2:num_pCells, :);
                % Def.-1
%                 ratio_     = event_p1p2;
%                 ratio_(isnan(ratio_)) = 1; ratio_(isinf(ratio_)) = 1;
                % Def.-2
                ratio_     = ( event_p1p2 ./ event_p1 ./ event_p2 );
                ratio_(isnan(ratio_)) = 1; ratio_(isinf(ratio_)) = 1;
                % Def.-3
%                 ratio_     = log2(event_p1p2)-log2(event_p1)-log2(event_p2);
%                 ratio_(isnan(ratio_)) = 0; ratio_(isinf(ratio_)) = 0;
                % Def.-4
%                 ratio_     = log2(event_p1p2)-log2(event_p1)-log2(event_p2);
%                 ratio_(isnan(ratio_)) = 0; ratio_(isinf(ratio_)) = 0;
%                 ratio_     = ratio_ .* event_p1p2;
                % synch_ratio
                data_ratio_ = nan(size(event_data_margn_));
                data_ratio_(1:2:(num_pCells-1), :) = ratio_;
                data_ratio_(2:2: num_pCells   , :) = ratio_;
                synch_ratio.(variable)(counter_tag).(event_type_name){counter_var, counter_ang} = data_ratio_;
            end
        end
    end
end
fprintf(' --> Completed. \n')
end
%% function synchrony_mutual_info
function synch_MI = synchrony_mutual_info(synch_joint_MI, synch_cross_MI, synch_margn_MI, variable_MI)
%% Handle inputs
if nargin < 3
    variable_MI = 'amp';
end
if nargin < 2
    synch_MI = [];
    return;
end

%% Calc necessary variables
JSP_global_variables();
global event_type_list inds_span
if isempty(event_type_list)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
length_trace = length(inds_span);
num_pCells  = size(synch_joint_MI.(variable_MI)(1).onset{1, 1}, 1);
num_tag_bin = length(synch_joint_MI.(variable_MI));
num_ang_bin = size(synch_joint_MI.(variable_MI)(1).onset, 2);
num_var_bin = size(synch_joint_MI.(variable_MI)(1).onset, 1);

%% Init population_avg_pCells population_sem_pCells
synch_MI = synch_joint_MI;

%% Compute population_avg_pCells
fprintf(['synchrony_mutual_info' ' ...'])
for counter_event_type = 1 : length(event_type_list)
    event_type_name = event_type_list{counter_event_type};
    for counter_tag = 1 : num_tag_bin
        for counter_ang = 1 : num_ang_bin
            for counter_var = 1 : num_var_bin
                event_data_joint_ = ...
                    synch_joint_MI.(variable_MI)(counter_tag).(event_type_name){counter_var, counter_ang};
                event_data_cross_ = ...
                    synch_cross_MI.(variable_MI)(counter_tag).(event_type_name){counter_var, counter_ang};
                event_data_margn_ = ...
                    synch_margn_MI.(variable_MI)(counter_tag).(event_type_name){counter_var, counter_ang};
                event_p1p2 = event_data_joint_(1:2:(num_pCells-1), :);
                event_n1n2 = event_data_joint_(2:2:num_pCells, :);
                event_p1n2 = event_data_cross_(1:2:(num_pCells-1), :);
                event_n1p2 = event_data_cross_(2:2:num_pCells, :);
                event_p1   = event_data_margn_(1:2:(num_pCells-1), :);
                event_p2   = event_data_margn_(2:2:num_pCells, :);
                event_n1   = 1 - event_p1;
                event_n2   = 1 - event_p2;
                % % % % % MI_
                term_p1p2 = event_p1p2 .* log2(event_p1p2 ./ event_p1 ./ event_p2);
                term_n1n2 = event_n1n2 .* log2(event_n1n2 ./ event_n1 ./ event_n2);
                term_p1n2 = event_p1n2 .* log2(event_p1n2 ./ event_p1 ./ event_n2);
                term_n1p2 = event_n1p2 .* log2(event_n1p2 ./ event_n1 ./ event_p2);
                term_p1p2(isnan(term_p1p2)) = 0; term_p1p2(isinf(term_p1p2)) = 0;
                term_n1n2(isnan(term_n1n2)) = 0; term_n1n2(isinf(term_n1n2)) = 0;
                term_p1n2(isnan(term_p1n2)) = 0; term_p1n2(isinf(term_p1n2)) = 0;
                term_n1p2(isnan(term_n1p2)) = 0; term_n1p2(isinf(term_n1p2)) = 0;
                MI_ = term_p1p2 + term_n1n2 + term_p1n2 + term_n1p2;
                % % % % %  CFI_MI_
                H_1 = -(event_p1.*log2(event_p1)) - (event_n1.*log2(event_n1));
                H_2 = -(event_p2.*log2(event_p2)) - (event_n2.*log2(event_n2));
                H_1(isnan(H_1)) = 0; H_1(isinf(H_1)) = 0;
                H_2(isnan(H_2)) = 0; H_2(isinf(H_2)) = 0;
                H_ = H_1;
                H_(:,:,2) = H_2;
                H_min = min(H_,[],3);
                MI_norm = MI_ ./ H_min;
                MI_norm(isnan(MI_norm)) = 0; MI_norm(isinf(MI_norm)) = 0;
                term_p_p1p2 = (event_p1p2./event_p2);
                term_p_n1n2 = (event_n1n2./event_n2);
                term_p_p1n2 = (event_p1n2./event_n2);
                term_p_n1p2 = (event_n1p2./event_p2);
                term_p_p1p2(isnan(term_p_p1p2)) = 0; term_p_p1p2(isinf(term_p_p1p2)) = 0;
                term_p_n1n2(isnan(term_p_n1n2)) = 0; term_p_n1n2(isinf(term_p_n1n2)) = 0;
                term_p_p1n2(isnan(term_p_p1n2)) = 0; term_p_p1n2(isinf(term_p_p1n2)) = 0;
                term_p_n1p2(isnan(term_p_n1p2)) = 0; term_p_n1p2(isinf(term_p_n1p2)) = 0;
                p_c  = 0.5 .* (term_p_p1p2 + term_p_n1n2);
                p_ac = 0.5 .* (term_p_p1n2 + term_p_n1p2);
                condition_ = ones(size(p_c));
                condition_(p_c < p_ac) = -1;
                CFI_MI_ = condition_ .* MI_norm;
                % % % % %  STTC
                P_1 = event_p1p2 ./ event_p1;
                P_2 = event_p1p2 ./ event_p2;
                T_1 = event_p1;
                T_2 = event_p2;
                P_1(isnan(P_1)) = 0; P_1(isinf(P_1)) = 0;
                P_2(isnan(P_2)) = 0; P_2(isinf(P_2)) = 0;
                term_STTC_1 = (P_1 - T_2) ./ (1 - P_1.*T_2);
                term_STTC_2 = (P_2 - T_1) ./ (1 - P_2.*T_1);
                term_STTC_1(isnan(term_STTC_1)) = 0; term_STTC_1(isinf(term_STTC_1)) = 0;
                term_STTC_2(isnan(term_STTC_2)) = 0; term_STTC_2(isinf(term_STTC_2)) = 0;
                STTC_ = 0.5 .* (term_STTC_1 + term_STTC_2);
                % % % % %  synch_MI
                data_MI_ = nan(size(event_data_margn_));
                data_MI_(1:2:(num_pCells-1), :) = MI_; % STTC_; % CFI_MI_; % 
                data_MI_(2:2: num_pCells   , :) = MI_; % STTC_; % CFI_MI_; % 
                synch_MI.(variable_MI)(counter_tag).(event_type_name){counter_var, counter_ang} = data_MI_;
            end
        end
    end
end
fprintf(' --> Completed. \n')
end

%% PLOT FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function plot_sac_sorter
function plot_sac_sorter(SACS_ALL_DATA, params)
%% Set parameters
amp_edges = -.25 : 0.5 : 15.25;
ang_edges = (-pi-(pi/16)) : (pi/8) : (pi-(pi/16));
react_edges = -12.5: 25 : 512.5;
num_row = 9;
num_col = 9;

%% Init plot
hFig = figure(1);
clf(hFig)
hold on

%% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for counter_tag = 1 : 9
if counter_tag == 8
    idx_tag = (SACS_ALL_DATA.tag==8) | (SACS_ALL_DATA.tag==9);
    title_ = [params.sac_tag_list{8} ' & ' params.sac_tag_list{9}];
elseif counter_tag == 9
    idx_tag = (SACS_ALL_DATA.tag==10);
    title_ = params.sac_tag_list{10};
else
    idx_tag = (SACS_ALL_DATA.tag==counter_tag);
    title_ = params.sac_tag_list{counter_tag};
end

axes_minor_nums = reshape(1:num_row*num_col, num_row, num_col)';
axes_main_row = floor((counter_tag - 1) / 3)+1;
axes_main_col = mod(counter_tag, 3); if (axes_main_col==0); axes_main_col=3; end
row1_ = ((axes_main_row-1)*3)+1;
row2_ = ((axes_main_row-1)*3)+2;
row3_ = ((axes_main_row-1)*3)+3;
col1_ = ((axes_main_col-1)*3)+1;
col2_ = ((axes_main_col-1)*3)+2;
col3_ = ((axes_main_col-1)*3)+3;

axes_trace = [axes_minor_nums(row1_,col1_), axes_minor_nums(row1_,col2_), axes_minor_nums(row1_,col3_)...
              axes_minor_nums(row2_,col1_), axes_minor_nums(row2_,col2_), axes_minor_nums(row2_,col3_) ];
axes_amp_dis = axes_minor_nums(row3_,col1_);
axes_ang_dis = axes_minor_nums(row3_,col2_);
axes_react_dis = axes_minor_nums(row3_,col3_);

subplot(num_row,num_col,axes_trace)
hold on
plot(SACS_ALL_DATA.eye_r_px(:,idx_tag), ...
     SACS_ALL_DATA.eye_r_py(:,idx_tag), 'k')
plot(SACS_ALL_DATA.eye_r_px_offset(:,idx_tag), ...
     SACS_ALL_DATA.eye_r_py_offset(:,idx_tag), 'om')
title([title_ ': ' num2str(nansum(idx_tag)) ' sac'], 'interpret', 'none');
xlim([-17, 17])
ylim([-15, 15])
% axis equal;

subplot(num_row,num_col,axes_amp_dis)
hold on
histogram(SACS_ALL_DATA.eye_r_amp_m(:,idx_tag), amp_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(SACS_ALL_DATA.eye_r_amp_m(:,idx_tag), amp_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 15])
set(gca, 'XTick', 0:3:15)
ylabel('Amplitude')

subplot(num_row,num_col,axes_ang_dis)
polarhistogram(deg2rad(SACS_ALL_DATA.eye_r_ang(:,idx_tag)), (ang_edges), 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
hold on
polarhistogram(deg2rad(SACS_ALL_DATA.eye_r_ang(:,idx_tag)), (ang_edges), 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
set(gca, 'ThetaTick', [])
set(gca, 'RTick', [])
set(gca, 'Title', [])

subplot(num_row,num_col,axes_react_dis)
hold on
histogram(SACS_ALL_DATA.reaction(:,idx_tag)/1000, react_edges/1000, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
histogram(SACS_ALL_DATA.reaction(:,idx_tag)/1000, react_edges/1000, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
xlim([0, 500]/1000)
set(gca, 'XTick', (0:200:500)/1000)
ylabel('Reaction')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add the title info
sgtitle([params.cell_name ', ' ...
    'trial: ' num2str(params.num_trials) ', ' ...
    'sac: ' num2str(length(SACS_ALL_DATA.validity)) ', ' ...
    'dur: ' num2str(params.duration/60,3) 'min' ...
    ], ...
    'interpret', 'none');
ESN_Beautify_Plot(hFig, [13 13], 8)

end
%% function plot_neural_properties
function plot_neural_properties(fig_num)
%% Load population_neural_properties
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
path_cell_data = [path_cell_data '..' filesep];
load([path_cell_data, 'population_neural_properties.mat'], 'population_neural_properties');

%% Init plot
hFig = figure(fig_num);
clf(hFig)
num_row_fig = 1;
num_col_fig = 5;
SS_firing_edges = 5:10:135;
CS_firing_edges = 0.25:0.1:1.55;
CS_suppression_edges = 5.5: 1 : 25.5;

%% Calc variables
JSP_global_variables();
global waveform_inds_span
if isempty(waveform_inds_span)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
num_pCells          = size(population_neural_properties.SS_firing_rate, 1);
% XProb
SS_firing_pCells    = population_neural_properties.SS_firing_rate;
CS_firing_pCells    = population_neural_properties.CS_firing_rate;
SS_waveform_pCells  = population_neural_properties.SS_waveform;
CS_waveform_pCells  = population_neural_properties.CS_waveform;
SS_xprob_pCells     = population_neural_properties.Corr_data_SS_SSxSS_AUTO;
CS_xprob_pCells     = population_neural_properties.Corr_data_CS_CSxSS_AUTO;
CS_suppression_time = population_neural_properties.CS_suppression_time;
SS_waveform_pCells_norm = SS_waveform_pCells ./ repmat(max(abs(SS_waveform_pCells), [],2), 1, size(SS_waveform_pCells, 2));
CS_waveform_pCells_norm = CS_waveform_pCells ./ repmat(max(abs(SS_waveform_pCells), [],2), 1, size(CS_waveform_pCells, 2)); % normalize to SS max and not CS max
SS_xprob_pCells_norm    = SS_xprob_pCells ./ repmat(SS_firing_pCells, 1, size(SS_xprob_pCells, 2)) .* 1000;
CS_xprob_pCells_norm    = CS_xprob_pCells ./ repmat(SS_firing_pCells, 1, size(CS_xprob_pCells, 2)) .* 1000; % normalize to SS firing and not CS firing
time_waveform = (waveform_inds_span ./ 30e3) * 1000;
time_xprob = nanmean(population_neural_properties.Corr_data_SS_inds_span .* ...
    repmat(population_neural_properties.Corr_data_SS_bin_size_time, 1, size(population_neural_properties.Corr_data_SS_inds_span, 2))) * 1000;

% Firing rate
SS_firing_pCells_mean = nanmean(SS_firing_pCells);
SS_firing_pCells_stdv = nanstd( SS_firing_pCells);
SS_firing_pCells_sem  = nanstd( SS_firing_pCells)./sqrt(num_pCells);
CS_firing_pCells_mean = nanmean(CS_firing_pCells);
CS_firing_pCells_stdv = nanstd( CS_firing_pCells);
CS_firing_pCells_sem  = nanstd( CS_firing_pCells)./sqrt(num_pCells);

% Waveform
SS_waveform_mean = nanmean(SS_waveform_pCells_norm);
SS_waveform_stdv = nanstd( SS_waveform_pCells_norm);%./sqrt(num_pCells);
SS_waveform_stdv_p = SS_waveform_mean + SS_waveform_stdv;
SS_waveform_stdv_m = SS_waveform_mean - SS_waveform_stdv;
SS_waveform_stdv_y_axes = [(SS_waveform_stdv_p) flip(SS_waveform_stdv_m)];
SS_waveform_stdv_x_axes = [(time_waveform) flip(time_waveform)];

CS_waveform_mean = nanmean(CS_waveform_pCells_norm);
CS_waveform_stdv = nanstd( CS_waveform_pCells_norm);%./sqrt(num_pCells);
CS_waveform_stdv_p = CS_waveform_mean + CS_waveform_stdv;
CS_waveform_stdv_m = CS_waveform_mean - CS_waveform_stdv;
CS_waveform_stdv_y_axes = [(CS_waveform_stdv_p) flip(CS_waveform_stdv_m)];
CS_waveform_stdv_x_axes = [(time_waveform) flip(time_waveform)];

SS_xprob_pCells_norm(:,round(size(SS_xprob_pCells_norm,2)/2)) = nan;
SS_xprob_mean = nanmean(SS_xprob_pCells_norm);
SS_xprob_stdv = nanstd( SS_xprob_pCells_norm);%./sqrt(num_pCells);
SS_xprob_stdv_p = SS_xprob_mean + SS_xprob_stdv;
SS_xprob_stdv_m = SS_xprob_mean - SS_xprob_stdv;
SS_xprob_stdv_y_axes = [(SS_xprob_stdv_p) flip(SS_xprob_stdv_m) (SS_xprob_stdv_p(1))];
SS_xprob_stdv_x_axes = [(time_xprob) flip(time_xprob) (time_xprob(1))];

CS_xprob_mean = nanmean(CS_xprob_pCells_norm);
CS_xprob_stdv = nanstd( CS_xprob_pCells_norm);%./sqrt(num_pCells);
CS_xprob_stdv_p = CS_xprob_mean + CS_xprob_stdv;
CS_xprob_stdv_m = CS_xprob_mean - CS_xprob_stdv;
CS_xprob_stdv_y_axes = [(CS_xprob_stdv_p) flip(CS_xprob_stdv_m)];
CS_xprob_stdv_x_axes = [(time_xprob) flip(time_xprob)];

% CS_suppression_time
CS_suppression_mean =  nanmean(CS_suppression_time);
CS_suppression_stdv =  nanstd( CS_suppression_time);

%% Firing rate
subplot(num_row_fig,num_col_fig, 1)
hold on
histogram(SS_firing_pCells, SS_firing_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'b')
histogram(SS_firing_pCells, SS_firing_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'b', 'FaceColor', 'none', 'linewidth', 1)
xline((SS_firing_pCells_mean+SS_firing_pCells_stdv), '-b', 'LineWidth', 0.5);
xline((SS_firing_pCells_mean-SS_firing_pCells_stdv), '-b', 'LineWidth', 0.5);
xline(SS_firing_pCells_mean,                         '-b', 'LineWidth', 1.0);
xlabel('SS Firing Rate (Hz)')
ylabel('Count (#)')
% title(stat_SS, 'interpreter', 'none')

subplot(num_row_fig,num_col_fig, 2)
hold on
histogram(CS_firing_pCells, CS_firing_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
histogram(CS_firing_pCells, CS_firing_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 1)
xline((CS_firing_pCells_mean+CS_firing_pCells_stdv), '-r', 'LineWidth', 0.5);
xline((CS_firing_pCells_mean-CS_firing_pCells_stdv), '-r', 'LineWidth', 0.5);
xline((CS_firing_pCells_mean),                       '-r', 'LineWidth', 1.0);
xlabel('CS Firing Rate (Hz)')
ylabel('Count (#)')
% title(stat_CS,  'interpreter', 'none')

%% Waveform
subplot(num_row_fig,num_col_fig, 3)
hold on
% plot(time_waveform, SS_waveform_stdv_m, '-b', 'LineWidth', 0.5)
% plot(time_waveform, SS_waveform_stdv_p, '-b', 'LineWidth', 0.5)
plot(SS_waveform_stdv_x_axes, SS_waveform_stdv_y_axes, '-b', 'LineWidth', 0.5)
plot(time_waveform, SS_waveform_mean, '-b', 'LineWidth', 1.0)

% plot(time_waveform, CS_waveform_stdv_m, '-r', 'LineWidth', 0.5)
% plot(time_waveform, CS_waveform_stdv_p, '-r', 'LineWidth', 0.5)
plot(CS_waveform_stdv_x_axes, CS_waveform_stdv_y_axes, '-r', 'LineWidth', 0.5)
plot(time_waveform, CS_waveform_mean, '-r', 'LineWidth', 1.0)
ylabel('waveform')
xlabel('Time (ms)')
ylim([-1.3 +1.2])
xlim([-2 4])

subplot(num_row_fig,num_col_fig, 4)
hold on
% plot(time_xprob, SS_xprob_stdv_p, '-b', 'LineWidth', 0.5)
% plot(time_xprob, SS_xprob_stdv_m, '-b', 'LineWidth', 0.5)
plot(SS_xprob_stdv_x_axes, SS_xprob_stdv_y_axes, '-b', 'LineWidth', 0.5)
plot(time_xprob, SS_xprob_mean, '-b', 'LineWidth', 1.0)

% plot(time_xprob, CS_xprob_stdv_p, '-r', 'LineWidth', 0.5)
% plot(time_xprob, CS_xprob_stdv_m, '-r', 'LineWidth', 0.5)
plot(CS_xprob_stdv_x_axes, CS_xprob_stdv_y_axes, '-r', 'LineWidth', 0.5)
plot(time_xprob, CS_xprob_mean, '-r', 'LineWidth', 1.0)
ylabel('prob')
xlabel('Time (ms)')
ylim([-0.2 +1.6])
xlim([-50 50])

%% Suppression
subplot(num_row_fig,num_col_fig, 5)
hold on
histogram(CS_suppression_time, CS_suppression_edges,  'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5])
histogram(CS_suppression_time, CS_suppression_edges,  'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 1)
xline((CS_suppression_mean+CS_suppression_stdv), '-r', 'LineWidth', 0.5);
xline((CS_suppression_mean-CS_suppression_stdv), '-r', 'LineWidth', 0.5);
xline((CS_suppression_mean),                     '-r', 'LineWidth', 1.0);
ylabel('Count')
xlabel('CS suppression (ms)')
%% ESN_Beautify_Plot
stat_SS = ['SS_firing, ', 'mean: ', num2str(SS_firing_pCells_mean,3), ', std: ', num2str(SS_firing_pCells_stdv,3)];
stat_CS = ['CS_firing, ', 'mean: ', num2str(CS_firing_pCells_mean,3), ', std: ', num2str(CS_firing_pCells_stdv,3)];
stat_suppression = ['Suppression, ', 'mean: ', num2str(CS_suppression_mean,3), ', std: ', num2str(CS_suppression_stdv,3)];
sgtitle({[stat_SS, ' | ', stat_CS], stat_suppression}, ...
    'interpret', 'none', 'FontSize', 8);
ESN_Beautify_Plot(hFig, [8, 2], 8)

%% Save figs
path_fig_ = [path_cell_data 'population_figs'];
if ~exist(path_fig_, 'dir')
    mkdir(path_fig_);
end
file_name_fig_ = 'neural_properties';
saveas(hFig,[path_fig_ filesep file_name_fig_], 'pdf');

end
%% function plot_CS_on_properties
function plot_CS_on_properties(fig_num)
%% Load population_neural_properties
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
path_cell_data = [path_cell_data '..' filesep];
load([path_cell_data, 'population_neural_properties.mat'], 'population_neural_properties');

%% Init plot
hFig = figure(fig_num);
clf(hFig)
num_row_fig = 1;
num_col_fig = 4;
num_pCells = size(population_neural_properties.CS_ang_avg, 1);
step_size_ = 22.5;
ang_edges = 0-(step_size_/2):step_size_:360-(step_size_/2);
sig_edges = 35:2:61;

CS_ang_avg = population_neural_properties.CS_ang_avg;
vonMises_std = population_neural_properties.vonMises_std;
CS_prob_avg_tuned = population_neural_properties.CS_prob_avg_tuned;
overall_prob_TUNED_mean = nanmean(CS_prob_avg_tuned, 1);
overall_prob_TUNED_stdv = nanstd(CS_prob_avg_tuned, 0, 1) ./ sqrt(num_pCells);
overall_prob_TUNED_stdv_p = overall_prob_TUNED_mean + overall_prob_TUNED_stdv;
overall_prob_TUNED_stdv_m = overall_prob_TUNED_mean - overall_prob_TUNED_stdv;

%% Plot CS-on distribution
subplot(num_row_fig, num_col_fig, 1);
idx_pCells = 1:num_pCells; % All
polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), deg2rad(ang_edges), 'DisplayStyle', 'bar','FaceColor',[0.6 0.6 0.6], 'EdgeColor', 'none')
hold on
polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), deg2rad(ang_edges), 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'r', 'linewidth', 1)
rlim([0 25])
set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:5:25,...
    'RTickLabel', {'', '', '10', '', '20', ''}, 'ThetaTickLabel', {'0','','90','', '180','','270', ''})
title('CS-on mean Dist.')

%% Plot std distribution
subplot(num_row_fig, num_col_fig, 2);
histogram(vonMises_std, sig_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', [0.6 0.6 0.6])
hold on
histogram(vonMises_std, sig_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 1)
xline(nanmean(vonMises_std),'Color', 'r', 'linewidth', 1);
ylabel('count')
xlabel('Circular std. (deg)')
title('CS-on stdv Dist.')

%% Plot CS tuning circular
JSP_global_variables();
global ang_values
if isempty(ang_values)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end

plot_data_amp_mean = [overall_prob_TUNED_mean, overall_prob_TUNED_mean(1), nan]';
plot_data_deg_mean = [ang_values, ang_values(1), nan]';

plot_data_amp_stdv_p = [overall_prob_TUNED_stdv_p, overall_prob_TUNED_stdv_p(1), nan]';
plot_data_deg_stdv_p = [ang_values, ang_values(1), nan]';

plot_data_amp_stdv_m = [overall_prob_TUNED_stdv_m, overall_prob_TUNED_stdv_m(1), nan]';
plot_data_deg_stdv_m = [ang_values, ang_values(1), nan]';

subplot(num_row_fig, num_col_fig, 3);
polarplot(deg2rad(plot_data_deg_stdv_p),plot_data_amp_stdv_m, '-k', 'LineWidth', 0.5)
hold on
polarplot(deg2rad(plot_data_deg_stdv_m),plot_data_amp_stdv_p, '-k', 'LineWidth', 0.5)
polarplot(deg2rad(plot_data_deg_mean),plot_data_amp_mean, '-k', 'LineWidth', 1)
% rlim([0 0.25])
set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:0.05:0.25, ...
    'RTickLabel', {'', '', '0.1', '', '0.2', ''}, 'ThetaTickLabel', {'0','','90','', '180','','270', ''})
title('CS Tuning', 'Interpreter', 'none');

%% Plot CS tuning linear, CS-on at center
subplot(num_row_fig, num_col_fig, 4);
hold on
% plot_order_ = [6 7 8 1 2 3 4 5 6];
plot_order_ = [5 6 7 8 1 2 3 4 5];
plot(overall_prob_TUNED_stdv_p(plot_order_), '-k', 'LineWidth', 0.5)
plot(overall_prob_TUNED_stdv_m(plot_order_), '-k', 'LineWidth', 0.5)
plot(overall_prob_TUNED_mean(plot_order_), '-k', 'LineWidth', 1)
ylabel('CS probability');
xlabel('Direction')
% set(gca, 'XTick', 1:1:8, 'XTickLabel', {'', '-90','','ON','','90','','180',''})
set(gca, 'XTick', 1:1:9, 'XTickLabel', {'-180', '', '-90','','ON','','90','','180'})

avg_prob_TUNED_mean = nanmean(nanmean(CS_prob_avg_tuned, 2));
avg_prob_TUNED_stdv = nanstd(nanmean(CS_prob_avg_tuned, 2), 0, 1) ./ sqrt(num_pCells);
avg_prob_TUNED_mean = repmat(avg_prob_TUNED_mean, 1, size(CS_prob_avg_tuned,2));
avg_prob_TUNED_stdv = repmat(avg_prob_TUNED_stdv, 1, size(CS_prob_avg_tuned,2));
avg_prob_TUNED_stdv_p = avg_prob_TUNED_mean + avg_prob_TUNED_stdv;
avg_prob_TUNED_stdv_m = avg_prob_TUNED_mean - avg_prob_TUNED_stdv;

plot(avg_prob_TUNED_stdv_p(plot_order_), '-k', 'LineWidth', 0.5)
plot(avg_prob_TUNED_stdv_m(plot_order_), '-k', 'LineWidth', 0.5)
plot(avg_prob_TUNED_mean(plot_order_), '-k', 'LineWidth', 1)
% ylim([0.05 0.25])
title('CS Tuning', 'Interpreter', 'none');

%% ESN_Beautify_Plot
ESN_Beautify_Plot(hFig, [8, 2], 8)

%% Save figs
path_fig_ = [path_cell_data 'population_figs'];
if ~exist(path_fig_, 'dir')
    mkdir(path_fig_);
end
file_name_fig_ = 'CS_on_properties';
saveas(hFig,[path_fig_ filesep file_name_fig_], 'pdf');

%% CS-on distribution for different vermal areas
% if num_pCells ~= 143
%     error('plot_CS_on_properties: number of P-cells is not 143. Please modify the code.')
% end
% 
% hFig = figure(fig_num+1);
% clf(hFig)
% num_row_fig = 2;
% num_col_fig = 3;
% 
% subplot(num_row_fig, num_col_fig, 2);
% idx_pCells = [1:18, 52:68];
% polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), deg2rad(ang_edges), 'DisplayStyle', 'bar','FaceColor',[0.6 0.6 0.6], 'EdgeColor', 'none')
% hold on
% polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), deg2rad(ang_edges), 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'r', 'linewidth', 1)
% rlim([0 10])
% set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:5:10,...
%     'RTickLabel', {'', '', '10'}, 'ThetaTickLabel', {'0','','90','', '180','','270', ''})
% title(['(M) center, n=', num2str(length(idx_pCells))])
% 
% subplot(num_row_fig, num_col_fig, 3);
% idx_pCells = 19:51;
% polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), deg2rad(ang_edges), 'DisplayStyle', 'bar','FaceColor',[0.6 0.6 0.6], 'EdgeColor', 'none')
% hold on
% polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), deg2rad(ang_edges), 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'r', 'linewidth', 1)
% rlim([0 10])
% set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:5:10,...
%     'RTickLabel', {'', '', '10'}, 'ThetaTickLabel', {'0','','90','', '180','','270', ''})
% title(['(M) right, n=', num2str(length(idx_pCells))])
% 
% subplot(num_row_fig, num_col_fig, 1);
% idx_pCells = 69:90;
% polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), deg2rad(ang_edges), 'DisplayStyle', 'bar','FaceColor',[0.6 0.6 0.6], 'EdgeColor', 'none')
% hold on
% polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), deg2rad(ang_edges), 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'r', 'linewidth', 1)
% rlim([0 10])
% set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:5:10,...
%     'RTickLabel', {'', '', '10'}, 'ThetaTickLabel', {'0','','90','', '180','','270', ''})
% title(['(M) left, n=', num2str(length(idx_pCells))])
% 
% subplot(num_row_fig, num_col_fig, 6);
% idx_pCells = 91:143;
% polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), deg2rad(ang_edges), 'DisplayStyle', 'bar','FaceColor',[0.6 0.6 0.6], 'EdgeColor', 'none')
% hold on
% polarhistogram(deg2rad(CS_ang_avg(idx_pCells)), deg2rad(ang_edges), 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'r', 'linewidth', 1)
% rlim([0 15])
% set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:5:15,...
%     'RTickLabel', {'', '', '10', ''}, 'ThetaTickLabel', {'0','','90','', '180','','270', ''})
% title(['(R) right, n=', num2str(length(idx_pCells))])
% 
% ESN_Beautify_Plot(hFig, [6, 4], 8)
% file_name_fig_ = 'CS_on_properties_2';
% saveas(hFig,[path_fig_ filesep file_name_fig_], 'pdf');

%% PAIR CS-on difference
% Load CS_on_population_pairs
% path_pair_data = uigetdir;
% if ~strcmp(path_pair_data(end), filesep);path_pair_data = [path_pair_data filesep];end
% path_pair_data = [path_pair_data '..' filesep];
% load([path_pair_data, 'CS_on_population_pairs.mat'],'CS_on_population');
% num_pairs = length(CS_on_population)/2;
% CS_ang_avg_pairs = nan(num_pairs,2);
% for counter_pair = 1 : num_pairs
%     CS_ang_avg_pairs(counter_pair, 1) = CS_on_population(2*counter_pair-1).CS_ang_avg;
%     CS_ang_avg_pairs(counter_pair, 2) = CS_on_population(2*counter_pair  ).CS_ang_avg;
% end
% x_values_ = cosd(CS_ang_avg_pairs); 
% y_values_ = sind(CS_ang_avg_pairs);
% diff_ang = acosd( (x_values_(:,1) .* x_values_(:,2)) + (y_values_(:,1) .* y_values_(:,2)) );
% vec_1 = [x_values_(:,1)'; y_values_(:,1)'; zeros(1,num_pairs)];
% vec_2 = [x_values_(:,2)'; y_values_(:,2)'; zeros(1,num_pairs)];
% cross_ang = cross(vec_1, vec_2);
% diff_CS_ang_avg_pairs = diff_ang .* sign(cross_ang(3,:)');
% 
% step_size_ = 22.5;
% ang_edges = -135-(step_size_/2):step_size_:135+(step_size_/2);
% 
% pair_distance = estimate_pair_distance();
% pair_distance = pair_distance(1:2:end);
% 
% hFig = figure(fig_num+2);
% clf(hFig)
% subplot(1,2,1)
% hold on
% histogram(diff_CS_ang_avg_pairs, ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
% histogram(diff_CS_ang_avg_pairs, ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 1)
% xline(mean(diff_CS_ang_avg_pairs),'Color', 'r', 'linewidth', 1);
% set(gca, 'XTick', -135:45:135)
% ylabel('count')
% xlabel('diff in CS-on dir (deg)')
% title('CS-on Pairs')
% 
% subplot(1,2,2)
% hold on
% x_axis_data = pair_distance;
% y_axis_data = diff_ang;
% plot(x_axis_data, y_axis_data, 'o','MarkerSize',3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None')
% 
% P_ = polyfit(x_axis_data, y_axis_data, 1);
% y_axis_hat = polyval(P_,x_axis_data);
% plot(x_axis_data, y_axis_hat, '-', 'LineWidth', 1)
% 
% [b,~,~,~,stats] = regress(y_axis_data,[ones(size(x_axis_data)) x_axis_data]);
% % p_value = stats(3)
% 
% xlabel('Contact distance (um)')
% ylabel('Absolute CS-on diff. (deg)')
% 
% ESN_Beautify_Plot(hFig, [3 1.5], 8)
% file_name_fig_ = 'CS_on_properties_3';
% saveas(hFig,[path_fig_ filesep file_name_fig_], 'pdf');
end
%% function plot_modulation_z_score()
function plot_modulation_z_score(fig_num)
%% Load sac_modulation
clc; close all;
flag_pair_list = false;
JSP_global_variables(flag_pair_list);
global length_trace
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
pCell_ids = build_pCell_ids(flag_pair_list);
num_pCells = size(pCell_ids, 1);

modulation_trace_mean    = nan(length_trace, num_pCells);
modulation_trace_stdv    = nan(length_trace, num_pCells);
modulation_range_mean    = nan(1, num_pCells);
modulation_range_stdv    = nan(1, num_pCells);
modulation_z_score  = nan(1, num_pCells);
modulation_baseline = nan(1, num_pCells);
is_modulated        = false(1, num_pCells);
% Loop over pCells
fprintf('\nLoading sac_modulation ')
for counter_pCell = 1 : num_pCells
    fprintf('.')
    %% load SACS_ALL_DATA
    cell_file_name = pCell_ids{counter_pCell, 1};
    load([path_cell_data cell_file_name], 'sac_modulation');
    %% Stack data
    modulation_trace_mean(:,counter_pCell) = sac_modulation.modulation_trace_mean;
    modulation_trace_stdv(:,counter_pCell) = sac_modulation.modulation_trace_stdv;
    modulation_range_mean(:,counter_pCell) = sac_modulation.modulation_range_mean;
    modulation_range_stdv(:,counter_pCell) = sac_modulation.modulation_range_stdv;
    modulation_z_score(:,counter_pCell)    = sac_modulation.modulation_z_score;
    modulation_baseline(:,counter_pCell)   = sac_modulation.modulation_baseline;
    is_modulated(:,counter_pCell)          = sac_modulation.is_modulated;
end
if length_trace ~= 500
    error('sac_modulation_index: length_trace is not 500. Please modify the code.')
end
data_baseline = repmat(modulation_baseline, length_trace, 1);
data_trace = (modulation_trace_mean - data_baseline)' * 1000;
data_smooth = data_trace(:,201:350);

fprintf(' --> Completed. \n')

%% Cluster Using UMAP
[~, pca_mat, ~] = pca(data_smooth);
[reduction, umap, clusterIdentifiers, extras]=run_umap(data_smooth);

% Cluster Using UMAP
data_gmm_ = [reduction(:, 1), reduction(:, 2)];
n_component_gmm_ = 2;
gmm_model_ = fitgmdist(data_gmm_,n_component_gmm_);
idx = cluster(gmm_model_,data_gmm_);

% Set pauser as idx_1, burster as idx_2, not_modulated as idx_3
idx(~is_modulated) = 3;
trace_idx_1 = nanmean(data_smooth(idx==1,:));
trace_idx_2 = nanmean(data_smooth(idx==2,:));
if max(trace_idx_1) > max(trace_idx_2)
    idx(idx==1) = 4;
    idx(idx==2) = 1;
    idx(idx==4) = 2;
end
[sum(idx==1) sum(idx==2) sum(idx==3)]

%% Save umap_data
idx_pauser = (idx==1);
idx_burster = (idx==2);
idx_modulated = ~(idx==3);
idx_not_modulated = (idx==3);
save([path_cell_data '..' filesep 'umap_data.mat'], ...
    'idx', 'idx_pauser','idx_burster','idx_modulated','idx_not_modulated', ...
    'pca_mat', 'reduction', 'umap', 'clusterIdentifiers', ...
    'data_gmm_', 'n_component_gmm_', 'gmm_model_', ...
    '-v7.3');

%% Init plot
close all;
hFig = figure(fig_num);
clf(hFig)
num_row_fig = 3;
num_col_fig = 3;

%% Plot modulation_z_score distribution
subplot(num_row_fig, num_col_fig, [1 4]);
z_score_edges = 0 : 1 : 21; 
histogram(modulation_z_score, z_score_edges,'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', [0.6 0.6 0.6])
hold on
histogram(modulation_z_score, z_score_edges,'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 1)
xline(3,'Color', 'k', 'linewidth', 1);
ylabel('count')
xlabel('Z-score')
title('SS Modulation Z-score')

%% Plot UMAP
subplot(num_row_fig, num_col_fig, [7]);
hold on
plot((-49:100)',data_smooth(idx==2,:)', 'r')
plot((-49:100)',nanmean(data_smooth(idx==2,:))', 'k', 'LineWidth', 2)
% plot((-49:100)',data_smooth(9,:)', 'm')
% plot((-49:100)',data_smooth(27,:)', 'c')
ylim([-50 150])
ylabel('SS firing rate')

subplot(num_row_fig, num_col_fig, [9]);
hold on
plot((-49:100)',data_smooth(idx==3,:)', 'k')
plot((-49:100)',nanmean(data_smooth(idx==3,:))', 'k', 'LineWidth', 2)
% plot((-49:100)',data_smooth(110,:)', 'm')
% plot((-49:100)',data_smooth(23,:)', 'c')
ylim([-50 50])

subplot(num_row_fig, num_col_fig, [8]);
hold on
plot((-49:100)',data_smooth(idx==1,:)', 'b')
plot((-49:100)',nanmean(data_smooth(idx==1,:))', 'k', 'LineWidth', 2)
% plot((-49:100)',data_smooth(110,:)', 'm')
% plot((-49:100)',data_smooth(23,:)', 'c')
ylim([-100 100])
xlabel('Saccade onset (ms)')

% ylim([0 inf])

MarkerSize_ = 3;
subplot(num_row_fig, num_col_fig, [2 5]);
hold on
plot(reduction(idx==1, 1), reduction(idx==1, 2), 'ob',...
    'MarkerFaceColor', 'b','MarkerEdgeColor', 'None', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
% plot(reduction(9, 1), reduction(9, 2), 'ob',...
%     'MarkerFaceColor', 'm', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
% plot(reduction(27, 1), reduction(27, 2), 'ob',...
%     'MarkerFaceColor', 'c', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(reduction(idx==2, 1), reduction(idx==2, 2), 'or', ...
    'MarkerFaceColor', 'r','MarkerEdgeColor', 'None', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
% plot(reduction(110, 1), reduction(110, 2), 'or',...
%     'MarkerFaceColor', 'm', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
% plot(reduction(23, 1), reduction(23, 2), 'or',...
%     'MarkerFaceColor', 'c', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(reduction(idx==3, 1), reduction(idx==3, 2), 'ok', ...
    'MarkerFaceColor', 'k','MarkerEdgeColor', 'None', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
xlabel('umap 1')
ylabel('umap 2')
title('umap')

subplot(num_row_fig, num_col_fig, [3 6]);
hold on
plot(pca_mat(idx==1, 1), pca_mat(idx==1, 2), 'ob',...
    'MarkerFaceColor', 'b','MarkerEdgeColor', 'None', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
% plot(pca_mat(9, 1), pca_mat(9, 2), 'ob',...
%     'MarkerFaceColor', 'm', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
% plot(pca_mat(27, 1), pca_mat(27, 2), 'ob',...
%     'MarkerFaceColor', 'c', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(pca_mat(idx==2, 1), pca_mat(idx==2, 2), 'or', ...
    'MarkerFaceColor', 'r','MarkerEdgeColor', 'None', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
% plot(pca_mat(110, 1), pca_mat(110, 2), 'or',...
%     'MarkerFaceColor', 'm', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
% plot(pca_mat(23, 1), pca_mat(23, 2), 'or',...
%     'MarkerFaceColor', 'c', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
plot(pca_mat(idx==3, 1), pca_mat(idx==3, 2), 'ok', ...
    'MarkerFaceColor', 'k','MarkerEdgeColor', 'None', 'MarkerSize', MarkerSize_, 'linewidth', 0.5)
xlabel('pca 1')
ylabel('pca 2')
title('pca')

ESN_Beautify_Plot(hFig, [6,3], 8)

%% Plot sample trace
hFig = figure;
hold on
id_pCell = 42;
y_mean = data_trace(id_pCell,201:350)';
y_stdv = modulation_trace_stdv(201:350,id_pCell)*1000;
y_stdv_p = y_mean + y_stdv;
y_stdv_m = y_mean - y_stdv;
x_mean = (-49:100)';
plot(x_mean, y_stdv_p, 'k', 'LineWidth', 0.5)
plot(x_mean, y_stdv_m, 'k', 'LineWidth', 0.5)
plot(x_mean, y_mean, 'k', 'LineWidth', 1)
yline(max(y_mean));
yline(min(y_mean));
yline(max(y_mean)-(modulation_range_stdv(1,id_pCell)*1000));
xlabel('Saccade onset (ms)')
ylabel('SS firing rate')

fprintf(['range mean: ' num2str(modulation_range_mean(1,id_pCell)*1000) '\n'])
fprintf(['range stdv: ' num2str(modulation_range_stdv(1,id_pCell)*1000) '\n'])
fprintf(['z-score: ' num2str(modulation_z_score(1,id_pCell)) '\n'])

ESN_Beautify_Plot(hFig, [1.7,1.4], 8)

end
%% function plot_population_data_iteratively
function plot_population_data_iteratively(fig_num)
%% Set variables
JSP_global_variables();
global event_type_list tag_name_list 
% event_type_list = {'onset'} ;
% tag_name_list = {'prim_success'}; 

data_type_list = {'SS', 'VM'};
CSYS_type_list = {'absol', 'tuned'};
% CSYS_type_list = {'absol'};

num_tag = length(tag_name_list);


%% Load population_neural_properties
path_cell_data = uigetdir;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
path_cell_data = [path_cell_data '..' filesep];

if ~exist([path_cell_data 'population_figs'], 'dir')
    mkdir([path_cell_data 'population_figs']);
end
% load('umap_data.mat', 'idx_pauser', 'idx_burster','idx_modulated','idx_not_modulated');
% flag_pair_list = false; [idx_mirza, idx_ramon] = idx_mirza_ramon(flag_pair_list);

%% Loop over conditions
params.variable        = 'amp';
% idx_pCells: is a boolean array. 1 for including a pCell, and 0 for exluding a pCell
% idx_mirza; % idx_ramon; % idx_pauser; % idx_burster; % idx_modulated; % idx_not_modulated;
% params.idx_pCells = idx_mirza | idx_ramon;
% params.idx_pCells = idx_ramon;

if ~exist('population_neural_properties', 'var')
    load([path_cell_data 'population_neural_properties' '.mat'], 'population_neural_properties')
end
for counter_CSYS_type = 1 : length(CSYS_type_list)
    params.CSYS_type       = CSYS_type_list{counter_CSYS_type};
    if ~exist(['num_sac_' params.CSYS_type], 'var')
        load([path_cell_data 'num_sac_' params.CSYS_type '.mat'], ['num_sac_' params.CSYS_type])
    end
    eval(['num_sac_data = ' 'num_sac_' params.CSYS_type ';']);
    clearvars(['num_sac_' params.CSYS_type]);
    for counter_data_type = 1 : length(data_type_list)
        params.data_type       = data_type_list{counter_data_type};
        if ~exist([params.data_type '_population_' params.CSYS_type], 'var')
            load([path_cell_data params.data_type '_population_' params.CSYS_type '.mat'], [params.data_type '_population_' params.CSYS_type])
        end
        eval(['population_data = ' params.data_type '_population_' params.CSYS_type ';']);
        clearvars([params.data_type '_population_' params.CSYS_type]);
        if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
            eval(['firing_rate = population_neural_properties.' params.data_type '_firing_rate'  ';']);
            population_data = population_data_subtract_baseline(population_data, firing_rate, params.variable);
        end
        [population_dir, num_sac_dir] = population_data_avg_over_levels(population_data, num_sac_data, params.variable, 1);
%         [population_dir, num_sac_dir] = population_data_combine_levels(population_dir, num_sac_dir, params.variable, 2, [3 7]);
%         [population_dir, num_sac_dir] = population_data_idx_pCells(population_dir, num_sac_dir, params.variable, params.idx_pCells);
        if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
            population_dir = population_data_smooth_pCells(population_dir, num_sac_dir, params.variable);
        end
        [population_avg_dir, population_sem_dir] = population_data_avg_over_pCells(population_dir, num_sac_dir, params.variable);
        
        [population_allDir, num_sac_allDir] = population_data_avg_over_levels(population_dir, num_sac_dir, params.variable, 2);
        [population_avg_allDir, population_sem_allDir] = population_data_avg_over_pCells(population_allDir, num_sac_allDir, params.variable);
        for counter_event_type = 1 : length(event_type_list)
            params.event_type_name = event_type_list{counter_event_type};
            for counter_tag = 1 : num_tag
                params.tag_id          = counter_tag;
                params.fig_num = fig_num;
                fprintf(['### Plotting: ' params.CSYS_type ', ' params.data_type ', ' params.event_type_name ', ' tag_name_list{counter_tag} '. ###\n'])
                
                data_avg_dir     = population_avg_dir.(params.variable)(params.tag_id).(params.event_type_name);
                data_sem_dir     = population_sem_dir.(params.variable)(params.tag_id).(params.event_type_name);
                data_avg_allDir = population_avg_allDir.(params.variable)(params.tag_id).(params.event_type_name);
                data_sem_allDir = population_sem_allDir.(params.variable)(params.tag_id).(params.event_type_name);
                
                plot_population_data(params.fig_num, data_avg_dir, data_avg_allDir, params, data_sem_dir, data_sem_allDir);
                %% Save figs
                path_fig_ = [path_cell_data 'population_figs' filesep params.CSYS_type filesep params.data_type filesep num2str(counter_tag) '_' tag_name_list{counter_tag}];
                if ~exist(path_fig_, 'dir')
                    mkdir(path_fig_);
                end
                file_name_fig_ = [params.CSYS_type '_' params.data_type '_' num2str(counter_tag) '_' tag_name_list{counter_tag} '_' num2str(counter_event_type) '_' params.event_type_name];
                hFig_ = figure(params.fig_num);
                saveas(hFig_,[path_fig_ filesep file_name_fig_], 'pdf');
                close(hFig_)
            end
        end % counter_event_type
        clearvars([params.data_type '_population_' params.CSYS_type], 'population_data');
    end % counter_data_type
    clearvars(['num_sac_' params.CSYS_type], 'num_sac_data');
end % counter_CSYS_type
fprintf('### ALL DONE. ###\n')

end
%% function plot_population_data
function plot_population_data(fig_num, data_avg_dir, data_avg_allDir, params, data_sem_dir, data_sem_allDir)
%% Handle inputs
if nargin < 5
    flag_std = false;
    data_sem_dir = data_avg_dir;
    data_sem_allDir = data_avg_allDir;
elseif (nargin >= 5) && (nargin <= 6)
    flag_std = true;
end

%% Init plot
hFig = figure(fig_num);
clf(hFig)
num_row_fig = 3;
num_col_fig = 3;
ax_ang_id = [6, 3, 2, 1, 4, 7, 8, 9, 5];
num_ang_bin = size(data_avg_dir, 2);
num_var_bin = size(data_avg_dir, 1);
line_colors_ = [0,0,0; pink(round(1.5*num_var_bin))];
flag_empty_data = false(num_ang_bin+1, 1);

%% Plot
JSP_global_variables();
global inds_span ang_values tag_name_list length_trace
if isempty(inds_span)
    fprintf('><ERROR><: Global variables are empty.\n');
    return;
end
ang_values_ = [ang_values nan];
clearvars h_ax
for counter_ang = 1 : num_ang_bin+1
    h_ax(counter_ang) = subplot(num_row_fig, num_col_fig, ax_ang_id(counter_ang));
    hold on;
    for counter_var = 1 : num_var_bin
        if counter_ang == (num_ang_bin+1)
            data_mean_ = data_avg_allDir{counter_var, 1};
            data_sem_ = data_sem_allDir{counter_var, 1};
        else
            data_mean_ = data_avg_dir{counter_var, counter_ang};
            data_sem_ = data_sem_dir{counter_var, counter_ang};
        end
        if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
            if isfield(params,'scale')
                scale = params.scale;
            else
                scale = 1000.0;
            end
            data_mean_ = data_mean_ .* scale; % convert Pr to Firing/Hz
            data_sem_ = data_sem_ .* scale;
        end
        if strcmp(params.data_type, 'CS')
            length_data_plot = 400;
        else
            length_data_plot = 200;
        end
        data_mean_x_axis = reshape(inds_span, 1, []);
        ind_diff_ = round(length_trace - length_data_plot) / 2;
        idx_plot = (ind_diff_+1) : 1 : (length_trace-ind_diff_);
        data_mean_ = data_mean_(1, idx_plot);
        data_sem_  = data_sem_(1, idx_plot);
        data_mean_x_axis = data_mean_x_axis(1, idx_plot);
        data_sem_p_ = data_mean_ + data_sem_;
        data_sem_m_ = data_mean_ - data_sem_;
        data_sem_y_axis_ = [(data_sem_p_) flip(data_sem_m_)];
        data_sem_x_axis_ = [(data_mean_x_axis) flip(data_mean_x_axis)];
        yline(0);
        xline(50);
        xline(-50);
        xline(0);
        if flag_std
            plot(data_sem_x_axis_, data_sem_y_axis_, 'LineWidth', 0.25, 'Color', line_colors_(counter_var, :))
        end
        plot(data_mean_x_axis, data_mean_, 'LineWidth', 1, 'Color', line_colors_(counter_var, :))
        ylim_ = ylim;
        if ( ylim_(1)==-1 ) && ( ylim_(2)==1 )
            flag_empty_data(counter_ang) = true;
        end
        if counter_ang == (num_ang_bin+1)
            title('all dir.')
        else
            title([num2str(ang_values_(counter_ang)) ' dir.'])
        end
        if ang_values_(counter_ang) == 270
            xlabel(['Time from ' params.event_type_name ' (ms)']);
        end
        if ang_values_(counter_ang) == 180
            if strcmp(params.data_type, 'SS')
%                 ylabel('SS firing (change, Hz)');
                 ylabel('SS firing (population average, Hz)');

            elseif strcmp(params.data_type, 'CS')
                ylabel('CS firing (change, Hz)');
            elseif strcmp(params.data_type, 'VT')
                ylabel('Tangent velocity (deg/s)');
            elseif strcmp(params.data_type, 'VM')
                ylabel('Velocity magnitude (deg/s)');
            end
        end
    end
end

y_lim_ = zeros(length(h_ax), 2);
for counter_ax = 1 : length(h_ax)
    y_lim_(counter_ax, :) = ylim(h_ax(counter_ax));
end
y_lim_(flag_empty_data, :) = [];
if strcmp(params.data_type, 'SS')
    y_lim__ = [0 +100];
elseif strcmp(params.data_type, 'CS')
    y_lim__ = [-1 +2];
elseif strcmp(params.data_type, 'VT')
    y_lim__ = [-25 +650];
elseif strcmp(params.data_type, 'VM')
    y_lim__ = [-25 +650];
else
    y_lim__ = [min(y_lim_(:,1)) max(y_lim_(:,2))];
end
if isfield(params,'ylim')
    if isempty(params.ylim)
        y_lim__ = [min(y_lim_(:,1)) max(y_lim_(:,2))];
    else
        y_lim__ = params.ylim;
    end
end
for counter_ax = 1 : length(h_ax)
%     set(h_ax(counter_ax), 'ylim', y_lim__ + 15);
% y_lim_modified = [0 100];
%      set(h_ax(counter_ax), 'ylim', y_lim_modified);

end

%% ESN_Beautify_Plot
sgtitle([...
    tag_name_list{params.tag_id} ', ' ...
    params.data_type ', ' ...
    params.CSYS_type ', ' ...
    params.event_type_name ', ' ...
    params.variable ...
    ], ...
    'interpret', 'none', 'FontSize', 8);

ESN_Beautify_Plot(hFig, [4, 4], 8)

end
%% function plot_population_data_single_condition
function plot_population_data_single_condition()
%% clear
clc; clear;
%% close all
close all;
%% set params
params.data_type       = 'SS'; % SS % VM
params.CSYS_type       = 'absol'; % absol
params.event_type_name = 'onset'; % visual % onset % vmax % offset
params.variable        = 'amp'; % amp % vel
params.tag_id          =1;
params.flag_smooth_plot = true; % false; % 
params.fig_num = 3;
params.plot_mode = 2; % mode=1 collapse the amp/vel, mode=2 plots the amp/vel
fprintf('params --> Completed. \n')
%% Load data
% load('umap_data.mat', 'idx_pauser', 'idx_burster','idx_modulated','idx_not_modulated');
% flag_pair_list = false; [idx_mirza, idx_ramon] = idx_mirza_ramon(flag_pair_list);
if ~exist([params.data_type '_population_' params.CSYS_type], 'var')
    load([params.data_type '_population_' params.CSYS_type '.mat'], [params.data_type '_population_' params.CSYS_type])
end
if ~exist(['num_sac_' params.CSYS_type], 'var')
    load(['num_sac_' params.CSYS_type '.mat'], ['num_sac_' params.CSYS_type])
end 
if ~exist('population_neural_properties', 'var')
    load(['population_neural_properties' '.mat'], 'population_neural_properties')
end
eval(['population_data = ' params.data_type '_population_' params.CSYS_type ';']);
eval(['num_sac_data = ' 'num_sac_' params.CSYS_type ';']);
if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
    eval(['firing_rate = population_neural_properties.' params.data_type '_firing_rate'  ';']);
    population_data = population_data_subtract_baseline(population_data, firing_rate, params.variable);
end
fprintf('Load data --> Completed. \n')

%% MODE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.plot_mode == 1
% [population_dir, num_sac_dir] = population_data_avg_over_levels(population_data, num_sac_data, params.variable, 1, [3 4 5 6 7]); % based on vel, 250-750
[population_dir, num_sac_dir] = population_data_avg_over_level(population_data, num_sac_data, params.variable, 1); % based on all vel

% [population_dir, num_sac_dir] = population_data_combine_levels(population_dir, num_sac_dir, params.variable, 2, [3 7]);
% [population_dir, num_sac_dir] = population_data_combine_tags(population_dir, num_sac_dir, params.variable, [1]);

% idx_pCells: is a boolean array. 1 for including a pCell, and 0 for exluding a pCell
% idx_mirza; % idx_ramon; % idx_pauser; % idx_burster; % idx_modulated; % idx_not_modulated;
% idx_pCells = idx_mirza | idx_ramon;
% [population_dir, num_sac_dir] = population_data_idx_pCells(population_dir, num_sac_dir, params.variable, idx_pCells);

%
if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
    population_dir = population_data_smooth_pCells(population_dir, num_sac_dir, params.variable);
end
[population_avg_dir, population_sem_dir] = population_data_avg_over_pCells(population_dir, num_sac_dir, params.variable);
[population_allDir, num_sac_allDir] = population_data_avg_over_levels(population_dir, num_sac_dir, params.variable, 2);
[population_avg_allDir, population_sem_allDir] = population_data_avg_over_pCells(population_allDir, num_sac_allDir, params.variable);

data_avg_dir     = population_avg_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_dir     = population_sem_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_avg_allDir = population_avg_allDir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_allDir = population_sem_allDir.(params.variable)(params.tag_id).(params.event_type_name);

plot_population_data(params.fig_num, data_avg_dir, data_avg_allDir, params, data_sem_dir, data_sem_allDir);
end

%% MODE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.plot_mode == 2
    
% [population_data, num_sac_data] = population_data_combine_levels(population_data, num_sac_data, params.variable, 2, [3 7]);
[population_data, num_sac_data] = population_data_combine_tags(  population_data, num_sac_data, params.variable, [1]);

% idx_pCells: is a boolean array. 1 for including a pCell, and 0 for exluding a pCell
% idx_mirza; % idx_ramon; % idx_pauser; % idx_burster; % idx_modulated; % idx_not_modulated;
% idx_pCells = idx_mirza | idx_ramon;
% [population_data, num_sac_data] = population_data_idx_pCells(population_data, num_sac_data, params.variable, idx_pCells);

if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
    population_data = population_data_smooth_pCells(population_data, num_sac_data, params.variable);
end
[population_avg_dir, population_sem_dir] = population_data_avg_over_pCells(population_data, num_sac_data, params.variable);
[population_allDir, num_sac_allDir] = population_data_avg_over_levels(population_data, num_sac_data, params.variable, 2);
[population_avg_allDir, population_sem_allDir] = population_data_avg_over_pCells(population_allDir, num_sac_allDir, params.variable);

data_avg_dir     = population_avg_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_dir     = population_sem_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_avg_allDir = population_avg_allDir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_allDir = population_sem_allDir.(params.variable)(params.tag_id).(params.event_type_name);

plot_population_data(params.fig_num, data_avg_dir, data_avg_allDir, params, data_sem_dir, data_sem_allDir);
end

end
%% function plot_synchronous_data_single_condition
function plot_synchronous_data_single_condition()
%% clear
clc; clear;
%% close all
close all;
%% set params
params.data_type       = 'SS';
params.CSYS_type       = 'absol';
params.event_type_name = 'onset';
params.variable        = 'amp';
params.tag_id          = 1;
params.flag_smooth_plot = true;
params.fig_num = 3;
params.plot_mode = 1; % mode=1 collapse the amp/vel, mode=2 plots the amp/vel

%% Load data
if ~exist([params.data_type '_synch_joint_' params.CSYS_type], 'var')
    load([params.data_type '_synch_joint_' params.CSYS_type '.mat'], [params.data_type '_synch_joint_' params.CSYS_type])
end
if ~exist(['num_synch_' params.CSYS_type], 'var')
    load(['num_synch_' params.CSYS_type '.mat'], ['num_synch_' params.CSYS_type])
end 
% if ~exist('population_neural_properties', 'var')
%     load(['population_neural_properties' '.mat'], 'population_neural_properties')
% end
eval(['population_data = ' params.data_type '_synch_joint_' params.CSYS_type ';']);
eval(['num_sac_data = ' 'num_synch_' params.CSYS_type ';']);
% if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
%     eval(['firing_rate = population_neural_properties.' params.data_type '_firing_rate'  ';']);
%     population_data = population_data_subtract_baseline(population_data, firing_rate, params.variable);
% end

% idx_pCells: is a boolean array. 1 for including a pCell, and 0 for exluding a pCell
num_pCells = size(population_data.amp(1).visual{1,1},1);
idx_pCells = true(num_pCells, 1);
idx_pCells(2:2:num_pCells, 1) = false;
[population_data, num_sac_data] = population_data_idx_pCells(population_data, num_sac_data, params.variable, idx_pCells);

%% MODE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.plot_mode == 1
[population_dir, num_sac_dir] = population_data_avg_over_levels(population_data, num_sac_data, params.variable, 1);

[population_dir, num_sac_dir] = population_data_combine_levels(population_dir, num_sac_dir, params.variable, 2, [3 7]);
[population_dir, num_sac_dir] = population_data_combine_levels(population_dir, num_sac_dir, params.variable, 2, [2 8]);
[population_dir, num_sac_dir] = population_data_combine_levels(population_dir, num_sac_dir, params.variable, 2, [4 6]);
% [population_dir, num_sac_dir] = population_data_combine_tags(  population_dir, num_sac_dir, params.variable, [1]);

if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
    population_dir = population_data_smooth_pCells(population_dir, num_sac_dir, params.variable);
end
[population_avg_dir, population_sem_dir] = population_data_avg_over_pCells(population_dir, num_sac_dir, params.variable);
[population_allDir, num_sac_allDir] = population_data_avg_over_levels(population_dir, num_sac_dir, params.variable, 2);
[population_avg_allDir, population_sem_allDir] = population_data_avg_over_pCells(population_allDir, num_sac_allDir, params.variable);

data_avg_dir     = population_avg_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_dir     = population_sem_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_avg_allDir = population_avg_allDir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_allDir = population_sem_allDir.(params.variable)(params.tag_id).(params.event_type_name);

plot_population_data(params.fig_num, data_avg_dir, data_avg_allDir, params, data_sem_dir, data_sem_allDir);
end
%% MODE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.plot_mode == 2
[population_data, num_sac_data] = population_data_combine_levels(population_data, num_sac_data, params.variable, 2, [3 7]);
[population_data, num_sac_data] = population_data_combine_levels(population_data, num_sac_data, params.variable, 2, [2 8]);
[population_data, num_sac_data] = population_data_combine_levels(population_data, num_sac_data, params.variable, 2, [4 6]);
% [population_data, num_sac_data] = population_data_combine_tags(  population_data, num_sac_data, params.variable, [1]);

if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
    population_data = population_data_smooth_pCells(population_data, num_sac_data, params.variable);
end
[population_avg_dir, population_sem_dir] = population_data_avg_over_pCells(population_data, num_sac_data, params.variable);
[population_allDir, num_sac_allDir] = population_data_avg_over_levels(population_data, num_sac_data, params.variable, 2);
[population_avg_allDir, population_sem_allDir] = population_data_avg_over_pCells(population_allDir, num_sac_allDir, params.variable);

data_avg_dir     = population_avg_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_dir     = population_sem_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_avg_allDir = population_avg_allDir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_allDir = population_sem_allDir.(params.variable)(params.tag_id).(params.event_type_name);

plot_population_data(params.fig_num, data_avg_dir, data_avg_allDir, params, data_sem_dir, data_sem_allDir);
end

end
%% function plot_synch_ratio()
function plot_synch_ratio()
%% clear
clc; clear;
%% close all
close all;
%% set params
params.data_type       = 'SS';
params.CSYS_type       = 'absol'; % should be fixed to 'tuned', do not change to 'absol'
params.event_type_name = 'visual';
params.variable        = 'amp';
params.tag_id          = 1;
params.tags            = [1 4 6 7 8 11];
params.flag_smooth_plot = true;
params.fig_num = 4;
params.plot_mode = 1; % mode=1 collapse the amp/vel, mode=2 plots the amp/vel
params.ylim = [];
params.scale = 1;
params.win_len = '5ms';

%% Load data
if ~exist([params.data_type '_synch_joint_' params.CSYS_type], 'var')
    load([params.data_type '_synch_joint_' params.win_len '_' params.CSYS_type '.mat'], [params.data_type '_synch_joint_' params.CSYS_type])
end
if ~exist([params.data_type '_synch_margn_' params.CSYS_type], 'var')
    load([params.data_type '_synch_margn_' params.win_len '_' params.CSYS_type '.mat'], [params.data_type '_synch_margn_' params.CSYS_type])
end
if ~exist(['num_synch_' params.CSYS_type], 'var')
    load(['num_synch_' params.win_len '_' params.CSYS_type '.mat'], ['num_synch_' params.CSYS_type])
end

eval(['synch_joint = ' params.data_type '_synch_joint_' params.CSYS_type ';']);
eval(['synch_margn = ' params.data_type '_synch_margn_' params.CSYS_type ';']);
eval(['num_synch = ' 'num_synch_' params.CSYS_type ';']);

%% calc synch_ratio
[synch_joint_dir, num_joint_dir] = population_data_avg_over_levels(synch_joint, num_synch, params.variable, 1);
% [synch_joint_dir, num_joint_dir] = population_data_combine_levels(synch_joint_dir, num_joint_dir, params.variable, 2, [3 7]);
% [synch_joint_dir, num_joint_dir] = population_data_combine_levels(synch_joint_dir, num_joint_dir, params.variable, 2, [2 8]);
% [synch_joint_dir, num_joint_dir] = population_data_combine_levels(synch_joint_dir, num_joint_dir, params.variable, 2, [4 6]);
 [synch_joint_dir, num_joint_dir] = population_data_combine_tags(  synch_joint_dir, num_joint_dir, params.variable, params.tags);

[synch_margn_dir, num_margn_dir] = population_data_avg_over_levels(synch_margn, num_synch, params.variable, 1);
% [synch_margn_dir, num_margn_dir] = population_data_combine_levels(synch_margn_dir, num_margn_dir, params.variable, 2, [3 7]);
% [synch_margn_dir, num_margn_dir] = population_data_combine_levels(synch_margn_dir, num_margn_dir, params.variable, 2, [2 8]);
% [synch_margn_dir, num_margn_dir] = population_data_combine_levels(synch_margn_dir, num_margn_dir, params.variable, 2, [4 6]);
 [synch_margn_dir, num_margn_dir] = population_data_combine_tags(  synch_margn_dir, num_margn_dir, params.variable, params.tags);

[synch_joint_allDir, num_joint_allDir] = population_data_avg_over_levels(synch_joint_dir, num_joint_dir, params.variable, 2);
[synch_margn_allDir, num_margn_allDir] = population_data_avg_over_levels(synch_margn_dir, num_margn_dir, params.variable, 2);

synch_ratio_dir = synchrony_data_ratio(synch_joint_dir, synch_margn_dir, params.variable);
synch_ratio_allDir = synchrony_data_ratio(synch_joint_allDir, synch_margn_allDir, params.variable);

%% remove douplicates
% idx_pCells: is a boolean array. 1 for including a pCell, and 0 for exluding a pCell
num_pCells = size(synch_joint.amp(1).visual{1,1},1);
idx_pCells = true(num_pCells, 1);
idx_pCells(2:2:num_pCells, 1) = false;
[synch_ratio_dir, num_joint_dir] = population_data_idx_pCells(synch_ratio_dir, num_joint_dir, params.variable, idx_pCells);
[synch_ratio_allDir, num_joint_allDir] = population_data_idx_pCells(synch_ratio_allDir, num_joint_allDir, params.variable, idx_pCells);

%% smooth data
if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
    synch_ratio_dir = population_data_smooth_pCells(synch_ratio_dir, num_joint_dir, params.variable);
    synch_ratio_allDir = population_data_smooth_pCells(synch_ratio_allDir, num_joint_allDir, params.variable);
end

%% avg over pCells
[synch_ratio_avg_dir, synch_ratio_sem_dir] = population_data_avg_over_pCells(synch_ratio_dir, num_joint_dir, params.variable);
[synch_ratio_avg_allDir, synch_ratio_sem_allDir] = population_data_avg_over_pCells(synch_ratio_allDir, num_joint_allDir, params.variable);

%% Plot results PH
id = {'02x03', '02x04_2', '03x04_2'};

for counter = 1:3
figure(counter)
subplot(3,3,1)
hold on
eval(['plot(-249:250,synch_ratio_dir.amp(1).' params.event_type_name '{1, 4} (counter,:));']);
xline(0);xline(50);xline(-50);
title('135 dir.');
subplot(3,3,2)
hold on
eval(['plot(-249:250,synch_ratio_dir.amp(1).' params.event_type_name '{1, 3} (counter,:));']);
xline(0);xline(50);xline(-50);
title('90 dir.');
subplot(3,3,3)
hold on
eval(['plot(-249:250,synch_ratio_dir.amp(1).' params.event_type_name '{1, 2} (counter,:));']);
xline(0);xline(50);xline(-50);
title('45 dir.');
subplot(3,3,4)
hold on
eval(['plot(-249:250,synch_ratio_dir.amp(1).' params.event_type_name '{1, 5} (counter,:));']);
xline(0);xline(50);xline(-50);
ylabel('Synchrony index')
title('180 dir.');
subplot(3,3,5)
hold on
eval(['plot(-249:250,synch_ratio_allDir.amp(1).' params.event_type_name '{1, 1} (counter,:));']);
xline(0);xline(50);xline(-50);
title('all dir.');
subplot(3,3,6)
hold on
eval(['plot(-249:250,synch_ratio_dir.amp(1).' params.event_type_name '{1, 1} (counter,:));']);
xline(0);xline(50);xline(-50);
title('0 dir.');
subplot(3,3,7)
hold on
eval(['plot(-249:250,synch_ratio_dir.amp(1).' params.event_type_name '{1, 6} (counter,:));']);
xline(0);xline(50);xline(-50);
title('225 dir.');
subplot(3,3,8)
hold on
eval(['plot(-249:250,synch_ratio_dir.amp(1).' params.event_type_name '{1, 7} (counter,:));']);
xline(0);xline(50);xline(-50);
xlabel(['Time from ' params.event_type_name ' (ms)'])
title('270 dir.');
subplot(3,3,9)
hold on
eval(['plot(-249:250,synch_ratio_dir.amp(1).' params.event_type_name '{1, 8} (counter,:));']);
xline(0);xline(50);xline(-50);
title('315 dir.');
sgtitle(['Synchrony: ' char(id(counter)) '| Tags: ' num2str(params.tags)])
ESN_Beautify_Plot
saveas(gcf, ['Synchrony_ ' params.data_type '_' params.event_type_name '_' char(id(counter)) '_tags' num2str(params.tags)], 'pdf' )
end



%% plot results
data_avg_dir     = synch_ratio_avg_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_dir     = synch_ratio_sem_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_avg_allDir = synch_ratio_avg_allDir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_allDir = synch_ratio_sem_allDir.(params.variable)(params.tag_id).(params.event_type_name);
% plot_population_data(params.fig_num, data_avg_dir, data_avg_allDir, params, data_sem_dir, data_sem_allDir);

figure(params.fig_num)
subplot(3,3,1)
hold on
plot(-249:250, data_avg_dir{1, 4}  , 'k', 'LineWidth', 1)
plot(-249:250, data_avg_dir{1, 4} + data_sem_dir{1, 4}  , 'k', 'LineWidth', 0.5)
plot(-249:250, data_avg_dir{1, 4} - data_sem_dir{1, 4}  , 'k', 'LineWidth', 0.5)
xline(0);xline(50);xline(-50);
title('135 dir.');
subplot(3,3,2)
hold on
plot(-249:250, data_avg_dir{1, 3}  , 'k', 'LineWidth', 1)
plot(-249:250, data_avg_dir{1, 3} + data_sem_dir{1, 3}  , 'k', 'LineWidth', 0.5)
plot(-249:250, data_avg_dir{1, 3} - data_sem_dir{1, 3}  , 'k', 'LineWidth', 0.5)
xline(0);xline(50);xline(-50);
title('90 dir.');
subplot(3,3,3)
hold on
plot(-249:250, data_avg_dir{1, 2}  , 'k', 'LineWidth', 1)
plot(-249:250, data_avg_dir{1, 2} + data_sem_dir{1, 2}  , 'k', 'LineWidth', 0.5)
plot(-249:250, data_avg_dir{1, 2} - data_sem_dir{1, 2}  , 'k', 'LineWidth', 0.5)
xline(0);xline(50);xline(-50);
title('45 dir.');
subplot(3,3,4)
hold on
plot(-249:250, data_avg_dir{1, 5}  , 'k', 'LineWidth', 1)
plot(-249:250, data_avg_dir{1, 5} + data_sem_dir{1, 5}  , 'k', 'LineWidth', 0.5)
plot(-249:250, data_avg_dir{1, 5} - data_sem_dir{1, 5}  , 'k', 'LineWidth', 0.5)
xline(0);xline(50);xline(-50);
ylabel('Synchrony index')
title('180 dir.');
subplot(3,3,5)
hold on
plot(-249:250, data_avg_dir{1, 1}  , 'k', 'LineWidth', 1)
plot(-249:250, data_avg_dir{1, 1} + data_sem_dir{1, 5}  , 'k', 'LineWidth', 0.5)
plot(-249:250, data_avg_dir{1, 1} - data_sem_dir{1, 5}  , 'k', 'LineWidth', 0.5)
xline(0);xline(50);xline(-50);
title('all dir.');
subplot(3,3,6)
hold on
plot(-249:250, data_avg_allDir{1, 1}  , 'k', 'LineWidth', 1)
plot(-249:250, data_avg_allDir{1, 1} + data_sem_allDir{1, 1}  , 'k', 'LineWidth', 0.5)
plot(-249:250, data_avg_allDir{1, 1} - data_sem_allDir{1, 1}  , 'k', 'LineWidth', 0.5)
xline(0);xline(50);xline(-50);
title('0 dir.');
subplot(3,3,7)
hold on
plot(-249:250, data_avg_dir{1, 6}  , 'k', 'LineWidth', 1)
plot(-249:250, data_avg_dir{1, 6} + data_sem_dir{1, 6}  , 'k', 'LineWidth', 0.5)
plot(-249:250, data_avg_dir{1, 6} - data_sem_dir{1, 6}  , 'k', 'LineWidth', 0.5)
xline(0);xline(50);xline(-50);
title('225 dir.');
subplot(3,3,8)
hold on
plot(-249:250, data_avg_dir{1, 7}  , 'k', 'LineWidth', 1)
plot(-249:250, data_avg_dir{1, 7} + data_sem_dir{1, 7}  , 'k', 'LineWidth', 0.5)
plot(-249:250, data_avg_dir{1, 7} - data_sem_dir{1, 7}  , 'k', 'LineWidth', 0.5)
xline(0);xline(50);xline(-50);
xlabel(['Time from ' params.event_type_name ' (ms)'])
title('270 dir.');
subplot(3,3,9)
hold on
plot(-249:250, data_avg_dir{1, 8}  , 'k', 'LineWidth', 1)
plot(-249:250, data_avg_dir{1, 8} + data_sem_dir{1, 8}  , 'k', 'LineWidth', 0.5)
plot(-249:250, data_avg_dir{1, 8} - data_sem_dir{1, 8}  , 'k', 'LineWidth', 0.5)
xline(0);xline(50);xline(-50);
title('315 dir.');
sgtitle(['Synchrony (averaged) | Tags: ' num2str(params.tags)])
ESN_Beautify_Plot


saveas(gcf, ['Synchrony_Averaged_ ' params.data_type '_' params.event_type_name '_tags' num2str(params.tags)], 'pdf' )
 
end
%% function plot_synch_MI()
function plot_synch_MI()
%% clear
clc; clear;
%% close all
close all;
%% set params
params.data_type       = 'SS';
params.CSYS_type       = 'absol';
params.event_type_name = 'onset';
params.variable        = 'amp';
params.tag_id          = 1;
params.flag_smooth_plot = true;
params.fig_num = 3;
params.plot_mode = 1; % mode=1 collapse the amp/vel, mode=2 plots the amp/vel
params.ylim = [];
params.scale = 1;
params.win_len = '1ms';

%% Load data
if ~exist([params.data_type '_synch_joint_' params.CSYS_type], 'var')
    load([params.data_type '_synch_joint_' params.win_len '_'  params.CSYS_type '.mat'], [params.data_type '_synch_joint_' params.CSYS_type])
end
if ~exist([params.data_type '_synch_cross_' params.CSYS_type], 'var')
    load([params.data_type '_synch_cross_' params.win_len '_'  params.CSYS_type '.mat'], [params.data_type '_synch_cross_' params.CSYS_type])
end
if ~exist([params.data_type '_synch_margn_' params.CSYS_type], 'var')
    load([params.data_type '_synch_margn_' params.win_len '_'  params.CSYS_type '.mat'], [params.data_type '_synch_margn_' params.CSYS_type])
end
if ~exist(['num_synch_' params.CSYS_type], 'var')
    load(['num_synch_' params.win_len '_'  params.CSYS_type '.mat'], ['num_synch_' params.CSYS_type])
end 

eval(['synch_joint = ' params.data_type '_synch_joint_' params.CSYS_type ';']);
eval(['synch_cross = ' params.data_type '_synch_cross_' params.CSYS_type ';']);
eval(['synch_margn = ' params.data_type '_synch_margn_' params.CSYS_type ';']);
eval(['num_synch = ' 'num_synch_' params.CSYS_type ';']);

%% calc synch_ratio
[synch_joint_dir, num_joint_dir] = population_data_avg_over_levels(synch_joint, num_synch, params.variable, 1);
[synch_joint_dir, num_joint_dir] = population_data_combine_levels(synch_joint_dir, num_joint_dir, params.variable, 2, [3 7]);
[synch_joint_dir, num_joint_dir] = population_data_combine_levels(synch_joint_dir, num_joint_dir, params.variable, 2, [2 8]);
[synch_joint_dir, num_joint_dir] = population_data_combine_levels(synch_joint_dir, num_joint_dir, params.variable, 2, [4 6]);
[synch_joint_dir, num_joint_dir] = population_data_combine_tags(  synch_joint_dir, num_joint_dir, params.variable, [1 4 6 7 8]);

[synch_cross_dir, num_cross_dir] = population_data_avg_over_levels(synch_cross, num_synch, params.variable, 1);
[synch_cross_dir, num_cross_dir] = population_data_combine_levels(synch_cross_dir, num_cross_dir, params.variable, 2, [3 7]);
[synch_cross_dir, num_cross_dir] = population_data_combine_levels(synch_cross_dir, num_cross_dir, params.variable, 2, [2 8]);
[synch_cross_dir, num_cross_dir] = population_data_combine_levels(synch_cross_dir, num_cross_dir, params.variable, 2, [4 6]);
[synch_cross_dir, num_cross_dir] = population_data_combine_tags(  synch_cross_dir, num_cross_dir, params.variable, [1 4 6 7 8]);

[synch_margn_dir, num_margn_dir] = population_data_avg_over_levels(synch_margn, num_synch, params.variable, 1);
[synch_margn_dir, num_margn_dir] = population_data_combine_levels(synch_margn_dir, num_margn_dir, params.variable, 2, [3 7]);
[synch_margn_dir, num_margn_dir] = population_data_combine_levels(synch_margn_dir, num_margn_dir, params.variable, 2, [2 8]);
[synch_margn_dir, num_margn_dir] = population_data_combine_levels(synch_margn_dir, num_margn_dir, params.variable, 2, [4 6]);
[synch_margn_dir, num_margn_dir] = population_data_combine_tags(  synch_margn_dir, num_margn_dir, params.variable, [1 4 6 7 8]);

[synch_joint_allDir, num_joint_allDir] = population_data_avg_over_levels(synch_joint_dir, num_joint_dir, params.variable, 2);
[synch_cross_allDir, num_cross_allDir] = population_data_avg_over_levels(synch_cross_dir, num_cross_dir, params.variable, 2);
[synch_margn_allDir, num_margn_allDir] = population_data_avg_over_levels(synch_margn_dir, num_margn_dir, params.variable, 2);

synch_MI_dir = synchrony_mutual_info(synch_joint_dir, synch_cross_dir, synch_margn_dir, params.variable);
synch_MI_allDir = synchrony_mutual_info(synch_joint_allDir, synch_cross_allDir, synch_margn_allDir, params.variable);

%% remove douplicates
% idx_pCells: is a boolean array. 1 for including a pCell, and 0 for exluding a pCell
num_pCells = size(synch_joint.amp(1).visual{1,1},1);
idx_pCells = true(num_pCells, 1);
idx_pCells(2:2:num_pCells, 1) = false;
[synch_MI_dir, num_joint_dir] = population_data_idx_pCells(synch_MI_dir, num_joint_dir, params.variable, idx_pCells);
[synch_MI_allDir, num_joint_allDir] = population_data_idx_pCells(synch_MI_allDir, num_joint_allDir, params.variable, idx_pCells);

%% smooth data
if strcmp(params.data_type, 'SS') || strcmp(params.data_type, 'CS')
    synch_MI_dir = population_data_smooth_pCells(synch_MI_dir, num_joint_dir, params.variable);
    synch_MI_allDir = population_data_smooth_pCells(synch_MI_allDir, num_joint_allDir, params.variable);
end

%% avg over pCells
[synch_MI_avg_dir, synch_MI_sem_dir] = population_data_avg_over_pCells(synch_MI_dir, num_joint_dir, params.variable);
[synch_MI_avg_allDir, synch_MI_sem_allDir] = population_data_avg_over_pCells(synch_MI_allDir, num_joint_allDir, params.variable);

%% plot results
data_avg_dir     = synch_MI_avg_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_dir     = synch_MI_sem_dir.(params.variable)(params.tag_id).(params.event_type_name);
data_avg_allDir = synch_MI_avg_allDir.(params.variable)(params.tag_id).(params.event_type_name);
data_sem_allDir = synch_MI_sem_allDir.(params.variable)(params.tag_id).(params.event_type_name);

plot_population_data(params.fig_num, data_avg_dir, data_avg_allDir, params, data_sem_dir, data_sem_allDir);

end
%% function plot_single_session_modulation
function plot_single_session_modulation(pCell_id)
%% clear
% clc; clear;
%% close all
% close all;
%% set params
flag_pair_list = false;
JSP_global_variables(flag_pair_list);
global ang_edges ang_values length_trace 
% params.pCell_id = '210507_174229_03_combine_1';
params.pCell_id = pCell_id;
params.tag_ids  = [1 6];
params.fig_num = 5;
fprintf([pCell_id '\n'])
fprintf('params --> Completed. \n')

%% Load data
% path_cell_data = uigetdir;
path_cell_data = pwd;
if ~strcmp(path_cell_data(end), filesep);path_cell_data = [path_cell_data filesep];end
load([path_cell_data params.pCell_id '.mat'], 'SACS_ALL_DATA', 'CS_on_data')
pCell_ids = build_pCell_ids(flag_pair_list);
counter_pCell = find(strcmp(pCell_ids, params.pCell_id));
fprintf('Load data --> Completed. \n')
%% Build sac_data
idx_sacs    = false(size(SACS_ALL_DATA.validity));
for counter_tag = 1 : length(params.tag_ids)
    idx_sacs    = idx_sacs | (SACS_ALL_DATA.tag == params.tag_ids(counter_tag));
end
idx_sacs    = idx_sacs & SACS_ALL_DATA.validity;

sac_data.time_onset  = SACS_ALL_DATA.time_onset(     :,idx_sacs);
sac_data.time_visual = SACS_ALL_DATA.time_visual(    :,idx_sacs);
sac_data.time_offset = SACS_ALL_DATA.time_offset(    :,idx_sacs);
sac_data.SS_onset    = SACS_ALL_DATA.neuro_SS_onset( :,idx_sacs);
sac_data.CS_onset    = SACS_ALL_DATA.neuro_CS_onset( :,idx_sacs);
sac_data.SS_visual   = SACS_ALL_DATA.neuro_SS_visual(:,idx_sacs);
sac_data.CS_visual   = SACS_ALL_DATA.neuro_CS_visual(:,idx_sacs);
sac_data.eye_vx_onset  = SACS_ALL_DATA.eye_vx_onset( :,idx_sacs);
sac_data.eye_vy_onset  = SACS_ALL_DATA.eye_vy_onset( :,idx_sacs);
sac_data.eye_vx_visual = SACS_ALL_DATA.eye_vx_visual(:,idx_sacs);
sac_data.eye_vy_visual = SACS_ALL_DATA.eye_vy_visual(:,idx_sacs);
sac_data.visual_px_offset = SACS_ALL_DATA.visual_px_offset(:,idx_sacs);
sac_data.visual_py_offset = SACS_ALL_DATA.visual_py_offset(:,idx_sacs);
sac_data.eye_r_px_onset   = SACS_ALL_DATA.eye_r_px_onset(  :,idx_sacs);
sac_data.eye_r_py_onset   = SACS_ALL_DATA.eye_r_py_onset(  :,idx_sacs);

sac_data.eye_vm_onset  = sqrt(sac_data.eye_vx_onset.^2  + sac_data.eye_vy_onset.^2 );
sac_data.eye_vm_visual = sqrt(sac_data.eye_vx_visual.^2 + sac_data.eye_vy_visual.^2);
sac_data.time_diff_visual_onset = sac_data.time_onset - sac_data.time_visual;
sac_data.time_diff_onset_offset = sac_data.time_offset - sac_data.time_onset;
%% Build sac_data_dir
sac_data.delta_x = sac_data.visual_px_offset - sac_data.eye_r_px_onset;
sac_data.delta_y = sac_data.visual_py_offset - sac_data.eye_r_py_onset;
sac_data.visual_ang = wrapTo360(atan2d(sac_data.delta_y, sac_data.delta_x));
sac_data.visual_ang_bin = discretize(sac_data.visual_ang, ang_edges);
last_bin_id = length(ang_edges) - 1;
sac_data.visual_ang_bin(sac_data.visual_ang_bin == last_bin_id) = 1; % wrap the circle around
% 1: 0deg % 2: 45deg % 3: 90deg % 4: 135deg % 5: 180deg % 6: 225deg % 7: 270deg % 8: 315deg

if length(ang_values) ~= 8
    error('plot_single_session_modulation: length ang_values is not 8. Please modify the code.')
end
if length_trace ~= 500
    error('sac_modulation_index: length_trace is not 500. Please modify the code.')
end

field_names_sac_data = fieldnames(sac_data);
sac_data_dir = struct;
for counter_dir = 1 : 8
for counter_field = 1 : length(field_names_sac_data)
    field_name = field_names_sac_data{counter_field};
    idx_ang = (sac_data.visual_ang_bin == counter_dir);
    sac_data_dir(counter_dir).(field_name) = sac_data.(field_name)(:,idx_ang);
end
end
%% Init Plot
hFig = figure(params.fig_num);
clf(hFig)
num_row_fig = 9;
num_col_fig = 9;
ax_ang_id = [6, 3, 2, 1, 4, 7, 8, 9];

% Figure parameters
Line_Color = lines(7);
color_SS     = Line_Color(1,:);
color_CS     = Line_Color(7,:);
color_visual = Line_Color(2,:);
color_onset  = Line_Color(2,:);
color_offset = Line_Color(2,:);
color_SS_firing = [0    0.3    0.5];
color_vm     = Line_Color(5,:);
clearvars h_ax

num_trial_dir = zeros(1,8);
for counter_dir = 1 : 8
    num_trial_dir(counter_dir) = size(sac_data_dir(counter_dir).time_onset,2);
end
num_trial_dir_max = round(median(num_trial_dir));
%% Plot data, Loop over dirs
for counter_dir = 1 : 8
    if isempty(sac_data_dir(counter_dir).SS_onset)
        continue;
    end
    ax_ang_id_ = ax_ang_id(counter_dir);
    row_ = ceil(ax_ang_id_ / 3);
    col_ = mod(ax_ang_id_ , 3);
    if col_==0
        col_ = 3;
    end
    rows_ = ((row_-1)*3) + [1 2 3]';
    cols_ = ((col_-1)*3) + [1 2 3];
    rows_ = repmat(rows_,1,3);
    cols_ = repmat(cols_,3,1);
    ax_ids = ( (rows_-1)*9 ) + cols_;
    num_trial_dir_max_ = min([num_trial_dir_max size(sac_data_dir(counter_dir).time_onset,2)]);
    
    %% Plot visual raster
    h_ax(ax_ang_id_, 1) = subplot(num_row_fig, num_col_fig, [ax_ids(1,1) ax_ids(2,1)]);
    hold on
    train_data_logic_SS_ = sac_data_dir(counter_dir).SS_visual(250:400,1:num_trial_dir_max_)';
    firing_SS_ = mean(sac_data_dir(counter_dir).SS_visual(:,1:num_trial_dir_max_)') * 1000;
    firing_SS_ = ESN_smooth(firing_SS_); % smooth(firing_SS_(:), 21, 'sgolay', 2)'; % 
    firing_SS_(firing_SS_<0)=0;
    if length(firing_SS_) == 1
        firing_SS_ = zeros(1,length(sac_data_dir(counter_dir).SS_visual(:,1:num_trial_dir_max_)));
    end
    firing_SS_ = firing_SS_(1, 250:400);
    train_data_logic_CS_ = sac_data_dir(counter_dir).CS_visual(250:400,1:num_trial_dir_max_)';
    train_data_logic_visual = false(size(train_data_logic_SS_));
    train_data_logic_visual(:,1) = true;
    train_data_logic_onset = false(size(train_data_logic_SS_));
    ind_onset = round(sac_data_dir(counter_dir).time_diff_visual_onset(:,1:num_trial_dir_max_) * 1000); % convert sec to ms
    ind_onset_out_of_range = false(size(ind_onset));
    ind_onset_out_of_range(ind_onset>150) = true; ind_onset_out_of_range(ind_onset<1) = true;
    onset_row_number_ = 1:size(train_data_logic_onset,1);
    onset_row_number_ = onset_row_number_(~ind_onset_out_of_range);
    ind_onset = ind_onset(~ind_onset_out_of_range);
    for counter_ind = 1 : sum(~ind_onset_out_of_range)
        onset_row_number__ = onset_row_number_(counter_ind);
        onset_col_number__ = ind_onset(counter_ind);
        train_data_logic_onset(onset_row_number__,onset_col_number__) = true;
    end

    inds_span = 0 : 150;
    [x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
    plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
%     [x_axis_visual_, y_axis_visual_] = ESN_raster_plot_axes(train_data_logic_visual, inds_span, 0.1);
%     plot(x_axis_visual_(:), y_axis_visual_(:), 'LineWidth', 2, 'Color', color_visual);
%     [x_axis_onset_, y_axis_onset_] = ESN_raster_plot_axes(train_data_logic_onset, inds_span, 0.1);
%     plot(x_axis_onset_(:), y_axis_onset_(:), 'LineWidth', 2, 'Color', color_onset);
    xlim([0 150])
    if size(train_data_logic_SS_,1) > 1
    ylim([1 size(train_data_logic_SS_,1)])
    else
    ylim([0 1])
    end
    set(gca, 'XTick', [0 50 100 150])
    set(gca, 'XTickLabel', {'',''})
    if (ax_ang_id_==1) || (ax_ang_id_==4) || (ax_ang_id_==7)
    ylabel('Saccade')
    end
    yyaxis right;
    plot(inds_span, firing_SS_, 'LineWidth', 1, 'Color', color_SS_firing)
    ylim([0 200])
    xlim([0 150])
    set(gca, 'YTick', [0 50 100 150 200])
    set(gca, 'yTickLabel', {'','','',''})
    set(gca, 'YColor', color_SS_firing)
%     if (ax_ang_id_==1) || (ax_ang_id_==4) || (ax_ang_id_==7)
%          xlabel('Time from visual (ms)')
%     end
    
    %% Plot onset raster
    h_ax(ax_ang_id_, 2) = subplot(num_row_fig, num_col_fig, [ax_ids(1,2) ax_ids(1,3) ax_ids(2,2) ax_ids(2,3)]);
    hold on
    train_data_logic_SS_ = sac_data_dir(counter_dir).SS_onset(200:350,1:num_trial_dir_max_)';
    firing_SS_ = mean(sac_data_dir(counter_dir).SS_onset(:,1:num_trial_dir_max_)') * 1000;
    firing_SS_ = ESN_smooth(firing_SS_); % smooth(firing_SS_(:), 21, 'sgolay', 2)'; % 
    firing_SS_(firing_SS_<0)=0;
    if length(firing_SS_) == 1
        firing_SS_ = zeros(1,length(sac_data_dir(counter_dir).SS_onset(:,1:num_trial_dir_max_)));
    end
    firing_SS_ = firing_SS_(1, 200:350);
    train_data_logic_CS_ = sac_data_dir(counter_dir).CS_onset(200:350,1:num_trial_dir_max_)';
    train_data_logic_onset = false(size(train_data_logic_SS_));
    train_data_logic_onset(:,51) = true;
    train_data_logic_offset = false(size(train_data_logic_SS_));
    ind_offset = round(sac_data_dir(counter_dir).time_diff_onset_offset(:,1:num_trial_dir_max_) * 1000)+51; % convert sec to ms
    ind_offset_out_of_range = false(size(ind_offset));
    ind_offset_out_of_range(ind_offset>100) = true; ind_offset_out_of_range(ind_offset<1) = true;
    offset_row_number_ = 1:size(train_data_logic_offset,1);
    offset_row_number_ = offset_row_number_(~ind_offset_out_of_range);
    ind_offset = ind_offset(~ind_offset_out_of_range);
    for counter_ind = 1 : sum(~ind_offset_out_of_range)
        offset_row_number__ = offset_row_number_(counter_ind);
        offset_col_number__ = ind_offset(counter_ind);
        train_data_logic_offset(offset_row_number__,offset_col_number__) = true;
    end
    
    inds_span = -50 : 100;
    [x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
    plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
    plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
    [x_axis_onset_, y_axis_onset_] = ESN_raster_plot_axes(train_data_logic_onset, inds_span, 0.5);
    plot(x_axis_onset_(:), y_axis_onset_(:), 'LineWidth', 2, 'Color', color_onset);
    [x_axis_offset_, y_axis_offset_] = ESN_raster_plot_axes(train_data_logic_offset, inds_span, 0.5);
    plot(x_axis_offset_(:), y_axis_offset_(:), 'LineWidth', 2, 'Color', color_offset);
    xlim([-50 100])
    if size(train_data_logic_SS_,1) > 1
    ylim([1 size(train_data_logic_SS_,1)])
    else
    ylim([0 1])
    end
    set(gca, 'XTick', [-50 0 50 100])
    set(gca, 'XTickLabel', {'','','', ''})
    set(gca, 'YTickLabel', [])
    
    yyaxis right;
    
    plot(inds_span, firing_SS_, 'LineWidth', 1, 'Color', color_SS_firing)
    ylim([0 200])
    xlim([-50 100])
    set(gca, 'YColor', color_SS_firing)
    if (ax_ang_id_==3) || (ax_ang_id_==6) || (ax_ang_id_==9)
        ylabel('SS Firing (spk/s)')
        set(gca, 'YTick', [0 50 100 150 200])
        set(gca, 'YTickLabel', {'0','','100', '', '200'})
    else
        set(gca, 'YTick', [0 50 100 150 200])
        set(gca, 'YTickLabel', {'','','', '', ''})
    end
       
    %% Plot cue vm
    h_ax(ax_ang_id_, 3) = subplot(num_row_fig, num_col_fig, ax_ids(3,1));
    hold on
% %     train_data_logic_CS_ = sac_data_dir(counter_dir).CS_visual(250:400,1:num_trial_dir_max_)';
%     firing_CS_ = mean(sac_data_dir(counter_dir).CS_visual(:,1:num_trial_dir_max_)') * 1000;
%     firing_CS_ = ESN_smooth(firing_CS_); % smooth(firing_SS_(:), 21, 'sgolay', 2)'; % 
%     firing_CS_(firing_CS_<0)=0;
%     firing_CS_ = firing_CS_(1, 250:400);
%     inds_span = 0 : 150;
%     plot(inds_span, firing_CS_, '-', 'LineWidth', 0.5, 'color', color_CS)
%     ylim([0 3])
%     xlim([0 150])
    vm_ = sac_data_dir(counter_dir).eye_vm_visual(250:400,1:num_trial_dir_max_)';
    inds_span = 0 : 150;
    vm_ = nanmean(vm_);
    plot(inds_span, vm_, '-k', 'LineWidth', 0.5)
    ylim([0 600])
    xlim([0 150])
    set(gca, 'XTick', [0 50 100 150])
    set(gca, 'YTick', [0 300 600])
    set(gca, 'XTickLabel', {'0', '', '', '150'})
    
    if (ax_ang_id_==1) || (ax_ang_id_==4) || (ax_ang_id_==7)
        ylabel('Sac. vel.')
    end
    
    if (ax_ang_id_==7)
         xlabel('Time from visual (ms)')
    end
    
    %% Plot onset vm
    h_ax(ax_ang_id_, 4) = subplot(num_row_fig, num_col_fig, [ax_ids(3,2) ax_ids(3,3)]);
    hold on
    vm_ = sac_data_dir(counter_dir).eye_vm_onset(200:350,1:num_trial_dir_max_)';
    inds_span = -50 : 100;
    vm_ = nanmean(vm_);
%     firing_CS_ = mean(sac_data_dir(counter_dir).CS_visual(:,1:num_trial_dir_max_)') * 1000;
%     firing_CS_ = ESN_smooth(firing_CS_); % smooth(firing_SS_(:), 21, 'sgolay', 2)'; % 
%     firing_CS_(firing_CS_<0)=0;
%     firing_CS_ = firing_CS_(1, 200:350);
    
    plot(inds_span, vm_, '-k', 'LineWidth', 0.5)
    xlim([-50 100])
    ylim([0 600])
    set(gca, 'XTick', [-50 0 50 100])
    set(gca, 'YTick', [0 300 600])
    set(gca, 'XTickLabel', {'-50' ,'0','', '100'})
    set(gca, 'YTickLabel', {'' ,'',''})
    
%     yyaxis right;
%     plot(inds_span, firing_CS_, '-', 'LineWidth', 0.5, 'color', color_CS)
%     ylim([0 3])
%     xlim([-50 100])
%     set(gca, 'YColor', color_CS)
        if (ax_ang_id_==9)
         xlabel('Time from onset (ms)')
    end
end
% % Plot CS-Tuning
% ax_ang_id_ = 5;
% row_ = ceil(ax_ang_id_ / 3);
% col_ = mod(ax_ang_id_ , 3);
% if col_==0
%     col_ = 3;
% end
% rows_ = ((row_-1)*3) + [1 2 3]';
% cols_ = ((col_-1)*3) + [1 2 3];
% rows_ = repmat(rows_,1,3);
% cols_ = repmat(cols_,3,1);
% ax_ids = ( (rows_-1)*9 ) + cols_;
% h_ax(ax_ang_id_, 1) = subplot(num_row_fig, num_col_fig, ax_ids(:));
% prob_amplitude = CS_on_data.CS_prob_avg;
% 
% vonMises_std = CS_on_data.vonMises_std;
% CS_ang_avg = CS_on_data.CS_ang_avg;
% CS_rho_avg = CS_on_data.CS_rho_avg;
% std_curv_ang = (CS_on_data.CS_ang_avg-vonMises_std) : 2 : (CS_on_data.CS_ang_avg+vonMises_std);
% std_curv_amp = repmat(CS_rho_avg, length(std_curv_ang), 1);
% 
% plot_data_amp_mean = [prob_amplitude, prob_amplitude(1), nan]';
% plot_data_deg_mean = [ang_values, ang_values(1), nan]';
% polarplot(deg2rad(plot_data_deg_mean),plot_data_amp_mean, '-', 'LineWidth', 1, 'Color', color_CS)
% hold on
% polarplot(deg2rad(std_curv_ang), std_curv_amp, '-', 'LineWidth', 1.5, 'Color', color_CS)
% polarplot([0 deg2rad(CS_ang_avg)],[0 CS_rho_avg], '-', 'LineWidth', 1.5, 'Color', color_CS)
% rlim([0 0.25])
% set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:0.05:0.25, ...
%     'RTickLabel', {'', '', '0.1', '', '0.2', ''}, 'ThetaTickLabel', {'0','','90','', '180','','270', ''})
% title('CS Tuning', 'Interpreter', 'none', 'Color',color_CS);

sgtitle([params.pCell_id ' | Tags: ' num2str(params.tag_ids) ]  , 'Interpreter', 'none')
ESN_Beautify_Plot(hFig, [11 8], 8)
% ESN_Beautify_Plot(hFig, [8 8], 12)
saveas(hFig,[params.pCell_id '_tags' num2str(params.tag_ids) ] , 'pdf');

end

%% SCRATCH AREA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function scratch_permutation_ss_population()
function scratch_permutation_ss_population()
%%
clc;
JSP_global_variables();
global length_trace
num_iterations = 1000;
path_cell_data = [pwd filesep ];
params.CSYS_type = 'tuned';
params.data_type = 'SS';
params.variable  = 'amp';
variable = params.variable;
counter_tag = 1;
event_type_name = 'onset';
counter_var = 1;
counter_ang = 5;
if ~exist(['num_sac_' params.CSYS_type], 'var')
    load([path_cell_data 'num_sac_' params.CSYS_type '.mat'], ['num_sac_' params.CSYS_type])
end
if ~exist([params.data_type '_population_' params.CSYS_type], 'var')
    load([path_cell_data params.data_type '_population_' params.CSYS_type '.mat'], [params.data_type '_population_' params.CSYS_type])
end
if ~exist('population_neural_properties', 'var')
    load([path_cell_data 'population_neural_properties' '.mat'], 'population_neural_properties')
end
eval(['firing_rate = population_neural_properties.' params.data_type '_firing_rate'  ';']);
eval(['num_sac_data = ' 'num_sac_' params.CSYS_type ';']);
eval(['population_data = ' params.data_type '_population_' params.CSYS_type ';']);
eval(['firing_rate = population_neural_properties.' params.data_type '_firing_rate'  ';']);
num_pCells = size(population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang},1);
population_data = population_data_subtract_baseline(population_data, firing_rate, params.variable);
[population_avg_levels, num_sac_data_avg] = population_data_avg_over_levels(population_data, num_sac_data, params.variable, 1);

%%
population_data = population_avg_levels;
num_sac_data = num_sac_data_avg;

event_data_raw = ...
    population_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};

baseline_ = nanmean(event_data_raw(:,50:150),2);
event_data_raw = event_data_raw - repmat(baseline_, 1, length_trace);

event_data_raw(isnan(event_data_raw)) = 0;
event_data_ = ESN_smooth(event_data_raw, 2);

num_sac_ = ...
    num_sac_data.(variable)(counter_tag).(event_type_name){counter_var, counter_ang};
num_sac_matrix_ = repmat(num_sac_, 1, length_trace);

% SEM based on bootstrapping on pCells
%
num_pCells_perm_1 = num_pCells;
event_data_perm_ = nan(num_iterations, length_trace);
for counter_iteration = 1 : num_iterations
    idx_iteration_ = randi(num_pCells, 1, num_pCells_perm_1);
    event_data_iteration_ = event_data_(idx_iteration_, :);
    num_sac_iteration_    = num_sac_matrix_(   idx_iteration_, :);
    avg_data_iteration_ = nansum(event_data_iteration_ .* num_sac_iteration_) ./ nansum(num_sac_iteration_);
    event_data_perm_(counter_iteration, :) = avg_data_iteration_;
end
%}
avg_perm_pCell_1 = nanmean(event_data_perm_);
std_perm_pCell_1 = nanstd(event_data_perm_);

trace_avg_weights = num_sac_matrix_ ./ repmat(nansum(num_sac_matrix_), num_pCells, 1);
avg_traces_sac_weight = nansum(event_data_ .* trace_avg_weights);
% sem_traces = nanstd(event_data_)./ sqrt(num_pCells);
sem_traces = sqrt(var(event_data_,trace_avg_weights(:,1),1,'omitnan'))./ sqrt(num_pCells);

%
hFig_ = figure(1);
clf(hFig_)
subplot(2,1,1)
hold on
% plot(avg_raw, 'linewidth', 1)
plot(avg_traces_sac_weight, 'linewidth', 1)
plot(avg_perm_pCell_1, 'linewidth', 1)
legend({...
    'avg_traces',...
    ['avg_perm_pCell_' num2str(num_pCells_perm_1)],...
    }, ...
    'interpreter', 'none', 'location', 'eastoutside')

subplot(2,1,2)
hold on
% plot(sem_raw, 'linewidth', 1)
plot(sem_traces, 'linewidth', 1)
plot(std_perm_pCell_1, 'linewidth', 1)

legend({...
    'sem_traces',...
    ['std_perm_pCell_' num2str(num_pCells_perm_1)],...
    },...
    'interpreter', 'none', 'location', 'eastoutside')

ESN_Beautify_Plot(hFig_, [8, 5], 12)
%%
end
%% function scratch_plot_vmax_vs_ss_modulation()
function scratch_plot_vmax_vs_ss_modulation()
%% for ploting vmax vs. ss_max_modulation: set ss_max
ss_max_avg = nan(8,8);
ss_max_sem = nan(8,8);
ss_min_avg = nan(8,8);
ss_min_sem = nan(8,8);

for counter_ang = 1 : 8
    for counter_var = 1 : 8
        [max_, idx_] = max(data_avg_dir{counter_var, counter_ang});
        ss_max_avg(counter_var, counter_ang) = max_;
        ss_max_sem(counter_var, counter_ang) = data_sem_dir{counter_var, counter_ang}(1,idx_);
        
        [min_, idx_] = min(data_avg_dir{counter_var, counter_ang});
        ss_min_avg(counter_var, counter_ang) = min_;
        ss_min_sem(counter_var, counter_ang) = data_sem_dir{counter_var, counter_ang}(1,idx_);
    end
end

%% for ploting vmax vs. ss_crossing: set ss_crossing
ss_crossing_avg   = nan(8,8);
ss_crossing_sem_p = nan(8,8);
ss_crossing_sem_m = nan(8,8);

for counter_ang = 5
    for counter_var = 1 : 8
        trace_ = data_avg_dir{counter_var, counter_ang};
        [~, idx_max_] = max(trace_);
        [~, idx_min_] = min(trace_(idx_max_:end));
        trace_max_min_ = trace_(idx_max_:idx_max_+idx_min_);
        ss_crossing_avg(counter_var, counter_ang) = find(trace_max_min_(1:end-1) >= 0 & trace_max_min_(2:end) < 0, 1, 'first') + idx_max_ - 250;
        
        trace_ = data_avg_dir{counter_var, counter_ang}+data_sem_dir{counter_var, counter_ang};
        [~, idx_max_] = max(trace_);
        [~, idx_min_] = min(trace_(idx_max_:end));
        trace_max_min_ = trace_(idx_max_:idx_max_+idx_min_);
        ss_crossing_sem_p(counter_var, counter_ang) = find(trace_max_min_(1:end-1) >= 0 & trace_max_min_(2:end) < 0, 1, 'first') + idx_max_ - 250;
        
        trace_ = data_avg_dir{counter_var, counter_ang}-data_sem_dir{counter_var, counter_ang};
        [~, idx_max_] = max(trace_);
        [~, idx_min_] = min(trace_(idx_max_:end));
        trace_max_min_ = trace_(idx_max_:idx_max_+idx_min_);
        ss_crossing_sem_m(counter_var, counter_ang) = find(trace_max_min_(1:end-1) >= 0 & trace_max_min_(2:end) < 0, 1, 'first') + idx_max_ - 250;
        
    end
end

%% for ploting vmax vs. ss_crossing_mid_point: set ss_crossing_mid_point
ss_crossing_mid_point_avg   = nan(8,8);
ss_crossing_mid_point_sem_p = nan(8,8);
ss_crossing_mid_point_sem_m = nan(8,8);

for counter_ang = 5
    for counter_var = 1 : 8
        trace_ = data_avg_dir{counter_var, counter_ang};
        [~, idx_max_] = max(trace_);
        [~, idx_min_] = min(trace_(idx_max_:end));
        ss_crossing_mid_point_avg(counter_var, counter_ang) = idx_max_ + (idx_min_/2) - 250;
        
        trace_ = data_avg_dir{counter_var, counter_ang}+data_sem_dir{counter_var, counter_ang};
        [~, idx_max_] = max(trace_);
        [~, idx_min_] = min(trace_(idx_max_:end));
        ss_crossing_mid_point_sem_p(counter_var, counter_ang) = idx_max_ + (idx_min_/2) - 250;
        
        trace_ = data_avg_dir{counter_var, counter_ang}-data_sem_dir{counter_var, counter_ang};
        [~, idx_max_] = max(trace_);
        [~, idx_min_] = min(trace_(idx_max_:end));
        ss_crossing_mid_point_sem_m(counter_var, counter_ang) = idx_max_ + (idx_min_/2) - 250;
        
    end
end

%% for ploting vmax vs. ss_max_modulation: set v_max
vm_max_avg = nan(8,8);
vm_max_sem = nan(8,8);
vm_min_avg = nan(8,8);
vm_min_sem = nan(8,8);

for counter_ang = 1 : 8
    for counter_var = 1 : 8
        [max_, idx_] = max(data_avg_dir{counter_var, counter_ang});
        vm_max_avg(counter_var, counter_ang) = max_;
        vm_max_sem(counter_var, counter_ang) = data_sem_dir{counter_var, counter_ang}(1,idx_);
        
        [min_, idx_] = min(data_avg_dir{counter_var, counter_ang});
        vm_min_avg(counter_var, counter_ang) = min_;
        vm_min_sem(counter_var, counter_ang) = data_sem_dir{counter_var, counter_ang}(1,idx_);
    end
end

%% for ploting vmax vs. ss_max_modulation: plot v_max vs ss_max
clf(figure(11))
hold on
x_axis_data = vm_max_avg(:,5);
y_axis_data = ss_max_avg(:,5)*1000;
err_x_axis = vm_max_sem(:,5);
err_y_axis = ss_max_sem(:,5)*1000;

errorbar(x_axis_data,y_axis_data,err_y_axis,'vertical')
plot(x_axis_data, y_axis_data)

P_ = polyfit(x_axis_data, y_axis_data, 1);
y_axis_hat = polyval(P_,x_axis_data);
plot(x_axis_data, y_axis_hat, '-', 'LineWidth', 1)

xlabel('Max eye velocity (deg/s)')
ylabel('Max SS modulation (spk/s)')

ESN_Beautify_Plot(gcf, [2 2], 8)

[b,~,~,~,stats] = regress(y_axis_data,[ones(size(x_axis_data)) x_axis_data]);
p_value = stats(3)

%% for ploting vmax vs. ss_crossing: plot v_max vs ss_crossing
clf(figure(11))
hold on
x_axis_data = vm_max_avg(:,5);
y_axis_data = ss_crossing_avg(:,5);
err_p_y_axis = ss_crossing_sem_p(:,5)-ss_crossing_avg(:,5);
err_m_y_axis = ss_crossing_avg(:,5)-ss_crossing_sem_m(:,5);
% x_axis_data = vm_max_avg(:,5);
% y_axis_data = ss_crossing_mid_point_avg(:,5);
% err_p_y_axis = ss_crossing_mid_point_sem_p(:,5)-ss_crossing_mid_point_avg(:,5);
% err_m_y_axis = ss_crossing_mid_point_avg(:,5)-ss_crossing_mid_point_sem_m(:,5);

errorbar(x_axis_data,y_axis_data,err_m_y_axis,err_p_y_axis)
plot(x_axis_data, y_axis_data)

P_ = polyfit(x_axis_data, y_axis_data, 1);
y_axis_hat = polyval(P_,x_axis_data);
plot(x_axis_data, y_axis_hat, '-', 'LineWidth', 1)

xlabel('Max eye velocity (deg/s)')
ylabel('SS crossing (ms)')

ESN_Beautify_Plot(gcf, [1.5 1.5], 8)

[b,~,~,~,stats] = regress(y_axis_data,[ones(size(x_axis_data)) x_axis_data]);
p_value = stats(3)

end
%% function scratch_compare_cs_on_001_004_008_148()
function scratch_compare_cs_on_001_004_008_148()
%% Compare CS_on_population_001, CS_on_population_004, CS_on_population_008
CS_on_population_001 = load('CS_on_population_001.mat');
CS_on_population_004 = load('CS_on_population_004.mat');
CS_on_population_008 = load('CS_on_population_008.mat');
CS_on_population_148 = load('CS_on_population_148.mat');
num_pCells = length(CS_on_population_148.CS_on_population);
for counter_pCell = 1 : num_pCells
    CS_on_population_001.CS_ang_avg(counter_pCell,1) = CS_on_population_001.CS_on_population(counter_pCell).CS_ang_avg;
    CS_on_population_004.CS_ang_avg(counter_pCell,1) = CS_on_population_004.CS_on_population(counter_pCell).CS_ang_avg;
    CS_on_population_008.CS_ang_avg(counter_pCell,1) = CS_on_population_008.CS_on_population(counter_pCell).CS_ang_avg;
    CS_on_population_148.CS_ang_avg(counter_pCell,1) = CS_on_population_148.CS_on_population(counter_pCell).CS_ang_avg;
end
CS_ang_avg_001 = CS_on_population_001.CS_ang_avg;
CS_ang_avg_004 = CS_on_population_004.CS_ang_avg;
% CS_ang_avg_008 = CS_on_population_008.CS_ang_avg;
CS_ang_avg_148 = CS_on_population_148.CS_ang_avg;
%
CS_ang_pairs = [CS_ang_avg_001 CS_ang_avg_004];
x_values_ = cosd(CS_ang_pairs); 
y_values_ = sind(CS_ang_pairs);
diff_ang = acosd( (x_values_(:,1) .* x_values_(:,2)) + (y_values_(:,1) .* y_values_(:,2)) );
vec_1 = [x_values_(:,1)'; y_values_(:,1)'; zeros(1,num_pCells)];
vec_2 = [x_values_(:,2)'; y_values_(:,2)'; zeros(1,num_pCells)];
cross_ang = cross(vec_1, vec_2);
diff_ang_001_004 = diff_ang .* sign(cross_ang(3,:)');

CS_ang_pairs = [CS_ang_avg_001 CS_ang_avg_148];
x_values_ = cosd(CS_ang_pairs); 
y_values_ = sind(CS_ang_pairs);
diff_ang = acosd( (x_values_(:,1) .* x_values_(:,2)) + (y_values_(:,1) .* y_values_(:,2)) );
vec_1 = [x_values_(:,1)'; y_values_(:,1)'; zeros(1,num_pCells)];
vec_2 = [x_values_(:,2)'; y_values_(:,2)'; zeros(1,num_pCells)];
cross_ang = cross(vec_1, vec_2);
diff_ang_001_148 = diff_ang .* sign(cross_ang(3,:)');

CS_ang_pairs = [CS_ang_avg_004 CS_ang_avg_148];
x_values_ = cosd(CS_ang_pairs); 
y_values_ = sind(CS_ang_pairs);
diff_ang = acosd( (x_values_(:,1) .* x_values_(:,2)) + (y_values_(:,1) .* y_values_(:,2)) );
vec_1 = [x_values_(:,1)'; y_values_(:,1)'; zeros(1,num_pCells)];
vec_2 = [x_values_(:,2)'; y_values_(:,2)'; zeros(1,num_pCells)];
cross_ang = cross(vec_1, vec_2);
diff_ang_004_148 = diff_ang .* sign(cross_ang(3,:)');

step_size_ = 22.5;
ang_edges = -135-(step_size_/2):step_size_:135+(step_size_/2);

hFig = figure(7);
clf(hFig)

subplot(1,3,1)
hold on
histogram(diff_ang_001_004, ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
histogram(diff_ang_001_004, ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 1)
xline(mean(diff_ang_001_004),'Color', 'r', 'linewidth', 1);
set(gca, 'XTick', -135:45:135)
ylim([0 70])
ylabel('Count')
title('tag1 - tag4')

subplot(1,3,2)
hold on
histogram(diff_ang_001_148, ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
histogram(diff_ang_001_148, ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 1)
xline(mean(diff_ang_001_148),'Color', 'r', 'linewidth', 1);
set(gca, 'XTick', -135:45:135)
ylim([0 70])
title('tag1 - tag148')

subplot(1,3,3)
hold on
histogram(diff_ang_004_148, ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'r')
histogram(diff_ang_004_148, ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 1)
xline(mean(diff_ang_004_148),'Color', 'r', 'linewidth', 1);
set(gca, 'XTick', -135:45:135)
ylim([0 70])
title('tag4 - tag148')

ESN_Beautify_Plot(hFig, [4 1.5], 8)

fprintf(['diff_ang_001_004, mean: ' num2str(mean(diff_ang_001_004)) '\n'])
fprintf(['diff_ang_001_004, sem : ' num2str(std( diff_ang_001_004)./sqrt(num_pCells)) '\n'])

fprintf(['diff_ang_001_148, mean: ' num2str(mean(diff_ang_001_148)) '\n'])
fprintf(['diff_ang_001_148, sem : ' num2str(std( diff_ang_001_148)./sqrt(num_pCells)) '\n'])

fprintf(['diff_ang_004_148, mean: ' num2str(mean(diff_ang_004_148)) '\n'])
fprintf(['diff_ang_004_148, sem : ' num2str(std( diff_ang_004_148)./sqrt(num_pCells)) '\n'])

% num_perm = 5000;
% diff_ang_001_004_perm = nan(num_perm, 1);
% diff_ang_001_148_perm = nan(num_perm, 1);
% diff_ang_004_148_perm = nan(num_perm, 1);
% for counter_perm = 1 : num_perm
%     diff_ang_001_004_perm(counter_perm, 1) = median(diff_ang_001_004(randi(num_pCells, num_pCells, 1)));
%     diff_ang_001_148_perm(counter_perm, 1) = median(diff_ang_001_148(randi(num_pCells, num_pCells, 1)));
%     diff_ang_004_148_perm(counter_perm, 1) = median(diff_ang_004_148(randi(num_pCells, num_pCells, 1)));
% end
% 
% fprintf(['diff_ang_001_004, median   : ' num2str(mean(diff_ang_001_004_perm)) '\n'])
% fprintf(['diff_ang_001_004, SE-median: ' num2str(std( diff_ang_001_004_perm)) '\n'])
% 
% fprintf(['diff_ang_001_148, median   : ' num2str(mean(diff_ang_001_148_perm)) '\n'])
% fprintf(['diff_ang_001_148, SE-median: '   num2str(std( diff_ang_001_148_perm)) '\n'])
% 
% fprintf(['diff_ang_004_148, median   : ' num2str(mean(diff_ang_004_148_perm)) '\n'])
% fprintf(['diff_ang_004_148, SE-median: '   num2str(std( diff_ang_004_148_perm)) '\n'])

%
% step_size_ = pi/8;
% ang_edges = -pi-(step_size_/2):step_size_:pi-(step_size_/2);
% 
% hFig = figure(7);
% clf(hFig)
% 
% subplot(1,3,1)
% polarhistogram(deg2rad(diff_ang_001_004), ang_edges, 'DisplayStyle', 'bar','FaceColor','red', 'EdgeColor', 'none')
% hold on
% polarhistogram(deg2rad(diff_ang_001_004), ang_edges, 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'r', 'linewidth', 1)
% set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:10:60)
% rlim([0 70])
% title('tag1 - tag4')
% 
% subplot(1,3,2)
% polarhistogram(deg2rad(diff_ang_001_148), ang_edges, 'DisplayStyle', 'bar','FaceColor','red', 'EdgeColor', 'none')
% hold on
% polarhistogram(deg2rad(diff_ang_001_148), ang_edges, 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'r', 'linewidth', 1)
% set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:10:60)
% rlim([0 70])
% title('tag1 - tag148')
% 
% subplot(1,3,3)
% polarhistogram(deg2rad(diff_ang_004_148), ang_edges, 'DisplayStyle', 'bar','FaceColor','red', 'EdgeColor', 'none')
% hold on
% polarhistogram(deg2rad(diff_ang_004_148), ang_edges, 'DisplayStyle', 'stairs','FaceColor','none', 'EdgeColor', 'r', 'linewidth', 1)
% set(gca, 'ThetaTick', 0:45:315, 'RTick', 0:10:60)
% rlim([0 70])
% title('tag4 - tag148')
% 
% ESN_Beautify_Plot(hFig, [5 3], 8)

end
%% function scratch_rename_meta_data()
function scratch_rename_meta_data()
%%
flag_pair_list = false; % This should be false. DO NOT change it to true
pCell_list = PH_build_pCell_list(flag_pair_list);
path_data_monkey_sorted = uigetdir;

if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
pCell_list_isstr = arrayfun(@iscellstr,pCell_list);
num_pCells = size(pCell_list, 1);

% Loop over pCells
for counter_pCell = 1 : 1 : num_pCells
    fprintf(['### ' 'Analyzing pCell no. ', num2str(counter_pCell), ' / ' num2str(num_pCells) ' ###' '\n']);
    num_recording = nansum(pCell_list_isstr(counter_pCell, :));
    for counter_recording = 1 : 1 : num_recording
        %% build plot_data address
        file_name_cell = pCell_list{counter_pCell, counter_recording}; % '190423_142023_01_sorted_ESN_plot_data';
        file_name_cell = file_name_cell(1:13);

        year_ = file_name_cell(1:2);
        month_ = file_name_cell(3:4);
        day_ = file_name_cell(5:6);
        hour_ = file_name_cell(8:9);
        minute_ = file_name_cell(10:11);
        second_ = file_name_cell(12:13);
        subFolder_month = ['20' year_ '-' month_ filesep];
        subFolder_day = ['20' year_ '-' month_ '-' day_ filesep];
        subFolder_recording = ['20' year_ '-' month_ '-' day_ '_' hour_ '-' minute_ '-' second_ filesep];
        subFolder_data = ['analyzed_data' filesep];
        %% load BEHAVE data
        file_path = [path_data_monkey_sorted subFolder_month subFolder_day subFolder_recording subFolder_data];
        cd(file_path)
        meta_data_file = dir('*_meta_data.json');
        if ~strcmp([meta_data_file.name], [file_name_cell '_meta_data.json'])
            movefile([meta_data_file.name], [file_name_cell '_meta_data.json']);
        end
    end
end
end
