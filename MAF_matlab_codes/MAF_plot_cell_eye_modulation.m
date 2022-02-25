function MAF_plot_cell_eye_modulation(SACS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS, CS_on_data, tag_id)

hFig = figure(1);

DB.file_name = EXPERIMENT_PARAMS.name;
DB.sample_rate = 3e4;
DB.waveform = Neural_Properties.waveform;
DB.type = Neural_Properties.type;
DB.SS_duration = Neural_Properties.SS_duration;
DB.ss_isi = diff(Neural_Properties.SS_time);
DB.cs_isi = diff(Neural_Properties.CS_time);
DB.ss_xprob = Neural_Properties.Corr_data_SS_SSxSS_AUTO;
DB.ss_xprob_span = Neural_Properties.Corr_data_SS_inds_span;
DB.cs_xprob = Neural_Properties.Corr_data_CS_CSxSS_AUTO;
DB.cs_xprob_span = Neural_Properties.Corr_data_CS_inds_span;
DB.CS_firing_rate = Neural_Properties.CS_firing_rate;
DB.SS_firing_rate = Neural_Properties.SS_firing_rate;
DB.numSS = Neural_Properties.SS_num;
DB.numCS = Neural_Properties.CS_num;

% Plot Properties

subplot(3, 4, 1);
hold on;
plot_waveform(DB);

subplot(3, 4, 5);
hold on;
plot_ss_isi(DB);

subplot(3, 4, 9);
hold on;
plot_xprob(DB);

plot_title(DB, hFig)

% Plot Modulation
time_ind = 1:500;
length_trace = length(time_ind);
inds_span = time_ind-250;

AxMod000 = subplot(3,4,8);
AxMod045 = subplot(3,4,4);
AxMod090 = subplot(3,4,3);
AxMod135 = subplot(3,4,2);
AxMod180 = subplot(3,4,6);
AxMod225 = subplot(3,4,10);
AxMod270 = subplot(3,4,11);
AxMod315 = subplot(3,4,12);

ang_step = 45;
ang_values   = (0) : ang_step : (360 - ang_step);

% Figure parameters
Line_Color = lines(7);
color_SS = Line_Color(1,:);
color_offset = Line_Color(2,:);
color_visual = Line_Color(3,:);
color_onset = Line_Color(4,:);
color_vm = Line_Color(5,:);
color_SS_firing = [0    0.3    0.5];
color_CS = Line_Color(7,:);

range_SS_Firing = [0 200];

tag_name_list = EXPERIMENT_PARAMS.sac_tag_list;
event_type_list       = {'visual', 'onset', 'vmax', 'offset', 'auditory'};
% Which alignment to plot, depending on the tag
which_event_for_tag = cell(length(tag_name_list),length(event_type_list));
which_event_for_tag(1,1:2) = {'visual','onset'}; % tag 1
which_event_for_tag(2,1) = {'onset'}; % tag 2
which_event_for_tag(3,1) = {'onset'}; % tag 3
which_event_for_tag(4,1:2) = {'visual','onset'}; % tag 4
which_event_for_tag(5,1) = {'onset'}; % tag 5
which_event_for_tag(6,1) = {'onset'}; % tag 6
which_event_for_tag(7,1) = {'onset'}; % tag 7
which_event_for_tag(8,1) = {'onset'}; % tag 8
which_event_for_tag(9,1) = {'onset'}; % tag 9
which_event_for_tag(10,1) = {'onset'}; % tag 10
which_event_for_tag(11,1) = {'visual'}; % tag 11
which_event_for_tag(12,1:2) = {'visual','onset'}; % tag 12
which_event_for_tag(13,1) = {'visual'}; % tag 13

which_event_for_tag_isstr = arrayfun(@iscellstr,which_event_for_tag);

alignment = 'onset';

axs = [AxMod000,...
    AxMod045,...
    AxMod090,...
    AxMod135,...
    AxMod180,...
    AxMod225,...
    AxMod270,...
    AxMod315];

sac_data_dir = buildSacData(SACS_ALL_DATA,tag_id,alignment);

num_trial_dir = zeros(1,8);
for counter_dir = 1 : 8
    num_trial_dir(counter_dir) = size(sac_data_dir(counter_dir).time_onset,2);
end
num_trial_dir_max = round(median(num_trial_dir));
% Plot data, Loop over dirs
for counter_dir = 1 : 8

    hold(axs(counter_dir),'on')

    if isempty(sac_data_dir(counter_dir).SS)
        continue;
    end
    num_trial_dir_max_ = min([num_trial_dir_max size(sac_data_dir(counter_dir).time,2)]);

    % Plot onset raster
    train_data_logic_SS_ = sac_data_dir(counter_dir).SS(time_ind,:)';
    firing_SS_ = mean(sac_data_dir(counter_dir).SS(:,:),2) * 1000;
    firing_SS_ = ESN_smooth(firing_SS_);
    firing_SS_(firing_SS_<0)=0;
    if length(firing_SS_) == 1
        firing_SS_ = zeros(1,length(sac_data_dir(counter_dir).SS(:,:)));
    end
    firing_SS_ = firing_SS_(time_ind);
    train_data_logic_CS_ = sac_data_dir(counter_dir).CS(time_ind,:)';
    
    % plot SS raster
    [x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes( ...
        train_data_logic_SS_, inds_span, 0.5);
    plot(axs(counter_dir), x_axis_SS_(:), ...
        y_axis_SS_(:),'Color',color_SS,'LineWidth',1);

    % plot CS raster
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes( ...
        train_data_logic_CS_, inds_span, 1);
    plot(axs(counter_dir), x_axis_CS_(:), y_axis_CS_(:),...
        'Marker','.','MarkerSize',10,'LineWidth',5,'Color',color_CS);

    num_sac = length(sac_data_dir(counter_dir).time_visual);

    % Visual raster
    time_visual_shifted = (sac_data_dir(counter_dir).time_visual - sac_data_dir(counter_dir).(sprintf('time_%s',alignment))).*1000;
    time_visual_shifted = ESN_Round(time_visual_shifted);
    idx_visual_shifted = time_visual_shifted + length_trace/2;
    train_data_logic_visual = repmat(1:length_trace,num_sac,1) == idx_visual_shifted';
    [x_axis_visual_, y_axis_visual_] = ESN_raster_plot_axes(train_data_logic_visual, inds_span, 1);
    plot(axs(counter_dir), x_axis_visual_(:), y_axis_visual_(:),...
        'LineWidth', 2, 'Color', color_visual);

    % Onset raster
    time_onset_shifted = (sac_data_dir(counter_dir).time_onset - sac_data_dir(counter_dir).(sprintf('time_%s',alignment))).*1000;
    time_onset_shifted = ESN_Round(time_onset_shifted);
    idx_onset_shifted = time_onset_shifted + length_trace/2;
    train_data_logic_onset = repmat(1:length_trace,num_sac,1) == idx_onset_shifted';
    [x_axis_onset_, y_axis_onset_] = ESN_raster_plot_axes(train_data_logic_onset, inds_span, 1);
    plot(axs(counter_dir), x_axis_onset_(:), y_axis_onset_(:),...
        'LineWidth', .1, 'Color', color_onset);

    % Offset raster
    time_offset_shifted = (sac_data_dir(counter_dir).time_offset - sac_data_dir(counter_dir).(sprintf('time_%s',alignment))).*1000;
    time_offset_shifted = ESN_Round(time_offset_shifted);
    idx_offset_shifted = time_offset_shifted + length_trace/2;
    train_data_logic_offset = repmat(1:length_trace,num_sac,1) == idx_offset_shifted';
    [x_axis_offset_, y_axis_offset_] = ESN_raster_plot_axes(train_data_logic_offset, inds_span, 1);
    plot(axs(counter_dir), x_axis_offset_(:), y_axis_offset_(:),...
        'LineWidth', 1, 'Color', color_offset);

    ylim(axs(counter_dir),[0,num_sac+.5]);
    xlim(axs(counter_dir),[inds_span(1),inds_span(end)]);

    % plot SS rate
    yyaxis(axs(counter_dir),'right');
    ylim(axs(counter_dir),range_SS_Firing);
    set(axs(counter_dir), 'YColor', color_SS_firing);

    % Plot vm
    vm_ = sac_data_dir(counter_dir).eye_vm(time_ind,:);
    vm_ = nanmean(vm_,2);
    area(axs(counter_dir),inds_span, vm_/3,...
        'FaceColor', color_vm,'FaceAlpha',.5);
    
    hold on;

    plot(axs(counter_dir), inds_span, ...
        firing_SS_,'-', 'LineWidth', 2, 'Color', color_SS_firing);

end

subplot(3,4,7);

if not(isempty(CS_on_data))
    % Plot CS-Tuning
    prob_amplitude = CS_on_data.CS_prob_avg;

    vonMises_std = CS_on_data.vonMises_std;
    CS_ang_avg = CS_on_data.CS_ang_avg;
    CS_rho_avg = CS_on_data.CS_rho_avg;
    std_curv_ang = (CS_on_data.CS_ang_avg-vonMises_std) : 2 : (CS_on_data.CS_ang_avg+vonMises_std);
    std_curv_amp = repmat(CS_rho_avg, length(std_curv_ang), 1);

    plot_data_amp_mean = [prob_amplitude,...
        prob_amplitude(1), nan]';
    plot_data_deg_mean = [ang_values, ang_values(1), nan]';

    polarplot(deg2rad(plot_data_deg_mean), ...
        plot_data_amp_mean, 'Color', color_CS);
    hold on;
    polarplot(deg2rad(std_curv_ang), std_curv_amp, 'Color', color_CS);
    polarplot([0 deg2rad(CS_ang_avg)], [0 CS_rho_avg],'Color',color_CS);
end
ESN_Beautify_Plot(hFig, [12, 8], 8)
end

%% function plot_title
function plot_title(DB, fig_handle)
file_name = DB.file_name;
type = DB.type;
duration = round(DB.SS_duration/60,1);
numCS = DB.numCS;
freqCS = DB.CS_firing_rate;
numSS = DB.numSS;
freqSS = DB.SS_firing_rate;
text = sprintf('%s (%s) :: Duration: %.1f min, numCS: %.0f, freqCS: %.2f Hz, numSS: %.0f, freqSS: %.2f Hz',...
    file_name, type, (duration), numCS, freqCS, numSS, freqSS);
sgtitle(fig_handle, text, 'Interpreter', 'none');
end

%% function plot_waveform
function plot_waveform(DB)
hold on
ss_wave = DB.waveform.ss_wave;
cs_wave = DB.waveform.cs_wave;
sample_rate = DB.sample_rate;

x_ = DB.waveform.ch_map.x * 4;
y_ = DB.waveform.ch_map.y * 100;

if not(isnan(cs_wave))
    n_ch = size(cs_wave,1);
    n_sig = length(cs_wave);
    wave_ch = DB.waveform.cs_wave_ch;
elseif not(isnan(ss_wave))
    n_ch = size(ss_wave,1);
    n_sig = length(ss_wave);
    wave_ch = DB.waveform.ss_wave_ch;
end

x = x_(wave_ch+1);
y = y_(wave_ch+1);
ch_num = DB.waveform.ch_map.map(wave_ch+1);

span_ind = (0:n_sig-1)/sample_rate;
span_group_ = repmat([span_ind,nan],n_ch,1);
span_group = reshape((span_group_*1e3+x)',1,n_ch*(n_sig+1));

if not(isnan(ss_wave))
    ss_wave_ = [ss_wave, nan(n_ch,1)];
    ss_wave_group = reshape((ss_wave_+y)',1,n_ch*(n_sig+1));
    plot(span_group, ss_wave_group, '-b', 'linewidth', 2)
end

if not(isnan(cs_wave))
    cs_wave_ = [cs_wave, nan(n_ch,1)];
    cs_wave_group = reshape((cs_wave_+y)',1,n_ch*(n_sig+1));
    plot(span_group, cs_wave_group, '-r', 'linewidth', 2);
end

ch_map = arrayfun(@num2str, ch_num+1,'UniformOutput', 0);

text(x-1,y,ch_map)

axis off

plot([0,1],[0,0]+min(y)-100,'black','LineWidth',1);

end

%% function plot_xprob
function plot_xprob(DB)
ss_xprob = DB.ss_xprob;
cs_xprob = DB.cs_xprob;
ss_xprob_span = DB.ss_xprob_span;
ss_xprob(end/2) = nan;
cs_xprob_span = DB.cs_xprob_span;
plot(ss_xprob_span, ss_xprob, '-b', 'linewidth', 2)
plot(cs_xprob_span, cs_xprob, '-r', 'linewidth', 2)
xlabel('Time (ms)')
ylabel('Cross-Probability')
end

%% function plot_ss_isi
function plot_ss_isi(DB)
hold on
ss_isi = DB.ss_isi;
ss_isi_min = 0;
ss_isi_max = .1;
ss_isi_binNum = 50;
ss_isi_edges = linspace(ss_isi_min, ss_isi_max, ss_isi_binNum);
histogram(ss_isi, ss_isi_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'b')
histogram(ss_isi, ss_isi_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'b', 'FaceColor', 'none', 'linewidth', 2)
xlabel('SS ISI (S)')
ylabel('Count (#)')
end

%% Function buildSacData
function sac_data_dir = buildSacData(SACS_ALL_DATA,tag_id,alignment)

ang_step     = 45;
ang_edges    = (0 - (ang_step/2)) : ang_step : (360 + (ang_step/2));
ang_values   = (0) : ang_step : (360 - ang_step);
length_trace = 500;

% Build sac_data
idx_sacs    = ismember(SACS_ALL_DATA.tag, tag_id);
idx_sacs    = idx_sacs & SACS_ALL_DATA.validity;

sac_data.time    = SACS_ALL_DATA.(['time_',alignment])(:,idx_sacs);
sac_data.SS      = SACS_ALL_DATA.(['neuro_SS_',alignment])(:,idx_sacs);
sac_data.CS      = SACS_ALL_DATA.(['neuro_CS_',alignment])(:,idx_sacs);
sac_data.eye_vx  = SACS_ALL_DATA.(['eye_vx_',alignment])(:,idx_sacs);
sac_data.eye_vy  = SACS_ALL_DATA.(['eye_vy_',alignment])(:,idx_sacs);
    
sac_data.eye_vm = sqrt(sac_data.eye_vx.^2 + sac_data.eye_vy.^2);

sac_data.eye_r_px_onset   = SACS_ALL_DATA.eye_r_px_onset(:,idx_sacs);
sac_data.eye_r_py_onset   = SACS_ALL_DATA.eye_r_py_onset(:,idx_sacs);

sac_data.visual_px_offset = SACS_ALL_DATA.visual_px_offset(:,idx_sacs);
sac_data.visual_py_offset = SACS_ALL_DATA.visual_py_offset(:,idx_sacs);

sac_data.time_onset   = SACS_ALL_DATA.time_onset (     :,idx_sacs);
sac_data.time_offset  = SACS_ALL_DATA.time_offset(     :,idx_sacs);
sac_data.time_visual  = SACS_ALL_DATA.time_visual(     :,idx_sacs);

sac_data.time_diff_visual_onset = sac_data.time_onset - sac_data.time_visual;
sac_data.time_diff_onset_offset = sac_data.time_offset - sac_data.time_onset;

% Build sac_data_dir
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
end