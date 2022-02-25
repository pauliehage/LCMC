function PGH_plot_cell_tongue_modulation(LICKS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS,params , funcs,  tag_id, alignment)
PGH_global_variables();
global ang_values length_trace inds_span

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
time_ind = 1:length_trace;

AxModneg090 = subplot(3,4,6);
AxModneg045 = subplot(3,4,2);
AxModpos000 = subplot(3,4,3);
AxModpos45 = subplot(3,4,4);
AxModpos90 = subplot(3,4,8);

AxBoutStr_r = subplot(3,4,12);
AxBoutStr_l = subplot(3,4,10);
AxBoutStr_g = subplot(3,4,11);


% Figure parameters
Line_Color = lines(7);
color_SS = Line_Color(1,:);
color_offset = Line_Color(2,:);
color_onset = Line_Color(3,:);
color_dmax = Line_Color(4,:);
color_dm = Line_Color(5,:);
color_SS_firing = [0    0.3    0.5];
color_CS = Line_Color(7,:);
color_CS_firing = [0.8    0.2   0];

range_SS_Firing = [0 200];

axs = [AxModneg090,...
    AxModneg045,...
    AxModpos000,...
    AxModpos45,...
    AxModpos90];

axs_bout = [AxBoutStr_r,...
    AxBoutStr_l,...
    AxBoutStr_g];

lick_data_dir = buildLickData(LICKS_ALL_DATA,params, funcs, tag_id,alignment);
bout_data_dir = buildBoutData(LICKS_ALL_DATA,params, funcs, 1:3,alignment);

num_lick_dir = zeros(1,size(lick_data_dir,2));
for counter_dir = 1 : size(lick_data_dir,2)
    num_lick_dir(counter_dir) = size(lick_data_dir(counter_dir).time_onset,2);
end

% Plot lick data, Loop over dirs
for counter_dir = 1 : size(lick_data_dir,2)

    hold(axs(counter_dir),'on')

    if isempty(lick_data_dir(counter_dir).SS)
        continue;
    end

    % Plot onset raster
    train_data_logic_SS_ = lick_data_dir(counter_dir).SS(time_ind,:)';
    firing_SS_ = mean(lick_data_dir(counter_dir).SS(:,:),2) * 1000;
    firing_SS_ = ESN_smooth(firing_SS_);
    firing_SS_(firing_SS_<0)=0;
    if length(firing_SS_) == 1
        firing_SS_ = zeros(1,length(lick_data_dir(counter_dir).SS(:,:)));
    end
    firing_SS_ = firing_SS_(time_ind);

    train_data_logic_CS_ = lick_data_dir(counter_dir).CS(time_ind,:)';
    firing_CS_ = mean(lick_data_dir(counter_dir).CS(:,:),2) * 1000;
    firing_CS_ = ESN_smooth(firing_CS_,50);
    firing_CS_(firing_CS_<0)=0;
    if length(firing_CS_) == 1
        firing_CS_ = zeros(1,length(lick_data_dir(counter_dir).CS(:,:)));
    end
    firing_CS_ = firing_CS_(time_ind)*10;

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

    num_lick = length(lick_data_dir(counter_dir).time_onset);

    % Onset raster
    time_onset_shifted = (lick_data_dir(counter_dir).time_onset - lick_data_dir(counter_dir).(sprintf('time_%s',alignment))).*1000;
    time_onset_shifted = ESN_Round(time_onset_shifted);
    idx_onset_shifted = time_onset_shifted + length_trace/2;
    train_data_logic_onset = repmat(1:length_trace,num_lick,1) == idx_onset_shifted';
    [x_axis_onset_, y_axis_onset_] = ESN_raster_plot_axes(train_data_logic_onset, inds_span, 1);
    plot(axs(counter_dir), x_axis_onset_(:), y_axis_onset_(:),...
        'LineWidth', .1, 'Color', color_onset);

    % Dmax raster
    time_dmax_shifted = (lick_data_dir(counter_dir).time_dmax - lick_data_dir(counter_dir).(sprintf('time_%s',alignment))).*1000;
    time_dmax_shifted = ESN_Round(time_dmax_shifted);
    idx_dmax_shifted = time_dmax_shifted + length_trace/2;
    train_data_logic_dmax = repmat(1:length_trace,num_lick,1) == idx_dmax_shifted';
    [x_axis_dmax_, y_axis_dmax_] = ESN_raster_plot_axes(train_data_logic_dmax, inds_span, 1);
    plot(axs(counter_dir), x_axis_dmax_(:), y_axis_dmax_(:),...
        'LineWidth', .1, 'Color', color_dmax);

    % Offset raster
    time_offset_shifted = (lick_data_dir(counter_dir).time_offset - lick_data_dir(counter_dir).(sprintf('time_%s',alignment))).*1000;
    time_offset_shifted = ESN_Round(time_offset_shifted);
    idx_offset_shifted = time_offset_shifted + length_trace/2;
    train_data_logic_offset = repmat(1:length_trace,num_lick,1) == idx_offset_shifted';
    [x_axis_offset_, y_axis_offset_] = ESN_raster_plot_axes(train_data_logic_offset, inds_span, 1);
    plot(axs(counter_dir), x_axis_offset_(:), y_axis_offset_(:),...
        'LineWidth', 1, 'Color', color_offset);

    ylim(axs(counter_dir),[0,num_lick+.5]);
    xlim(axs(counter_dir),[inds_span(1),inds_span(end)]);

    % plot SS rate
    yyaxis(axs(counter_dir),'right');
    ylim(axs(counter_dir),range_SS_Firing);
    set(axs(counter_dir), 'YColor', color_SS_firing);

    % Plot dm
    dm_ = lick_data_dir(counter_dir).tongue_dm(time_ind,:);
    dm_ = nanmean(dm_,2);
    area(axs(counter_dir),inds_span, dm_*10,...
        'FaceColor', color_dm,'FaceAlpha',.5);

    hold on;
    plot(axs(counter_dir), inds_span, ...
        firing_SS_,'-', 'LineWidth', 2, 'Color', color_SS_firing);
    %      plot(axs(counter_dir), inds_span, ...
    %         firing_CS_,'-', 'LineWidth', 2, 'Color', color_CS_firing);

end

% Plot bout data, Loop over dirs
for counter_dir = 1 : size(bout_data_dir,2)

    hold(axs_bout(counter_dir),'on')

    if isempty(bout_data_dir(counter_dir).SS)
        continue;
    end

    % Plot onset raster
    train_data_logic_SS_ = bout_data_dir(counter_dir).SS(time_ind,:)';
    firing_SS_ = mean(bout_data_dir(counter_dir).SS(:,:),2) * 1000;
    firing_SS_ = ESN_smooth(firing_SS_);
    firing_SS_(firing_SS_<0)=0;
    if length(firing_SS_) == 1
        firing_SS_ = zeros(1,length(bout_data_dir(counter_dir).SS(:,:)));
    end
    firing_SS_ = firing_SS_(time_ind);

    train_data_logic_CS_ = bout_data_dir(counter_dir).CS(time_ind,:)';
    firing_CS_ = mean(bout_data_dir(counter_dir).CS(:,:),2) * 1000;
    firing_CS_ = ESN_smooth(firing_CS_,50);
    firing_CS_(firing_CS_<0)=0;
    if length(firing_CS_) == 1
        firing_CS_ = zeros(1,length(bout_data_dir(counter_dir).CS(:,:)));
    end
    firing_CS_ = firing_CS_(time_ind)*10;

    % plot SS raster
    [x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes( ...
        train_data_logic_SS_, inds_span, 0.5);
    plot(axs_bout(counter_dir), x_axis_SS_(:), ...
        y_axis_SS_(:),'Color',color_SS,'LineWidth',1);

    % plot CS raster
    [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes( ...
        train_data_logic_CS_, inds_span, 1);
    plot(axs_bout(counter_dir), x_axis_CS_(:), y_axis_CS_(:),...
        'Marker','.','MarkerSize',10,'LineWidth',5,'Color',color_CS);

    num_lick = length(bout_data_dir(counter_dir).time_onset);

    % Onset raster
    time_onset_shifted = (bout_data_dir(counter_dir).time_onset - bout_data_dir(counter_dir).(sprintf('time_%s',alignment))).*1000;
    time_onset_shifted = ESN_Round(time_onset_shifted);
    idx_onset_shifted = time_onset_shifted + length_trace/2;
    train_data_logic_onset = repmat(1:length_trace,num_lick,1) == idx_onset_shifted';
    [x_axis_onset_, y_axis_onset_] = ESN_raster_plot_axes(train_data_logic_onset, inds_span, 1);
    plot(axs_bout(counter_dir), x_axis_onset_(:), y_axis_onset_(:),...
        'LineWidth', .1, 'Color', color_onset);

    % Dmax raster
    time_dmax_shifted = (bout_data_dir(counter_dir).time_dmax - bout_data_dir(counter_dir).(sprintf('time_%s',alignment))).*1000;
    time_dmax_shifted = ESN_Round(time_dmax_shifted);
    idx_dmax_shifted = time_dmax_shifted + length_trace/2;
    train_data_logic_dmax = repmat(1:length_trace,num_lick,1) == idx_dmax_shifted';
    [x_axis_dmax_, y_axis_dmax_] = ESN_raster_plot_axes(train_data_logic_dmax, inds_span, 1);
    plot(axs_bout(counter_dir), x_axis_dmax_(:), y_axis_dmax_(:),...
        'LineWidth', .1, 'Color', color_dmax);

    % Offset raster
    time_offset_shifted = (bout_data_dir(counter_dir).time_offset - bout_data_dir(counter_dir).(sprintf('time_%s',alignment))).*1000;
    time_offset_shifted = ESN_Round(time_offset_shifted);
    idx_offset_shifted = time_offset_shifted + length_trace/2;
    train_data_logic_offset = repmat(1:length_trace,num_lick,1) == idx_offset_shifted';
    [x_axis_offset_, y_axis_offset_] = ESN_raster_plot_axes(train_data_logic_offset, inds_span, 1);
    plot(axs_bout(counter_dir), x_axis_offset_(:), y_axis_offset_(:),...
        'LineWidth', 1, 'Color', color_offset);

    ylim(axs_bout(counter_dir),[0,num_lick+.5]);
    xlim(axs_bout(counter_dir),[inds_span(1),inds_span(end)]);

    % plot SS rate
    yyaxis(axs_bout(counter_dir),'right');
    ylim(axs_bout(counter_dir),range_SS_Firing);
    set(axs_bout(counter_dir), 'YColor', color_SS_firing);

    % Plot dm
    dm_ = bout_data_dir(counter_dir).tongue_dm(time_ind,:);
    dm_ = nanmean(dm_,2);
    area(axs_bout(counter_dir),inds_span, dm_*10,...
        'FaceColor', color_dm,'FaceAlpha',.5);

    hold on;
    plot(axs_bout(counter_dir), inds_span, ...
        firing_SS_,'-', 'LineWidth', 2, 'Color', color_SS_firing);
    %      plot(axs(counter_dir), inds_span, ...
    %         firing_CS_,'-', 'LineWidth', 2, 'Color', color_CS_firing);

end
subplot(3,4,7);
CS_on_data = PGH_CS_on_analysis(LICKS_ALL_DATA, params, funcs, tag_id);
if not(isempty(CS_on_data))
    % Plot CS-Tuning
    prob_amplitude = CS_on_data.CS_fr_avg;

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
    ax = gca;
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    hold on;
    polarplot(deg2rad(std_curv_ang), std_curv_amp, 'Color', color_CS);
    ax = gca;
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    polarplot([0 deg2rad(CS_ang_avg)], [0 CS_rho_avg],'Color',color_CS);
    ax = gca;
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';

    thetaticks([-90 -45 0 45 90 ])
    thetaticklabels({ '-90', '-45', '0','45', '90',})
    thetalim([-90 90]);

end

ESN_Beautify_Plot(hFig, [20, 10], 8)

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

%% function buildLickData
function lick_data_dir = buildLickData(LICKS_ALL_DATA, tag_id,alignment, params, funcs)
PGH_global_variables();
global length_trace  ang_edges ang_values

% Build lick_data
idx_licks    = ismember(LICKS_ALL_DATA.tag, tag_id);
idx_licks    = idx_licks & LICKS_ALL_DATA.validity;

lick_data.time    = LICKS_ALL_DATA.(['time_',alignment])(:,idx_licks);
lick_data.SS      = LICKS_ALL_DATA.(['neuro_SS_',alignment])(:,idx_licks);
lick_data.CS      = LICKS_ALL_DATA.(['neuro_CS_',alignment])(:,idx_licks);
lick_data.tongue_dm  = LICKS_ALL_DATA.(['tongue_dm_',alignment])(:,idx_licks);
lick_data.tongue_vm  = LICKS_ALL_DATA.(['tongue_vm_',alignment])(:,idx_licks);
lick_data.tongue_ang = LICKS_ALL_DATA.(['tongue_ang_',alignment])(:,idx_licks);

lick_data.tongue_dm_max =  LICKS_ALL_DATA.tongue_dm_max(    :,idx_licks);
lick_data.tongue_vm_max =  LICKS_ALL_DATA.tongue_vm_max(    :,idx_licks);
lick_data.tongue_vm_min =  LICKS_ALL_DATA.tongue_vm_min(    :,idx_licks);
lick_data.tongue_ang_max =  LICKS_ALL_DATA.tongue_ang_max(    :,idx_licks);

lick_data.time_onset   = LICKS_ALL_DATA.time_onset (     :,idx_licks);
lick_data.time_vmax   = LICKS_ALL_DATA.time_vmax (     :,idx_licks);
lick_data.time_dmax  = LICKS_ALL_DATA.time_dmax   (     :,idx_licks);
lick_data.time_vmin   = LICKS_ALL_DATA.time_vmin (     :,idx_licks);
lick_data.time_offset  = LICKS_ALL_DATA.time_offset(     :,idx_licks);

lick_data.time_diff_onset_offset = lick_data.time_offset - lick_data.time_onset;
lick_data.time_diff_onset_dmax= lick_data.time_dmax - lick_data.time_onset;
lick_data.time_diff_offset_dmax= lick_data.time_offset - lick_data.time_onset;

% Build lick_data_dir
lick_data.tongue_ang_bin = discretize(lick_data.tongue_ang_max, ang_edges);

% 1: -90deg % 2: -45deg % 3: 0deg % 4: 45deg % 5: 90deg

if length(ang_values) ~= 5
    error('plot_single_session_modulation: length ang_values is not 5. Please modify the code.')
end
if length_trace ~= 500
    error('lick_modulation_index: length_trace is not 500. Please modify the code.')
end

field_names_lick_data = fieldnames(lick_data);
lick_data_dir = struct;
for counter_dir = 1 : length(ang_values)
    for counter_field = 1 : length(field_names_lick_data)
        field_name = field_names_lick_data{counter_field};
        idx_ang = (lick_data.tongue_ang_bin == counter_dir);
        lick_data_dir(counter_dir).(field_name) = lick_data.(field_name)(:,idx_ang);
    end
end
end
%% function buildBoutData
function bout_data_dir = buildBoutData(LICKS_ALL_DATA, params, funcs, tag_id,alignment)
% Build lick_data
tag_id = tag_id;
idx_bout    = ismember(LICKS_ALL_DATA.tag_bout, tag_id);
idx_bout    = idx_bout & LICKS_ALL_DATA.validity;
% idx_bout = LICKS_ALL_DATA.tag_bout(LICKS_ALL_DATA.tag_bout ~= 0);

bout_data.time    = LICKS_ALL_DATA.(['time_',alignment])(:,idx_bout);
bout_data.SS      = LICKS_ALL_DATA.(['neuro_SS_',alignment])(:,idx_bout);
bout_data.CS      = LICKS_ALL_DATA.(['neuro_CS_',alignment])(:,idx_bout);
bout_data.tongue_dm  = LICKS_ALL_DATA.(['tongue_dm_',alignment])(:,idx_bout);
bout_data.tongue_vm  = LICKS_ALL_DATA.(['tongue_vm_',alignment])(:,idx_bout);
bout_data.tongue_ang = LICKS_ALL_DATA.(['tongue_ang_',alignment])(:,idx_bout);

bout_data.tongue_dm_max =  LICKS_ALL_DATA.tongue_dm_max(    :,idx_bout);
bout_data.tongue_vm_max =  LICKS_ALL_DATA.tongue_vm_max(    :,idx_bout);
bout_data.tongue_vm_min =  LICKS_ALL_DATA.tongue_vm_min(    :,idx_bout);
bout_data.tongue_ang_max =  LICKS_ALL_DATA.tongue_ang_max(    :,idx_bout);

bout_data.time_onset   = LICKS_ALL_DATA.time_onset (     :,idx_bout);
bout_data.time_vmax   = LICKS_ALL_DATA.time_vmax (     :,idx_bout);
bout_data.time_dmax  = LICKS_ALL_DATA.time_dmax   (     :,idx_bout);
bout_data.time_vmin   = LICKS_ALL_DATA.time_vmin (     :,idx_bout);
bout_data.time_offset  = LICKS_ALL_DATA.time_offset(     :,idx_bout);

bout_data.time_diff_onset_offset = bout_data.time_offset - bout_data.time_onset;
bout_data.time_diff_onset_dmax= bout_data.time_dmax - bout_data.time_onset;
bout_data.time_diff_offset_dmax= bout_data.time_offset - bout_data.time_onset;

idx_bout_ = LICKS_ALL_DATA.tag_bout(LICKS_ALL_DATA.tag_bout ~= 0);
idx_bout_ = idx_bout_(ismember(idx_bout_, tag_id)) ;

field_names_bout_data = fieldnames(bout_data);
bout_data_dir = struct;
for counter_dir = 1 : length(tag_id)
    for counter_field = 1 : length(field_names_bout_data)
        field_name = field_names_bout_data{counter_field};
        idx_tag = (idx_bout_ == counter_dir);
        bout_data_dir(counter_dir).(field_name) = bout_data.(field_name)(:,idx_tag);
    end
end

end

