function MAF_plot_NeuralProperties(Neural_Properties, params)

DB.duration = Neural_Properties.SS_duration;
DB.sample_rate = 3e4;
DB.file_name = params.cell_name;
DB.waveform = Neural_Properties.waveform;
DB.type = Neural_Properties.type;
DB.numCS = length(Neural_Properties.CS_time);
DB.numSS = length(Neural_Properties.SS_time);

DB.ss_isi = diff(Neural_Properties.SS_time);
DB.cs_isi = diff(Neural_Properties.CS_time);
DB.ss_xprob = Neural_Properties.Corr_data_SS_SSxSS_AUTO;
DB.ss_xprob_span = Neural_Properties.Corr_data_SS_inds_span;
DB.cs_xprob = Neural_Properties.Corr_data_CS_CSxSS_AUTO;
DB.cs_xprob_span = Neural_Properties.Corr_data_CS_inds_span;

%% Plot data
fig_handle = figure('Units','normalized','Position',[0,0,1,1]);
clf(fig_handle);
num_row = 2;
num_col = 2;

subplot(num_row, num_col, 1);
plot_waveform(DB);

subplot(num_row, num_col, 2);
plot_ss_isi(DB);

subplot(num_row, num_col, 3);
plot_xprob(DB);

subplot(num_row, num_col, 4);
plot_cs_isi(DB);

plot_title(DB, fig_handle)

ESN_Beautify_Plot(fig_handle, [15,9])


end

%% function plot_title
function plot_title(DB, fig_handle)
file_name = DB.file_name;
type = DB.type;
duration = DB.duration;
numCS = DB.numCS;
freqCS = numCS / duration;
numSS = DB.numSS;
freqSS = numSS / duration;
text = sprintf('%s (%s) :: Duration: %.1f min, numCS: %.0f, freqCS: %.2f Hz, numSS: %.0f, freqSS: %.2f Hz',...
    file_name, type, (duration / 60.), numCS, freqCS, numSS, freqSS);
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
hold on
ss_xprob = DB.ss_xprob;
cs_xprob = DB.cs_xprob;
ss_xprob_span = DB.ss_xprob_span;
cs_xprob_span = DB.cs_xprob_span;
ss_xprob(ss_xprob_span == 0) = nan;
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

%% function plot_cs_isi
function plot_cs_isi(DB)
hold on
cs_isi = DB.cs_isi;
cs_isi_min = 0;
cs_isi_max = 2;
cs_isi_binNum = 25;
cs_isi_edges = linspace(cs_isi_min, cs_isi_max, cs_isi_binNum);
histogram(cs_isi, cs_isi_edges, 'DisplayStyle', 'bar',  'EdgeColor', 'none', 'FaceColor', 'r')
histogram(cs_isi, cs_isi_edges, 'DisplayStyle', 'stairs',  'EdgeColor', 'r', 'FaceColor', 'none', 'linewidth', 2)
xlabel('CS ISI (S)')
ylabel('Count (#)')
end