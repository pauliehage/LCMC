%% function plot_sac_sorter
function PGH_plot_lick_sorter(LICKS_ALL_DATA, params)
% Set parameters
amp_edges  = params.lick.amp_edges;
dur_edges  = params.lick.dur_edges;
ang_edges  = params.lick.ang_edges;
% Set parameters
num_row = 9;
num_col = 9;

% Init plot
hFig = figure;
clf(hFig)
hold on

% plot
for counter_tag = 1 : 9

    idx_tag = (LICKS_ALL_DATA.tag==counter_tag);
    title_ = params.lick.tag_name_list{counter_tag};


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
    plot(LICKS_ALL_DATA.tongue_tip_px(:,idx_tag), ...
        LICKS_ALL_DATA.tongue_tip_py(:,idx_tag), 'k')
    plot(LICKS_ALL_DATA.rew_r_px_dmax(:,idx_tag), ...
        LICKS_ALL_DATA.rew_r_py_dmax(:,idx_tag), 'oy')
    plot(LICKS_ALL_DATA.tongue_tip_px_dmax(:,idx_tag), ...
        LICKS_ALL_DATA.tongue_tip_py_dmax(:,idx_tag), 'om')
    plot(LICKS_ALL_DATA.rew_r_px_dmax(:,idx_tag), ...
        LICKS_ALL_DATA.rew_r_py_dmax(:,idx_tag), 'oy')
    plot(LICKS_ALL_DATA.rew_l_px_dmax(:,idx_tag), ...
        LICKS_ALL_DATA.rew_l_py_dmax(:,idx_tag), 'oy')
    plot(LICKS_ALL_DATA.nose_r_px_dmax(:,idx_tag), ...
        LICKS_ALL_DATA.nose_r_py_dmax(:,idx_tag), 'ok')
    plot(LICKS_ALL_DATA.nose_l_px_dmax(:,idx_tag), ...
        LICKS_ALL_DATA.nose_l_py_dmax(:,idx_tag), 'ok')
    plot(LICKS_ALL_DATA.rtube_r_px_dmax(:,idx_tag), ...
        LICKS_ALL_DATA.rtube_r_py_dmax(:,idx_tag), 'sr')
    plot(LICKS_ALL_DATA.rtube_l_px_dmax(:,idx_tag), ...
        LICKS_ALL_DATA.rtube_l_py_dmax(:,idx_tag), 'sr')
    plot(LICKS_ALL_DATA.ltube_r_px_dmax(:,idx_tag), ...
        LICKS_ALL_DATA.ltube_r_py_dmax(:,idx_tag), 'sb')
    plot(LICKS_ALL_DATA.ltube_l_px_dmax(:,idx_tag), ...
        LICKS_ALL_DATA.ltube_l_py_dmax(:,idx_tag), 'sb')
    title([title_ ': ' num2str(sum(idx_tag)) ' lick'], 'interpret', 'none');
    xlim([-20, 20]);
    ylim([-20, 20]);
    set(gca, 'YDir','reverse')

    subplot(num_row,num_col,axes_amp_dis)
    hold on
    histogram(LICKS_ALL_DATA.tongue_dm_max(:,idx_tag), amp_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
    histogram(LICKS_ALL_DATA.tongue_dm_max(:,idx_tag), amp_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
    xlim([amp_edges(1), amp_edges(end)]);
    set(gca, 'XTick', 0:5:25)
    ylabel('Amplitude')

    subplot(num_row,num_col,axes_ang_dis)
    hold on
    histogram(LICKS_ALL_DATA.tongue_ang_max(:,idx_tag), ang_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
    histogram(LICKS_ALL_DATA.tongue_ang_max(:,idx_tag), ang_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
    xlim([ang_edges(1), ang_edges(end)]);
    set(gca, 'XTick', -90:45:90)
    ylabel('Angle')

    subplot(num_row,num_col,axes_react_dis)
    hold on
    histogram(LICKS_ALL_DATA.duration_lick(:,idx_tag)*1000, dur_edges, 'DisplayStyle', 'bar', 'EdgeColor', 'none', 'FaceColor', 'k')
    histogram(LICKS_ALL_DATA.duration_lick(:,idx_tag)*1000, dur_edges, 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'FaceColor', 'none', 'linewidth', 2)
    xlim([0, 500]);
    set(gca, 'XTick', (0:100:500))
    ylabel('Duration')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add the title info
sgtitle([params.cell_name ', ' ...
    'lick: ' num2str(length(LICKS_ALL_DATA.duration_lick)) ', ' ...
    'bout: ' num2str(length(LICKS_ALL_DATA.duration_bout)) ', ' ...
    'harvest: ' num2str(length(LICKS_ALL_DATA.duration_harvest)) ', ' ...
    'dur: ' num2str(params.duration/60,3) 'min' ...
    ], ...
    'interpret', 'none');
ESN_Beautify_Plot(hFig, [20 10], 7)
end
