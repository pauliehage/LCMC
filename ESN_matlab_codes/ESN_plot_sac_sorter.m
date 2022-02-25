%% function plot_sac_sorter
function ESN_plot_sac_sorter(SACS_ALL_DATA, params)
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
ESN_Beautify_Plot(gcf, [20 10])
end
