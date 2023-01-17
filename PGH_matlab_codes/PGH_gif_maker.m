function PGH_gif_maker(DLC, n_lick,format)
%% Build variables
fprintf('Building variables  ... ');
% load([file_path file_name]);
d_tip = DLC.KINEMATIC.d_tip;
v_tip = DLC.KINEMATIC.v_tip;
angle_midtip = DLC.KINEMATIC.angle_midtip;
tip_tongue_x = DLC.POINTS.tip_tongue_x;
tip_tongue_y = DLC.POINTS.tip_tongue_y;
r_tongue_x = DLC.POINTS.r_tongue_x;
r_tongue_y = DLC.POINTS.r_tongue_y;
l_tongue_x = DLC.POINTS.l_tongue_x;
l_tongue_y = DLC.POINTS.l_tongue_y;
mid_tongue_x = DLC.POINTS.mid_tongue_x;
mid_tongue_y = DLC.POINTS.mid_tongue_y;
r_nose_x =  DLC.POINTS.r_nose_x;
r_nose_y = DLC.POINTS.r_nose_y;
l_nose_x =  DLC.POINTS.l_nose_x;
l_nose_y = DLC.POINTS.l_nose_y;
r_tube_r_x = DLC.POINTS.r_tube_r_x;
r_tube_r_y = DLC.POINTS.r_tube_r_y;
r_tube_l_x = DLC.POINTS.r_tube_l_x;
r_tube_l_y = DLC.POINTS.r_tube_l_y;
l_tube_r_x = DLC.POINTS.l_tube_r_x;
l_tube_r_y = DLC.POINTS.l_tube_r_y;
l_tube_l_x = DLC.POINTS.l_tube_l_x;
l_tube_l_y = DLC.POINTS.l_tube_l_y;
r_food_x = DLC.POINTS.r_food_x;
r_food_y = DLC.POINTS.r_food_y;
l_food_x = DLC.POINTS.l_food_x;
l_food_y = DLC.POINTS.l_food_y;
d_tip_r_food = sqrt((tip_tongue_x - r_food_x).^2 + (tip_tongue_y - r_food_y).^2);
d_tip_l_food = sqrt((tip_tongue_x - l_food_x).^2 + (tip_tongue_y - l_food_y).^2);
% num_lick = DLC.IND.num_lick;
ind_lick_onset = DLC.IND.ind_lick_onset;
ind_lick_offset = DLC.IND.ind_lick_offset;
time_1K = DLC.TIME.time_1K;

if length(n_lick) > 1
    inds_ = ind_lick_onset(n_lick(1))-1:ind_lick_offset(n_lick(end)) +1;
    inds_v = ind_lick_onset(n_lick(1))-2:ind_lick_offset(n_lick(end)) +2;

else
    inds_ = ind_lick_onset(n_lick)-1:ind_lick_offset(n_lick) +1;
    inds_v = ind_lick_onset(n_lick)-2:ind_lick_offset(n_lick) +2;

end

x_lim = (time_1K(inds_(end)) - time_1K(inds_(1)))*1000;
fprintf(' --> Completed. \n')

%% Build gif images
fprintf('Building gif images  ... ');
path_to_analyzed = DLC.FILE.path_to_analyzed;
path_to_raw = 'D:\data_59d\2021-02\2021-02-03\2021-02-03_15-11-06\raw_data\';

vid = dir([path_to_analyzed '*.mp4']);
% vid_name = vid.name;
vid_name = vid.name(1:13);
fig = figure;
for counter_frame = 1:length(inds_)
    subplot(3,3,1) % video
%      frame = imshow(read(VideoReader([path_to_analyzed vid_name]), inds_(counter_frame)));
    frame = imshow(read(VideoReader([path_to_raw vid_name '.mp4']), inds_(counter_frame)));

    title('Video frame')

    subplot(3,3,2) % trace
    hold on
    line([r_tube_r_x(inds_(counter_frame)) r_tube_r_x(inds_(counter_frame))] ...
        ,[r_tube_r_y(inds_(counter_frame)) 20], 'Color', 'black');
    line([r_tube_l_x(inds_(counter_frame)) r_tube_l_x(inds_(counter_frame))] ...
        ,[r_tube_l_y(inds_(counter_frame)) 20], 'Color', 'black');
    line([r_tube_r_x(inds_(counter_frame)) r_tube_l_x(inds_(counter_frame))] ...
        , [r_tube_r_y(inds_(counter_frame)) r_tube_l_y(inds_(counter_frame))], 'Color', 'black');
    line([l_tube_r_x(inds_(counter_frame)) l_tube_r_x(inds_(counter_frame))] ...
        ,[-20 l_tube_r_y(inds_(counter_frame))], 'Color', 'black');
    line([l_tube_l_x(inds_(counter_frame)) l_tube_l_x(inds_(counter_frame))] ...
        ,[-20 l_tube_l_y(inds_(counter_frame))], 'Color', 'black');
    line([l_tube_r_x(inds_(counter_frame)) l_tube_l_x(inds_(counter_frame))] ...
        , [l_tube_r_y(inds_(counter_frame)) l_tube_l_y(inds_(counter_frame))], 'Color', 'black');

    plot(tip_tongue_x(inds_(1):inds_(counter_frame)),tip_tongue_y(inds_(1):inds_(counter_frame)),'.-k', 'Linewidth', 1);
%     plot(mid_tongue_x(inds_(1):inds_(counter_frame)),mid_tongue_y(inds_(1):inds_(counter_frame)),'.-g', 'Linewidth', 1);
%     plot(r_tongue_x(inds_(1):inds_(counter_frame)),r_tongue_y(inds_(1):inds_(counter_frame)),'.-r', 'Linewidth', 1);
%     plot(l_tongue_x(inds_(1):inds_(counter_frame)),l_tongue_y(inds_(1):inds_(counter_frame)),'.-b', 'Linewidth', 1);
    plot(r_tube_r_x(inds_(counter_frame)),r_tube_r_y(inds_(counter_frame)),'sk','markersize', 5);
    plot(r_tube_l_x(inds_(counter_frame)),r_tube_l_y(inds_(counter_frame)),'sk','markersize', 5);
    plot(l_tube_r_x(inds_(counter_frame)),l_tube_r_y(inds_(counter_frame)),'sk','markersize', 5);
    plot(l_tube_l_x(inds_(counter_frame)),l_tube_l_y(inds_(counter_frame)),'sk','markersize', 5);
    plot(r_nose_x(inds_(counter_frame)),r_nose_y(inds_(counter_frame)),'ok','markersize', 5);
    plot(l_nose_x(inds_(counter_frame)),l_nose_y(inds_(counter_frame)),'ok','markersize', 5);
    plot(r_food_x(inds_(counter_frame)),r_food_y(inds_(counter_frame)),'oy','markersize', 5);
    plot(l_food_x(inds_(counter_frame)),l_food_y(inds_(counter_frame)),'oy','markersize', 5);

    xlabel('x (mm)');
    ylabel('y (mm)');
    set(gca, 'YDir','reverse')
    xlim([-10 20])
    ylim([-25 25])
    title('Trajectory trace')

    subplot(3,3,3) % geo
    geo_tongue = polyshape([tip_tongue_x(inds_(counter_frame)),l_tongue_x(inds_(counter_frame)),...
        mid_tongue_x(inds_(counter_frame)),r_tongue_x(inds_(counter_frame))],[tip_tongue_y(inds_(counter_frame)),...
        l_tongue_y(inds_(counter_frame)),mid_tongue_y(inds_(counter_frame)),r_tongue_y(inds_(counter_frame))]);

    geo_r_tube_empty = polyshape([r_food_x(inds_(counter_frame)),r_tube_r_x(inds_(counter_frame)),...
        r_tube_r_x(inds_(counter_frame)),r_tube_l_x(inds_(counter_frame)), r_tube_l_x(inds_(counter_frame))],[r_food_y(inds_(counter_frame)),...
        r_food_y(inds_(counter_frame)), r_tube_r_y(inds_(counter_frame)),r_tube_l_y(inds_(counter_frame)), r_food_y(inds_(counter_frame))]);

    geo_r_tube_full = polyshape([r_tube_r_x(inds_(counter_frame)), r_tube_r_x(inds_(counter_frame)),...
        r_tube_l_x(inds_(counter_frame)), r_tube_l_x(inds_(counter_frame)), r_food_x(inds_(counter_frame))],...
        [r_food_y(inds_(counter_frame)), r_tube_r_y(inds_(counter_frame)) + abs(max(r_food_y)-min(r_food_y)), ...
        r_tube_l_y(inds_(counter_frame)) + abs(max(r_food_y)-min(r_food_y)), r_food_y(inds_(counter_frame)), r_food_y(inds_(counter_frame))]);

    geo_l_tube_empty = polyshape([l_food_x(inds_(counter_frame)),l_tube_r_x(inds_(counter_frame)),...
        l_tube_r_x(inds_(counter_frame)),l_tube_l_x(inds_(counter_frame)), l_tube_l_x(inds_(counter_frame))],[l_food_y(inds_(counter_frame)),...
        l_food_y(inds_(counter_frame)), l_tube_r_y(inds_(counter_frame)),l_tube_l_y(inds_(counter_frame)), l_food_y(inds_(counter_frame))]);

    geo_l_tube_full = polyshape([l_tube_r_x(inds_(counter_frame)), l_tube_r_x(inds_(counter_frame)),...
        l_tube_l_x(inds_(counter_frame)), l_tube_l_x(inds_(counter_frame)), l_food_x(inds_(counter_frame))],...
        [l_food_y(inds_(counter_frame)), l_tube_r_y(inds_(counter_frame)) - abs(max(l_food_y)-min(l_food_y)), ...
        l_tube_l_y(inds_(counter_frame)) - abs(max(l_food_y)-min(l_food_y)), l_food_y(inds_(counter_frame)), l_food_y(inds_(counter_frame))]);
    geo_inter_tongue_r_tube_empty = intersect(geo_tongue, geo_r_tube_empty);

    geo_inter_tongue_r_tube_full = intersect(geo_tongue, geo_r_tube_full);

    geo_inter_tongue_l_tube_empty = intersect(geo_tongue, geo_l_tube_empty);

    geo_inter_tongue_l_tube_full = intersect(geo_tongue, geo_l_tube_full);

    geo_all = [geo_tongue geo_r_tube_empty geo_r_tube_full geo_l_tube_empty geo_l_tube_full...
        geo_inter_tongue_r_tube_empty geo_inter_tongue_r_tube_full geo_inter_tongue_l_tube_empty...
        geo_inter_tongue_l_tube_full ];

    [cent_tongue_r_tube_empty_x, cent_tongue_r_tube_empty_y ] = centroid(geo_inter_tongue_r_tube_empty);
    [cent_tongue_r_tube_full_x, cent_tongue_r_tube_full_y] = centroid(geo_inter_tongue_r_tube_full);
    [cent_tongue_l_tube_empty_x,cent_tongue_l_tube_empty_y]  = centroid(geo_inter_tongue_l_tube_empty);
    [cent_tongue_l_tube_full_x,cent_tongue_l_tube_full_y] = centroid(geo_inter_tongue_l_tube_full);

    bool_overlaps_all = overlaps(geo_all);
    bool_tongue_r_tube_empty = bool_overlaps_all(1,2);
    bool_tongue_r_tube_full = bool_overlaps_all(1,3);
    bool_tongue_l_tube_empty= bool_overlaps_all(1,4);
    bool_tongue_l_tube_full = bool_overlaps_all(1,5);

    area_tongue(counter_frame) = area(geo_tongue);
    area_r_tube_empty = area(geo_r_tube_empty);
    area_r_tube_full = area(geo_r_tube_full);
    area_l_tube_empty = area(geo_l_tube_empty);
    area_l_tube_full = area(geo_l_tube_full);
    area_inter_tongue_r_tube_empty = area(geo_inter_tongue_r_tube_empty);
    area_inter_tongue_r_tube_full = area(geo_inter_tongue_r_tube_full);
    area_inter_tongue_l_tube_empty = area(geo_inter_tongue_l_tube_empty);
    area_inter_tongue_l_tube_full = area(geo_inter_tongue_l_tube_full);

    hold on
    plot(geo_tongue,'LineWidth', 1, 'FaceColor','red');
    plot(geo_r_tube_empty,'LineWidth', 1, 'FaceColor','white');
    plot(geo_l_tube_empty,'LineWidth', 1, 'FaceColor','white');
    plot(geo_r_tube_full,'LineWidth', 1, 'FaceColor','yellow');
    plot(geo_l_tube_full,'LineWidth', 1, 'FaceColor','yellow');
    plot(geo_inter_tongue_r_tube_empty,'LineWidth', 1, 'FaceColor','black');
    plot(geo_inter_tongue_r_tube_full,'LineWidth', 1, 'FaceColor','blue');
    plot(geo_inter_tongue_l_tube_empty,'LineWidth', 1, 'FaceColor','black');
    plot(geo_inter_tongue_l_tube_full,'LineWidth', 1, 'FaceColor','blue');
    plot(cent_tongue_r_tube_empty_x, cent_tongue_r_tube_empty_y,'*m','markersize', 7 );
    plot(cent_tongue_r_tube_full_x, cent_tongue_r_tube_full_y,'*m','markersize', 7);
    plot(cent_tongue_l_tube_empty_x, cent_tongue_l_tube_empty_y,'*m','markersize', 7);
    plot(cent_tongue_l_tube_full_x, cent_tongue_l_tube_full_y,'*m','markersize', 7);
    plot(tip_tongue_x(inds_(counter_frame)),tip_tongue_y(inds_(counter_frame)),'ok','markersize', 7);
    plot(r_tongue_x(inds_(counter_frame)),r_tongue_y(inds_(counter_frame)),'or','markersize', 7);
    plot(l_tongue_x(inds_(counter_frame)),l_tongue_y(inds_(counter_frame)),'ob','markersize', 7);
    plot(mid_tongue_x(inds_(counter_frame)),mid_tongue_y(inds_(counter_frame)),'og','markersize', 7);
    plot(r_nose_x(inds_(counter_frame)),r_nose_y(inds_(counter_frame)),'ok','markersize', 5);
    plot(l_nose_x(inds_(counter_frame)),l_nose_y(inds_(counter_frame)),'ok','markersize', 5);
    xlabel('x (mm)');
    ylabel('y (mm)');
    xlim([-10 20])
    ylim([-25 25])
    set(gca, 'YDir','reverse');
    title('Geometrized frame')

    subplot(3,3,4) % displacement
    plot((time_1K(inds_(1):inds_(counter_frame)) - time_1K(inds_(1)))*1000, d_tip(inds_(1):inds_(counter_frame)),'.-k', 'Linewidth', 1);
    xlabel('Time (ms)');
    ylabel('Displacement (mm)');
    ylim([0 25])
    xlim([0 x_lim])
    title('Displacement')

    subplot(3,3,5) % velocity
    plot((time_1K(inds_v(1):inds_v(counter_frame)) - time_1K(inds_v(1)))*1000, v_tip(inds_v(1):inds_v(counter_frame)),'.-k', 'Linewidth', 1);
    xlabel('Time (ms)');
    ylabel('Velocity (mms)');
    ylim([-650 650])
    xlim([0 x_lim])
    title('Velocity')

    subplot(3,3,6) % angle
    plot((time_1K(inds_(1):inds_(counter_frame)) - time_1K(inds_(1)))*1000, angle_midtip(inds_(1):inds_(counter_frame)),'.-k', 'Linewidth', 1);
    xlabel('Time (ms)');
    ylabel('Angle (deg)');
    ylim([-120 120])
    xlim([0 x_lim])
    title('Angle')

    subplot(3,3,7) % distance from left tube
    plot((time_1K(inds_(1):inds_(counter_frame)) - time_1K(inds_(1)))*1000, d_tip_l_food(inds_(1):inds_(counter_frame)),'.-k', 'Linewidth', 1);
    xlabel('Time (ms)');
    ylabel('Distance (mm)');
    ylim([0 50])
    xlim([0 x_lim])
    title('Distance to left food')

    subplot(3,3,8) % distance from right tube
    plot((time_1K(inds_(1):inds_(counter_frame)) - time_1K(inds_(1)))*1000, d_tip_r_food(inds_(1):inds_(counter_frame)),'.-k', 'Linewidth', 1);
    xlabel('Time (ms)');
    ylabel('Distance (mm)');
    ylim([0 50])
    xlim([0 x_lim])
    title('Distance to right food')

    subplot(3,3,9) % tongue area
    plot((time_1K(inds_(1):inds_(counter_frame)) - time_1K(inds_(1)))*1000, area_tongue(1 : counter_frame),'.-k', 'Linewidth', 1);
    xlabel('Time (ms)');
    ylabel('area (mm^s)');
    ylim([0 30])
    xlim([0 x_lim])
    title('Tongue area')

    sgtitle([vid_name(1:13) ' | lick: ' num2str(n_lick(1)) ' - ' num2str(n_lick(end)) ], 'interpret', 'none');

    ESN_Beautify_Plot(fig,[10 8])
    saveas(frame, ['Frame_' num2str(counter_frame)], format);
    %     close(gcf);
    clf
end
fprintf(' --> Completed. \n')


end
