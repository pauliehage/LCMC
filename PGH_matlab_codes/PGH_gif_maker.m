function PGH_gif_maker(DLC)
if nargin<1
[file_name,file_path] = uigetfile([pwd filesep '*.mat'], 'Select DLC file');

end
%% load DLC
fprintf(['Loading ', file_name, ' ... ']);
load([file_path file_name]);
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
num_lick = DLC.IND.num_lick;
ind_lick_onset = DLC.IND.ind_lick_onset;
ind_lick_offset = DLC.IND.ind_lick_offset;

fprintf(' --> Completed. \n')

%% Load video
[file_name,file_path] = uigetfile([pwd filesep '*.mp4'], 'Select MP4 file');
fprintf(['Loading ', file_name, ' ... ']);
fprintf(' --> Completed. \n')

%% Specify lick #
n_lick = input('Specify lick #: ');
if length(n_lick) > 1
inds_ = ind_lick_onset(n_lick(1))-1:ind_lick_offset(n_lick(end)) +1;
inds_v = ind_lick_onset(n_lick(1)):ind_lick_offset(n_lick(end)) +2;

else
inds_ = ind_lick_onset(n_lick)-1:ind_lick_offset(n_lick) +1;
inds_v = ind_lick_onset(n_lick):ind_lick_offset(n_lick) +2;

end
x_lim = (time_vid(inds_(end)) - time_vid(inds_(1)))*1000;
fprintf(' --> Completed. \n')

%% Grab frames from video
fprintf('Grabbing frames  ... ');
for counter_frame = 1:length(inds_)
    frame = imshow(read(VideoReader(file_name), inds_(counter_frame)));
    saveas(frame, ['Vid_' num2str(counter_frame)], 'png');
    close(gcf);
end
fprintf(' --> Completed. \n')

%% Grab frames from geo
fprintf('Geometrizing tongue ...');
for counter_frame = 1:1:length(inds_)   
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
    
    area_tongue = area(geo_tongue);
    area_r_tube_empty = area(geo_r_tube_empty);
    area_r_tube_full = area(geo_r_tube_full);
    area_l_tube_empty = area(geo_l_tube_empty);
    area_l_tube_full = area(geo_l_tube_full);
    area_inter_tongue_r_tube_empty = area(geo_inter_tongue_r_tube_empty);
    area_inter_tongue_r_tube_full = area(geo_inter_tongue_r_tube_full);
    area_inter_tongue_l_tube_empty = area(geo_inter_tongue_l_tube_empty);
    area_inter_tongue_l_tube_full = area(geo_inter_tongue_l_tube_full);
    
    figure;
    hold on;
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
    plot(tip_tongue_x(inds_(counter_frame)),tip_tongue_y(inds_(counter_frame)),'or','markersize', 7);
    plot(r_tongue_x(inds_(counter_frame)),r_tongue_y(inds_(counter_frame)),'oc','markersize', 7);
    plot(l_tongue_x(inds_(counter_frame)),l_tongue_y(inds_(counter_frame)),'og','markersize', 7);
    plot(mid_tongue_x(inds_(counter_frame)),mid_tongue_y(inds_(counter_frame)),'ob','markersize', 7);
    plot(r_nose_x(inds_(counter_frame)),r_nose_y(inds_(counter_frame)),'ok','markersize', 5);
    plot(l_nose_x(inds_(counter_frame)),l_nose_y(inds_(counter_frame)),'ok','markersize', 5);
    xlabel('x (mm)');
    ylabel('y (mm)');
    xlim([-10 20])
    ylim([-20 20])
    
    set(gca, 'YDir','reverse');
    saveas(gcf, ['Geo_' num2str(counter_frame)], 'png');
    close(gcf);
end
fprintf(' --> Completed. \n')

%% Grab traces from DLC
for counter_frame = 1:1:length(inds_)
    figure;
    hold on;
    plot(tip_tongue_x(inds_(1):inds_(counter_frame)),tip_tongue_y(inds_(1):inds_(counter_frame)),'.-', 'Linewidth', 0.8);
    plot(mid_tongue_x(inds_(1):inds_(counter_frame)),mid_tongue_y(inds_(1):inds_(counter_frame)),'.-', 'Linewidth', 0.8);
    plot(r_tongue_x(inds_(1):inds_(counter_frame)),r_tongue_y(inds_(1):inds_(counter_frame)),'.-', 'Linewidth', 0.8);
    plot(l_tongue_x(inds_(1):inds_(counter_frame)),l_tongue_y(inds_(1):inds_(counter_frame)),'.-', 'Linewidth', 0.8);
    plot(r_tube_r_x(inds_(counter_frame)),r_tube_r_y(inds_(counter_frame)),'sk','markersize', 5);
    plot(r_tube_l_x(inds_(counter_frame)),r_tube_l_y(inds_(counter_frame)),'sk','markersize', 5);
    plot(l_tube_r_x(inds_(counter_frame)),l_tube_r_y(inds_(counter_frame)),'sk','markersize', 5);
    plot(l_tube_l_x(inds_(counter_frame)),l_tube_l_y(inds_(counter_frame)),'sk','markersize', 5);
    plot(r_nose_x(inds_(counter_frame)),r_nose_y(inds_(counter_frame)),'ok','markersize', 5);
    plot(l_nose_x(inds_(counter_frame)),l_nose_y(inds_(counter_frame)),'ok','markersize', 5);
    xlabel('x (mm)');
    ylabel('y (mm)');
    set(gca, 'YDir','reverse')
    xlim([-10 20])
    ylim([-20 20])
    saveas(gcf, ['Trace_' num2str(counter_frame)], 'png');
    close(gcf);
end

%% Grab distance from DLC
for counter_frame = 1:1:length(inds_)  
    figure;
    hold on;
    plot((time_vid(inds_(1):inds_(counter_frame)) - time_vid(inds_(1)))*1000, d_tip(inds_(1):inds_(counter_frame)),'.-', 'Linewidth', 0.8);
    xlabel('Time (ms)');
    ylabel('Displacement (mm)');
    ylim([0 25])
    xlim([0 x_lim])
    saveas(gcf, ['Displacement_' num2str(counter_frame)], 'png');
        close(gcf);

end

%% Grab velocity from DLC

for counter_frame = 1:1:length(inds_v)   
    figure;
    hold on;
    plot((time_vid(inds_v(1):inds_v(counter_frame)) - time_vid(inds_v(1)))*1000, v_tip(inds_v(1):inds_v(counter_frame)),'.-', 'Linewidth', 0.8);
    xlabel('Time (ms)');
    ylabel('Velocity (mms)');
    ylim([-600 600])
    xlim([0 x_lim])
    saveas(gcf, ['Velocity_' num2str(counter_frame)], 'png');
        close(gcf);

end

%% Grab distance_food from DLC
for counter_frame = 1:1:length(inds_)
    figure;
    hold on;
    plot((time_vid(inds_(1):inds_(counter_frame)) - time_vid(inds_(1)))*1000, d_tip_l_food(inds_(1):inds_(counter_frame)),'.-', 'Linewidth', 0.8);
    xlabel('Time (ms)');
    ylabel('Distance (mm)');
    ylim([0 25])
    xlim([0 x_lim])
    saveas(gcf, ['Food_' num2str(counter_frame)], 'png');
        close(gcf);

end

%% Grab angle of lick
for counter_frame = 1:1:length(inds_)
    figure;
    hold on;
    plot((time_vid(inds_(1):inds_(counter_frame)) - time_vid(inds_(1)))*1000, angle_midtip(inds_(1):inds_(counter_frame)),'.-', 'Linewidth', 0.8);
    xlabel('Time (ms)');
    ylabel('Angle (deg)');
    ylim([-100 100])
    xlim([0 x_lim])
    saveas(gcf, ['Angle_' num2str(counter_frame)], 'png');
    close(gcf);

end
end
