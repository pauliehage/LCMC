%% function PGH_reorganize files
function PGH_reorganize_files(rec_path, params, funcs)
% local
path_to_analyzed_data = [rec_path , 'analyzed_data', filesep];
path_to_sorted_data = [path_to_analyzed_data,'sorted_data',filesep];
path_to_analyzed_eye = [path_to_analyzed_data, 'behavior_data', ...
    filesep, 'eye', filesep];
path_to_analyzed_tongue = [path_to_analyzed_data, 'behavior_data', ...
    filesep, 'tongue', filesep];
path_to_raw_data = [rec_path , 'raw_data', filesep];


%% reorganize old structure
is_reorganize_old_structure = 0;
if is_reorganize_old_structure == 1

    mkdir(path_to_analyzed_eye);
    mkdir(path_to_sorted_data);
    mkdir(path_to_analyzed_tongue);

    % delete folders from old structure
    delete_fig_folders = 1;
    if delete_fig_folders == 1
        try
            rmdir([rec_path, 'analyzed_figs'],'s');
            rmdir([rec_path, 'bundle_figs'],'s');
        end
    end

    % delete old mat files from old structure
    try
        delete([path_to_analyzed_data,'*.mat']);
    end

    % delete old pdf files from old structure
    try
        delete([path_to_analyzed_data,'*.pdf']);
    end

    % move sorted files
    try
        movefile([path_to_analyzed_data, '*.psort'],path_to_sorted_data);
    end

    % move h5 files
    try
        movefile([path_to_analyzed_data, '*.h5'],path_to_sorted_data);
    end

    % move labeled video files
    try
        labeled_video_file = dir([path_to_analyzed_tongue,'*_labeled.mp4']);
        labeled_video_file_name = labeled_video_file.name;
        labeled_video_file_name = labeled_video_file_name(1:13);
        movefile([path_to_analyzed_tongue, labeled_video_file_name '_labeled.mp4'],[path_to_analyzed_tongue,labeled_video_file_name '_DLC.mp4' ]);
    end
    try
        labeled_video_file = dir([path_to_raw_data,'*_labeled.mp4']);
        labeled_video_file_name = labeled_video_file.name;
        labeled_video_file_name = labeled_video_file_name(1:13);
        movefile([path_to_raw_data, labeled_video_file_name '_labeled.mp4'],[path_to_analyzed_tongue,labeled_video_file_name '_DLC.mp4' ]);
    end
    try
        labeled_video_file = dir([path_to_raw_data,'*_DLC.mp4']);
        labeled_video_file_name = labeled_video_file.name;
        labeled_video_file_name = labeled_video_file_name(1:13);
        movefile([path_to_raw_data, labeled_video_file_name '_DLC.mp4'],[path_to_analyzed_tongue,labeled_video_file_name '_DLC.mp4' ]);
    end
    try
        converted_video_file = dir([path_to_raw_data,'*_converted.mp4']);
        converted_video_file_name = converted_video_file.name;
        converted_video_file_name = converted_video_file_name(1:13);
        movefile([path_to_raw_data, converted_video_file_name '_converted.mp4'],[path_to_raw_data,converted_video_file_name '_.mp4' ]);
    end

    % move DLC mp4 and csv file from recording folder to analyzed tongue
    try
        DLC_video_file = dir([rec_path,'*_DLC.mp4']);
        DLC_video_file_name = DLC_video_file.name;
        DLC_video_file_name = DLC_video_file_name(1:13);
        movefile([rec_path, DLC_video_file_name '_DLC.mp4'],[path_to_analyzed_tongue,DLC_video_file_name '_DLC.mp4' ]);
    end
    try
        DLC_csv_file = dir([rec_path,'*_DLC.csv']);
        DLC_csv_file_name = DLC_csv_file.name;
        DLC_csv_file_name = DLC_csv_file_name(1:13);
        movefile([rec_path, DLC_csv_file_name '_DLC.csv'],[path_to_analyzed_tongue,DLC_csv_file_name '_DLC.csv' ]);
    end

    % move labeled csv files
    try
        labeled_csv_file = dir([path_to_analyzed_tongue,'*_labeled.csv']);
        labeled_csv_file_name = labeled_csv_file.name;
        labeled_csv_file_name = labeled_csv_file_name(1:13);
        movefile([path_to_analyzed_tongue, labeled_csv_file_name '_labeled.csv'],[path_to_analyzed_tongue,labeled_csv_file_name '_DLC.csv' ]);
    end
    try
        labeled_csv_file = dir([path_to_raw_data,'*_labeled.csv']);
        labeled_csv_file_name = labeled_csv_file.name;
        labeled_csv_file_name = labeled_csv_file_name(1:13);
        movefile([path_to_raw_data, labeled_csv_file_name '_labeled.csv'],[path_to_analyzed_tongue,labeled_csv_file_name '_DLC.csv' ]);
    end
    try
        labeled_csv_file = dir([path_to_raw_data,'*_DLC.csv']);
        labeled_csv_file_name = labeled_csv_file.name;
        labeled_csv_file_name = labeled_csv_file_name(1:13);
        movefile([path_to_raw_data, labeled_csv_file_name '_DLC.csv'],[path_to_analyzed_tongue,labeled_csv_file_name '_DLC.csv' ]);
    end

    % fix raw video names
    try
        temp_avi_file = dir([path_to_raw_data,'*.avi']);
        temp_avi_file_name = temp_avi_file.name;
        if length(temp_avi_file_name) >= 24
            temp_avi_file_name_new = [temp_avi_file_name(12:13) temp_avi_file_name(6:9) '_' temp_avi_file_name(14:19) ];
            movefile([path_to_raw_data, temp_avi_file_name],[path_to_raw_data,temp_avi_file_name_new '.avi' ]);
        else
            avi_file = dir([path_to_raw_data,'*.avi']);
            avi_file_name = avi_file.name;
            avi_file_name_new = avi_file_name(1:13);
            movefile([path_to_raw_data, avi_file_name],[path_to_raw_data,avi_file_name_new '.avi' ]);
        end
    end
    try
        temp_mp4_file = dir([path_to_raw_data,'*.mp4']);
        temp_mp4_file_name = temp_mp4_file.name;
        if length(temp_mp4_file_name) >= 24
            temp_mp4_file_name_new = [temp_mp4_file_name(12:13) temp_mp4_file_name(6:9) '_' temp_mp4_file_name(14:19) ];
            movefile([path_to_raw_data, temp_mp4_file_name],[path_to_raw_data,temp_mp4_file_name_new '.mp4' ]);
        else
            mp4_file = dir([path_to_raw_data,'*.mp4']);
            mp4_file_name = mp4_file.name;
            mp4_file_name_new = mp4_file_name(1:13);
            movefile([path_to_raw_data, mp4_file_name],[path_to_raw_data, mp4_file_name_new, '.mp4']);
        end
    end

    % delete reevaluated .csv from raw files
    try
        delete([path_to_raw_data,'*_Revaluated*']);
    end

    % delete existing converted DLC mat file from analyzed tongue
    delete_DLC_mat = 0;
    if delete_DLC_mat == 1
        try
            delete([path_to_analyzed_tongue,'*DLC.mat']);
        end
    end

    % delete existing _video metadata file
    delete_video_metadata = 0;
    if delete_video_metadata == 1
        try
            delete([path_to_analyzed_tongue,'*_video.mat']);
        end
        try
            delete([path_to_analyzed_tongue,'*_aligned.mat']);
        end
    end

    % delete existing files in analyzed eye
    delete_eye_mat = 0;
    if delete_eye_mat == 1
        try
            delete([path_to_analyzed_eye,'*.mat']);
        end
    end

end
%% move lick data to server
is_move_lick_to_server = 1;
if is_move_lick_to_server == 1
    % server
    rec_path_split = regexp(rec_path,'\','split');
    monkey_id = char(rec_path_split(2));
    month = char(rec_path_split(3));
    day = char(rec_path_split(4));
    rec = char(rec_path_split(5));

    rec_path_server = ['X:' filesep 'Ephys' filesep monkey_id filesep month filesep day filesep rec filesep];

    path_to_analyzed_data_server = [rec_path_server , 'analyzed_data', filesep];
    path_to_analyzed_tongue_server = [path_to_analyzed_data_server, 'behavior_data', ...
        filesep, 'tongue', filesep];
    path_to_raw_data_server = [rec_path_server , 'raw_data', filesep];

    % delete temp video
    delete_temp_vid = 1;
    if delete_temp_vid == 1
        try
            delete([path_to_raw_data_server,'*.avi']);
        end
    end

    % copy raw videos: avi and mp4
    move_raw_video_file = 1 ;
    if move_raw_video_file == 1
        try
            avi_file = dir([path_to_raw_data,'*.avi']);
            avi_file_name = avi_file.name(1:13);
            copyfile([path_to_raw_data, avi_file_name '.avi'],[path_to_raw_data_server,avi_file_name '.avi' ]);

            mp4_file = dir([path_to_raw_data,'*.mp4']);
            mp4_file_name = mp4_file.name(1:13);
            copyfile([path_to_raw_data, mp4_file_name '.mp4'],[path_to_raw_data_server,mp4_file_name '.mp4' ]);
        end

    end

    % copy tongue analyzed data
    copy_analyzed_tongue = 1 ;
    if copy_analyzed_tongue == 1
        try
            copyfile(path_to_analyzed_tongue(1:end-1),path_to_analyzed_tongue_server(1:end-1));
        end
    end


end

end


