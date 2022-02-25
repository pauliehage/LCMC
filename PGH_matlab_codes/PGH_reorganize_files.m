%% function PGH_reorganize files
function PGH_reorganize_files(rec_path, params, funcs)

path_to_analyzed_data = [rec_path , 'analyzed_data', filesep];
path_to_sorted_data = [path_to_analyzed_data,'sorted_data',filesep];
path_to_analyzed_eye = [path_to_analyzed_data, 'behavior_data', ...
    filesep, 'eye', filesep];
path_to_analyzed_tongue = [path_to_analyzed_data, 'behavior_data', ...
    filesep, 'tongue', filesep];
path_to_raw_data = [rec_path , 'raw_data', filesep];

mkdir(path_to_analyzed_eye);
mkdir(path_to_sorted_data);
mkdir(path_to_analyzed_tongue);

% delete folders from old structure
delete_fig_folders = 0;
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

% delete reevaluated .csv from raw files
try
       delete([path_to_raw_data,'*_Revaluated*']);

end

% delete existing converted DLC mat file from analyzed tongue
try
    delete([path_to_analyzed_tongue,'*DLC.mat']);
end

% delete existing _video metadata file
delete_video_metadata = 1;
if delete_video_metadata == 1
    try
        delete([path_to_analyzed_tongue,'*_video.mat']);
    end
    try
        delete([path_to_analyzed_tongue,'*_aligned.mat']);
    end
end

end

