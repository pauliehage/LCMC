%% Process videos with DLC
clear all
conda.setenv('dlc')

% set parameters to be passed to python script
% indicate which monkey (where the config file is located)
subj_folder = "/home/ijang/2_Subjects_Tongue/59d_Post_COVID-ijang-2020-09-30";
% path to raw video
vid_path = subj_folder; % just as e.g.
% video name (optional; will process all videdos in vid_path if not
% provided)
vid_name = "210503_160305_converted.mp4";
if isempty(vid_name)
    command = sprintf('ipython script.py %s %s', subj_folder, vid_path);
else
    command = sprintf('ipython script.py %s %s %s', subj_folder, vid_path, vid_name);
end
[~, cmdout] = system(command, '-echo');

%% rename files produced
addpath(vid_path)
files = dir(strcat(vid_path,'/*filtered.csv'));
for id = 1:length(files)
    [~, f] = fileparts(files(id).name);
    p2 = f(1:13);
    newname = p2 + "_labeled.csv"; 
    movefile(strcat(vid_path, '/', files(id).name), strcat(vid_path, '/', newname));
end
%%
files2 = dir(strcat(vid_path,'/*labeled.mp4'));
for id = 1:length(files2)
    [~, f] = fileparts(files2(id).name);
    p2 = f(1:13);
    newname = p2 + "_labeled.mp4";
    movefile(strcat(vid_path, '/', files2(id).name), strcat(vid_path, '/', newname));
end