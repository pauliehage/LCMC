function PGH_run_python_script(path_to_DLC_network, path_to_video, video_name)
conda.setenv('dlc')
if isempty(vid_name)
    command = sprintf('ipython run_DLC.py %s %s', path_to_DLC_network, path_to_video);
else
    command = sprintf('ipython run_DLC.py %s %s %s', path_to_DLC_network, path_to_video, video_name);
end
[~, cmdout] = system(command, '-echo');
end

