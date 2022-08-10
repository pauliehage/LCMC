import os
import sys
import deeplabcut

## Helper function to get subfolders
def getsubfolders(folder):
	return [os.path.join(folder,p) for p in os.listdir(folder) if os.path.isdir(os.path.join(folder,p))]

## Helper function to get video
def getvideos(path):
	return [os.path.join(path, p) for p in os.listdir(path) if p.endswith('.mp4')]

## Main function
def main():
	path_to_DLC_network= sys.argv[1]
	path_to_video= sys.argv[2]
	shuffle = 8

	config = os.path.join(path_to_DLC_network, 'config.yaml')

	# run (to be updated with new folder structures)
	if len(sys.argv) == 3:
		# only process specified video
		video_name = sys.argv[3]
		file = os.path.join(path_to_video, video_name)
		deeplabcut.analyze_videos(config, [file], shuffle=shuffle, videotype='mp4', save_as_csv=True)
		deeplabcut.filterpredictions(config, [file], shuffle=shuffle, videotype='mp4')
		deeplabcut.create_labeled_video(config, [file], shuffle=shuffle, videotype='mp4', filtered=True)

	else:
		files = getvideos(path_to_video)
		for file in files:
			deeplabcut.analyze_videos(config, [file], shuffle=shuffle, videotype='mp4', save_as_csv=True)
			deeplabcut.filterpredictions(config, [file], shuffle=shuffle, videotype='mp4')
			deeplabcut.create_labeled_video(config, [file], shuffle=shuffle, videotype='mp4', filtered=True)

if __name__ == "__main__": main()
