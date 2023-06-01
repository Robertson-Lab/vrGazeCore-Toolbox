import os, sys
import argparse

def get_args_parser():

	parser = argparse.ArgumentParser('vrGazeCore - Parameters', add_help=False)

	# Set top-level directory structure
	parser.add_argument('--project_dir', type=str, help='Absolute path to directory of current project')
	parser.add_argument('--raw_data_folder', type=str, help='Name or path of raw data folder')
	parser.add_argument('--stim_folder', type=str, help='Name or path of stimuli folder')

	# Processing Options
	parser.add_argument('--run_find_fix', action='store_true', default=False, help="Run find fixations for all individuals")
	parser.add_argument('--run_heatmapping_individ', action='store_true', default=False, help="Run individual heatmapping for each subject")
	parser.add_argument('--run_heatmapping_group', action='store_true', default=False, help="Run group heatmapping for the cohort")
	parser.add_argument('--run_timecourse_gif_individ', action='store_true', default=False, help="Make individual heatmapping timecourse gifs for each subject")
	parser.add_argument('--run_timecourse_gif_group', action='store_true', default=False, help="Make individual heatmapping timecourse gifs for each subject")

	# Headset Type
	parser.add_argument('--headset_type', type=int, default=0, choices=[0, 1, 2, 3], help="Headset type (0=DK2, 1=Vive, 2=Vive Eye, 3=Oculus Go)")

	# Select Subjects
	parser.add_argument('--cohort_name', type=str, default='cohortName', help="Name for the group of subjects")
	parser.add_argument('--list_subject_names', action='store_true', default=False, help="List subject names manually")

	# Heatmapping Options
	parser.add_argument('--heatmap_timesteps', type=int, default=1, help="Number of time segments to divide scene into")

	# Additional Scene Parameters
	parser.add_argument('--min_samples', type=int, default=100, help="Minimum number of samples to consider a scene")
	parser.add_argument('--min_time_in_scene', type=int, default=5, help="Minimum time in seconds for a scene")
	parser.add_argument('--scanned_filter', action='store_true', default=False, help="Remove scenes where they were not explored")
	parser.add_argument('--scanned_thresh', type=int, default=66, help="Percent threshold to exclude scenes not explored based on head direction")

	# Manually Exclude Scenes
	parser.add_argument('--exclude_scenes', type=str, nargs='+', default=[], help="List of scenes to exclude")

	# Ignore List
	parser.add_argument('--ignore_list', type=str, nargs='+', default=['fixate*'], help="List of scenes to ignore")

	# Load Scenes
	parser.add_argument('--scene_dir', type=str, default='path_to_scenes_directory', help="Directory containing scene files")

	# Eyetracking options
	parser.add_argument('--gaze_type', type=int, default=0, choices=[0, 1], help="Gaze type (0=2D tracking, 1=3D)")

	# Fixations Options
	parser.add_argument('--exclude_first_n_sec', type=int, default=0, help="Exclude n seconds from the start of the trial")
	parser.add_argument('--min_mad', type=int, default=50, help="Minimum mean absolute deviation of velocity for potential fixations")
	parser.add_argument('--max_mad', type=int, default=100, help="Maximum mean absolute deviation of velocity for potential fixations")
	parser.add_argument('--exclude_fix_durs_less_than', type=float, default=0.1, help="Exclude fixations with a duration less than a certain threshold")

	# Eye Selection
	parser.add_argument('--use_eye', type=int, default=3, choices=[0, 1, 2, 3], help="Eye selection (0=eye0, 1=eye1, 2=choose best eye, 3=average eyes)")

	# Locked Head to Center
	parser.add_argument('--head_locked', action='store_true', default=False, help="Lock head direction to center (0=No, 1=Yes)")

	# Confidence Filter
	parser.add_argument('--min_conf_thresh', type=float, default=0.25, help="Minimum confidence threshold")
	parser.add_argument('--max_conf_percent', type=int, default=75, help="Maximum percentage of discarded data due to confidence")

	# Eccentricity Filter
	parser.add_argument('--ecc_filt_x', type=int, default=100, help="Maximum eccentricity in the x dimension")
	parser.add_argument('--ecc_filt_y', type=int, default=100, help="Maximum eccentricity in the y dimension")

	# Smoothing
	parser.add_argument('--use_smoothing', action='store_true', default=False, help="Use smoothing (0=No, 1=Yes)")
	parser.add_argument('--use_interpolation', action='store_true', default=False,  help="Use linear interpolation (0=No, 1=Yes)")
	parser.add_argument('--duration_for_interpolation', type=float, default=0.10, help="Duration for interpolation")

	# Fixation Calculation
	parser.add_argument('--fix_type', type=int, default=1, choices=[1, 2], help="Fixation type (1=gaze fixation, 2=head fixation)")
	parser.add_argument('--fix_spatial_dist', type=int, default=2, help="Spatial distance for fixation grouping")
	parser.add_argument('--fix_temp_dist', type=float, default=0.15, help="Temporal distance for fixation grouping")
	parser.add_argument('--fix_validation', action='store_true', default=False, help="Fixation validation (0=No, 1=Yes)")
	parser.add_argument('--val_window', type=int, default=2, help="Number of MAD samples for fixation validation")
	parser.add_argument('--exclude_last_n_sec', type=int, default=0, help="Exclude n seconds from the end of the trial")
	parser.add_argument('--exclude_first_n_fix', type=int, default=1, help="Exclude n fixations from the start of the trial")
	parser.add_argument('--exclude_last_n_fix', type=int, default=0, help="Exclude n fixations from the end of the trial")

	# Drift Correction
	parser.add_argument('--drift_correction', action='store_true', default=False, help="Perform drift correction (0=No, 1=Yes)")

	# Pre-Trial Fixations
	parser.add_argument('--concat_sanity', action='store_true', default=False, help="Set Equator scenes to the correct name for plotting")
	parser.add_argument('--avg_pre_trial', action='store_true', default=False, help="Calculate pretrial fixation distance and threshold")
	parser.add_argument('--avg_pre_trial_thresh', type=float, default=5, help="Tolerance for determining if a scene should be drift corrected or skipped")
	parser.add_argument('--exclude_by_pre_trial',action='store_true', default=False, help="Exclude the next scene if the previous pre-trial was bad")
	parser.add_argument('--avg_pre_trial_x', type=float, default=-1.26, help="True x coordinate of stimulus")
	parser.add_argument('--avg_pre_trial_y', type=float, default=-0.54, help="True y coordinate of stimulus")

	return parser


def set_paths(args):

	# gaze_core_dir = os.path.join(args.project_dir, 'vrGazeCore')

	# Set stimuli and data folder to use
	project_raw_data_dir = os.path.join(args.project_dir, args.raw_data_folder)
	project_stim_dir = os.path.join(args.project_dir, args.stim_folder)

	# Based on the top-level folders, populate path name variables
	# Set analysis results folders
	project_data_dir = os.path.join(args.project_dir, 'eyeTrackResults')
	project_fix_data_dir = os.path.join(project_data_dir, 'fixations')
	project_fix_mat_dir = os.path.join(project_fix_data_dir, 'npy')
	project_fix_plots_dir = os.path.join(project_fix_data_dir, 'plots')
	project_heat_dir = os.path.join(project_data_dir, 'heatMaps')
	project_heat_mat_dir = os.path.join(project_heat_dir, 'npy')
	project_heat_plots_dir = os.path.join(project_heat_dir, 'plots')
	project_timecourse_plot_dir = os.path.join(project_data_dir, 'timecourseHeat')

	# Set meta and log folders
	project_anal_logs_dir = os.path.join(args.project_dir, 'eyeTrackLogs')
	project_meta_data_dir = os.path.join(project_anal_logs_dir, 'meta')
	project_logs_dir = os.path.join(project_anal_logs_dir, 'logs')

	# Print check about paths that need to be changed
	print("Check that the following paths are correct:")
	# print(f"gaze_core_dir = {gaze_core_dir}")
	print(f"project_raw_data_dir = {project_raw_data_dir}")
	print(f"project_stim_dir = {project_stim_dir}")

	# If correct, input 1
	check_path = int(input("Are the paths correct?\n1 = Yes\n2 = No\nEnter:"))
	if check_path == 2:
		# break script, prompt to go back
		print("batch_process_pars will stop running. Go to set_paths to change the incorrect paths.")
		raise ValueError("Paths are not correct.")

	# Dictionary to store paths
	paths = {
		'project_dir': args.project_dir,
		# 'gaze_core_dir': gaze_core_dir,
		'project_raw_data_dir': project_raw_data_dir,
		'project_stim_dir': project_stim_dir,
		'project_data_dir': project_data_dir,
		'project_fix_data_dir': project_fix_data_dir,
		'project_fix_mat_dir': project_fix_mat_dir,
		'project_fix_plots_dir': project_fix_plots_dir,
		'project_heat_dir': project_heat_dir,
		'project_heat_mat_dir': project_heat_mat_dir,
		'project_heat_plots_dir': project_heat_plots_dir,
		'project_timecourse_plot_dir': project_timecourse_plot_dir,
		'project_anal_logs_dir': project_anal_logs_dir,
		'project_meta_data_dir': project_meta_data_dir,
		'project_logs_dir': project_logs_dir,
	}

	# Check if necessary directories exist; if not, make them.
	for path in paths.values():
		if not os.path.exists(path):
			os.makedirs(path)

	return paths
