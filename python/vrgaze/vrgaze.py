import os
import glob
import pickle
import pandas as pd
import numpy as np
from datetime import datetime as dt
import cv2

from matplotlib import pyplot as plt

from PIL import Image

from sklearn.metrics.pairwise import haversine_distances
from vrgaze.utils import sliding_window_mad, get_fixation_indices, sphere_centroid, calculate_fixation_centroids, concatenate_fixations, degrees_to_pixels, scale_durations, apply_gaussian_smoothing

class vrGazeData(object):
    
    def __init__(self, subject, trial_name, trial_number, df, params, paths):
        
        self.subject = subject
        self.trial_name = trial_name
        self.trial_number = trial_number
        self.raw_data = df
        self.confidence_filter = np.zeros(df.shape[0]) if df is not None else None
        self.eccentricity_filter = np.zeros(df.shape[0]) if df is not None else None
        self.preprocessed_data = None
        self.filtered_data = None
        self.fixations = None
        self.density_maps = None
        
        self.subject_list = None
        self.params = params
        self.paths = paths
        
        self.image_path = self.set_image_path()
    
    def set_image_path(self):
        
        image_path = glob.glob(os.path.join(self.paths['project_stim_dir'], f"{self.trial_name}*"))
        
        if image_path:
            return image_path[0]
        else:
            return None

    def get_image_path(self):
        return self.image_path

    def get_raw_data(self):
        if self.raw_data is None:
            print (f'Raw dataframe has not been created!')
            return
        return self.raw_data.copy()

    def set_confidence_filter(self, indices):
        self.confidence_filter = indices

    def get_confidence_filter(self):
        return self.confidence_filter.copy()

    def set_eccentricity_filter(self, indices):
        self.eccentricity_filter = indices

    def get_eccentricity_filter(self):
        return self.eccentricity_filter.copy()

    def set_preprocessed_data(self, df):
        self.preprocessed_data = df

    def get_preprocessed_data(self):
        if self.preprocessed_data is None:
            print (f'Preprocessed dataframe has not been created!')
            return
        return self.preprocessed_data.copy()

    def set_fixations(self, df):
        self.fixations = df

    def get_fixations(self):
        if self.fixations is None:
            print (f'Fixations dataframe has not been created!')
            return
        return self.fixations.copy()

    def set_density_maps(self, array):
        self.density_map = array

    def get_density_map(self):
        if self.density_map is None:
            print (f'Density map has not been created!')
            return
        return self.density_map.copy()

    def set_subject_list(self, subject_list):
        self.subject_list = subject_list

    def get_subject_list(self):
        if self.density_map is None:
            print (f'Subject list has not been set! This usually means you only have one subject')
            return
        return self.subject_list

    def get_filters(self):
        # a more flexible implementation would be to add any given filter together then turn to bool
        filters = ~np.logical_or(self.confidence_filter, self.eccentricity_filter)
        return filters

    def write_data(self, out_dir, time_step = 0):
        
        if time_step > 0:
            out_dir = os.path.join(out_dir, self.subject, f'{time_step}-timeSegments')
        else:
            out_dir = os.path.join(out_dir, self.subject)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        out_fn = os.path.join(out_dir, f'{self.subject}_{self.trial_name}_{str(self.trial_number).zfill(5)}.pkl')

        with open(out_fn, 'wb') as f:
            pickle.dump(self, f)

        return out_fn

    def load_data(path):
        with open(path, 'rb') as f:
            self = pickle.load(f)
        return self

class vrGazeCore:

    def __init__(self, params, paths):
        # initialize the data class using 
        self.params = self.setFixedParams(params)
        self.paths = paths
        
    def setFixedParams(self, params):
        """
        Params is the result of an argparser. We need to modify as a dictionary to add 
        a few things.
        
        Equivalent to fixed params setting of setParams.m
        """
        
        d = vars(params)
        
        # DK2 + Pupil Labs
        if params.headset_type == 0:
            d['fov_x'] = 90
            d['fov_y'] = 100
        # Vive + Pupil Labs
        elif params.headset_type == 1:
            d['fov_x'] = 145
            d['fov_y'] = 160
        # Vive Eye
        elif params.headset_type == 2:
            d['fov_x'] = 75
            d['fov_y'] = 90
        # Oculus Go / Quest --> Headtracking only
        elif params.headset_type == 3:
            params.fix_type = 2
            params.use_smoothing = False
            params.use_interpolation = False
            params.avg_pre_trial = False
            d['fov_x'] = 100
            d['fov_y'] = 100
            
        if d['fov_x'] > d['fov_y']:
            d['max_fov'] = d['fov_x'] 
        else:
            d['max_fov'] = d['fov_y'] 
            
        return params

    def loadRawData(self, filename, header=None):
        
        # Set the path to the file based on our paths dictionary and load using pandas
        file_path = os.path.join(self.paths['project_raw_data_dir'], filename)
        raw_data = pd.read_csv(file_path, sep = ",", header=header)
        
        # Add default header labels #TODO: parameterize
        if header is None:
            column_labels = [
                'trial', 
                'data', 
                'core_time', 
                'exp_time', 
                'pitch', 
                'yaw', 
                'roll', 
                'right_x', 
                'right_y', 
                'left_x', 
                'left_y', 
                'right_conf', 
                'left_conf', 
                # 'rotation'
            ]
            
            raw_data.columns = column_labels
            
        # TODO: parameterize
        # make all sanity trials the same 
        sanity_trials = raw_data['trial'].str.contains('sanityTarget360')
        raw_data.loc[sanity_trials, 'trial'] = '_sanityTarget360_0000'
        
        # Set the core
        # raw_data['CoreTime'] = raw_data['CoreTime'].apply(lambda x: dt.strptime(x.strip(), '%H:%M:%S.%f').microsecond)
        
        return raw_data
    
    def processRawData(self, raw_data):
        
        # these are the basic columns we need for head fixations
        columns = ['trial', 'exp_time', 'yaw', 'pitch', 'roll']
        
        #### Head-Tracking Data ####
        if self.params.headset_type == 3:
            return raw_data[columns]
        
        #### Eye-Tracking Data ####
        
        # Right: eye coordinates + confidence
        if self.params.use_eye == 0:
            print ('Using right eye coordinates!')
            columns.extend(['right_x', 'right_y', 'conf'])
            raw_data = raw_data[columns].rename({'right_x': 'eye_x', 'right_y': 'eye_y', 'right_conf': 'conf'})
            return raw_data
        
        # Left: eye coordinates + confidence
        elif self.params.use_eye == 1:
            print ('Using left eye coordinates!')
            columns.extend(['left_x', 'left_y', 'left_conf'])
            raw_data = raw_data[columns].rename({'left_x': 'eye_x', 'left_y': 'eye_y', 'left_conf': 'conf'})
            return raw_data
        
        # Best: left or right eye coordinates based on best confidence
        elif self.params.use_eye == 2:
            raise NotImplementedError('Haven\'t implemented best confidence eye-choice') 
#             raw_data['Conf'] = raw_data[['RightConf', 'LeftConf']].agg(np.max, axis=1)
        # Average: average of left/right coordinates + minimum of confidence
        else:
            print ('Averaging coordinates from both eyes!')
            raw_data['eye_x'] = raw_data[['right_x', 'left_x']].agg(np.mean, axis=1)
            raw_data['eye_y'] = raw_data[['right_y', 'left_y']].agg(np.mean, axis=1)
            raw_data['conf'] = raw_data[['right_conf', 'left_conf']].agg(np.min, axis=1)
            columns.extend(['eye_x', 'eye_y', 'conf'])
            return raw_data[columns]
            
    def parseTrials(self, data, subject, use_dataframe=False):
        '''
        Takes a large dataframe and separates into a dataframe for each trial. 
        
        '''
        
        #create a row comparison for the dataframe, then split into trials based on row comparison
        def comp_prev(trial):
            return np.concatenate(([False], trial[1:] == trial[:-1]))

        trial_changes = comp_prev(data['trial'].values)

        #split into multiple dataframes based on trial
        parsed_data = np.split(data, np.where(trial_changes == 0)[0])[1:]

        if use_dataframe:
            parsed_data = [df.reset_index(drop=True) for df in parsed_data] 
            return parsed_data
        else:
            parsed_data = [vrGazeData(
                subject=subject,
                trial_name=sorted(df['trial'].unique())[0],
                trial_number=i+1, 
                df=df.reset_index(drop=True), 
                params=self.params, 
                paths=self.paths) for i, df in enumerate(parsed_data)]
            return parsed_data

    def parsedDataKey(self, parsed_data):
        parsed_data_key = pd.DataFrame()
        scene_idx = 0
        for trials in parsed_data:
            if hasattr(trials, '__iter__'):  # for group data
                scene_data_key = pd.DataFrame({
                            'scene_idx': scene_idx,
                            'subject_idx': [scene[0] for scene in enumerate(trials)],
                            'subject': [scene[1].subject for scene in enumerate(trials)],
                            'trial_name': [scene[1].trial_name for scene in enumerate(trials)],
                            'trial_number': [scene[1].trial_number for scene in enumerate(trials)]
                        })
                parsed_data_key = pd.concat([parsed_data_key, scene_data_key], ignore_index = True)
                scene_idx = scene_idx + 1
            else:  # for a single subject
                parsed_data_key = pd.DataFrame({
                    'scene_idx': [i for i in range(len(parsed_data))],
                    'subject': [trial.subject for trial in parsed_data],
                    'trial_name': [trial.trial_name for trial in parsed_data],
                    'trial_number': [trial.trial_number for trial in parsed_data]
                })
  
        return parsed_data_key
    
    def runFindFixations(self, trial, use_dataframe=False):
        """
        Runs vrGazeCore pipeline over a vrGazeData class
        """
        # retrieve raw data from the class
        data = trial.get_raw_data()

        ### RUN FILTERS ####
        # First run confidence filter
        confidence_filter, percent_removed = self.confidenceFilter(data)
        trial.set_confidence_filter(confidence_filter)

        # Then run eccentricity filter
        eccentricity_filter, percent_removed = self.eccentricityFilter(data)
        trial.set_eccentricity_filter(eccentricity_filter)

        # Map the coordinates on screen to FOV and roll to radians
        trial_filters = trial.get_filters()
        preprocessed_data = self.preprocessTrialData(data, filters=trial_filters)
        trial.set_preprocessed_data(preprocessed_data)

        trial_fixations = self.calculateFixations(preprocessed_data)
        trial.set_fixations(trial_fixations)

        if self.params.plot_fixations:
            out_dir = os.path.join(self.paths['project_fix_plots_dir'], trial.subject)

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            
            self.plotFixations(
                df_fixations=trial.get_fixations(), 
                image_path=trial.get_image_path(),
                out_path=os.path.join(out_dir, f'{trial.subject}_{trial.trial_name}_{str(trial.trial_number).zfill(5)}.png')
            )

        trial.write_data(out_dir=self.paths['project_fix_pkl_dir'])
        return trial

    def loadGroupFixations(self, subjects): #TODO: make it be able to run scenes not shared across all subjects
        # grab all filenames that are shared across subjects

        data = []
        all_stimuli = sorted(glob.glob(os.path.join(self.paths['project_stim_dir'], '*')))

        for stim in all_stimuli:
            stim_name = os.path.splitext(os.path.basename(stim))[0]

            subject_fns = []
            for subject in subjects:
                subject_fn = glob.glob(os.path.join(self.paths['project_fix_pkl_dir'], subject, f'*{stim_name}*'))
                if subject_fn:
                    subject_fns.append(subject_fn[0])
            if len(subject_fns) == len(subjects):
                data.append([vrGazeData.load_data(fn)for fn in subject_fns])
            # else:
                # print (f'Subject fixations missing for {os.path.basename(stim)}')
        # now load fixations
        return data

    def createGroupFixations(self, trials):
        """
        Given a list of data, create a trial for the group
        """

        # create a dummy trial for the group --> this will be used to save the data
        trial_group = vrGazeData(
            subject=f'group-{self.params.cohort_name}',
            trial_name=trials[0].trial_name,
            trial_number=None,
            df=None,
            params=self.params,
            paths=self.paths,
        )

        df_fixations = pd.concat([t.get_fixations() for t in trials])
        subject_list = [t.subject for t in trials]
        
        trial_group.set_fixations(df_fixations)
        trial_group.set_subject_list(subject_list)
        return trial_group

    def runHeatmapping(self, trial):
        """
        If trials is a list, then we concatenate the dataframes
        """
        
        ## TLB --> GROUP HASN'T BEEN TESTED YET
        if isinstance(trial, list):
            # create an aggregated fixations list for the group
            trial = self.createGroupFixations(trial)
        
        # if pre-trial, skip fixation density mapping
        if any([item in trial.trial_name for item in self.params.pretrial_list]):
            print (f'Skipping fixation density mapping of {trial.trial_name} - not a scene!')
            return trial
        
        df_fixations = trial.get_fixations()
        image_path = trial.get_image_path()

        df_timesteps = self.splitFixationTimesteps(df_fixations)
        
        # if there are not timesteps with items in it
        if df_timesteps is None or not any(list(map(np.size, df_timesteps))):
            return trial
        
        density_maps = np.stack([self.createDensityMap(df) for df in df_timesteps])

        trial.set_density_maps(density_maps)
        
        # save data
        # trial.write_data(out_dir=self.paths['project_timecourse_pkl_dir'] if self.params.heatmap_timesteps > 1 else self.paths['project_heat_pkl_dir'], time_step=self.params.heatmap_timesteps)
        trial.write_data(out_dir=self.paths['project_heat_pkl_dir'], time_step=self.params.heatmap_timesteps)

        ## TLB FIX --> TRIAL NUMBER IS NOT SET AT GROUP LEVEL
        out_fn = f'{trial.subject}_{trial.trial_name}_{str(trial.trial_number).zfill(5)}'

        if self.params.plot_density_maps:

            # create the output directory for the subject or group
            # out_dir = self.paths['project_timecourse_plots_dir'] if self.params.heatmap_timesteps > 1 else self.paths['project_heat_plots_dir']
            if self.params.heatmap_timesteps > 1:
                out_dir = os.path.join(self.paths['project_heat_plots_dir'], trial.subject,f'{self.params.heatmap_timesteps}-timeSegments',trial.trial_name)
            elif self.params.heatmap_timesteps == 1:
                out_dir = os.path.join(self.paths['project_heat_plots_dir'], trial.subject,f'{self.params.heatmap_timesteps}-timeSegments')

            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
                
            timestep_bounds = np.linspace(0, self.params.scene_length, self.params.heatmap_timesteps + 1)
            
            for i, density_map in enumerate(density_maps):
                start_dur = timestep_bounds[i]
                end_dur = timestep_bounds[i+1]
                self.plotFixationDensity(
                    density_map=density_map, 
                    image_path=image_path,
                    start_dur=start_dur,
                    end_dur=end_dur,
                    out_path=os.path.join(out_dir, f'{out_fn}_{self.params.heatmap_timesteps}-timeSegments_step-{str(i+1).zfill(3)}.png') if len(density_maps) > 1 else os.path.join(out_dir, f'{out_fn}.png')
                )
    
        if self.params.make_density_map_gif and self.params.heatmap_timesteps > 1:
            density_map_dir = os.path.join(self.paths['project_heat_plots_dir'], trial.subject,f'{self.params.heatmap_timesteps}-timeSegments',trial.trial_name) 
            self.makeDensityMapGIF(
                plot_dir=density_map_dir,
                out_path=os.path.join(density_map_dir, f'{out_fn}_{self.params.heatmap_timesteps}-timeSegments.gif')
            )
            
        return trial

    def confidenceFilter(self, trial):
        """

        This function is used for filtering scene data by the confidence threshold. The default
        requires a confidence of 0.6 to keep a data point, otherwise the point will be removed.

        By default, the threshold will use the average of the two eyes confidence (eyeType=2). The 
        following options are supported:
            - eyeType=0: right-eye only
            - eyeType=1: left-eye only

        """


        df = trial.copy()

        # Find datapoints falling below threshold
        confidence_filter = df['conf'] <= self.params.min_conf_thresh # remove datapoint in the trial that fell below confidence threshold
        
        percent_removed = 100*sum(confidence_filter) / len(confidence_filter)
        print (f"CONFIDENCE FILTER - percent removed: {percent_removed:.2f}%")
            
        return confidence_filter, percent_removed


    def eccentricityFilter(self, trial):

        """

        Remove any datapoints with eye-tracking coords beyond a given degree.

        tlb ---> figure this out

        """

        df = trial.copy()

        # go from 100 to +/- 50, since 0,0 is in center of screen
        filter_x, filter_y = self.params.ecc_filt_x / 2 , self.params.ecc_filt_y / 2
        
        # Convert eye coordinates to degrees
        
        eye_x_degrees = df['eye_x'] - 0.5 * self.params.fov_x
        eye_y_degrees = df['eye_y'] - 0.5 * self.params.fov_y
        
        # find indices within threshold
        within_x = np.logical_and(eye_x_degrees < filter_x, eye_x_degrees > -1*filter_x)
        within_y = np.logical_and(eye_y_degrees < filter_y, eye_y_degrees > -1*filter_y)
        
        # find the intersection of the two and invert to make a filter
        eccentricity_filter = ~np.logical_and(within_x, within_y)
        
        # then calculate percent removed
        percent_removed = 100*sum(eccentricity_filter) / len(eccentricity_filter)
        
        print (f"ECCENTRICITY FILTER - percent removed: {percent_removed:.2f}%")
        
        return eccentricity_filter, percent_removed

    def preprocessTrialData(self, trial, filters=None):
        """
        Functional equivalent of the internals of parseDS.m
        
        Push yaw/pitch into positive domain (e.g., -90/90 to 0/180)
        Convert roll to radians
        Map eye coordinates in norm screen space to degrees
        
        """
        df = trial.copy()
        
        # No need to preprocess because no eye-tracking data involved
        if self.params.headset_type == 3:
            return df
        
        # Vive Eye 2D tracking needs to be normalized to same field as pupil
        if self.params.headset_type == 2 and self.params.gaze_type == 0:
            # bring positive then half
            df[['eye_x', 'eye_y']] = (df[['eye_x', 'eye_y']] + 1) / 2
        
        # Convert roll to radians
        df['roll'] = df['roll'].apply(np.deg2rad)
        
        ##### Shift to positive domain #####
        # Pitch: shift from -90 to 90 to 0 to 180
        # Yaw: shift from -180 to 180 to 0 to 360
        df['pitch'] = np.mod(df['pitch'] + 90, 180)
        df['yaw'] = np.mod(df['yaw'] + 180, 360)
        
        df = self.mapScreenToFOV(df)
        df = self.rectifyGaze(df)
        df = self.applyFilters(df, filters=filters)
        
        return df
        

    def mapScreenToFOV(self, trial):
        """
        
        Maps eye-coordinates on screen to FOV. Incorporates scaling for screen warp.
        
        Within MATLAB code, we begin referring to this as gaze here --> not entirely correct 
        
        tlb 9-02-19: useful line explaining the conversion below
        
        Ehringer (2019): we converted the x (and y) gaze points of the raw samples from screen coordinates in pixels, 
        to spherical angles in degree (with a reference system centered on the subject):
        
            Bx = 2 * atan2(px * m, d)
        
        where Bx denotes the azimuth angle (equivalent to the horizontal position) of the gaze points 
        in visual degrees from the monitor center
            - px: horizontal position relative to the center of the monitor in pixel, 
            - m: unit conversion of pixel to mm of the monitor, 
            - d: the distance to the monitor in mm.    
        
        """
        
        if self.params.headset_type == 0:
            # normalize center of the screen (eye data) to 0,0. because y
            # calibration coord is slightly below center of screen, shift a
            # bit more than X.
            
            # X,Y coordinates are presented on a plane ~0.8 units from the camera
            map_eye_x = lambda x: np.rad2deg((2*self.params.fov_x / self.params.max_fov) * np.arctan2(x - 0.5, 0.8))
            map_eye_y = lambda x: np.rad2deg((2*self.params.fov_y / self.params.max_fov) * np.arctan2(x - 0.505, 0.8))
        
        if self.params.headset_type == 2:
            # Vive Eye - 2D tracking
            if self.params.gaze_type == 1:
                map_eye_x = lambda x: np.rad2deg((self.params.max_fov / self.params.fov_x) * np.arcsin(-1*x))
                map_eye_y = lambda x: np.rad2deg((self.params.max_fov / self.params.fov_y) * np.arcsin(x + 0.005))
            
            # Vive Eye - 3D tracking --> still need to correct for distortion
            else:
                map_eye_x = lambda x: np.rad2deg((2*self.params.fov_x / self.params.max_fov) * np.arctan2(x - 0.5, 0.55))
                map_eye_y = lambda x: np.rad2deg((2*self.params.fov_y / self.params.max_fov) * np.arctan2(x - 0.5, 0.55))
        
        
        trial['eye_x_fov'] = map_eye_x(trial['eye_x'])
        trial['eye_y_fov'] = map_eye_y(trial['eye_y'])
        
        return trial
    
    def rectifyGaze(self, trial):
        """
        Rotates all gaze points based on roll, rectifying gaze data with head data
        """
        df = trial.copy()
        
        # move back -90 to 90 / -180 to 180
        pitch, yaw = df['pitch'] - 90, df['yaw'] - 180
        
        # Calculate change in pitch and yaw accounting for head rotation (roll)
        # dx = change in yaw coordinate
        # dy = change in pitch coordinate
        dx = df['eye_x_fov'] * np.cos(df['roll']) - df['eye_y_fov'] * np.sin(df['roll'])
        dy = df['eye_y_fov'] * np.cos(df['roll']) + df['eye_x_fov'] * np.sin(df['roll'])
        
        # Create gaze coordinates by shifting pitch and yaw
        # Scale yaw by pitch (operating in lat long)
        df['gaze_yaw'] = np.mod(yaw + (dx * np.cos(np.deg2rad(pitch))) + 180, 360)
        df['gaze_pitch'] = np.mod(pitch - dy + 90, 180)

        return df

    def applyFilters(self, trial, filters=None):
        if filters is not None:
            trial_filtered = trial.loc[filters]
            return trial_filtered.iloc[1:]
        else:
            return trial.iloc[1:]
    
    def calculateFixations(self, trial):
    
        """
        The meat of the code
        """
        df = trial.copy()

        # get gaze x and gaze y on sphere --> move back into radians
        coordinates = df[['gaze_pitch','gaze_yaw']].values - [90, 180]

        # get all pairwise distances, then take immediate off-diagonal to get 
        # distances between a point and subsequent points
        pairwise_distances = np.rad2deg(haversine_distances(np.deg2rad(coordinates)))
        distances = np.diagonal(pairwise_distances, offset=1)

        # calculate velocity by dividing distance by time
        velocity = distances / np.diff(df['exp_time'])

        # calcualte MAD veloctiy --> chop off first time that was removed in distance calculation
        times = df['exp_time'].values
        mad_velocity, trimmed_time = sliding_window_mad(velocity, times[1:], window_size=0.1)

        # get indices of each fixation window and the length of each fixation
        fixation_indices, length_fixations = get_fixation_indices(
            stat = mad_velocity if self.params.fix_type == 1 else velocity,
            threshold = self.params.min_mad
        )

        # Now calculate fixation centroids --> first level derivative before concatenating fixations
        lat, lon = coordinates.T

        df_fixations = calculate_fixation_centroids(
            lat=lat,
            lon=lon,
            time=trimmed_time,
            fixation_indices = fixation_indices
        )

        df_fixations = concatenate_fixations(
            df = df_fixations,
            spatial_distance = self.params.fix_spatial_dist,
            temporal_distance = self.params.fix_temp_dist
        )

        # Apply filtering to exclude first and last fixations, if not pretrial
        if any([item in trial.trial for item in self.params.pretrial_list]):
                print('FIXATION TRIM FILTER - not applied to non-scene trials!')
        else:
            total_exclusions = self.params.exclude_first_n_fix + self.params.exclude_last_n_fix

            if len(df_fixations) > self.params.exclude_first_n_fix and self.params.fix_type == 1:
                print (f'FIXATION TRIM FILTER - removed {total_exclusions} out of {len(df_fixations)} fixations')
                df_fixations = df_fixations.loc[self.params.exclude_first_n_fix:len(df_fixations) - self.params.exclude_last_n_fix].reset_index(drop=True)

        # Filter out short fixations
        duration_filter = df_fixations['duration'] < self.params.exclude_fix_durs_less_than

        print (f'FIXATION DURATION FILTER - removed {sum(duration_filter)} out of {len(duration_filter)} fixations')
        df_fixations = df_fixations[~duration_filter].reset_index(drop=True)

        ## REPORT SOME FIXATION STATS
        durations = df_fixations['duration']
        print (f'\nFIXATION STATISTICS')
        print (f'Num fixations: {len(durations)}')
        print (f'Avg fixation duration: {durations.mean():.3f} seconds')
        print (f'STD fixation duration: {durations.std():.3f} seconds')
        print (f'Min fixation duration: {durations.min():.3f} seconds')
        print (f'Max fixation duration: {durations.max():.3f} seconds')

        df_fixations['fix_yaw'] = np.mod(df_fixations['fix_yaw'] + 180, 360)
        df_fixations['fix_pitch'] = np.mod(df_fixations['fix_pitch'] + 90, 180)

        # add normalized start and end times (useful for heatmapping)
        df_fixations['norm_start_time'] = df_fixations['start_time'] - df['exp_time'].iloc[0]
        df_fixations['norm_end_time'] = df_fixations['end_time'] - df['exp_time'].iloc[0]

        return df_fixations

    def plotFixations(self, df_fixations, image_path, out_path=None, fig_size=(20,10)):
        
        if image_path is None:
            print(f'No image path provided!')
            return 

        image_name = os.path.basename(image_path)

        if not os.path.exists(image_path):
            print (f'No image found at provided path for {image_name}')
            
        img = cv2.imread(image_path)
        res = cv2.resize(img, dsize=(self.params.plotting_image_width, self.params.plotting_image_height), interpolation=cv2.INTER_CUBIC)[..., ::-1]
        
        x_img, y_img = degrees_to_pixels(df_fixations['fix_yaw'], df_fixations['fix_pitch'], self.params.plotting_image_width, self.params.plotting_image_height)

        # grab the spread of points
        point_spread = df_fixations['spread']
        x_spread, y_spread = degrees_to_pixels(point_spread, point_spread, self.params.plotting_image_width, self.params.plotting_image_height)

        fig = plt.figure(figsize=fig_size)
        plt.axis('off')
        plt.imshow(res)
        
        plt.errorbar(x_img, y_img, xerr=x_spread, yerr=y_spread, linewidth=2.75, linestyle="None", color='r', zorder=1)
        plt.scatter(x_img, y_img, s=2.75*sum(fig_size), c=df_fixations['duration'], cmap='hot')
        
        plt.title(image_name)
        
        if out_path:
            plt.savefig(out_path)
            plt.close(fig)

    def splitFixationTimesteps(self, df_fixations, overlapping=False):
        """

        Splits a dataframe of fixations into chunks of non-overlapping time 
        based on a flag of how many chunks
        """

        if not len(df_fixations):
            print ('No fixations! Cannot split into timesteps.')
            return None

        timestep_bounds = np.linspace(0, self.params.scene_length, self.params.heatmap_timesteps + 1);
        
        timesteps = []

        for i in range(len(timestep_bounds)-1):
            # find the indices of fixations for this timestep
            indices = np.where(
                np.logical_and(
                    df_fixations['norm_start_time'] >= timestep_bounds[i],
                    df_fixations['norm_start_time'] < timestep_bounds[i+1])
            )[0]
            
            # if there are any indices for this timestep
            if indices.any():
                timesteps.append(
                    df_fixations.iloc[indices.min():indices.max()]
                )
            else:
                timesteps.append(
                    np.array([])
                )
                
        return timesteps

        # # calculate timestep indices as points when fixation falls between subsequent steps of bounds
        # timestep_indices = [np.where(
        #     np.logical_and(
        #         df_fixations['norm_start_time'] >= timestep_bounds[i], 
        #         df_fixations['norm_start_time'] < timestep_bounds[i+1])
        # )[0].max() for i in range(len(timestep_bounds)-1)]

        # timesteps = np.split(df_fixations, np.array(timestep_indices) + 1)[:-1]
        # return timesteps

    def createDensityMap(self, df_fixations, base_width=200):
        
        # if no fixations during duration, create density map of zeros
        if np.any(df_fixations) == False:
           density_map = np.zeros((self.params.plotting_image_height,self.params.plotting_image_width))
           return density_map
 
        # map location of fixations to equirectangular space
        x_pix, y_pix = degrees_to_pixels(df_fixations['fix_yaw'], df_fixations['fix_pitch'], self.params.plotting_image_width, self.params.plotting_image_height)

        # normalize durations and trim outliers
        norm_durations = scale_durations(df_fixations['duration'], bound_filtering=self.params.bound_filtering)
        
        # # save normalized durations
        # duration_idx = df_fixations.columns.get_loc('duration')
        # df_fixations.insert(duration_idx,'norm_duration', norm_durations)  # TODO: not working
        
        # create a base image to use 
        density_map = np.zeros((self.params.plotting_image_height,self.params.plotting_image_width))
        
        # set points in the image as fixation durations
        for i, (col, row) in enumerate(zip(x_pix, y_pix)):
            density_map[row, col] += norm_durations.iloc[i]

        # pad the image --> to account for wraparound, we append the right edge onto to
        # the left part of the image, and the left edge to the right part of the image
        density_map = np.hstack((
            density_map[:,-self.params.plotting_image_width//4:], # right edge
            density_map, 
            density_map[:, :self.params.plotting_image_width//4], # left edge
        ))

        # apply variable width smoothing over pitch (axis=0), but regular gaussian over columns
        density_map = apply_gaussian_smoothing(density_map, axis=0, variable_width=True)
        density_map = apply_gaussian_smoothing(density_map, axis=1, variable_width=False)

        # unpad image
        density_map = density_map[:, self.params.plotting_image_width//4:-self.params.plotting_image_width//4]
            
        return density_map
    
    def plotFixationDensity(self, density_map, image_path, start_dur, end_dur, alpha=0.6, out_path=None, fig_size=(20,10), vmin = None, vmax = None):
        
        if image_path is None:
            print(f'No image path provided!')
            return 
        
        image_name = os.path.basename(image_path)

        if not os.path.exists(image_path):
            print (f'Image not found for {image_name}')

        img = cv2.imread(image_path)
        res = cv2.resize(img, dsize=(self.params.plotting_image_width, self.params.plotting_image_height), interpolation=cv2.INTER_CUBIC)[..., ::-1]
        
        #overlay the heatmap on top of the panoramic image

        # define colormap range
        if vmin is None:
            vmin = density_map.min()
        if vmax is None:
            vmax = density_map.max()  

        fig = plt.figure(figsize=fig_size)
        plt.axis('off')
        plt.imshow(res)
        plt.imshow(density_map, alpha=alpha, vmin=vmin, vmax=vmax)
        
        plt.suptitle(image_name)
        plt.title(f'Duration: {start_dur} to {end_dur} seconds')
        

        if out_path:
            plt.savefig(out_path)
            plt.close(fig)
        
    def makeDensityMapGIF(self,plot_dir=None,out_path=None):
        
        if plot_dir is None:
            print(f'No density map directory provided!')
            return
        
        if len(os.listdir(plot_dir)) == 0:
            print (f'Plot directory is empty. Skipping!')
            return
            
        frames = [Image.open(image) for image in glob.glob(f'{plot_dir}/*.png')]
        frames[0].save(out_path, format='GIF', save_all=True, append_images=frames[1:], duration=self.params.gif_frame_duration, loop=0) 

    def trimTrialLength(self, trial):
        """
        
        Excludes the first and last N seconds of trial based on parameters
        """
        raise NotImplementedError()
    
    def calculateTrialFPS(self, trial):
        """
        Calculates each trial's frames per second based on Unity Engine time
        """
        
        trial_fps = np.mean(1.0/trial['exp_time'].diff())

        return trial_fps