import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import haversine_distances

from scipy.signal.windows import gaussian
from scipy.signal import convolve

def get_window_indices(x, window_size):
        
    # Given a x (a series of points/times) and y (a start index)
    def get_forward_index(x, index):
        window_forward = abs(x - (x[index] + window_size/2))
        return np.argwhere(window_forward == min(window_forward))[0]
    
    def get_backward_index(x, index):
        window_backward = abs(x - (x[index] - window_size/2))
        return np.argwhere(window_backward == min(window_backward))[0]
    
    # Calculate start and end of the mad sliding window
    start_index = get_forward_index(x, 0)
    end_index = get_backward_index(x, len(x)-1)
    
    indices = [np.arange(get_backward_index(x, i), get_forward_index(x,i)+1) for i in np.arange(start_index, end_index+1)]
    
    return indices, start_index, end_index

def mad(data, axis=None):
    return np.mean(abs(data - np.mean(data, axis)), axis)

def sliding_window_mad(x, time, window_size):
    
    """
    Rolling sliding window of Mean Absolute Deviation
    
    window_size is expected to be time
    """
    
    indices, start_index, end_index = get_window_indices(time, window_size)
    mean_absolute_deviation = np.asarray([mad(x[i]) for i in indices])

    trimmed_time = time[np.arange(start_index,end_index)]
    
    return mean_absolute_deviation, trimmed_time

######### FIXATION FUNCTIONS #########

def get_fixation_indices(stat, threshold):
    
    # Create a list of potential fixation indices based on moments where it fell below
    # MAD threshold
    potential_fixation_indices = np.where(stat < threshold)[0]

    # Find when indices were consecutive --> 1 when consecutive, >1 when not consecutive
    potential_consecutive_indices = np.diff(potential_fixation_indices)

    # Find  all non-consecutive potential fixes, will be potential_fix_diff will be > 1 
    # when changing b/w two discrete gaze points
    
    potential_fixation_start = np.concatenate([
        np.zeros(1).astype(int),
        np.where(potential_consecutive_indices > 1)[0] + 1 # adjust for missing the first index
    ])

    # Fixation end is right before the next fixation start
    potential_fixation_end = np.concatenate([
        np.where(potential_consecutive_indices > 1)[0],
        potential_consecutive_indices.shape
    ])

    # Get the indices for the original velocity/MAD array
    begin_fixation_indices = potential_fixation_indices[potential_fixation_start] - 1
    end_fixation_indices = potential_fixation_indices[potential_fixation_end]

    # gets sets of indices of potential fixations
    fixation_indices = list(map(lambda x: np.arange(*x), zip(begin_fixation_indices, end_fixation_indices)))

    # Get length of each fixation
    length_fixations = list(map(len, fixation_indices))
    
    return fixation_indices, length_fixations

def calculate_fixation_centroids(lat, lon, time, fixation_indices):
    """
    
    """
    
    df = pd.DataFrame(columns=['fix_yaw', 'fix_pitch', 'start_time', 'end_time', 'spread'])
    
    for i, indices in enumerate(fixation_indices):

        lat_sample, lon_sample = lat[indices], lon[indices]
        latbar, lonbar = sphere_centroid(lat_sample, lon_sample)

        point_spread = np.rad2deg(haversine_distances(
            np.deg2rad(np.stack([latbar, lonbar])[np.newaxis]),
            np.deg2rad(np.stack([lat_sample, lon_sample]).T)
        ))

        df.loc[len(df)] = {
            'fix_yaw': lonbar,
            'fix_pitch': latbar,
            'start_time': time[indices[0]],
            'end_time': time[indices[-1]],
            'spread': point_spread.mean()
        }
        
    df['duration'] = df['end_time'] - df['start_time']
    
    return df

def concatenate_fixations(df, spatial_distance, temporal_distance):
    """
    
    Concatenates fixations that fall within a range
    
    """

    n_fixations = len(df)
    
    fixation_spatial_distance = np.rad2deg(haversine_distances(np.deg2rad(df[['fix_pitch', 'fix_yaw']].values)))
    fixation_spatial_distance = np.diagonal(fixation_spatial_distance, offset=1)

    # Find difference in time between subsequent fixations
    start_time, end_time = df[['start_time', 'end_time']].values.T
    fixation_time_distance = start_time[1:] - end_time[:-1]

    concatenate_indices = np.where(np.logical_and(
        fixation_spatial_distance < spatial_distance, 
        fixation_time_distance < temporal_distance
    ))[0]
    
    n_concatenations = len(concatenate_indices)

    for index in concatenate_indices:
        df.loc[index] = df.loc[index:index+1].agg(np.mean)
        df = df.drop(index=[index+1]).reset_index(drop=True)
        
        concatenate_indices -= 1

    print (f'FIXATION CONCATENATION: concatenated {n_concatenations} out of {n_fixations} fixations')
        
    return df

def scale_durations(durations, bound_filtering=False):
    norm_durations = durations.copy()
    if bound_filtering:
        norm_durations[durations > np.percentile(durations,95)] = np.percentile(durations,95)
        norm_durations[durations < np.percentile(durations,0.1)] = np.percentile(durations,0.1)

    # make sure all values aren't the same --> otherwise normalization produces nan
    if len(set(norm_durations)) > 1:
        norm_durations = (norm_durations - min(norm_durations)) / (max(norm_durations) - min(norm_durations))
    norm_durations = 0.1 + norm_durations * 0.9
    
    return norm_durations

### PLOTTING FUNCTIONS ####

def degrees_to_pixels(x, y, width, height):
    """
    Scales x,y coordinates in degrees to the size of an image
    """
    
    x_img = np.mod(np.round((x*width) / 360).astype(int), width)
    y_img = np.mod(np.round((y*height) / 180).astype(int), height)
    
    return x_img, y_img

def pixels_to_degrees(x, y, width, height):
    """
    Scales x,y coordinates in degrees to the size of an image
    """
    
    yaw = np.mod(x*360/width, 360)
    pitch = np.mod(y*180/height, 180)
    
    return yaw, pitch


### GAUSSIAN SMOOTHING ###

def get_gaussian_window(n_points, width_factor=2.5):
    """
    width_factor matches how matlab calculates the sigma of the gaussian
    """
    
    # calculate sigma based on a width factor (this matches matlab's implementation)
    sigma = (n_points - 1) / (2*width_factor)
    window = gaussian(n_points, sigma)
    
    return window

def apply_gaussian_smoothing(image, axis, gaussian_base_width=200, variable_width=False):
    
    n_items = image.shape[axis]
    
    for i in range(n_items):
        indices = np.arange(i, i+1)[..., np.newaxis]
        
        if variable_width:
            _, pitch = pixels_to_degrees(i, i, 2000, 1000)
            pitch = pitch - 90
            
            # scale the number of points in the distribution by pitch
            n_points = np.round(gaussian_base_width*(1/np.cos(np.deg2rad(pitch)))) - 1
            n_points = n_points if n_points < 1e6 else 1e6
        else:
            n_points = gaussian_base_width
                    
        # grab gaussian window for the current pitch
        gauss_window = get_gaussian_window(n_points)
        
        # find the values along the current axis --> convolve with the window
        vals = np.take_along_axis(image, indices=indices, axis=axis).squeeze()
        smoothed = convolve(vals, gauss_window, mode='same')
        
        # insert values back into the axis
        np.put_along_axis(image, indices=indices, values=np.expand_dims(smoothed, axis=axis), axis=axis)
    
    return image

### FOR CARTESIAN/SPHERE MAPPING

from typing import Tuple, Union
from math import sin, cos, atan2, sqrt


Number = Union[int, float]
Vector = Tuple[Number, Number, Number]


def distance(a: Vector, b: Vector) -> Number:
    """Returns the distance between two cartesian points."""
    x = (b[0] - a[0]) ** 2
    y = (b[1] - a[1]) ** 2
    z = (b[2] - a[2]) ** 2
    return (x + y + z) ** 0.5

  
def magnitude(x: Number, y: Number, z: Number) -> Number:
    """Returns the magnitude of the vector."""
    return np.sqrt(x * x + y * y + z * z)


def to_spherical(x: Number, y: Number, z: Number) -> Vector:
    """Converts a cartesian coordinate (x, y, z) into a spherical one (radius, theta, phi)."""
    radius = magnitude(x, y, z)
    phi = np.arctan2(z, np.sqrt(x * x + y * y))
    theta = np.arctan2(y, x)
    return (radius, theta, phi)


def to_cartesian(radius: Number, theta: Number, phi: Number) -> Vector:
    """Converts a spherical coordinate (radius, theta, phi) into a cartesian one (x, y, z)."""
    x = radius * np.cos(phi) * np.cos(theta)
    y = radius * np.cos(phi) * np.sin(theta)
    z = radius * np.sin(phi)
    return (x, y, z)

def sphere_centroid(lat, lon):
    """
    Designed to match MATLAB meanm function w/ assumed radius of 1

    theta=azimuth (longitude)
    phi=elevation (latitude)
    """
    
    # convert to radians
    lat, lon = map(np.deg2rad, [lat, lon])

    # convert to x,y,z coordinates in cartesian space
    points_cartesian = to_cartesian(radius=1, theta=lon, phi=lat)

    # vector sum of cartesian points
    x, y, z = map(np.sum, points_cartesian)

    # map back to spherical space
    r, lon, lat = to_spherical(x,y,z)

    # # convert back to degrees
    latbar, lonbar = map(np.rad2deg, [lat, lon])
    
    return latbar, lonbar