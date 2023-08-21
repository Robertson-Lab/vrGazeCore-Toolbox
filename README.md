# vrGazeCore

Still in development, use at your own discretion!

When using the toolbox, please cite our [CCN Abstract](https://2023.ccneuro.org/view_paper.php?PaperNum=1555).

vrGazeCore was authored by Deepasri Prasad, Amanda J. Haskins, Thomas L. Botch, Jeff Mentch, and Caroline E. Robertson.

vrGazeCore is a toolbox that processes raw eye-tracking data from head mounted virtual reality displays and does basic gaze analysis, including determining fixations and creating fixation density maps. It is available in both MATLAB and Python.

This README will be a simple overview of the dependencies needed to run vrGazeCore, the results it outputs, the parameters that can be adjusted by the end user, and the scripts contained.

## Overview

This codebase is has two main functions:  

1. Processes raw gaze data into fixations for each participant & each trial. (`findFixations`)
2. Creating duration-weighted and time-segmented fixation density maps at the participant and group level (`fix2Heat`)

## Requirements & Dependencies

1) Unity Version 2020.3.32 or later
2) MATLAB 2020 or later
3) [MATLAB Mapping Toolbox](https://www.mathworks.com/products/mapping.html)
4) [MATLAB Signal Processing Toolbox](https://www.mathworks.com/products/signal.html)
5) [MATLAB Image Processing Toolbox](https://www.mathworks.com/products/image.html)
6) Project Directory which includes:
   1) Raw eyetracking data folder
   2) Stimulus files folder
   3) A copy of vrGazeCore toolbox

### To Use

1. Download vrGazeCore to your project folder
   1. If you have multiple projects that use vrGazeCore, it's recommended to have separate `setPaths.m` and `setParams.m` for each project
2. Change your directory paths in `setPaths.m`
   1. set `projectDir` to your project directory
   2. set `gazeCoreDir` to where you downloaded vrGazeCore toolbox
   3. set `projectRawDataDir` to where your raw eye-tracking data was saved from the experiment
   4. set `projectStimDir` to where your stimuli (e.g., scene image) files are located
3. Set the necessary parameters in `setParams.m`
   1. Processing options:
      - Do you want to *find fixations* in your data? Necessary for all other processing steps!
        >set `runFindFix = 1`
      - Do you want to make *fixation density maps*?
        >**for each subject**: set `runHeatmappingIndivid = 1`  
        >**for your cohort of subjects**: set `runHeatmappingGroup = 1`
      - Do you want to make *time-segmented* fixation density maps? Necessary for making GIFs!
        > set `heatmapTimesteps` as an array of time segmentations (e.g., `heatmapTimesteps = [1 4];` will generate one fixation density map spanning the whole scene duration and four fixation density maps, each spanning consecutive quarters of the scene duration)
      - Do you want to make a *GIF of the heatmap timecourse*?
        >**for each subject**: set `runTimecourseGifIndivid* = 1`  
        >**for your cohort of subjects**: `runTimecourseGifGroup* = 1`
   2. VR equipment parameters
      - Which kind of headset did you use?
        >**DK2** `headsetType = 0`  
        >**Vive** `headsetType = 1`  
        >**ViveEye** `headsetType = 2`  
        >**Oculus Go** `headsettype = 3`  
   3. Selecting subjects
      - Set `cohortName` to the name of your group of subjects! Used when doing group level analysis!
      - Do you want to manually list subjects or run a directory of subjects?
        > if **manually**: list subjects under `subjectNames`, with each subject separated by a `;`  
        > if **directory**: will use the subjects in the raw data directory `projectRawDataDir`
   4. Scene options
       - Set `sceneLength` to the experiment scene viewing duration in seconds
4. Run `batchProcessPars.m` in MATLAB

### The Results

**Results Directory Structure**  
Results are located in the 'eyeTrackResults' subfolder of the project directory. Each processing step has its own directory within 'eyeTrackResults'. Results are further separated into 'mat' folders, which contain the MATLAB results and 'plot' folders, which contain any image files generated. Each subject and/or cohort has their own folder within the 'mat' and 'plot' folders.

1. **'eyeTrackResults/fixations'** (results of  `findFixations`)
   1. '/mat' - .mat for each scene and sanity target  
      a. **rawData**:  structure containing raw head and eye movement data as well as the raw confidence  
      b. **fixData**:  contains the data for the fixations  
      c. **rawFixData**: ??
   2. '/plot' - .jpg for each scene  
      a. **'b_fixations_Durations.jpg'** - plot of raw fixation points coded by color onto each image
      b. **'_combinedPlot.jpg'** - plot of combined head, raw, fix data on one image; head in white, gaze black, fix cross
2. **'eyeTrackResults/heatMaps'** (results of `fix2Heat`)
   1. '/mat' - .mat for each scene
      a. **histogram**: histogram of heatmap (`imshow(histogram)` will return a black & white version of the heatmap overlay)
   2. '/plot' - **heatmap of duration-weighted fixations** for each scene

## Parameters

### Need to Adjust Parameters

#### Processing Options

Set these parameters to **1 to run the function**

`runFindFix`: find fixations in data

`runHeatmappingIndivid`: create fixation density maps for each subject

`runHeatmappingGroup`: create fixation density maps for the entire cohort

`runTimecourseGifIndivid`: create a .gif file of the fixation density timecourse for each subject

`runTimecourseGifGroup`: create a .gif file of the fixation density timecourse for the entire cohort

#### Headset Options

`headsetType`: type of headset used

    0: DK2
    1: Vive
    2: ViveEye    
    3: Oculus Go

#### Subject Options

`cohortName`: name of cohort, used to label group processing files

`listSubjectNames`: if `0`, use raw data directory; if `1`, list subjects manually under `subjectNames`

#### Heatmapping Options

`heatmapTimesteps`: array of how many time segments to divide scene into

#### Scene Parameter

`sceneLength`: scene duration in seconds

### Additional Scene Parameters

`minSamples`: skip scenes with less than set number of samples; default is 100

`minTimeinScene`: skip scenes shorter than set duration; default is 5 seconds

`scannedFilter`: if `1`, removes scenes that were not explored; default is 1

`scannedThresh`: set to a percent threshold to exclude scenes not explored by that much based on head direction; default is 66

`excludeScenes`: list of scenes to exclude by subject

`ignoreList`: anything in scene directory that should be ignored for all subjects (e.g., subfolders)

`sceneFiles`: stimulus files pulled from stimulus directory whose names are added to `sceneList`

### Eyetracking Options

`gazeType`: type of tracking

    0: 2D tracking (default)
    1: 3D tracking

`excludeFirstNSec`: exclude set number of seconds at the start of the trial; default is 0 seconds

`minMad`: windows with a mean absolute deviation less than set degrees per second are classified as potential fixations; default is 50 degrees per second

`maxMad`: windows with a mean absolute deviation greater than set degrees per second are classified as likely saccades; default is 100 degrees per second

`excludeFixDursLessThan`: exclude any fixation with less than set duration; default is 0.1 seconds

`useEye`: if stereo eye-tracking, which eye to use

    0: eye 0
    1: eye 1
    2: choose best eye at each point
    3: average eyes and set confidence to lowest confidence of 2 eyes (default)

`headLocked`: set to `1` if used a  head-locked paradigm, sets all head direction to straight ahead

`minConfThresh`: data lower than set pupil confidence value is discarded when applying confidence filtering; default is 0.25

`maxConfPercent`: if more than set percentage of data is discarded due to confidence, the scene is skipped; default is 75

### Plotting & Heatmap Options

`plotVisiblity`: if `off`, figures are not visible when making them; default is `off`

`plotHeadFixationInfo`: if `1`, plot head yaw path (with fixation points), velocity, acceleration, and mean absolute deviation on image

`plotCombinedGazeFlag`: if `1`, plot combined head, gaze, and fixation points on image

`plotHeadRawFlag`: if `1`, plot the raw HMD direction center point on image

`plotGazeRawFlag`: if `1`, plot the raw gaze points on image

`plotFixFlag`: if `1`, plot the fixations on image

`plotMADFlag`: if `1`, plot viewport position & mean absolute deviation for each trial over time

`trimFactor`: if `1`, trim top and bottom by trim factor; defaul is 0

`padAmount`: set amount of padding in pixels around the border of image to account for smoothing at edges; default is 100 pixels

`boundFiltering`: if `1`, apply upper and lower bound filtering to density map; default is 0

`delayTime`: delay time in seconds between images in GIF

`deleteTimecourseJPEG`: if `1`, delete timecourse .jpeg images after making gif

`imDimX`: horizontal dimension (in pixels) of equirectangular images

`imDimY`: vertical dimension (in pixels) of equirectangular images

### Fixed Parameters (Change at Your Own Risk!)

`fovX`: horizontal field of view of the headset in degrees (depends on headset type)

`fovY`: vertical field of view of the headset in degrees (depends on headset type)

`maxFOV`: maximal field of view in degrees

`norm2DVA.csv`: values needed to convert normalized screen coordinates to degrees of visual angle
