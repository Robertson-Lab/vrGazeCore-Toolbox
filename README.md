# vrGazeCore

This is the core code for gaze analysis including fixation calculation and heatmapping. JSM 2019

## Overview

This codebase is able to do four main functions:  

1. Processes raw gaze data into fixations for each participant & each trial. (`findFixations`)
2. Averages fixations into heatmaps for each participant & across participants for each scene. (`fix2Heat`)
3. Averages fixations over a series of timepoints for a scene, to create a time course. You can do this for each participant and at the group level. (`fix2HeatTimecourse`)
4. Makes GIFs of the time courses. (`plotTimecourseGif`)

## Requirements

-1) MATLAB 20--
2) [MATLAB Mapping Toolbox](https://www.mathworks.com/products/mapping.html)
3) [MATLAB Signal Processing Toolbox]
4) [MATLAB Image Processing Toolbox]
5) Project Directory which includes:
   1) Raw eyetracking data
   2) Stimulus files
   3) A copy of vrGazeCore toolbox

## Plug-and-Play

### To Use

1. Download vrGazeCore to your project folder
2. Change your directory paths in `setPaths.m`
   1. set `projectDir` to your project directory
   2. set `gazeCoreDir` to where you downloaded vrGazeCore toolbox
   3. set `projectRawDataDir` to where your raw eye-tracking data was saved from the experiment
   4. set `projectStimDir` to where your stimuli (e.g., scene image) files are located
3. Set the necessary parameters in `setParams.m`
   1. Processing options:
      - Do you want to *find fixations* in your data? Necessary for all other processing steps!
        >set `runFindFix = 1`
      - Do you want to make *heatmaps of the fixations*?
        >**for each subject**: set `runHeatmappingIndivid = 1`  
        >**for your cohort of subjects**: set `runHeatmappingGroup = 1`
      - Do you want to make *heatmaps of the timecourse* of fixations? Necessary for making GIFs!
        >**for each subject**: set `runTimecourseIndivid = 1`  
        >**for your cohort of subjects**: set `runTimecourseGroup = 1`
      - Do you want to make a *GIF of the heatmap timecourse*?
        >**for each subject**: set `runTimecourseGifIndivid* = 1`  
        >**for your cohort of subjects**: `runTimecourseGifGroup* = 1`
      - Do you want to *delete the timecourse images* after making the GIFs?
        >set `deleteTimecourseJPEG = 1`
   2. VR equipment parameters
      - Which version of Unity did you use?
        >`unityProjectVersion = 1`  
        >`unityProjectVersion = 0`  
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
       - Set `sceneLength` to the scene viewing duration in the experiment
4. Run `batchProcessPars.m` in MATLAB

### The Results

**Results Directory Structure**  
Results are located in the 'analResults' subfolder of the project directory. Each processing step has its own directory within 'analResults'. Results are further separated into 'mat' folders, which contain the MATLAB results and 'plot' folders, which contain any image files generated. Each subject and/or cohort has their own folder within the 'mat' and 'plot' folders.

1. **'analResults/fixations'** (results of  `findFixations`)
   1. '/mat' - .mat for each scene and sanity target  
      a. **rawData**:  structure containing raw head and eye movement data as well as the raw confidence  
      b. **fixData**:  contains the data for the fixations  
      c. **rawFixData**: ??
   2. '/plot' - .jpg for each scene  
      a. **'b_fixations_Durations.jpg'** - plot of raw fixation points coded by color onto each image
      b. **'_combinedPlot.jpg'** - plot of combined head, raw, fix data on one image; head in white, gaze black, fix cross
2. **'analResults/heatMaps'** (results of `fix2Heat`)
   1. '/mat' - .mat for each scene
      a. **histogram**: histogram of heatmap (`imshow(histogram)` will return a black & white version of the heatmap overlay)
   2. '/plot' - **heatmap of average fixations** for each scene
3. **'analResults/timecourseHeat'**
   1. '/mat' - .mat for each scene
      a. **histogram**: histogram of heatmap (`imshow(histogram)` will return a black & white version of the heatmap overlay) --BROKEN
   2. '/plot' - **heatmaps of average fixations for each time point**
      a. if `deleteTimecourseJPEG = 1`: static time course heatmaps (JPG files) deleted after GIF for scene generated
      b. if `deleteTimecourseJPEG = 0`: static time course heatmaps remain after making GIFs

## The Inner Workings

### Main Caculations

### Parameters

### Processing Options

Set these parameters to **1 to run the function**

    *runFindFix* find fixations in data
    *runHeatmappingIndivid*: create heatmap of fixations for each subject
    *runHeatmappingGroup*: create heatmap of fixations for entire cohort
    *runTimecourseIndivid*: create a timecourse of fixations for each subject
    *runTimecourseGroup*: create a timecourse of fixations for entire cohort
    *runTimecourseGifIndivid*: create a .gif file of the fixation timecourse for each subject

### Unity & Headset Options

    unityProjectVersion
    headsetType

### Subject Options

    cohortName
    listSubjectNames
    subjectNames

### Scene Parameters

    sceneFiles
    sceneList
    excludeScenes
    minSamples
    minTimeInScene
    scannedFilter
    scannedThresh
    sceneLength

### Eyetracking Options

    gazeType
    excludeFirstNSec
    minMad
    maxMad
    excludeFixDursLessThan
    useEye
    headLocked
    minConfThresh
    maxConfPercent
    eccFiltX
    eccFiltY
    useSmoothing
    useInterpolation
    durationForInterpolation
    fixType
    fixValidation
    valWindow
    excludeLastNSec
    excludeFirstNFix
    excludeLastNFix
    driftCorrection
    concatSanity
    avgPreTrial
    avgPreTrialThresh
    excludeByPreTrial
    avgPreTrailX
    avgPreTrialY

### Plotting & Heatmap Options

    plotHeadFixationInfo
    plotCombinedGazeFlag
    plotHeadRawFlag
    plotGazeRawFlag
    plotFixFlag
    trimFactor
    heatPlot
    padAmount
    boundFiltering
    imDimX
    imDimY
    timeStep
    timeWindow

### Fixed Parameters

    fovX
    fovY
    maxFOV
    norm2DVA.csv

### Description of the scripts

## Heirarchy of Scripts

    coreParams
    batchProcessPars
        loadCoreParams
            coreParams
        findFixations
            processRawData
            parseRaw2struct
            confidenceFilter
            eccFilter
            parseDS
            rectifyGaze
            wrapPointsEquirect
            bilatFilt
            calculateFixations
                wrapPointsEquirect
                degreesToPixels
            viewportNorm2dva
                loadn2dcoef
                    norm2DVA.csv
            plotFix
                degreesToPixels
                pointSpread
                    degreesToPixels
            plotCombinedGaze
                wrapPointsEquirect
                degreesToPixels
                pointSpread
            plotHeadRaw
                degreesToPixels
            plotGazeRaw
                degreesToPixels
            plotFix
        fix2Heat
            gaussianFilterEquirect
                pixelsToDegrees
        fix2HeatTimecourse
            antonioGaussian
            zScoreScale

    Helper Functions
        batchPlotMeta
        getSphereCoords
        SpiralSampleSphere
        trimMap

    Orphaned Code
        calculateSaccades [commented out; does work according to Tommy]
        runPlotTimecourseGif [standalone]
            plotTimecourseGif 
        plotRingDVA [seems to be replaced by loadn2dcoef]
