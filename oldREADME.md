# vrGazeCore

This is the core code for gaze analysis including fixation calculation and heatmapping. JSM 2019

## Overview

This codebase basically does three things:

1) processes raw gaze data into fixations for each participant/trial (findFixations.m)
2) averages fixations into heatmaps for each participant/scene (fix2Heat.m)
3) creates average heatmaps across participants for each scene, or across scenes for each participant (fix2Heat.m)

## Getting started

-As a user, you mostly interface with the top level scripts:

    --coreParams (modify your directories and other parameters here)
    --loadCoreParams (loads these parameters for the analysis)
    --batchProcessPars (runs the three analyses described above)

-To get started (short version):

1) Make changes to coreParams to match your study needs (e.g. projectDir and sceneDir)
2) Run loadCoreParams
3) Run batchProcessPars

Parameters that must be changed!

    cohortName
    listSubjectNames
        - 0 = point to subject raw data directory
        - 1 = manually list subjects in section

## Short description of all scripts

Top level scripts:

    --coreParams (user changes this script to modify your directories and other parameters)
    --loadCoreParams (called by ....
    --batchProcessPars (...

Mid-level scripts:

    --findFixations
    --fix2Heat

Sub-routines:

## Legacy information

fixData had a number of fields which have been renamed

## Parameters

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

## Paths

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

    Orphaned Code
        calculateSaccades [commented out; does work according to Tommy]
        runPlotTimecourseGif [standalone]
            plotTimecourseGif 
        plotRingDVA [seems to be replaced by loadn2dcoef]
        batchPlotMeta
        trimMap
        getSphereCoords
        SpiralSampleSphere