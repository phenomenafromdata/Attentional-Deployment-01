These files are Matlab scripts that analyze the data from the task and produce figures showing these analyses.
Each file with name "Figure_XX_AttDeploym01" is intended to produce a single figure and the corresponding analysis. Most of these files call a script named "map_asc2mat_AttDeployment.m", that is also included.

The code analyzes three kinds of data:

* Behavioral data (response times, intensity and valence ratings)

* Eye-tracking data (gaze position during the image presentation periods)

* ANT data (application of Attention Network Test to a subset of the participants).

These data can be found at: [https://figshare.com/account/home#/projects/232625](https://figshare.com/account/home#/projects/232625)

Basic task structure:

- Task consists of 20 trials
- Before each trial, participant receives one out of two possible instructions: "observe freely" or "focus your gaze inside blue circle"
- Then 5 images are presented in sequence, with each one lasting 4 seconds
- Then participant has to rate the images in emotional intensity (1 to 9)
-  Then participant has to rate the images in emotional valence (1 to 9)
-  Participant views 5 screens of different grey levels (4 seconds each) containing only a fixation oval at the center
-  Next trial begins, process repeats

Images:
All images belong to IAPS data bank.
There are 5 possible image types (on each trial all the images presented belong to one of these types)

1. Image Neutral, no circle (no attentional focus)
2. Image Neutral, circle (attentional focus, non-arousing)
3. Image Unpleasant, no circle (no attentional focus)
4. Image Unpleasant, circle (attentional focus, located in an emotionally non-arousing area of the image)
5. Image Unpleasant, circle (attentional focus, located in an emotionally arousing area of the image)

Task design was based on:
Jamie Ferri, Joseph Schmidt, Greg Hajcak, Turhan Canli (2013). Neural correlates of attentional deployment within unpleasant pictures. Neuroimage 15:70:268-77.
doi: 10.1016/j.neuroimage.2012.12.030.

