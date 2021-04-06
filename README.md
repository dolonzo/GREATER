# GREATER Processing Pipeline for TMS-EEG Recordings
# 
# The functions contained within this repository allow for implementation within EEGLAB of the GREATER processing pipeline as described in [insert citation to paper]. The pipeline consists of the following steps: 
       1. Remove and interpolate the pulse artifact
       2. Downsample the data
       3. Remove interpolated data prior to spatial filtering
       4. Separate Intertrain Interval data from TMS pulse epochs
       5. Reject bad channels based on FASTER criteria applied to the
       Intertrain Interval data
       6. Reject bad portions of data from Intertrain Interval data based
       on FASTER criteria
       7. Compute the Generalized Singular Value Decomposition (GSVD)
       jointly on the Intetrain Interval data and the TMS Pulse epochs
       8. (Optionally) Plot the components
       9. (Optionally) Reject artifactual GSVD components from the
       dataset. This can be done based on amplitude in a post-pulse
       window or based on the singular values
       10. Recombine the TMS pulse epochs with the Intertrain Interval
       data for further processing and analysis
# Further processing used in the paper referenced after these steps included filtering (high- and low-pass), creating TMS pulse epochs, rejecting TMS pulse epochs contaminated by artifacts, ICA decomposition, and automated rejection of artifactual ICA components. These steps were not the focus of the GREATER pipeline, and thus not included within this repository. The script implementing these additional processing steps can be made available upon request.
      
# [Need to edit to actually include license language]. The GREATER pipeline includes dependencies from FASTER (insert citation) and TESA (insert citation). These are included within this repository as a convenience under the GNU General Public License which they are distributed under. Users of GREATER will likely find other functions from both plugins helpful, and installation of the full distributions are encouraged.

# Implementation details for each function can be found in the header to the function or by calling help function_name within MATLAB
