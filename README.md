This repository contains all code necessary to reproduce the results from “Loftus et al. 2022. Ecological and social pressures interfere with homeostatic sleep regulation in the wild. eLife. https://doi.org/10.7554/eLife.73695”. To reproduce the analysis, run the scripts in the order they are numbered. Please email me with any questions at jcloftus@ucdavis.edu.


Main analysis:

The script “00_processing_raw_acc.ipynb” takes the raw GPS and accelerometry data from 2012, downloaded with “all sensor types” (i.e. including both GPS and accelerometry data) from Movebank as an input – “Collective movement in wild baboons.csv”. The script then removes columns that are not needed, and downsamples and interpolates the daytime accelerometry data to match the nighttime accelerometry bursts. It outputs the file “all_burst_acc.csv”.

The script “01_acc_to_vedba.R” takes “all_burst_acc.csv” as an input, calculates the average VeDBA and log VeDBA for each accelerometry burst, and outputs the file “full_night_and_day_data.csv”

The script “02_acc_1min_vedba_2012.m” takes the raw data as an input (“Collective movement in wild baboons.csv”, downloaded from Movebank), calculates the average VeDBA for each minute of the day, using “02a_calc_vedba.m” and “02b_calc_stat_only_vedba.m”, and produces “vedba_mean_2012.csv” as an output.

The script “03_baboon_sleep_analysis.R” performs most of the analysis associated with the manuscript. It requires the following files inputs (with the locations of the inputs in parentheses): “full_night_and_day_data.csv” (produced by script above), “sleep trees.cpg/.dbf/.prj/.shp/shx” (Dyrad), “env_data.csv-6150899038464587825.csv” (Dryad), “vedba_mean_2012.csv” (both produced by script above, and also available on Dryad). The script runs the sleep classification algorithm using the log VeDBA data, adds all data that is used for predictor variables in the models (e.g. temperature, moon phase), runs most of the models used in the analysis, and plots several the results from these models.

The script “04_social_sleep_analysis.R” takes “full_night_and_day_data.csv” (produced by “01_acc_to_vedba.R”) and “final_sleep.csv" (produced by “03_baboon_sleep_analysis.R”) as inputs. The script performs permutations to test the sentinel hypothesis (whether at least one group member is awake more often than expected by chance) and whether the group exhibits synchronization in their sleep-wake patterns during the night. The script then runs a model to test whether baboons are more likely to synchronize their sleep-wake patterns when they are sleeping in the same tree.

The script “05_arousal_threshold_analysis.R” takes “sleep_analysis_pub_code_pre_mods.RData” (produced by “03_baboon_sleep_analysis.R”) as an input. The script further prepares the data for analysis and runs a model to test whether baboons are less likely to wake in response to the waking activity of their neighbors (i.e. if they have a higher arousal threshold) following nights of poor sleep.

The raw data is available for download from the project called  “Collective movement in wild baboons” on Movebank.


Validation study:

The script “00_processing_raw_acc_2019.ipynb” takes the raw GPS and accelerometry data from 2019, downloaded with “all sensor types” (i.e. including both GPS and accelerometry data) from Movebank as an input – “Papio Anubis Mpala 2019.csv”. The script then removes columns that are not needed as well as the data from individuals whose collars were programmed to sample on a different schedule than in 2012. With a different sampling schedule, we could not apply the same sleep algorithm to these individuals. The script then produces the output “2019_Papio_anubis_acc_Loftus_et_al_Dryad.csv", which is published on Dryad (note: the full 2019 dataset is not yet publicly available on Movebank), and this file becomes the input for the rest of the script. The rest of the script interpolates both the daytime and nighttime accelerometry data to match the sampling rate of the 2012 accelerometry bursts. It outputs the file “validation_burst_acc.csv”.

The script “01_acc_to_vedba_2019.R” takes “validation_burst_acc.csv” as an input, trims the daytime accelerometry to match the nighttime accelerometry bursts, calculates the average VeDBA and log VeDBA for each burst, and outputs this information in the dataframe “2019_full_night_and_day_data.csv”.
The script “02_sleep_validation.R” performs the actual validation of the sleep algorithm. It takes “2019_full_night_and_day_data.csv” and files within the folder “loopy_focal_follows_2021_09_17”, which contain the behavioral observations from the thermal imagery, as inputs. It also requires the input file “tag_metadata.csv”, and .txt files associated with the thermal videos that are saved within the archive of the MPI-AB EAS department storage (these raw inputs are not needed if downloading the data from Dryad). The script first runs the sleep classification algorithm on the log VeDBA data, and saves the sleep classification. The script then aggregates the behavioral observation data into one dataframe, adds the real (absolute) timestamps based on the timestamps of the frames of the videos, trims the dataframe to only the relevant observations, and produces the file “2019_Papio_anubis_behavioral_scoring_Loftus_et_al_Dryad.csv”, which is available for download on Dryad (the raw behavioral observations are not available for download). The rest of the script uses the sleep classification from above within this script, and the behavioral scoring csv that is available on Dryad as inputs, and applies a time correction to he behavioral scoring, so that it matches GPS time, then compares the sleep classification from the algorithm to the sleep classification from behavioral observations to determine the accuracy of, and produce a confusion matrix for, the accelerometry-based sleep classification algorithm.

