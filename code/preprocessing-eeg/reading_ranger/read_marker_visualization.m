clear
clc

%% Quick Epoch Visualization Script
% This script loads one participant's data and shows epochs around S255/S127 markers
% Created for visualizing the audio wave used in the READ study from
% StimTrak
% NDCLab, Charles Knowlton

%% Setup paths (from your script)
main_dir = '/home/data/NDClab/datasets/read-study2-dataset';
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b'));

% Remove octave conflict
rmpath(['/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b' filesep 'functions' filesep 'octavefunc' filesep 'signal'])

%% Start EEGLAB
eeglab;

%% Load your participant's raw data
% 1 random participant to use for visualizing
participant_file = '/home/data/NDClab/datasets/read-study2-dataset/sourcedata/raw/s1_r1/eeg/sub-3300008/sub-3300008_all_eeg_s1_r1_e1.vhdr';

% Split into path and filename separately
[filepath, filename, ext] = fileparts(participant_file);
filename_with_ext = [filename ext];

% Load the data
EEG = pop_loadbv(filepath, filename_with_ext);

%% Load channel locations
channel_locations = loadbvef('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeg_preprocessing/chan_locs_files/electrode_locs_files/CACS-128-X7-FIXED-64only.bvef');
EEG.chanlocs = channel_locations;
EEG = eeg_checkset(EEG);

%% Check what events exist in your data
% Should be the flanker markers and the 2 StimTrak markers (S127 and S255)
disp('Events in your data:');
unique({EEG.event.type})

%% Epoch around S255 and S127 markers
epoch_markers = {'S255', 'S127'};  % READ Stim Markers
epoch_window = [-1 2];  % 1 second before to 2 seconds after (in seconds)

EEG_epoched = pop_epoch(EEG, epoch_markers, epoch_window, 'epochinfo', 'yes');
EEG_epoched = eeg_checkset(EEG_epoched);

disp(['Number of epochs created: ' num2str(EEG_epoched.trials)]);

%% Visualize the epochs - Method 1: Scroll through data
% This opens the EEGLAB data scroll viewer
pop_eegplot(EEG_epoched, 1, 1, 1);  % Shows all channels, all epochs
