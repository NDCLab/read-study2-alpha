%Compute TF, ITPS, ICPS, and wPLI measures for EEG data
%Maureen Bowers 6/29/2021 based on scripts by Ranjan Debnath
%Modified to run on an HPC environment by Charles Knowlton April 2026

clear all;
clc;

%%
%%%%% PARALLEL POOL SETUP %%%%%
%%
 
% Create a local cluster object
cluster = parcluster('local');
 
% Start matlabpool with max workers set in the slurm file
if ~isempty(getenv('SLURM_CPUS_PER_TASK'))
    num_workers = str2num(getenv('SLURM_CPUS_PER_TASK'));
    parpool(cluster, num_workers)
    fprintf('Started parallel pool with %d workers from SLURM\n', num_workers);
else
    % parpool('local', 4);
    fprintf('No SLURM environment detected. Start parpool manually if desired.\n');
end
 
% Check pool status
pool = gcp('nocreate');
if isempty(pool)
    warning('No parallel pool running. Processing will be sequential.');
else
    fprintf('Parallel pool with %d workers is running.\n', pool.NumWorkers);
end

%%
%%%%% Setting paths %%%%%
%1. Data Location - specifically modified scripts for tf
data_location = '/home/data/NDClab/datasets/read-study2-dataset/derivatives/preprocessed/t-fModified';

%2. Save Data Location
save_location = '/home/data/NDClab/analyses/read-study2-alpha/derivatives/tf';

%3. Scripts Location
scripts_location = '/home/data/NDClab/analyses/read-study2-alpha/code/postprocessing/eeg/tf';

%4. Set EEGLab path
addpath(genpath('/home/data/NDClab/tools/lab-devOps/scripts/MADE_pipeline_standard/eeglab13_4_4b'));

%%
%%%%% Information about dataset and analysis procedures %%%%%%
%5. Number of channels
nbchan = 64;

%6. What kind of data? Resting State or Event-Related Data
RestorEvent = 0; %1 = rest, 0 = event

%7. What are your conditions of interest if using Event-Related Data? This
%naming convention should come from the Edit_events.m script provided.
%Note: not needed for resting state data.
Conds = {'resp_s_i_0','resp_s_i_1', 'resp_ns_i_1', 'stim_ns_i_0'};

%8. Minimum number of trials to analyze
mintrialnum =6; %If the participant does not have enough trials in a condition based on this cutoff, a "notenoughdata.mat" file will be saved into save_location.

%9. Would you like to baseline correct your data? NOTE: ICPS calculated over time will not be baseline corrected.
BaselineCorrect = 1; %1=Yes, 0=No

%10. What time period would you like to use to baseline correct
BaselineTime = [-400 -100];
%Put in time in ms for event-related data. For example, if you want to baseline correct from -100 to 0ms before the 
% event for evernt-related paradigms, put [-100 0].

%11. Would you like to downsample the output to 125Hz? This is done after the time-frequency computations and will have minor impacts on the resolution. 
%We recommend downsampling to reduce file size for ease of storage. All data will be initially downsampled to 250Hz. 
Downsample = 1; %0=No, 1=Yes

%12. Dataset Name - will be appended to the saved files
DatasetName = '_read_study2';

%%
%%%%% Settings for Time-Frequency Measures %%%%%

%13. Minimum and Maximum Frequency, Number of Frequency Bins, and range cycles to calculate complex Morlet wavelet decomposition
min_freq = 1;
max_freq = 30;
num_frex = 59; % number of frequency bins between minimum and maximum frequency
range_cycles = [3 10]; % wavelet cycles: min 3 max 10


%%
%%%%% Questions about phase-based measures %%%%%

%14. Would you like to calculate inter-trial phase synchrony (ITPS) in addition to TF?
ITPS_calc = 0; %(1=yes,0=no)

%15. Would you like to subsample trials? This is recommended for event-related paradigms, 
% especially when there are uneven numbers of trials in conditions. 
Subsample = 1; %1=Yes, 0=No
%How many trials to pull for each subsample? 
NumTrialsPulled = 6;
%How many times to do subsampling? We recommend doing at least 10 subsamples and to have the possiblility of using all your data
% (e.g., if you have 150 trials, do 15 subsamples of 10 trials).
NumSubsamples = 100;

%16. Would you like to calculate inter-channel phase synchrony (ICPS or wPLI) in addition to TF?
ICPS_calc = 0; %(1=yes,0=no) 

%17. Would you like to caluclate coherence or weighted phaselagidx?
ICPS_or_wPLI = 1; %(1=coherence, 0=wPLI)

%18. Inter-channel phase synchrony over trials or connectivity over time?
TimeOrTrials = 1; %0 = over time, 1 = over trials
% NOTE: over time calculations will not be able to be subsampled or
% baseline corrected

%19. Type of Connectivity to compute 
ConnectType = 1; %(0= all-to-all connectivity; 1=seed-based connectivity)
%If seed based, choose seed and which electrodes to compute connectivity:
Seed = '1';
Elecs4Connect = { '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '34' '35' '36' '37' '38' '39' '40' '41' '42' '43' '44' '45' '46' '47' '48' '49' '50' '51' '52' '53' '54' '55' '56' '57' '58' '59' '60' '61' '62' '63' '64' }; 

% 20. Create List of subjects to loop through
subnum = dir([data_location filesep '*.set']); % Use regex to find your files 
subject= {subnum.name};
for ii=1:length(subject)
    subject_list{ii}=subject{ii};   
end

% We recommend checking that your subject list looks like you would expect:
% disp('SUBJECTS TO BE PROCESSED:')
% subject_list

eeglab % Loading EEGLAB

% Add scripts location to path for all workers
addpath(scripts_location);

%%%%%%%%%%%%%%%%%%%% COMPUTATIONS BEGIN BELOW HERE %%%%%%%%%%%%%%%%
TrialNums_all = cell(length(subject_list), 1);
%% loop through all subject
parfor sub=1:length(subject_list)
    try
        subjStart = tic;
        % Initialize objects for this participant:
        timefreqs_data = [];
        phase_data=[];
        ITPS_all=[];
        ICPS_all=[];
        wPLI_all=[];
        EEG=[];
        TrialNums = struct();

        subject = subject_list{sub};
        fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', sub, subject);

        % Load data
        EEG=pop_loadset('filename', [subject], 'filepath', data_location);
        EEG = pop_selectevent( EEG, 'latency','-.1 <= .1','deleteevents','on');
    
        %Downsample to 250Hz
        %EEG = pop_resample(EEG, 250);
        EEG = eeg_checkset( EEG );
    
        %Keep only markers of interest
        EEG = pop_selectevent( EEG, 'Condition',Conds, 'deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = eeg_checkset( EEG );
    
        %save times
        time = EEG.times;
    
        %save channel locations
        channel_location = EEG.chanlocs;
    
        %% Compute complex morlet wavelet time frequency decomposition
        TrialNums = timefreq_script(EEG, sub, subject, save_location, ...
        scripts_location, Conds, mintrialnum, BaselineCorrect, BaselineTime, ...
        Downsample, DatasetName, min_freq, max_freq, num_frex, range_cycles, ...
        ITPS_calc, Subsample, NumTrialsPulled, NumSubsamples, ICPS_calc, ...
        ICPS_or_wPLI, TimeOrTrials, ConnectType, Seed, Elecs4Connect, ...
        nbchan, RestorEvent, time, channel_location);
    catch ME
        fprintf('\n\nERROR: Processing failed for subject %d (%s)\n', sub, subject_list{sub});
        fprintf('Error message: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('Error in: %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
        end
        diary off
        continue %Skip participant
    end
end

% Combine TrialNums from all subjects
TrialNums = struct();
counter = 1;
for sub = 1:length(subject_list)
    if ~isempty(TrialNums_all{sub})
        % Copy all fields from this subject's TrialNums
        subj_trials = TrialNums_all{sub};
        for idx = 1:length(subj_trials)
            TrialNums(counter) = subj_trials(idx);
            counter = counter + 1;
        end
    end
end

%Save out table with trial numbers
if ~isempty(fieldnames(TrialNums))
    TrialNums_table = struct2table(TrialNums);
    writetable(TrialNums_table,[save_location filesep 'TrialNums.csv']);
    fprintf('\n\nTrialNums table saved to: %s\n', [save_location filesep 'TrialNums.csv']);
else
    warning('No TrialNums data to save.');
end
 
fprintf('\n\n*** All subjects processed! ***\n\n');