% runMIBmeg_analysis
% run preproc and freqanalysis 

clear; close all

%%
% restoredefaultpath
basepath = '/path/'; % on the cluster
addpath(fullfile(basepath, 'MATLAB/tools/qsub_tardis'))
addpath(genpath(fullfile(basepath, 'MATLAB/critEEG_analysis')))    
rmpath(genpath(fullfile(basepath, 'MATLAB/critEEG_analysis/.git' )))
addpath(fullfile(basepath, 'MATLAB/tools/fieldtrip-20170611')) % fieldtrip-20170611
ft_defaults

%%
runcfg = [];
runcfg.batchlists = {
    'batch_critEEG_Erik_ses1'  ...
    % etc for each Subj/session
};

%% Peersetup preproc scripts
runcfg.overwrite = 1;
runcfg.sourceloc = 'no';

% % preproc runanalysis settings
runcfg.preproc.loaddata = 'no'; %load in data to do visual muscle rejection
runcfg.preproc.loaddatapath = '/home/mpib/kloosterman/projectdata/critEEG';

runcfg.preproc.artf_rejection = 'yes';
runcfg.preproc.artf_feedback = 'yes';
runcfg.preproc.loadartf = 'no';

if ismac
    runcfg.preproc.parallel = 'local'; %torque peer local qsublocal parfor
else
    runcfg.preproc.parallel = 'torque'; %torque peer local qsublocal parfor
end
runcfg.preproc.compile = 'no';

fprintf('Running critEEG preproc analysis . . .\n\n')
disp(runcfg.preproc.parallel); disp(runcfg.batchlists); disp(runcfg);

critEEG_preproc_peersetup

%% After preproc, we analyze per subject

SUBJ  = { 
    'Subj1'     ...
    'Subj2' ...
    'Subj3'    ...
    'Subj4'   ... only 2 sessions
    'Subj5'   ...
    'Subj6'    ...
    'Subj7'       'Subj8'  'Subj9' ...
    'Subj10'   ...
    'Subj11' ...
    'Subj12' ... % has no low crit for session 3
    'Subj13'    'Subj14'    'Subj15'   ...
    'Subj16'
};   

critEEG_runICA(SUBJ) %ICA for EOG. Run per session after concatenating runs
%
% % close all
critEEG_rejectICAcomps(SUBJ) % reject components
% %
critEEG_rejectEOGartf(SUBJ) % blinks during Figs, threshold artifacts, SCD computation


%% run freqanalysis for all bands, types etc.
triggers = {'stim' 'resp'};
freqtypes = {'totalpow' 'evoked' 'induced'};
freqbands =  {'low' 'high'};

for itrg=1:length(triggers)

    % timelock erp analysis
    critEEG_setup_timelockanalysis(SUBJ,  triggers{itrg})

    for ifrq=1:length(freqtypes)
        for iband=1:length(freqbands)
            disp( [triggers{itrg}, freqtypes{ifrq}, freqbands{iband}] )

            critEEG_freqanalysis_torque(SUBJ, triggers{itrg}, freqtypes{ifrq}, freqbands{iband})
            critEEG_concat_runs(SUBJ, triggers{itrg}, freqtypes{ifrq}, freqbands{iband})
        end
    end
end

%% Behavioral analysis
%
% critEEG_analysebehavioral(SUBJ) % run sorttrials
behav = critEEG_analyze_trialinfo(SUBJ) % collect trialinfo's and plot c, d' etc
critEEG_plot_behavior

critEEG_plot_hddmparas_JW % plot ddm behavior

%% LOAD in data for plotting
%% power and erp analysis
respavg = critEEG_load_respavg()

critEEG_plotTFR_poolings
critEEG_corr_gamma_vs_ddmDC

% plot TFR's of prestim period
critEEG_prestim_stats

%% link prestim alpha to poststim gamma, spectra all SDT 
critEEG_plot_psalpha_gain_gamma
