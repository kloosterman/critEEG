function [respavg] = critEEG_load_respavg()
% load freqanalysis and timelockanalysis data into memory, concat over
% subjects

basepath = '/path'; % on the cluster

%% Set up arrays
inputfolder = 'freq'; 
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
    'Subj13'    'Subj14'   ... Subj14 lowfreq noise
     'Subj15'   ...
    'Subj16'
};   

disp(length(SUBJ))
PREOUT = '/plots';
respavg = [];
respavg.SUBJ = SUBJ;
respavg.PREOUT = PREOUT;

%totalpow settings
% trigger_leg = {'stim' 'resp'}; %  'stim' resp
trigger_leg = {'stim' }; %  'stim' resp
freqtype = 'totalpow';  %evoked induced totalpow
% freqtype = 'induced';  %evoked induced totalpow
% freqtype = 'induced_robdetr';  %evoked induced totalpow
freq_leg = {'low' 'high'}; % low high
% TIMLO(1,1) = -0.8; %stim
% TIMHI(1,1) = 1.5; 
TIMLO(1,1) = -0.8; %stim
TIMHI(1,1) = 0.9; 
TIMLO(2,1) = -0.5;%resp
TIMHI(2,1) = 0.75; %dat{itrig}.TIMHI;

respavg.trigger_leg = trigger_leg;
respavg.freqtype = freqtype;
respavg.freqband = freq_leg;

FREQLO(1,1) = 3; %lowfreq
% FREQHI(1,1) = 30;
FREQHI(1,1) = 35;
FREQLO(2,1) = 36; % highfreq
FREQHI(2,1) = 100;
% FREQHI(2,1) = 127;
% FREQHI(2,1) = 110;

frind = {};
faxis = {};
tind = {};
taxis = {};
for itrig = 1:length(trigger_leg)
    for iband = 1:length(freq_leg)
        PREIN = fullfile(basepath, 'projectdata/critEEG', inputfolder, trigger_leg{itrig});
        
        subjdir = fullfile(PREIN, SUBJ{1});
        fileload = fullfile(subjdir, sprintf('%s_%s_%s_%sfreq_respavg.mat', SUBJ{1}, trigger_leg{itrig}, freqtype, freq_leg{iband}));

        dat = load(fileload);
        
        frind{iband} = find((dat.faxis >= FREQLO(iband,1)) & (dat.faxis <= FREQHI(iband,1))); % select frind of interest
        faxis{iband} = dat.faxis(frind{iband});
        
        if iband == 1
            tind{itrig} = find((dat.taxis >= TIMLO(itrig,1)) & (dat.taxis <= TIMHI(itrig,1)));
            taxis{itrig} = dat.taxis(tind{itrig});
            % for ERP:
            tind2{itrig} = find((dat.taxis2 >= TIMLO(itrig,1)) & (dat.taxis2 <= TIMHI(itrig,1)));
            taxis2{itrig} = dat.taxis2(tind2{itrig});
        end
    end    
end

taxis_lengths = cellfun(@length, taxis);
faxis_lengths = cellfun(@length, faxis);
respavg.freq = faxis; %[faxis{1} faxis{2}];
respavg.time = taxis; %[taxis{1} taxis{2}];

respavg.time_erp = taxis2; %[taxis{1} taxis{2}];

taxis2_lengths = cellfun(@length, taxis2);

chlabel = dat.chlabel;
respavg.label = chlabel;
respavg.sens = critEEG_sensorselection();

% Condition labels
respavg.stim_conds = {'Fig' 'Hom' 'allstim' 'Fig - Hom' };
respavg.resp_conds = {'Report' 'Noreport' 'Reportcomb' 'Report - Noreport' };
% respavg.behav_conds = {'Highcrit' 'Lowcrit' 'allbehavconds' 'High-Lowcrit'};
respavg.behav_conds = {'Conservative' 'Liberal' 'allbehavconds' 'Liberal-Cons'};
respavg.sdt_conds = {'Hit' 'Miss' 'Hit+Miss' 'Hit-Miss'; ...
    'FA'  'CR'  'FA+CR' 'FA-CR' ; ...
    'Hit+FA' 'Miss+CR' 'allSDT' 'Hit+FA-Miss+CR'; ...
    'Hit-FA' 'Miss-CR'  'Hit-FA+Miss-CR' 'Hit-FA-Miss-CR'};

behavdat = cell(size(SUBJ));
respavg.dat = nan(length(SUBJ),length(chlabel),max(faxis_lengths), max(taxis_lengths),2,2, ...
    4,4,4,4, 'single');
respavg.dimord = 'subj_chan_freq_time ... band_trig ... ses_cond_stim_resp';  
respavg.dimordsize = size(respavg.dat);

respavg.basespec = nan(length(SUBJ), 4, length(chlabel),max(faxis_lengths),2, 2,2,2, 'single'); %subj label freq band cond stim resp

respavg.pow = nan(length(SUBJ),length(chlabel),max(faxis_lengths), max(taxis_lengths),2,2, ...
    4,4,4,4, 'single');

respavg.dat_erp = nan(length(SUBJ),length(chlabel), max(taxis2_lengths),2, ...
    4,4,4,4, 'single');
respavg.dimord_erp = 'subj_chan_time ... trig ... ses_cond_stim_resp';  


for isub = 1:length(SUBJ)
    for iband = 1:length(freq_leg)
        for itrig = 1:length(trigger_leg)
            PREIN = fullfile(basepath, 'projectdata/critEEG', inputfolder, trigger_leg{itrig});
            
            subjdir = fullfile(PREIN, SUBJ{isub});
            fileload = fullfile(subjdir, sprintf('%s_%s_%s_%sfreq_respavg.mat', SUBJ{isub}, trigger_leg{itrig}, freqtype, freq_leg{iband}));
            fprintf('loading %s %s %s\n', fileload, trigger_leg{itrig}, freq_leg{iband})

            temp = load(fileload);

            respavg.dat(isub,:,1:faxis_lengths(iband),1:taxis_lengths(itrig), iband,itrig, :,1:3, 1:3, 1:3) = ...
                temp.subjrespavg(:,frind{iband},tind{itrig}, :,:,:,:);

            if iband == 1
                respavg.dat_erp(isub,:,1:taxis2_lengths(itrig), itrig, :,1:3, 1:3, 1:3) = ...
                    temp.subj_erpavg(:,tind2{itrig}, :,:,:,:);
            end

            if itrig == 1                
%                 respavg.basespec(isub,:,1:faxis_lengths(iband), iband, :,:,:) = ...
%                     squeeze(nanmean(temp.basespec_ses(:,:,frind{iband},1:2,1:2,1:2),1)); % chan freq band cond stim resp
                respavg.basespec(isub,1:3,:,1:faxis_lengths(iband), iband, :,:,:) = ...
                    squeeze(temp.basespec_ses(:,:,frind{iband},1:2,1:2,1:2)); % chan freq band cond stim resp
                % load single trial power, concat
                respavg.pow_singletrial{isub,iband} = temp.pow_singletrial;
                % raw power
                respavg.pow(isub,:,1:faxis_lengths(iband),1:taxis_lengths(itrig), iband,itrig, :,1:3, 1:3, 1:3) = ...
                    temp.subjpowavg(:,frind{iband},tind{itrig}, :,:,:,:);

            end
        end
    end
end

fprintf('Calculating averages and contrasts . . .\n')
respavg.dat = respavg.dat .* 100; % make psc
respavg.dat(:,:,:,:, :,:, 4,:,:,:) = nanmean(respavg.dat(:,:,:,:, :,:, 1:3,:,:,:), 7); % ses avg
respavg.dat(:,:,:,:, :,:, :,3,:,:) = nanmean(respavg.dat(:,:,:,:, :,:, :,1:2,:,:), 8); % cond avg
respavg.dat(:,:,:,:, :,:, :,:,3,:) = nanmean(respavg.dat(:,:,:,:, :,:, :,:,1:2,:), 9); % stim avg
respavg.dat(:,:,:,:, :,:, :,:,:,3) = nanmean(respavg.dat(:,:,:,:, :,:, :,:,:,1:2), 10); % resp avg

respavg.dat(:,:,:,:, :,:, :,4,:,:) = respavg.dat(:,:,:,:, :,:, :,2,:,:) - respavg.dat(:,:,:,:, :,:, :,1,:,:); %behavcond diff

respavg.pow(:,:,:,:, :,:, 4,:,:,:) = nanmean(respavg.pow(:,:,:,:, :,:, 1:3,:,:,:), 7); % ses avg

respavg.basespec(:,:,:,:,:, 3,:,:)  = nanmean(respavg.basespec(:,:,:,:,:, 1:2,:,:),6);
respavg.basespec(:,:,:,:,:, 4,:,:)  = respavg.basespec(:,:,:,:,:, 1,:,:) - respavg.basespec(:,:,:,:,:, 2,:,:); % cond

respavg.basespec(:,4,:,:,:, :,:,:)  = nanmean(respavg.basespec(:,:,:,:,:, :,:,:), 2);

respavg.dat_erp(:,:,:, :, 4,:,:,:) = nanmean(respavg.dat_erp(:,:,:, :, 1:3,:,:,:), 5); % ses avg

% put behavior in
behavin = fullfile(basepath, 'projectdata', 'critEEG', 'behavior');
load(fullfile(behavin, 'behavstruct.mat'));
respavg.behavior = behav;
respavg.behavior.dimord = 'subj_ses_cond';

% add DDM 10/5 bins results
respavg.ddmpars10bins = critEEG_load_binned_paras(10);
respavg.ddmpars5bins = critEEG_load_binned_paras(5);
respavg.ddmpars3bins = critEEG_load_binned_paras(3);
