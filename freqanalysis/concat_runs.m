function [ subjrespavg ] = concat_runs( cfg )
%Concatenate runs for each subject and save subjrespavg
%   V2: have separate baselines for NoMiss and NoFA
fprintf('%s; %s; %s\n', cfg.trigger, cfg.freqtype, cfg.freqband)

PREIN = cfg.PREIN;
disp(PREIN)
PREOUT = cfg.PREOUT;
cd(PREIN);
cd('ses2')
[~,SUBJ] = fileparts(PREIN);
examplefreq = dir(fullfile(PREIN, 'ses2', sprintf('%s_ses%d*_%s_%sfreq.mat', SUBJ,  2, cfg.freqtype, cfg.freqband )));
load(examplefreq.name)
exampletimelock = dir(fullfile(PREIN, 'ses2', sprintf('%s_ses%d*timelock.mat', SUBJ,  2 )));
load(exampletimelock.name)

% respavgout = fullfile(PREIN, 'respavg');
% mkdir(respavgout)

% set up arrays for analysis
%--------------------------------------------------------------------------
switch cfg.trigger
    case 'stim'
        basetind = find((freq.time>= -0.4) & (freq.time<= 0));
        TIMLO = -1; TIMHI = 2;
    case 'resp'
        %         basetind = find((freq.time>= -0.75) & (freq.time<= -0.5));
        TIMLO = -1.25; TIMHI = 1.25; % resp
end
freq.time=round(freq.time*100)/100; 
tind = find((freq.time>=TIMLO) & (freq.time<=TIMHI));
% tind = find(freq.time);
taxis = freq.time(tind);

CUTLO  = 0;  CUTHI  = 200;       % cutoff for reading in data
frind = find((freq.freq>=CUTLO) & (freq.freq<=CUTHI));
faxis = freq.freq(frind);

chlabel = freq.label;

%timelock time ind
tind2 = find((timelock.time>=TIMLO) & (timelock.time<=TIMHI));
taxis2 = timelock.time(tind2);

% 1     2    3      4        5          6       7
% subj nchan nfreq ntimebins cond(2)  stim(2) resp(2)
subjrespavg = nan( length(chlabel),length(frind), length(taxis), 4,3,3,3, 'single' ); % ses cond stim resp
subjpowavg = subjrespavg;

pow_singletrial = cell( 10, 4,3,3,3); % tfoi ses cond stim resp

switch cfg.freqband
    case 'low'
        subj_erpavg = nan( length(chlabel), length(taxis2), 4,3,3,3, 'single' ); % ses cond stim resp
end

powdat_ses = cell(4,1);
respdat_ses = cell(4,1);
trialinfo_ses = cell(4,1);
basespec_ses = nan([3 length(chlabel) length(frind) 2 2 2 ]); % ses chan freq cond stim resp
timelock_ses = {};

for ises = 1:3
    sesdir = fullfile(PREIN, sprintf('ses%d', ises));
    if ~exist(sesdir, 'dir')
        fprintf('%s not found\n', sesdir)
        continue
    end
    
    cd(sesdir)
    
    trialinfo = [];
    fprintf('\n\nSubject session directory: %s  . . . \n', sesdir)
    fprintf('Concatenating runs Subject %s Session %d: . . .\n', SUBJ, ises)
    
    % load freq
    runfile = dir( sprintf('%s_ses%d*_%s_%sfreq.mat', SUBJ,  ises, cfg.freqtype, cfg.freqband ));
    if isempty(runfile);     break;            end
    fprintf('Loading run %s . . . ', runfile.name)
    load(runfile.name); % subj
    
    fprintf('%d trials\n', size(freq.powspctrm,1))
    powdat = freq.powspctrm; % convert to single after normalization
    trialinfo = freq.trialinfo;
    ntrl = size(powdat,1);
    
    % get single trial prestim alpha power
    sens = critEEG_sensorselection();
    cfg2 = [];
    cfg2.avgoverchan = 'yes';
    cfg2.avgoverfreq = 'yes';
    cfg2.avgovertime = 'yes';
    switch cfg.freqband
        case 'low'
            for isoi = 1:length(sens.ind)
                cfg2.channel = sens.ind{isoi}; % append vectors for soi
                cfg2.frequency = [8 12];
                %                 cfg2.latency = [-0.4 -0.2];
                cfg2.latency = [-0.8 -0.2];
                cfg2.nanmean = 'yes';
                tempfreq = ft_selectdata(cfg2, freq);
                logdat = log(tempfreq.powspctrm);
                
                % figure; histogram(logdat)                
                zdat = (logdat - nanmean(logdat)) ./ nanstd(logdat);
                % figure; histogram(zdat)

                tempfreq.powspctrm(zdat > 3) = NaN; % set outliers to NaN
%                 freq_pre_a{ises}(:,isoi) = tempfreq.powspctrm / max(tempfreq.powspctrm); %normalize per session
                freq_pre_a{ises}(:,isoi) = tempfreq.powspctrm; % do not normalize
            end
            % load timelock
            runfile = dir( sprintf('%s_ses%d*timelock.mat', SUBJ, ises));
            if isempty(runfile);     break;            end
            fprintf('Loading run %s . . . ', runfile.name)
            load(runfile.name); % subj
            
            % do single trial baseline correction TODO get baseline from
            % stimlocked for resplocked FOR TIME DOMAIN data, ie erp's
            cfg2 = [];
            cfg2.baseline     = [-0.2 0];
            cfg2.parameter    = 'trial';
            % cfg.parameter    = {'trial', 'avg'};
            timelock_ses{ises} = ft_timelockbaseline(cfg2, timelock);
    end
    
    % compute or load basespec
        basepath = '/path';
    
    file = sprintf('%sfreq_basespec.mat', cfg.freqband);
    tok = tokenize(PREIN, '/'); % get input folder
    inputfolder = tok{end-2}; % before SUBJ and trigger folder
    filesave = fullfile(basepath, inputfolder, 'stim', SUBJ, sprintf('ses%d', ises), file);
    
    if strcmp(cfg.trigger, 'stim') && (strcmp(cfg.freqtype, 'totalpow') || strcmp(cfg.freqtype, 'induced_robdetr'))  % % %         compute basespec within session and condition
        fprintf('Computing baseline spectrum . . .\n')
        basespec_trials = nanmean(powdat(:,:,:,basetind),4);     %average over basetind timebins
        basespec = squeeze(nanmean(basespec_trials));
        if ~strcmp(cfg.freqtype, 'induced_robdetr') % use own baseline for  induced_robdetr
            fprintf('Saving %s . . . \n', filesave)
            save(filesave, 'basespec', 'basespec_trials') %
        end
    elseif strcmp(cfg.trigger, 'stim') && strcmp(cfg.freqtype, 'evoked')  % % %         compute basespec within session and condition
        fprintf('Loading %s . . . \n', filesave)
        load(filesave) % use for denominator basespec from totalpow
        basespec_trials = nanmean(powdat(:,:,:,basetind),4);     %use for subtraction of evoked data: average over basetind timebins
    else
        %         basespec_file = dir(sprintf('*totalpow*%sfreq_basespec.mat', cfg.freqband));
        fprintf('Loading %s . . . \n', filesave)
        load(filesave) % for basespec (denominator)
        %         basespec_trials = nanmean(powdat(:,:,:,basetind),4);     % subtraction (nominator) is taken per condition
    end
      
    %normalize
    global_baseline = 1;
    respdat = nan(ntrl,length(chlabel),length(frind),length(taxis));
    respdat_pertrial = nan(ntrl,length(chlabel),length(frind),length(taxis));
    for icond = 1:2
        for istim = 1:2
            for iresp = 1:2
                cond_ind = trialinfo(:,1) == icond & trialinfo(:,2) == istim & trialinfo(:,3) == iresp;  % TODO add RT stratification: ... &
                
                if global_baseline
                    basespec_within = squeeze(nanmean(basespec_trials, 1));
                else
                    basespec_within = squeeze(nanmean(basespec_trials(cond_ind,:,:), 1));
                end
                basespec_ses(ises, :,:, icond, istim, iresp) = basespec_within; % saved below
                
                trialind = find(cond_ind);
                if length(trialind) > 0
                    fprintf('Normalizing %d trials . . .\n', length(trialind))
                else
                    fprintf('No trials present for this condition. . .\n')
                    continue
                end
                ft_progress('init', 'etf',     'Please wait...');
                for ich = 1:length(chlabel)
                    ft_progress(ich/length(chlabel), 'Processing channel %d from %d', ich, length(chlabel));
                    basedat = squeeze(repmat(squeeze(basespec(ich,:)),[1,size(taxis)]));
                    basedat_within = squeeze(repmat(squeeze(basespec_within(ich,:)),[1,size(taxis)]));
                    for itrial = trialind'
                        respdat(itrial,ich,:,:) = single( (squeeze(powdat(itrial,ich,frind,tind)) - basedat_within) ./ basedat ); % take out baseline differences between NOFA and Nomiss
                        %                         respdat(itrial,ich,:,:) = single( (squeeze(powdat(itrial,ich,frind,tind)) - basedat) ./ basedat );
                        % single trial baseline correction
                        respdat_pertrial(itrial,ich,:,:) = single( squeeze(powdat(itrial,ich,frind,tind) - basespec_trials(itrial,ich,frind)) ./ basedat ); % take out baseline differences between NOFA and Nomiss
                    end
                end
                ft_progress('close')
                
                
                %     %             %% ff plotten
                %             close all
                %             NKsensorselection
                %             test = nanmean(respdat(:,occind,:,:),2);
                %             test = squeeze(nanmean(test));
                %             figure
                %             imagesc(taxis,faxis,test, [-0.05 0.05]);
                %             colorbar
                %             set(gca,'Box','off','XTick',[-1 -0.5 0 0.5 1],...    [-0.5,0,0.5,1,1.5,2,2.5]
                %                 'YDir','normal','YTick',[0:20:140],...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
                %                 'TickDir','out', 'FontSize', 12);
                
            end
        end
    end
    
    % get single trial stimrelated responses
    %     , + SSVEP (or poststim gamma, ERP, alpha
    %     % suppression?
    
%     freq.powspctrm = respdat_pertrial;
%     freq.powspctrm = respdat; % normalized by all trial average
    freq.powspctrm = powdat; % just raw power
%     clear respdat_pertrial
    
    cfg2 = [];
    cfg2.avgoverchan = 'yes';
    for isoi = 1:length(sens.ind)
        cfg2.channel = sens.ind{isoi}; % append vectors for soi
        switch cfg.freqband
            case 'low'
                cfg2.avgovertime = 'yes';
                cfg2.avgoverfreq = 'yes';
                
                % poststim low freq suppression
                cfg2.frequency = [8 15];
                cfg2.latency = [0.2 0.6];
                tempfreq = ft_selectdata(cfg2, freq);
                freq_post_ab{ises}(:,isoi) = tempfreq.powspctrm;
                
                % ssvep
                cfg2.frequency = [23 27];
                cfg2.latency = [0.2 0.6];
                tempfreq = ft_selectdata(cfg2, freq);
                freq_post_ssvep{ises}(:,isoi) = tempfreq.powspctrm / max(tempfreq.powspctrm);
                
                % take P1 of the ERP,
                cfg2.avgovertime = 'no';
                cfg2.avgoverfreq = 'no';
                cfg2.latency = [0.15 0.25];
                temptimelock = ft_selectdata(cfg2, timelock_ses{ises});
                timelock_P1{ises}(:,isoi) = max(temptimelock.trial, [], 3); % baseline to peak
                %                 timelock_P1{ises}(:,isoi) = timelock_P1{ises}(:,isoi) ./ max(timelock_P1{ises}(:,isoi)); % normalize
                
                % take N2 of the ERP
                cfg2.latency = [0.3 0.5];
                temptimelock = ft_selectdata(cfg2, timelock_ses{ises});
                timelock_N2{ises}(:,isoi) = min(temptimelock.trial, [], 3); % baseline to trough
                %                 timelock_N2{ises}(:,isoi) = timelock_N2{ises}(:,isoi) ./ min(timelock_N2{ises}(:,isoi));  % normalize
                
                % take CPP of the ERP
                cfg2.avgovertime = 'yes';
                cfg2.latency = [0.2 0.7];
                temptimelock = ft_selectdata(cfg2, timelock_ses{ises});
                timelock_CPP{ises}(:,isoi) = temptimelock.trial; % mean
                %                 timelock_CPP{ises}(:,isoi) = timelock_CPP{ises}(:,isoi) ./ max(timelock_CPP{ises}(:,isoi));
                
            case 'high'
                cfg2.avgovertime = 'yes';
                cfg2.avgoverfreq = 'yes';
                
                cfg2.frequency = [42 58];
                cfg2.latency = [0.2 0.4];
                tempfreq = ft_selectdata(cfg2, freq);
                freq_post_ssvep_high{ises}(:,isoi) = tempfreq.powspctrm;
                
                %                 cfg2.frequency = [60 100]; % old
                
                %                 cfg2.frequency = [59 110];  %works well, but keep under 100 Hz
                %                 cfg2.latency = [0.2 0.4];
                cfg2.frequency = [59 100];
                cfg2.latency = [0.2 0.4];
                tempfreq = ft_selectdata(cfg2, freq);
                freq_post_gamma{ises}(:,isoi) = tempfreq.powspctrm;
                
                cfg2.frequency = [36 100];
                cfg2.latency = [0.2 0.4];
                tempfreq = ft_selectdata(cfg2, freq);
                freq_post_gamma_all{ises}(:,isoi) = tempfreq.powspctrm;
                
        end
    end
    
    % save raw powdat also
    powdat_ses{ises} = powdat;
    %     powdat_ses{ises} = (powdat - nanmean(powdat)) ./ nanstd(powdat);    % zscore powdat
    
    % normalization: subtract and divide by session average
        respdat_ses{ises} = respdat; % concat sessions; each session has its own basespec
        respdat_ses{4} = [respdat_ses{4}; respdat]; % collapse over sessions
    
    trialinfo_ses{ises} = trialinfo; % concat sessions; each session has its own basespec
    trialinfo_ses{4} = [trialinfo_ses{4}; trialinfo];
    
    %%     % export csv's for HDDM and single trial analyses
    % 1 condition: Nomiss of Nofa
    % 2 signal (fig or hom)
    % 3 resp (yes or no)
    % 4 RT
    % 5 trial index
    % 6 signal orientation
    % 7 ses nr
    % 8 run nr
    
    switch cfg.freqband
        case 'low'
            exp_path = fullfile(basepath, 'projectdata/critEEG/export');
            mkdir(exp_path)
            trl_outfile = sprintf('%s_ses%d.csv', lower(SUBJ), ises);

            % columns: subj_idx, response, RT, session_nr, run_nr, cond
            SUBJlist  = {     'Erik'     'Esther'     'Greke'        'Isabelle'     'Jesse'       'Lotte'        'Matthijs'       'Nicoletta'  'Niels'     'Rebecca'       'Rein'     'Rinske'    'Ruben'    'Stephen'    'Vera'       'Zimbo'};
%             outmat = zeros(size(trialinfo,1),17); % if we put binno's in boolean, 10 alpha bins
            outmat = zeros(size(trialinfo,1),9); %
            outmat(:,1) = find(strcmp(SUBJlist, SUBJ)) - 1;

%             % accuracy coding:
%             outmat(:,2) = (trialinfo(:,2) == 1 &  trialinfo(:,3) == 1) | ...
%                 (trialinfo(:,2) == 2 &  trialinfo(:,3) == 2); %0 incorrect, 1 correct
            % stimulus codin: response col should be choice: (1=yes, 2=no resp -> flip that so 0=no resp and 1=yes):
            outmat(:,2) = mod(trialinfo(:,3),2);

            outmat(:,3) = trialinfo(:,4) / 256; %RT
            outmat(trialinfo(:,4)==-1, 3) = NaN; % replace M and CR RT with NaN
            outmat(:,4) = trialinfo(:,7);
            outmat(:,5) = trialinfo(:,8);
            outmat(:,6) = trialinfo(:,1) - 1;
            outmat(:,7) = trialinfo(:,6);
            outmat(:,8) =  mod(trialinfo(:,2),2); %  % add stim column (1=fig, 2=hom -> flip that so 0=hom and 1=fig):

            % make bins based on occipital alpha
            freq_occalpha = log(freq_pre_a{ises}(:,1));
            %
            %     % plot alpha distribution
            %     %             figure; histogram(freq_alpha.powspctrm, 100); title(filename)
            %     %             ax=gca;
            %     %             ax.FontSize =18;
            %     %             xlabel('log alpha power')
            %     %             ylabel('Frequency of occurrence')
            %     %             print out -dpdf
            %     % 5 alpha bins: 30%, 15%, 10%, 15%, 30%
            range = linspace(min(freq_occalpha), max(freq_occalpha), 100);

            %             bin_ranges = [0 33 67 100];
            %             bin_ranges = 0:10:100;

            % actual binning done for gamma:
%             bin_ranges = [0    15;   5    25;   15    35;    25    45;    35    55;    45    65;   55 75; 65 85; 75 95;  85 100];

            % 5 bins, 50% overlap:
            bin_ranges = [0 33;  16 50; 33 67; 50 83; 66 100 ];

            % make boolean cols indicating if trial belongs to bin or not.
            % Forces freeing all parameters in HDDM, not good
            % %             outmat(:,8) = discretize(freq_occalpha, bin_edges, 0:length(bin_edges)-2 );
            %             for ibin = 1:length(bin_ranges)
            %                 bin_edges = prctile(range, bin_ranges(ibin,:));
            %                 temp = discretize(freq_occalpha, bin_edges, 1 );
            %                 temp(isnan(temp)) = 0; % set nans to 0
            %                 outmat(:,7+ibin) = temp;
            %             end

            % jw: Je zou ook de csv zou kunnen organiseren dat voor de trials die in twee bins voorkomen je
            % gewoon twee regels hebt, en één kolom voor alpha_bin (0-4).
            outmat2 = [];
            for ibin = 1:length(bin_ranges)
                bin_edges = prctile(range, bin_ranges(ibin,:));
                trl = discretize(freq_occalpha, bin_edges, 1 );
                trl(isnan(trl)) = 0; % set nans to 0
                outmat(find(trl),9) = ibin-1;
                outmat2 = [outmat2; outmat(find(trl),:)];
            end



            %             % remove zero RT trials, NK11, bug???
            %             figure; hist( outmat(:,3), 100)
            %             outmat = outmat(outmat(:,3) > 0.1 ,:);


            % write header and matrix to file
            fid = fopen(fullfile(exp_path, trl_outfile), 'w') ;
            fprintf(fid, 'subj_idx,response,rt,session_nr,run_nr,cond,orientation,stim_col,occ_alphabin\n');
%             fprintf(fid, 'subj_idx,response,rt,session_nr,run_nr,cond,orientation,occalphabin0,occalphabin1,occalphabin2,occalphabin3,occalphabin4,occalphabin5,occalphabin6,occalphabin7,occalphabin8,occalphabin9,stim_col\n');
%             dlmwrite(fullfile(exp_path, trl_outfile), outmat, '-append', 'precision', 10)
            dlmwrite(fullfile(exp_path, trl_outfile), outmat2, '-append', 'precision', 10); % outmat2 !!!
            fclose(fid);
    end
    clear powdat basespec respdat
end

%% Split up by condition and average over trials
% freq.trialinfo   % 1  condition% 2 stim (fig or hom) % 3 orientation % 4 RT (-1 for no bp) % 5 trial index

% 1     2    3      4        5          6       7
% subj nchan nfreq ntimebins cond(2)  stim(2) resp(2)
fprintf('\n\n\nSplit up by condition and average over trials\n')
ctr=0;
conds = 1:3;
stims = 1:3;
resps = 1:3;
ft_progress('init', 'etf',     'Please wait...');
nconds=length(conds) * length(stims) * length(resps);
for ises = 1:3
    for istim = stims % 1 = left, 2 = right
        for iresp = resps % left or right
            for icond = conds % NoMiss NoFA
                ctr=ctr+1;
                if isempty(trialinfo_ses{ises})
                    continue
                end
                cond_ind = trialinfo_ses{ises}(:,1) == icond;
                stim_ind = trialinfo_ses{ises}(:,2) == istim;
                resp_ind = trialinfo_ses{ises}(:,3) == iresp;
                if icond==3, cond_ind(:)=true; end
                if istim==3, stim_ind(:)=true; end
                if iresp==3, resp_ind(:)=true; end

                trial_inds = cond_ind & stim_ind & resp_ind;
                ft_progress(ctr/nconds, 'Processing cond %d stim %d resp %d: %d trials', icond, istim, iresp, length(find(trial_inds)));

                if length(find(trial_inds)) < 1
                    warning('No trials remain for this condition')
                    continue
                end
                % TODO get strongest 20 % of the power
                subjrespavg(:,:,:, ises, icond, istim, iresp) = nanmean(respdat_ses{ises}(trial_inds,:,:,:),1);

                subjpowavg(:,:,:, ises, icond, istim, iresp) = nanmean(powdat_ses{ises}(trial_inds,:,:,:),1);
%                 subjpowavg(:,:,:, ises, icond, istim, iresp) = (powdat_ses{ises}(trial_inds,:,:,:) - nanmean(powdat_ses{ises}(trial_inds,:,:,:))) ...
%                     ./ nanstd(powdat_ses{ises}(trial_inds,:,:,:));   % zscore

                switch cfg.freqband % put single trial data in cell TODO first cell can go into the matrix such that dimord rpt soi tfoi
                    case 'low'
                        tind_temp = timelock_ses{ises}.time >= TIMLO & timelock_ses{ises}.time <= TIMHI; % incase not all timepoints are present
                        subj_erpavg(:,:, ises, icond, istim, iresp) = nanmean(timelock_ses{ises}.trial(trial_inds,:,tind_temp),1);
                        % TODO compute var for quenching analysis

                        pow_singletrial{1, ises, icond, istim, iresp} = freq_pre_a{ises}(trial_inds,:);
                        pow_singletrial{2, ises, icond, istim, iresp} = freq_post_ab{ises}(trial_inds,:);
                        pow_singletrial{3, ises, icond, istim, iresp} = freq_post_ssvep{ises}(trial_inds,:);

                        pow_singletrial{4, ises, icond, istim, iresp} = timelock_P1{ises}(trial_inds,:);
                        pow_singletrial{5, ises, icond, istim, iresp} = timelock_N2{ises}(trial_inds,:);
                        pow_singletrial{6, ises, icond, istim, iresp} = timelock_CPP{ises}(trial_inds,:);

                    case 'high'
                        pow_singletrial{7, ises, icond, istim, iresp} = freq_post_ssvep_high{ises}(trial_inds,:);
                        pow_singletrial{8, ises, icond, istim, iresp} = freq_post_gamma{ises}(trial_inds,:);
                        pow_singletrial{9, ises, icond, istim, iresp} = freq_post_gamma_all{ises}(trial_inds,:);
                end
                % add per trial behavior (stim and resp given) to compute c and d' per alpha bin
                pow_singletrial{10, ises, icond, istim, iresp} = trialinfo_ses{ises}(trial_inds,2:3);
                
                if ises == 3 % concatenate sess
                    for im = 1:10
                        pow_singletrial{im, 4, icond, istim, iresp} = cat(1, pow_singletrial{im, :, icond, istim, iresp} );
                    end
                end

            end
        end
    end
    ft_progress('close')
end

clear powdat respdat_ses powdat_ses timelock_ses timelock

% save respavg per subj
filesave = fullfile(PREOUT, sprintf('%s_%s_%s_%sfreq_respavg.mat', SUBJ, cfg.trigger, cfg.freqtype, cfg.freqband ));

fprintf('Saving %s . . . \n', filesave)
save(filesave) % save all vars



