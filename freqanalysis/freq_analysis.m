function freq_analysis(cfg)
% Do freqanalysis for stim and resp, totalpow, evoked and induced, low,
% high and full freq

PREIN = cfg.PREIN;
PREOUT = cfg.PREOUT;
trigger = cfg.trigger;
freqtype = cfg.freqtype;
freqband = cfg.freqband;

mkdir(PREOUT)
cd(PREIN)
% inputfile = dir('*CSD.mat');
inputfile = dir('*costrap_CSD.mat');

fprintf('Loading %s from...\n %s\n', inputfile.name, PREIN)
load(inputfile.name);

% Robust detrend
for itrial=1:length(data.trial)
    data.trial{itrial} = transpose(nt_detrend(data.trial{itrial}', 2));
end

cfg = [];
cfg.output = 'pow';
cfg.channel = 'EEG';
cfg.keeptapers  = 'no';
cfg.pad = 7;
cfg.method = 'mtmconvol';
cfg.trials = 'all';
switch freqband
    case 'high'
        cfg.taper = 'dpss'; % high frequency-optimized analysis (smooth)
        cfg.keeptrials  = 'yes';
        cfg.foi = 36:2:127;
        cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.4;
        cfg.tapsmofrq = ones(length(cfg.foi),1) .* 8;
    case 'low'
        cfg.taper = 'hanning'; % low frequency-optimized analysis
        cfg.keeptrials  = 'yes'; % needed for fourier-output
        %             cfg.keeptapers = 'yes'; % idem
        cfg.foi = 3:35;
        cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.4;
        cfg.tapsmofrq = ones(length(cfg.foi),1) .* 4.5;
    case 'full'
        cfg.taper = 'dpss'; % high frequency-optimized analysis (smooth)
        cfg.keeptrials  = 'yes';
        cfg.foi = 5:2:127;
        cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.4;   % length of time window = 0.4 sec
        cfg.tapsmofrq = ones(length(cfg.foi),1) .* 5;
    otherwise
        error('Unexpected analysisname. abort.');
end %switch

data_stim = data; % keep the stim if needed in induced
switch trigger
    case 'resp'
        cfg_resp=[];
        cfg_resp.offset = -data.trialinfo(:,4);
        cfg_resp.trials = find(cfg_resp.offset < 1);
        data = ft_redefinetrial(cfg_resp, data);
        %         cfg.toi = -1.25:0.05:1;
        cfg.toi = -1.25:0.05:1.25;
    case 'stim'
        %         cfg.toi = -0.5:0.05:1;
        cfg.toi = -0.8:0.05:2;
end


switch freqtype
    case 'evoked'
        ctr = 0;                timelock = {};                trialinfo_ssvep = [];
        trialinfo_ses = data.trialinfo;
        conds = 1:3;                stims = 1:3;                resps = 1:3;
        for icond = conds % NoMiss NoFA
            for istim = stims % 1 = left, 2 = right
                for iresp = resps % left or right
                    cond_ind = trialinfo_ses(:,1) == icond;
                    if icond==3, cond_ind(:)=true; end
                    stim_ind = trialinfo_ses(:,2) == istim;
                    if istim==3, stim_ind(:)=true; end
                    resp_ind = trialinfo_ses(:,3) == iresp;
                    if iresp==3, resp_ind(:)=true; end
                    
                    cfg_time = [];
                    cfg_time.vartrllength = 2;
                    cfg_time.trials = find(cond_ind & stim_ind & resp_ind);
                    if cfg_time.trials > 0
                        ctr=ctr+1;
                        timelock{ctr} = ft_timelockanalysis(cfg_time, data);
                        trialinfo_ssvep = [trialinfo_ssvep; icond istim iresp];
                    else
                        warning('No trials for this condition present')
                    end
                end
            end
        end
        %         cfg2=[];  % timelock reduces time axis
        %         cfg2.appenddim = 'rpt';
        %         ft_appendtimelock(cfg2, timelock{:})
        
        data = ft_appenddata([], timelock{:});
        data.trialinfo = trialinfo_ssvep;
        data.cfg.previous = data.cfg.previous{1}; % otherwise 27 cfg's are saved, huge output files
        
        [~,file] = fileparts(inputfile.name);
        outputfile = [file '_' freqtype '.mat'];
        fprintf('Saving %s to...\n %s\n', outputfile, PREIN)
        save(fullfile(PREIN, outputfile), 'data');
        
    case 'induced' % load in evoked (ssvep) data, subtract from each trial, run freqanalysis TODO fix for resp-locked
        data = remove_ERP_fromdata(data_stim, 'subtract');

    case 'induced_robdetr' % load in evoked (ssvep) data, subtract from each trial, run freqanalysis TODO fix for resp-locked
        data = remove_ERP_fromdata(data_stim, 'nt_detrend');

%         inputfile = dir('*CSD_evoked.mat');
%         fprintf('Loading %s from...\n %s\n', inputfile.name, PREIN)
%         load(inputfile.name);
%         data_ssvep = data;
%         
%         %         ctr = 0; timelocks2append = {};
%         %         conds = 1:2;                stims = 1:2;                resps = 1:2;
%         %         for icond = conds % NoMiss NoFA
%         %             for istim = stims % 1 = left, 2 = right
%         %                 for iresp = resps % left or right
%         %                     cond_ind = data_ssvep.trialinfo(:,1) == icond; % trialinfo nicely coded for ssvep
%         %                     stim_ind = data_ssvep.trialinfo(:,2) == istim;
%         %                     resp_ind = data_ssvep.trialinfo(:,3) == iresp;
%         
%         cond_ind = data_ssvep.trialinfo(:,1) == 3; % take avg across all trials
%         stim_ind = data_ssvep.trialinfo(:,2) == 3;
%         resp_ind = data_ssvep.trialinfo(:,3) == 3;
%         
%         cfg_time = [];
%         cfg_time.vartrllength = 2;
%         cfg_time.trials = find(cond_ind & stim_ind & resp_ind); % 1 trial
%         %                     if isempty(cfg_time.trials)
%         %                         warning('No trials for this condition present, continuing')
%         %                         continue
%         %                     end
%         timelock_ssvep = ft_timelockanalysis(cfg_time, data_ssvep); % just to get the single trial in a nice matrix
%         
%         %                     % now get the single trials from the overall data
%         %                     cond_ind = data_stim.trialinfo(:,1) == icond;
%         %                     stim_ind = data_stim.trialinfo(:,2) == istim;
%         %                     resp_ind = data_stim.trialinfo(:,3) == iresp;
%         %                     cfg_time.trials = find(cond_ind & stim_ind & resp_ind);
%         
%         cfg_time.keeptrials = 'yes';
%         timelock = ft_timelockanalysis(cfg_time, data_stim); % take all trials
%         
%         % repmat timelock_ssvep and subtract from trials
%         ntrials=size(timelock.trial,1);
%         ssvepdat_repmat = shiftdim(timelock_ssvep.avg,-1); % make leading singleton
%         ssvepdat = repmat(ssvepdat_repmat, [ntrials 1 1]);
%         try % TODO fix time axis for SSVEP subtraction
%             timelock.trial = timelock.trial - ssvepdat; % subtract the ssvep!
%         catch
%             timelock;
%         end
%         %                     ctr=ctr+1;
%         %                     timelocks2append{ctr} = timelock;
%         %                 end
%         %             end
%         %         end
%         %         data = ft_appenddata([], timelocks2append{:});
%         %         data.cfg.previous = data.cfg.previous{1};
%         data = timelock;
%         
%         %                     save([inputfile '_removed'], 'data')  % save ssvep-removed data
%         
%         if strcmp(trigger, 'resp')
%             cfg_resp=[];
%             cfg_resp.offset = -data.trialinfo(:,4);
%             cfg_resp.trials = find(cfg_resp.offset < 1);
%             data = ft_redefinetrial(cfg_resp, data); %overwrite the (totalpow) resplocked data already in memory
%         end
        
    case 'totalpow' % Nothing to do, just enter all trials in ft_freqanalysis
        
    otherwise
        disp('Unknown freqtype, aborting . . .')
        return
end

freq = ft_freqanalysis(cfg, data); % run it

[~, datafilename]= fileparts(inputfile.name);
% outputfile = sprintf('%s_%s_%sfreq.mat', datafilename, freqtype, freqband);
outputfile = sprintf('%s_%s_%sfreq.mat', datafilename, freqtype, freqband);
fprintf('Saving %s to...\n %s\n', outputfile, PREOUT)
if exist(fullfile(PREOUT, outputfile))
    fprintf('Deleting old file . . .\n')
    delete(fullfile(PREOUT, outputfile)); % cluster does not overwrite? NO
end
save(fullfile(PREOUT, outputfile), 'freq');

%                         temp={}; cfg_sel = []; % workaround to get back to datatype data, for redefinetrial
%                         for itrial = 1:size(data.trialinfo,1)
%                             fprintf('Trial %d . . .  ', itrial)
%                             cfg_sel.trials = itrial;
%                             temp{itrial} = ft_selectdata(cfg_sel, data);
%                         end
%                         data = ft_appenddata([], temp{:});
%                         data.cfg.previous = data.cfg.previous{1};
%                         clear temp

%                     cfg=[];
%                     cfg.parameter = 'powspctrm';
%                     cfg.outputfile = sprintf('%s%s_%s_%s_%sfreq', PREOUT, SUBJ{isub}, session_name, freqtype, freqband);
%                     freq = ft_appendfreq(cfg, freq{:});
%                 FIXME:
%                 freq.cfg.previous = freq.cfg.previous{1}; % otherwise 8
%                 cfg's are saved, huge output files
%
%                 cfg_save=[];
%                 cfg_save.outputfile = sprintf('%s%s_%s_%s_%sfreq', PREOUT, SUBJ{isub}, session_name, freqtype, freqband);
%                 freq2 = ft_freqdescriptives(cfg_save, freq); % only to save fhe freq

%                         cfg = [];
%     cfg.component = [1:48];       % specify the component(s) that should be plotted
%     %             cfg.component = [13:24];       % specify the component(s) that should be plotted
%     cfg.layout = 'biosemi64.lay'; % specify the layout file that should be used for plotting
%     cfg.comment   = 'no';
%     figure
%     ft_topoplotIC(cfg, comp)
%
%         cfg = [];
%     cfg.layout = 'biosemi64.lay'; % specify the layout file that should be used for plotting
%     cfg.viewmode = 'component';
%     ft_databrowser(cfg, comp)

