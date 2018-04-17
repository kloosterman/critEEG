function [ behav ] = critEEG_analyze_trialinfo( SUBJ )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

addpath('/path/')
% close all
cfg=[]; behav=[];

load('all_events_unbalanced.mat')
cfg.eventall = event; % JJ's file
basepath = '/path/';

behavout = fullfile(basepath, 'behavior');
behav.PREOUT = '/plots/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RT_edges = 0:0.04:0.84; % count from target onset

behav.nhits = nan(length(SUBJ), 3, 2);
behav.hitrate = nan(length(SUBJ), 3, 2);
behav.farate = nan(length(SUBJ), 3, 2);
behav.RT =  nan(length(SUBJ), 3, 2, 2);
behav.RTsum =  nan(length(SUBJ), 3, 2, 2);
behav.RThist =  nan(length(SUBJ), 3, 2, 2, length(RT_edges));
for isub = 1:length(SUBJ)
    
    load(fullfile(behavout, [SUBJ{isub} '_trialinfo.mat'])) % trialinfo comes out
    
    nonzerorts = trialinfo(:,4) > 0;
    trialinfo(nonzerorts,4) = trialinfo(nonzerorts,4) - round(0.16*256); % MAKE RT's FROM TARGET ONSET!
    
    ses_nos = unique(trialinfo(:,7))';
    ses_nos(end+1) = 4; % fixed effects across sessions in 4
    for ises = ses_nos
        ses_ind = trialinfo(:,7) == ises;
        if ises==4, ses_ind(:) = true; end
        
        if isempty(ses_ind) && ises~=4 ; continue; end % % session not present
        for icond=1:3
            cond_ind = trialinfo(:,1) == icond;
            if icond==3, cond_ind(:) = true; end
            
            nhits = length(find(...
                ses_ind ...  % trialinfo(:,7) == ises ...
                & cond_ind  ... % trialinfo(:,1) == icond ...
                & trialinfo(:,2) == 1 ... % Fig
                & trialinfo(:,3) == 1)); % report
            ncr = length(find(...
                ses_ind ...  % trialinfo(:,7) == ises ...
                & cond_ind  ... % trialinfo(:,1) == icond ...
                & trialinfo(:,2) == 2 ... % Hom
                & trialinfo(:,3) == 2)); % No report
           
            behav.nhits(isub, ises, icond) = nhits;
            behav.hitrate(isub, ises, icond) = nhits ...
                / length(find(ses_ind & cond_ind & trialinfo(:,2) == 1));
            behav.propcorr(isub, ises, icond) = (nhits + ncr) ...
                / length(find(ses_ind & cond_ind));
            
            trialinds = ...
                ses_ind ... % trialinfo(:,7) == ises ...
                & cond_ind  ... %trialinfo(:,1) == icond ...
                & trialinfo(:,2) == 2 ... % hom
                & trialinfo(:,3) == 1; % Report
            trialsoi = trialinfo(trialinds,:);
            
            nlooseFAs = sum(unique(trialsoi(:,9))); % sum up over the 3 runs per cond
            behav.nlooseFAs(isub, ises, icond) = nlooseFAs;
                        
            % Or only include FA's within train
            behav.farate(isub, ises, icond) = (size(trialsoi,1) ) ...
                / length(find( ses_ind & cond_ind & trialinfo(:,2) == 2));
            
            for istim=1:3
                stim_ind = trialinfo(:,2) == istim;
                resp_ind = trialinfo(:,3) == 1;
                %                 if icond==3, cond_ind(:)=true; end
                if istim==3, stim_ind(:)=true; end
                %                 if iresp==3, resp_ind(:)=true; end
                
                trialinds = find(ses_ind & cond_ind & stim_ind & resp_ind);
                
                behav.RT(isub, ises, istim, icond) = mean(trialinfo(trialinds,4)/256);
                behav.RTsum(isub, ises, istim, icond) = sum(trialinfo(trialinds,4)/256);
                behav.RTstd(isub, ises, istim, icond) = std(trialinfo(trialinds,4)/256);
                
                % TODO rt variability
                behav.RThist(isub, ises, istim, icond,:) = histc(trialinfo(trialinds,4)/256, RT_edges) / length(trialinds);
                
                behav.meantime_to_next_resp(isub, ises, istim, icond) = nanmean(trialinfo(trialinds,10)/256);
                behav.mediantime_to_next_resp(isub, ises, istim, icond) = nanmedian(trialinfo(trialinds,10)/256);
            end
        end
        
    end
end

behav.subj = SUBJ;

behav.hitrate(behav.hitrate == 0) = 0.01; % 
behav.farate(behav.farate == 0) = 0.01; % 
behav.hitrate(:,:,4) = behav.hitrate(:,:,2) - behav.hitrate(:,:,1);
behav.farate(:,:,4) = behav.farate(:,:,2) - behav.farate(:,:,1);

behav.hitminusfarate = behav.hitrate - behav.farate;

% %
behav.dprime = norminv(behav.hitrate) - norminv(behav.farate);
% behav.dprime(:,4,:) = nanmean(behav.dprime,2);

behav.dprime(:,:,4) = behav.dprime(:,:,2) - behav.dprime(:,:,1); % NOFA - Nomiss    isub, ises, icond

behav.criterion = -0.5 * (norminv(behav.hitrate) + norminv(behav.farate));
% behav.criterion(:,4,:) = nanmean(behav.criterion,2); % sessions
behav.criterion(:,:,4) = behav.criterion(:,:,2) - behav.criterion(:,:,1); % Nomiss -  NOFA

behav.JJcriterion = normpdf(norminv(behav.hitrate))./normpdf(norminv(behav.farate)); %% Johannes implementation
% behav.JJcriterion(:,4,:) = nanmean(behav.JJcriterion,2); % sessions
behav.JJcriterion(:,:,4) = behav.JJcriterion(:,:,2) - behav.JJcriterion(:,:,1); % NOFA - Nomiss

% behav.RT(:,4,:,:) = nanmean(behav.RT,2);
behav.RT(:,:,3,:) = nanmean(behav.RT,3); % average over stim
behav.RT(:,:,:,4) = behav.RT(:,:,:,2) - behav.RT(:,:,:,1) ; % Nomiss -  NOFA  dims (isub, ises, istim, icond)

% behav.meantime_to_next_resp(:,4,:,:) = nanmean(behav.meantime_to_next_resp,2);
behav.meantime_to_next_resp(:,:,3,:) = nanmean(behav.meantime_to_next_resp,3); % average over stim
behav.meantime_to_next_resp(:,:,:,4) = behav.meantime_to_next_resp(:,:,:,2) - behav.meantime_to_next_resp(:,:,:,1) ; %  Nomiss-NOFA  dims (isub, ises, istim, icond)

% behav.mediantime_to_next_resp(:,4,:,:) = nanmean(behav.mediantime_to_next_resp,2);
behav.mediantime_to_next_resp(:,:,3,:) = nanmean(behav.mediantime_to_next_resp,3); % average over stim
behav.mediantime_to_next_resp(:,:,:,4) = behav.mediantime_to_next_resp(:,:,:,2) - behav.mediantime_to_next_resp(:,:,:,1) ; % Nomiss-NOFA  dims (isub, ises, istim, icond)


behav.RThist(:,4,:,:,:) = nanmean(behav.RThist,2);

behav.nlooseFAs(:, 4, :)  = nanmean(behav.nlooseFAs,2); % behav.nlooseFAs(isub, ises, icond)
behav.nlooseFAs(:, :, 4)  = behav.nlooseFAs(:,:,2) - behav.nlooseFAs(:,:,1) ; % behav.nlooseFAs(isub, ises, icond)


% nrand=10000; % 10000 is default
stat = [];
for icond = 1:2
    stat.dprime(icond) = permtest(behav.dprime(:,4,icond) ); % (isub, ises, icond)
    stat.criterion(icond) = permtest(behav.criterion(:,4,icond) ); % (isub, ises, icond)
    stat.JJcriterion(icond) = permtest(behav.JJcriterion(:,4,icond) ); % (isub, ises, icond)
end
icond = 3;
stat.dprime(icond) = permtest(behav.dprime(:,4,1) , behav.dprime(:,4,2) );
stat.criterion(icond) = permtest(behav.criterion(:,4,1) , behav.criterion(:,4,2) );
stat.JJcriterion(icond) = permtest(behav.JJcriterion(:,4,1) , behav.JJcriterion(:,4,2) );


RT_prob = [];
for istim=1:2
    stat.RT(istim) = permtest(squeeze(behav.RT(:,4,istim, 1)) , squeeze(behav.RT(:,4,istim, 2)));
end
stat.RT(3) = permtest(squeeze(behav.RT(:,4,1, 3)) , squeeze(behav.RT(:,4,2, 3))); % H vs FA both conds (isub, ises, istim, icond)

% p_c = randtest1d(behav.criterion(:,1) - behav.criterion(:,2), zeros(size(behav.criterion(:,1))), 0, 1000 )
% p = randtest1d(behav.criterion(:,1), zeros(size(behav.criterion(:,1))), 0, nrand )

behav.sesleg = {'Ses1', 'Ses2', 'Ses3', 'SesAvg'};
% condleg = {'NoFA', 'NoMiss'};
behav.condleg = {'Conservative', 'Liberal', 'critcomb', 'Liberal-Conservative'};

behav.behav_measures = {'dprime', 'criterion', 'JJcriterion'};
behav.stat = stat;
behav.RT_edges = RT_edges;

% RCS rate correct score: ncorrect / sum (allRTs)
% RTsum(isub, ises, istim, icond)    nhits(isub, ises, icond) 
behav.rcs = squeeze(behav.nhits(:,4,:)) ./ squeeze(behav.RTsum(:,4,3,:));
p = permtest(behav.rcs(:,1), behav.rcs(:,2) );

% RCS rate correct score: propcorr / meanRT
% RTsum(isub, ises, istim, icond)    nhits(isub, ises, icond);
% This score can be interpreted as the number of correct responses per second of activity.
behav.rcs = squeeze(behav.propcorr(:,4,:)) ./ squeeze(behav.RT(:,4,3,1:3))
p = permtest(behav.rcs(:,1), behav.rcs(:,2) )

% p.655: a new integrated measure is proposed here that is based on a linear combination of RT and PE. 
% This linear integrated speed-accuracy score (LISAS) is defined as: -- Highlighted Nov 28, 2017
% LISAS = RTj + Srt/Spe * PEj
% As an example, consider an average correct RT of 500 with a standard deviation of 100, 
% and an error rate of .10 with a standard deviation of .05, LISAS will be 500 + 100/.05 × .1 = 700.
% this measure yields an estimate of RT corrected for the number of errors.  -- Highlighted Nov 28, 2017
RTj = squeeze(behav.RT(:,4, 3, 1:3));
Srt = squeeze(behav.RTstd(:,4, 3, :));
PEj = squeeze(behav.propcorr(:,4,:));
Spe = squeeze(std(behav.propcorr(:,1:3,:), 0, 2));  % take std across sessions
LISAS = RTj + Srt./Spe .* PEj;

save(fullfile(behavout, 'behavstruct'), 'behav' )
