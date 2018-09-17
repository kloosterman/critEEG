%% PLot prestim alpha TFR + topo, then sorted alpha power, then gain in gamma + bars of quadratic coeff

% addpath(genpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/plotting'));
% addpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/cbrewer')

% respavg = critEEG_load_respavg()

%% plot binned by alpha: alpha itself and gamma, 10 bins
% close all

meas_leg = {'Alpha pre' 'Abeta' 'SSVEP' 'P1' ... 1-4
    'N1' 'CPP' 'SSVEP high' 'gamma'... 5-8
    'allhighfreq' 'sdt_crit' 'ddm_dc' 'sdt_dprime' 'RT'}; % 9-13

saveddm = 0;
nbins = 10;
binoverlap = 'yes';
binmethod = 'trlcount'; % datarange or trlcount
niteroutliers = 1; % niteroutliers = [];  1
outlierzthr = 3;
zscorevars = 0; % zscore single trials wrt subject mean and std across conds

gammafocus = 'transient';  % transient  overall
switch gammafocus
    case 'transient'
        logvars = 0;  % logvars = 0;
        normalizetype = 'nozscoring';  % zacrosscond nozscoring
        bl_correct = 'avgprestim_psc'; % avgprestim_psc  singletrialprestim avgprestim_log   avgprestim_psc_singletrial
    case 'overall'
        logvars = 1;
        normalizetype = 'zacrosscond';  % zacrosscond nozscoring
        bl_correct = 'none'; % avgprestim_psc  singletrialprestim avgprestim_log   avgprestim_psc_singletrial
end

% normalizetype = 'psc';  %'zwithincond_sharedsd' zwithincond_sharedmean
% normalizetype = 'psc_onlog';
% normalizetype = 'zwithincond_sharedsd';  %'zwithincond_sharedsd'
% normalizetype = 'zwithincond';  %'zwithincond_sharedsd'
% normalizetype = 'zwithincond_sharedmean';  %'zwithincond_sharedsd' ?
% normalizetype = 'zwithincondnodemeaning';
% normalizetype = 'zsdwithincond_sharedmean';

xaxisunit = 'alphavals';%     alphavals           binno logalphavals

if nbins == 10 % used for ddm binning!!!:
    %     bin_ranges = [0    20;   5    25;   15    35;    25    45;    35    55;    45    65;   55 75; 65 85; 75 95;  80 100];
    if strcmp(binoverlap, 'yes')
        bin_ranges = [0    15;   5    25;   15    35;    25    45;    35    55;    45    65;   55 75; 65 85; 75 95;  85 100];
    else
        bin_ranges = [0    10;   10    20;   20    30;    30    40;    40    50;    50    60;   60 70; 70 80; 80 90; 90 100];
    end
    %     bin_ranges = [0    9;   10    19;   20    29;    30    39;    40    49;    50    59;   60 69; 70 79; 80 89; 90 100];
elseif nbins == 9
    bin_ranges = [0 20;  10 30; 20 40; 30 50; 40 60; 50 70; 60 80; 70 90; 80 100 ];
elseif nbins == 5
    %     bin_ranges = [0 33;  16 50; 33 67; 50 83; 66 100 ];
    %     bin_ranges = [0 15;  15 40; 40 60; 60 85; 85 100 ];
    bin_ranges = [0 20;  20 40; 40 60; 60 80; 80 100 ];
    %     bin_ranges = [0 19;  20 39; 40 59; 60 79; 80 100 ];
elseif nbins == 3
    bin_ranges = [0 33; 33 67; 67 100];
elseif nbins == 20
    bin_ranges = [0:5:95; 5:5:100]';
%         bin_ranges = [bin_ranges(:,1)-2.5 bin_ranges(:,2)+2.5];
%         bin_ranges(1) = 0;    bin_ranges(end) = 100;
elseif nbins == 50
    bin_ranges = [0:2:98; 2:2:100]';
%     bin_ranges = [bin_ranges(:,1)-2 bin_ranges(:,2)+2];
%     bin_ranges(1) = 0;bin_ranges(end) = 100;
elseif nbins == 100
    bin_ranges = [0:1:99; 1:100]';
elseif nbins == 1
    bin_ranges = [0 100];
end

dat_bins = nan(length(respavg.SUBJ), 3, 13, nbins, 2); %dimord subj ses measureno binno cond
dat_bin_count  = nan(length(respavg.SUBJ), 3, 13, nbins, 2); %for histogram of alpha vals
dat_bin_edges  = nan(length(respavg.SUBJ), 3, 13, nbins, 2, 2); %for histogram of alpha vals

ises = 4; istim = 3; iresp = 3; isoi = 1; %     'occipital'    'motor'    'allsens'    'occpar'    'frontal'    'occlatr'    'posterior'    'Pz' 'POz'

alpha_soi = 1; alpha_var = []; avgalphapow = []; gammapow = []; trlcount =[]; outmat = []; histdat = []; histedgedat = [];
for isub = 1:15 % 1 :nsub %   % 1:nsub  %
    fprintf('SUBJ %s %d\n', respavg.SUBJ{isub}, ises)
    if strcmp(normalizetype, 'psc')  % load high freq basespec,
        basespecfolder = fullfile(respavg.inputfolder, 'stim', respavg.SUBJ{isub});
        file = dir(sprintf('%s/*highfreq*', basespecfolder));
        load(fullfile(basespecfolder, file.name), 'basespec_ses'); % % ses chan freq cond stim resp
        foi = respavg.freq{2} >= 59 & respavg.freq{2} <= 100; % band 2
        basedat = basespec_ses(:,respavg.sens.ind{isoi}, foi,:,:,:);
        basedat = squeeze(mean(basedat(:,:),2));
    end
    % get condition avg and sd of alpha and im for zscoring
    conddat_alpha = [];
    for ises = 1:3 % concat sessions
        if isub == 4 && ises == 3; disp('Ses not found'); continue; end
        if isub == 12 && ises == 3 && icond == 2; disp('Ses not found'); continue; end
%         temp = [log(respavg.pow_singletrial{isub,1}{1,ises,1,istim,iresp}(:,isoi)); ...
%             log(respavg.pow_singletrial{isub,1}{1,ises,2,istim,iresp}(:,isoi))];
        temp = [respavg.pow_singletrial{isub,1}{1,ises,1,istim,iresp}(:,isoi); ...
            respavg.pow_singletrial{isub,1}{1,ises,2,istim,iresp}(:,isoi)];
        conddat_alpha = [conddat_alpha; temp];
    end
    if logvars
        conddat_alpha = log(conddat_alpha);
    end
    edges = 0:max(conddat_alpha)/50:max(conddat_alpha);
    
    for im = [1 3 8 10] %[1 8 10] %1:length(meas_leg)  % 1:5 %2:5  % ssvep response %2:5
    
        if im < 7; iband = 1; else iband = 2; end
        % get condition avg and sd of im for zscoring 
        conddat_im = [];
        for ises = 1:3 % concat sessions
            if isub == 4 && ises == 3; disp('Ses not found'); continue; end
            if isub == 12 && ises == 3 && icond == 2; disp('Ses not found'); continue; end
            temp = [respavg.pow_singletrial{isub,iband}{im,ises,1,istim,iresp}(:,isoi); ...
                respavg.pow_singletrial{isub,iband}{im,ises,2,istim,iresp}(:,isoi)];
            conddat_im = [conddat_im; temp];
        end
        if logvars
            conddat_im = log(conddat_im);
        end

        for icond = 1:2
            dat=[];
            for ises = 1:3 % concat sessions
                if isub == 4 && ises == 3; disp('Ses not found'); continue; end
                if isub == 12 && ises == 3 && icond == 2; disp('Ses not found'); continue; end
                temp = [];
%                 temp(:,1) = log(respavg.pow_singletrial{isub,1}{1,ises,icond,istim,iresp}(:,1)); %prestim alpha, log works better even with trialcount? range?
                temp(:,1) = respavg.pow_singletrial{isub,1}{1,ises,icond,istim,iresp}(:,1); %prestim alpha
                if im < 10
                    if strcmp(normalizetype, 'psc') && im == 8 % take psc wrt basedat
                        temp(:,2) = respavg.pow_singletrial{isub,iband}{im,ises,icond,istim,iresp}(:,isoi);
                        temp(:,2) = (temp(:,2) - basedat(ises)) ./ basedat(ises) * 100;
                    elseif strcmp(normalizetype, 'psc_onlog') && im == 8  % take log
                        temp(:,2) = log(respavg.pow_singletrial{isub,iband}{im,ises,icond,istim,iresp}(:,isoi));
                    elseif strcmp(bl_correct, 'singletrialprestim') && im == 8  % subtract pre from post
%                         temp(:,2) = respavg.pow_singletrial{isub,iband}{8,ises,icond,istim,iresp}(:,isoi) - ...
%                             respavg.pow_singletrial{isub,iband}{11,ises,icond,istim,iresp}(:,isoi); % 11 = prestim gamma
                        temp(:,2) = log(respavg.pow_singletrial{isub,iband}{8,ises,icond,istim,iresp}(:,isoi)) - ...
                            log(respavg.pow_singletrial{isub,iband}{11,ises,icond,istim,iresp}(:,isoi)); % 11 = prestim gamma
                        
                    elseif strcmp(bl_correct, 'avgprestim_log') && im == 8
                        prestimavg = nanmean(respavg.pow_singletrial{isub,iband}{11,ises,icond,istim,iresp}(:,isoi)); % 11 = prestim gamma
                        temp(:,2) = log(respavg.pow_singletrial{isub,iband}{8,ises,icond,istim,iresp}(:,isoi)) - log(prestimavg);

                    elseif strcmp(bl_correct, 'avgprestim_psc') && im == 8
                        prestimdat = respavg.pow_singletrial{isub,iband}{11,ises,icond,istim,iresp}(:,isoi); % 11 = prestim gamma
                        poststimdat = respavg.pow_singletrial{isub,iband}{8,ises,icond,istim,iresp}(:,isoi);    

%                         % remove outliers prestimdat > or < 3
%                         for iter = niteroutliers % over sessions
%                             zdat = (log(prestimdat)  - nanmean(log(prestimdat) )) ./ nanstd(log(prestimdat )); % take log to make normal
%                             prestimdat = prestimdat(abs(zdat) < outlierzthr); % remove outliers, take abs to remove both sides
%                         end

                        % normalize to psc
                        temp(:,2) = (poststimdat - nanmean(prestimdat)) ./ nanmean(prestimdat) * 100;

                    elseif strcmp(bl_correct, 'avgprestim_psc') && im == 3
                        
                        prestimdat = respavg.pow_singletrial{isub,iband}{12,ises,icond,istim,iresp}(:,isoi); % 12 = prestim ssvep
                        poststimdat = respavg.pow_singletrial{isub,iband}{3,ises,icond,istim,iresp}(:,isoi); % 4 = ssvep                       % TODO reject outliers before taking avg for psc
                        
                        % normalize to psc
                        temp(:,2) = (poststimdat - nanmean(prestimdat)) ./ nanmean(prestimdat) * 100;

                    elseif strcmp(bl_correct, 'avgprestim_psc_singletrial') && im == 8
                        prestim = respavg.pow_singletrial{isub,iband}{11,ises,icond,istim,iresp}(:,isoi); % 11 = prestim gamma
                        temp(:,2) = (respavg.pow_singletrial{isub,iband}{8,ises,icond,istim,iresp}(:,isoi) - prestim) ./ nanmean(prestim) * 100;
                    else 
                        temp(:,2) = respavg.pow_singletrial{isub,iband}{im,ises,icond,istim,iresp}(:,isoi);
                    end
                else
                    temp = [temp respavg.pow_singletrial{isub,1}{im,ises,icond,istim,iresp}];
                end
                dat = [dat; temp]; % concat sessions                
            end
            dat = dat(~isnan(dat(:,1)),:); % remove nans

            alpha_var(isub, icond) = var(dat(:,1));
            
            if logvars
                dat = log(dat); %%%%%% take log
%                 dat(:,2) = log(dat(:,2)); %%%%%% take log
            end
            alpha_kurt(isub, icond) = kurtosis(dat(:,1));

            % zscore both vars with cond avg and sd
            if zscorevars
                dat(:,1) = (dat(:,1) - nanmean(conddat_alpha)) ./  nanstd(conddat_alpha);
                dat(:,2) = (dat(:,2) - nanmean(conddat_im)) ./  nanstd(conddat_im);
            end
            
            % remove outliers > or < 3
            for iter = niteroutliers % over sessions
                for i=1:2 % both vars
                    if zscorevars
                        dat = dat(abs(dat(:,i)) < outlierzthr,:); % remove outliers, take abs to remove both sides
                    else
                        if i == 1
                            zdat = (log(dat(:,i))  - nanmean(log(conddat_alpha) )) ./ nanstd(log(conddat_alpha )); % take log to make normal
                        else
                            %                             zdat = (log(dat(:,i))  - nanmean(log(conddat_im) )) ./ nanstd(log(conddat_im )); % take log to make normal
                            %                             zdat = (dat(:,i)  - nanmean(conddat_im) ) ./ nanstd(conddat_im ); % take log to make normal
                            zdat = (dat(:,i)  - nanmean(dat(:,i))) ./ nanstd(dat(:,i)); % take log to make normal

                        end
                        dat = dat(abs(zdat) < outlierzthr,:); % remove outliers, take abs to remove both sides
                    end
                end
            end
            if im == 8
                gammapow(isub, icond) = mean(dat(:,2));
            end

            if im == 1 % plot histogram later
                [histdat(:,icond,isub) histedgedat(:,icond,isub)] = histcounts(dat(:,2), 50, 'Normalization', 'probability');
%                 [histdat(:,icond,isub) histedgedat(:,icond,isub)] = histcounts(dat(:,2), edges, 'Normalization', 'probability');
%                 [histdat(:,icond,isub) histedgedat(:,icond,isub)] = histcounts(dat(:,2), 25);
            end

            ises = 4;
            if isempty(dat);
                fprintf('dat not found\n');
                continue;
            end
            if im == 1; avgalphapow(isub, icond) = mean(dat(:,1)); end % for corr with crit
                if icond == 2
%                     figure; scatter(dat(:,1), dat(:,2))
                end
         
            
            for ib = 1:size(bin_ranges,1)
                switch binmethod
                    case 'trlcount'
                        bin_edges = prctile(dat(:,1), bin_ranges(ib,:)); % use  percentiles itself
                    case 'datarange' % control bin_ranges, e.g. to get 30 % of the observations in each bin, cf Rajagovindan etal
                        data_range = linspace(min(dat(:,1)), max(dat(:,1)), 100);
                        bin_edges = prctile(data_range, bin_ranges(ib,:)); % use data range, not percentiles itself
                end
                bin_ind = discretize(dat(:,1), bin_edges, 1);
                
                if im == 10 % SDT
                    ydat = dat(bin_ind == 1, 2:4);
                    dat(bin_ind == 1, 5) = ib; % remember bin no
                    hitrate = length(find(ydat(:,1)==1 & ydat(:,2)==1)) / length(find(ydat(:,1)==1));
                    farate = length(find(ydat(:,1)==2 & ydat(:,2)==1)) / length(find(ydat(:,1)==1));
                    
                    % make correction dependent on N nontarget trials for FA and target trials for H
                    ntrials = size(ydat,1);
                    trlcount(isub, ib, icond) = ntrials;
                    if hitrate == 0;                        disp('H 0')
                        hitrate = 0.5*length(find(ydat(:,1) == 1)) / ntrials;
                    end
                    if hitrate == 1;                        disp('H 1')
                        hitrate = 1 - 0.5*length(find(ydat(:,1) == 1)) / ntrials;
                    end
                    if farate == 0;                        disp('FA 0')
                        farate = 0.5*length(find(ydat(:,1) == 2)) / ntrials;
                    end
                    if farate == 1;                        disp('FA 1')
                        farate = 1-0.5*length(find(ydat(:,1) == 2)) / ntrials;
                    end
                    
                    % criterion:
                    dat_bins(isub, ises, 10, ib, icond) =  -0.5 * (norminv(hitrate) + norminv(farate));
                    % dprime:
                    dat_bins(isub, ises, 12, ib, icond) =  norminv(hitrate) - norminv(farate); % dprime in 12
                    dat_bins(isub, ises, 13, ib, icond) =  median(ydat(ydat(:,3) > 0,3)) / 256;  % rt
                    if nbins == 10 % output for ddm, allow overlapping bins
                        ddmdat = zeros(size(ydat,1),6);
                        ddmdat(:,1) = isub-1;
                        ddmdat(:,2) = icond-1; % 0 is cons, 1=lib
                        ddmdat(:,3) = mod(ydat(:,1),2); % 1=fig, 2=hom -> flip that so 0=hom and 1=fig)
                        ddmdat(:,4) = mod(ydat(:,2),2); % (1=yes, 2=no resp -> flip that so 0=no resp and 1=yes):
                        ydat(ydat(:,3) < 0, 3) = NaN; % set missing (negative) RT's to NaN
                        ddmdat(:,5) = ydat(:,3) / 256; % RT in s
                        ddmdat(:,6) = ib; % 10 alpha bins
%                         ctr=0; % 5 bins: take 10 bins and put in 5 bins
%                         for ib = 1:2:10
%                             ctr = ctr+1;
%                             inds = ddmdat(:,6) == ib | ddmdat(:,6) == ib+1;
%                             ddmdat(inds,7) = ctr;
%                         end
                        ddmdat(:,6) = ddmdat(:,6)-1; % 10 alpha bins
%                         ddmdat(:,7) = ddmdat(:,7)-1; 
                        if isub < 8
                            outmat = [outmat; ddmdat];
                        elseif isub > 8
                            ddmdat(:,1) = ddmdat(:,1) - 1;
                            outmat = [outmat; ddmdat];
                        end
                    end
                else
                    ydat = dat(bin_ind == 1, 2);
                    dat_bins(isub, ises, im, ib, icond) = mean(ydat);
%                     dat_bins(isub, ises, im, ib, icond) = median(ydat);
                    dat_bin_count(isub, ises, im, ib, icond) = length(ydat);
                    dat_bin_edges(isub, ises, im, ib, icond,:) = bin_edges;
                end
            end
        end
    end
end
if im == 10 && nbins == 10 && saveddm    % write header and matrix to file
    exp_path = fullfile('/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/DDMexport');
    mkdir(exp_path)
    trl_outfile = sprintf('critEEG_data_binby%s_overlap%s.csv', binmethod, binoverlap);
    fid = fopen(fullfile(exp_path, trl_outfile), 'w') ;
%     fprintf(fid, 'subj_idx,cond,stimulus,response,rt,alpha10bins,alpha5bins\n');
    fprintf(fid, 'subj_idx,cond,stimulus,response,rt,alpha10bins\n');
    dlmwrite(fullfile(exp_path, trl_outfile), outmat, '-append', 'precision', 10)
    fclose(fid);
end

% % add ddm dc 10 bins
nsub = length(respavg.SUBJ);
if nsub == 15
%     if nbins == 10
% %         dat_bins(:, 4, 11, 1:10, 1:2) = respavg.ddmpars10bins.dc([1:13, 15:end],:,: ); % remove 8 and 14
% %         dat_bins(:, 4, 11, 1:10, 1:2) = respavg.ddmpars10bins.dc; % remove 8 and 14
%     elseif nbins == 5
%         % add ddm dc 5
%         dat_bins(:, 4, 11, 1:5, 1:2) = respavg.ddmpars5bins.dc([1:13, 15:end],:,: ); % remove 14
%     elseif nbins == 3
%         dat_bins(:, 4, 11, 1:3, 1:2) = respavg.ddmpars3bins.dc([1:13, 15:end],:,: ); % remove 14
%     end
end

if  nsub == 15
%     bin_fits = bin_fits([1:7, 9:15],:,:,:,:);
    alpha_var = alpha_var([1:7, 9:15],:);
    % remove subj 8, very low gamma
    dat_bins = dat_bins([1:7, 9:15],:,:,:,:);
end

if nsub == 15
    if nbins == 10
        dat_bins(:, 4, 11, 1:10, 1:2) = respavg.ddmpars10bins.dc;
%         % zscoren:
%         dat_bins(:, 4, 11, 1:10, 1:2) = (respavg.ddmpars10bins.dc - mean(respavg.ddmpars10bins.dc, 2)) ...
%             ./ std(respavg.ddmpars10bins.dc,[],2);

        %         dat_bins(:, 4, 11, 1:10, 1:2) = (respavg.ddmpars10bins.dc - mean(respavg.ddmpars10bins.dc, 2));
        
%         meandat = mean(dat_bins, 4); % mean across conditions
% %         meandat(:,:,:,:,3) = mean(meandat,5);
%         sddat = std( dat_bins, [], 4); % sd across conditions
%         
% %         meandat = mean(reshape( dat_bins, 14, 4, 13, [] ), 4); % mean across conditions
% %         sddat = std(reshape( dat_bins, 14, 4, 13, [] ), [], 4); % sd across conditions
%         %         dat_bins = (dat_bins - meandat) ./ sddat;
%         %         dat_bins(:,:,1,:,:) = (dat_bins(:,:,1,:,:) - meandat(:,:,1)) ./ sddat(:,:,1);
% %         dat_bins(:,:,11,:,:) = (dat_bins(:,:,11,:,:) - meandat(:,:,11)) ./ sddat(:,:,11);
% %         dat_bins(:,:,11,:,:) = (dat_bins(:,:,11,:,:) - meandat(:,:,11)) ./ meandat(:,:,11) *100;
% %         dat_bins(:,:,11,:,:) = (dat_bins(:,:,11,:,:) - meandat(:,:,11,:,:)) ./ meandat(:,:,11,:,1:2) *100;
        %         dat_bins = dat_bins - mean(dat_bins(:,:,:,:,1), [], 4) ;

    elseif nbins == 5
%         dat_bins(:, 4, 11, 1:5, 1:2) = respavg.ddmpars5bins.dc;
    end
end

switch normalizetype
    case 'zacrosscond'
        sddat = std(reshape( dat_bins, 14, 4, 13, [] ), [], 4); % sd across conditions
        meandat = mean(reshape( dat_bins, 14, 4, 13, [] ), 4); % mean across conditions
%         dat_bins = (dat_bins - meandat) ./ sddat;
%         dat_bins(:,:,1,:,:) = (dat_bins(:,:,1,:,:) - meandat(:,:,1)) ./ sddat(:,:,1);
        dat_bins(:,:,8,:,:) = (dat_bins(:,:,8,:,:) - meandat(:,:,8)) ./ sddat(:,:,8);
%         dat_bins = dat_bins - mean(dat_bins(:,:,:,:,1), [], 4) ;
    case 'zwithincond'
        sddat = std(dat_bins, [], 4);
        meandat = mean(dat_bins, 4);
        dat_bins = (dat_bins - meandat) ./ sddat;
    case 'zwithincond_sharedsd'
        sddat = std(reshape( dat_bins, 14, 4, 13, [] ), [], 4); % sd across conditions
        meandat = mean(dat_bins, 4); % mean within condition
        dat_bins = (dat_bins - meandat) ./ sddat;    
    case 'zsdwithincond_sharedmean'
        sddat = std(dat_bins, [], 4);
        meandat = mean(reshape( dat_bins, 14, 4, 13, [] ), 4); % mean across conditions
        dat_bins = (dat_bins - meandat) ./ sddat;
    case 'zwithincondnodemeaning'
        sddat = std(dat_bins, [], 4);
        dat_bins = dat_bins ./ sddat;
    case 'psc_onlog'
%         meandat = meandat + 100; % make logs positive
        meandat = mean(reshape( dat_bins, 14, 4, 13, [] ), 4); % mean across conditions
        dat_bins = (dat_bins - meandat) ./ meandat * 100;
end

dat_bins(:,:,:,:,4) = dat_bins(:,:,:,:,2) - dat_bins(:,:,:,:,1); % lib- cons in 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting bins
close all

SAV = 1;
nbins = size(bin_ranges,1);
linecol = cbrewer('qual', 'Set1',3);
linecol(4,:) = [0 0 0];
set(0, 'Defaultaxesfontsize', 8)
f = figure; hold on; iplot =0;
% f.Position = [ 1989         642         600         300];  %  567         314        1354         628
f.Position = [ 567         314        500         250];
keepdat = []; % for ratio lib vs cons
% for im = [1 3]% [1 11] % [1 8] % [1 3]% % % [1 8 10 12 13]  [1 8 11]% [1,8, 10, 12, 13] % [1,8, 11]  % [1 8] % %[1,8,9] %  % % % %  % [1:4, 6:9] 1:8% 1:length(meas_leg)  % 1:5 %2:5  % ssvep response %2:5
% for im = [1 11]% [1 11] % [1 8] % [1 3]% % % [1 8 10 12 13]  [1 8 11]% [1,8, 10, 12, 13] % [1,8, 11]  % [1 8] % %[1,8,9] %  % % % %  % [1:4, 6:9] 1:8% 1:length(meas_leg)  % 1:5 %2:5  % ssvep response %2:5
for im = [1 8]% [1 11] % [1 8] % [1 3]% % % [1 8 10 12 13]  [1 8 11]% [1,8, 10, 12, 13] % [1,8, 11]  % [1 8] % %[1,8,9] %  % % % %  % [1:4, 6:9] 1:8% 1:length(meas_leg)  % 1:5 %2:5  % ssvep response %2:5
    iplot=iplot+1;
    subplot(1,2,iplot); hold on; axis square; box on
    h = []; dat = []; clear s;
    for icond =[1 2]% [1,2,4] %1:2%:2%:2%:2
        
        %         if iplot == 1
        %         else
        % plot binned gamma for alpha levels
        %             % within subj error bars
        %             if im == 11 % drift bias
        %                 temp = squeeze(dat_bins(:,4,im,:,1:2));
        %                 dat = temp(:,:,icond);
        %             else
        %                 temp = squeeze(dat_bins(:,4,im,:,1:2)); % -  squeeze(mean(dat_bins(:,4,im,:,icond),4))
        temp = squeeze(dat_bins(:,4,im,:,1:2));
        condavg = squeeze(nanmean(temp,3));
        grandavg = nanmean(condavg);
        if icond< 4
            dat = temp(:,:,icond) - condavg + grandavg;
        else
            dat = squeeze(dat_bins(:,4,im,:,icond)); % - condavg + grandavg;
        end
        %             end
        %             dat(:,icond) = squeeze(mean(dat_bins(:,4,im,:,icond))); %
            %             actual data error bars
            %             sem = squeeze(std(dat_bins(:,4,im,1:nbins,icond))) / sqrt(size(dat_bins(:,4,im,:,icond),1));
            %             dat(:,icond) = flipud(dat(:,icond));
            %             s(icond) = shadedErrorBar(0:nbins-1, dat(1:nbins,icond), sem, {'Color', linecol(icond,:)}, 1);
            
            sem = nanstd(dat,0,1) / sqrt(size(dat,1));
            % plot bin no
            %             s(icond) = shadedErrorBar(0:nbins-1, nanmean(dat(:,1:nbins)), sem, {'Color', linecol(icond,:)}, 1);
            %             scatter(0:nbins-1, mean(dat(:,1:nbins)), 'filled')
            %     xlim([0 nbins-1])
            
            switch xaxisunit
                case 'alphavals'% plot actual alpha on x-axis
%                     xdat = squeeze(nanmean(dat_bins(:,4,1,:,icond)))';
%                     XLIM = squeeze(nanmean(dat_bins(:,4,1,:,1:2)))';
                    xdat = squeeze(nanmedian(dat_bins(:,4,1,:,icond)))'; % use median
                    XLIM = squeeze(nanmedian(dat_bins(:,4,1,:,1:2)))';
                    XLIM = [min(XLIM(:)) max(XLIM(:))];
                    XLAB = '8 - 12 Hz power (V/m2)';
                case 'logalphavals'% plot actual alpha on x-axis
%                     xdat = log(squeeze(nanmean(dat_bins(:,4,1,:,icond)))');
%                     XLIM = log(squeeze(nanmean(dat_bins(:,4,1,:,1:2)))');
                    xdat = log(squeeze(nanmedian(dat_bins(:,4,1,:,icond)))');
                    XLIM = log(squeeze(nanmedian(dat_bins(:,4,1,:,1:2)))');
                    XLIM = [min(XLIM(:)) max(XLIM(:))];
                    XLAB = '8 - 12 Hz power (log)';
                case 'binno'
                    xdat = 1:nbins;
                    XLIM = [0 11];
                    XLAB = '8 - 12 Hz power (bin #)';
            end
            ydat = nanmean(dat(:,1:nbins));
            s(icond) = shadedErrorBar(xdat, ydat, sem, {'Color', linecol(icond,:)}, 1);
            sc = scatter(xdat, ydat, 'filled');
            sc.CData = linecol(icond,:);
            keepdat(:,:,icond) = dat;
            %             title([meas_leg{im} ' ' normalizetype])
            %             if iplot == 2
            %                 if zscore; xlabel('8 - 12 Hz power (log)'); else; xlabel('8 - 12 Hz power (log)'); end
            %                 if zscore; ylabel('59 - 100 Hz power (zscore)'); else; ylabel('59 - 100 Hz power (log)'); end
            %                 xlabel('8 - 12 Hz power (bin #)');
            %             xlabel('8 - 12 Hz power (V/m^2)');
                 
            xlabel(XLAB);
            if im == 8
                 ylabel('59 - 100 Hz power modulation (%)')
            elseif im == 11
                ylabel('DDM drift bias (Z-score)')
            elseif im == 3
                ylabel('SSVEP strength (%)')
            end
            %                  ylabel('Post-prestim gamma power (log)')
            %                  ylabel('Gamma power (psc)')
            %                 ylabel('59 - 100 Hz power (psc)')
            %                 if zscore; xlabel('8 - 12 Hz power (bin #)'); else; xlabel('8 - 12 Hz power (bin #)'); end
            %                 if zscore; ylabel('59 - 100 Hz power (zscore)'); else; ylabel('59 - 100 Hz power (log)'); end
            if icond == 2
                legend([s.mainLine], respavg.behav_conds(1:2), 'Location', 'Northeast'); legend boxoff
            elseif icond == 4
                legend([s.mainLine], respavg.behav_conds(4), 'Location', 'Northeast'); legend boxoff
            end
%             end
            %             xlim([-17.2 -14.4])
            %             ylim([-19.775 -19.5])
            h = gca;
%             h.XLim = [-16.7 -14.9];
            h.XLim = XLIM;
%             h.YLim = [-4 4];
%             h.YTick = -4:2:4;
%             h.YTick = -0.8:0.4:0.8;
% ref = refline(0,0);
% ref.Color = 'k';
% ref.LineStyle = '--';
%         end
    end
    conddiff = squeeze(dat_bins(:,4,im,:,2)) - squeeze(dat_bins(:,4,im,:,1));
    im
    [h,p]= ttest(conddiff)
end

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prealpha');
    
    mkdir(outpath)
    %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    outfile = fullfile(outpath, sprintf('alphagroupfits_%s_%s_%dbins_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim,iresp}, nbins, binmethod, normalizetype ));
    disp(outfile)
    %     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
%         print(outfile, '-dpdf') %'-png',  '-pdf',
        print(outfile, '-dpng') %'-png',  '-pdf',
    f.Renderer = 'Painters';
    print(outfile, '-depsc2') %'-png',  '-pdf',
    print(outfile, '-dpng') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% plot alpha distributions lib and cons
SAV = 1;
f = figure; hold on; iplot =0;
% f.Position = [ 1989         642         600         300];  %  567         314        1354         628
f.Position = [ 567         314        700         350];
iplot=iplot+1;
subplot(1,2,iplot); hold on; axis square; box on
h = []; dat = []; clear s;
for icond = 1:2%:2%:2%:2
    
    b = bar( mean(histedgedat(1:size(histedgedat,1)-1,icond,:),3), mean(histdat(:,icond,:),3), 'histc');
    b.FaceColor = linecol(icond,:);
    b.EdgeColor = b.FaceColor;
    b.FaceAlpha = 0.5;
    b.EdgeAlpha = b.FaceAlpha;
    h = findobj(gca,'Type','line');   set(h,'Marker','none');
    %             xlabel('8 - 12 Hz power (log)')
    xlabel('8 - 12 Hz power (V/m^2)');
    ylabel('Frequency of occurrence')
    title('Alpha power distribution')
    h = gca;
    legend(respavg.behav_conds{1:2}); legend boxoff
end

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prealpha');
    
    mkdir(outpath)
    %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    outfile = fullfile(outpath, sprintf('alphadistributions_%s_%s_%dbins_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim,iresp}, nbins, binmethod, normalizetype ));
    disp(outfile)
    %     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
%         print(outfile, '-dpdf') %'-png',  '-pdf',
        print(outfile, '-dpng') %'-png',  '-pdf',
    f.Renderer = 'Painters';
    print(outfile, '-depsc2') %'-png',  '-pdf',
    print(outfile, '-dpng') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% plot 3 model predictions: alpha down, gamma up, more gaussian
plotdat = [];
% plot alpha cons vs lib, get data
iband = 1; itrig = 1; ises=4;    istim = 3; isoi=1;  iresp = 3; 
sub_ind = [1:7, 9:15];

freq = [];
freq.dimord = 'subj_chan_freq_time';
freq.label = respavg.label;
freq.freq = respavg.freq{iband};
freq.time = respavg.time{itrig};
cfg = [];
cfg.channel = respavg.sens.ind{isoi}; cfg.avgoverchan = 'yes';
cfg.latency = [-0.8 -0.2]; cfg.avgovertime = 'yes';
cfg.frequency = [8 12];  cfg.avgoverfreq = 'yes';
alpha = {};
for icond = 1:2
    freq.powspctrm  = squeeze(respavg.pow(sub_ind,:, 1:length(respavg.freq{iband}), 1:length(freq.time),   iband, itrig,     ises,icond,istim,iresp ));
    alpha{icond} = ft_selectdata(cfg, freq);
end
plotdat(:,:,1) = [alpha{1}.powspctrm alpha{2}.powspctrm];

% test bins 4:7 against 1:3 and 8:10, gamma, binned by alpha % plot gamma cons vs lib, get data
im = 8; 
% dathi = squeeze(mean(dat_bins(:,4,im, 4:7, 1:2), 4));
plotdat(:,:,2) = squeeze(mean(dat_bins(:,4,im, 5:6, 1:2), 4));
% % dathi = squeeze(mean(dat_bins(:,4,im,  7, 1:2), 4));
% pperm = permtest(dathi(:,1), dathi(:,2))
% [~,p]= ttest(dathi(:,1), dathi(:,2))


% plot gaussianity cons vs lib get data

im = 8;
testdat = squeeze(dat_bins(:,4,im, :, 1:2));
refdat = [-1000        -991        -825         295        2521        2521         295        -825        -991       -1000];
% refdat = repmat(refdat, 14,1);
fitdat=[];
for icond = 1:2
    for isub = 1:14
        fitdat(isub,icond) = goodnessOfFit(testdat(isub,:,icond)', refdat', 'NMSE'); % NRMSE MSE NMSE
    end
end
plotdat(:,:,3) = fitdat

% mean(fitdat)
% pval = [];
% pval(1) = permtest(fitdat(:,1));
% pval(2) = permtest(fitdat(:,2));
% pval(3) = permtest(fitdat(:,1), fitdat(:,2));
% fitdat = fitdat-mean(fitdat,2) + mean(fitdat,1) ;

% % plot bars
% figure; barweb(mean(fitdat), std(fitdat) / sqrt(14))
% ax=gca;
% ax.FontSize = 12;
% title(pval)

close all
f = figure;
f.Position = [  680   937   450   200];
titleg = {'Prestimulus alpha', 'Poststimulus gamma', 'Gaussianity'};
yleg = {'Scalp current density (V/m2)', 'Modulation (%)', 'Gaussian goodness of fit' };
for iplot=1:3
    pval = [];
    pval(1) = permtest(plotdat(:,1,iplot), plotdat(:,2,iplot));
    if iplot == 2
        pval(2) = permtest(plotdat(:,1,iplot));
        pval(3) = permtest(plotdat(:,2,iplot));
    end
    
    subplot(1,3,iplot)
    semdat = plotdat(:,:,iplot) - mean(plotdat(:,:,iplot),2) + mean(plotdat(:,:,iplot),1);
    b= barweb(mean(plotdat(:,:,iplot)), std(semdat)/sqrt(14), 0.5);
    b.bars(1).FaceColor = linecol(1,:);
    b.bars(2).FaceColor = linecol(2,:);
    b.bars(1).EdgeColor = linecol(1,:);
    b.bars(2).EdgeColor = linecol(2,:);
    ax=gca;
    ax.FontSize = 8;
    title(sprintf('%s', titleg{iplot}))
    ylabel(yleg{iplot})
    if iplot == 2; legend(respavg.behav_conds{1:2}); legend boxoff; end
    xlabel(pval)
    if iplot == 1
       ylim([1.8e-7 2.2e-7])
    end
end
if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prealpha');
    
    mkdir(outpath)
    %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    outfile = fullfile(outpath, sprintf('threepredictions') );
    disp(outfile)
    
    print(outfile, '-dpng') %'-png',  '-pdf',
    f.Renderer = 'Painters';

        export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
%     print(outfile, '-dpdf', '-fillpage') %'-png',  '-pdf',
    print(outfile, '-depsc2') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% JJ idea: plot average dc and gamma over subj, corr bins

% close all
SAV=1;
f = figure;
% f.Position = [2080         438 300 150];
f.Position = [2080         438 600 150];
corrkeep = [];
% corrtype = 'Spearman';
corrtype = 'Pearson';
iplot=0;
correlationdata = [];
% linestyle = {'-k', '-g'};
linestyle = {'-', '--'};
linecol = cbrewer('qual', 'Set1',3);
linecol(4,:) = [0 0 0];
dotsymbol = {'s', '^'};
dotsize = 200;
% cmap = colormap(parula(nbins)); % colormap(parula(nbins));  %     cmap = colormap(cool(nbins));    %     cmap = colormap(jet(nbins));
cmap = zeros(nbins,3); % colormap(parula(nbins));  %     cmap = colormap(cool(nbins));    %     cmap = colormap(jet(nbins));
% for icond = [2,1,4]
for icond = [2,1]
%     iplot = iplot+1;
%     subplot(2,2,iplot); hold on;  box on; axis square;
%     colormap(cmap)
    ctr=0;    corrdat = []; clear s; clear sc
%     yyaxis left
    for im = [8 11]   %[1 11] %   %  %
%     for im = [ 3 11]   %[1 11] %   %  %
        ctr=ctr+1;
        
        temp = squeeze(dat_bins(:,4,im,:,1:2)); % -  squeeze(mean(dat_bins(:,4,im,:,icond),4)) 
%         temp(:,:,4) = squeeze(dat_bins(:,4,im,:,2)) - squeeze(dat_bins(:,4,im,:,1));
        temp(:,:,4) = squeeze(dat_bins(:,4,im,:,4));
        if im == 1
            corrdat(:,:,ctr) = log(temp(:,:,icond));
        else
            corrdat(:,:,ctr) = temp(:,:,icond);
        end
    end
    
    % plot single subj after demeaning both vars (rmcorr)
    corrdat = corrdat - mean(corrdat,2);
    corrdat = reshape(corrdat, [],2);
    
    mdl = fitlm( corrdat(:,1), corrdat(:,2) );
    rho_incoutliers = sqrt(mdl.Rsquared.Ordinary);

    outliers = find((mdl.Diagnostics.CooksDistance) > 5 *mean(mdl.Diagnostics.CooksDistance));
%     outliers = find((mdl.Diagnostics.CooksDistance) > (4/140));
    %     figure;    plotDiagnostics(mdl, 'cookd')
    corrdat = reshape(corrdat, 14,10,2);

    iplot = iplot+1;
    subplot(1,3,iplot); hold on; box on; axis square;
    colormap(jet(14))
    l=[]; colors=[]; yvals = [];
    
    % plot avg corr 
    corrdat = reshape(corrdat, [],2);
    x = corrdat(:,1); y = corrdat(:,2);
    [fitresult, gof] = fit(x,y, 'poly1');
%     [fitresult, gof] = fit(x,y,'poly2');
    confint = predint(fitresult,x,0.95,'functional','on');
    pl = plot(fitresult,x,y); hold on, plot(x,confint,'k--');
    pl(1).Marker = 'o';
    pl(1).MarkerSize = 4;
    pl(1).Color = linecol(icond,:);
    pl(2).Color = 'k';    % pl(2).Color = linecol(icond,:);
    pl(1).LineWidth = 0.5;
    pl(2).LineWidth = 0.5;
    legend off
%     text(0.25,1,sprintf('R^2 = %g', gof.rsquare))
%     text(0.25,1,sprintf('r = %g', sqrt(gof.rsquare)))

    corrdat = reshape(corrdat, [],2);    
    sc = scatter(corrdat(outliers,1), corrdat(outliers,2), 'x', 'MarkerEdgeColor', [0 0 0], 'Sizedata', 24); % [0.5 0.5 0.5] 
    
%     corrdat(outliers,:) = NaN;

    mdl = fitlm( corrdat(:,1), corrdat(:,2) );
    corrdat = reshape(corrdat, 14,10,2);
 
    rho_outliersremoved = sqrt(mdl.Rsquared.Ordinary);
    pvals = [9e-6 5e-9 nan 0.24 ]; % from rmcorr
    text(-23, 1.25, sprintf('linear r = %1.2f, p = %g', rho_incoutliers, pvals(icond)), 'FontSize', 8, 'FontWeight', 'bold' )
    title(respavg.behav_conds{icond} )
%     title(sprintf('linear r = %1.2f (%d removed: r = %1.2f)', rho_incoutliers, length(outliers), rho_outliersremoved))
%     title(sprintf('%s, linear r = %g', respavg.behav_conds{icond}, rho_incoutliers))
%     title(sprintf('%s, linear r = %g', respavg.behav_conds{icond}, rho_incoutliers))

%     xlabel('Residual 59-100 Hz modulation (%)')
    xlabel('SSVEP strength (%)')
%     xlabel('Residual 8-12 Hz power (log)')
    ylabel('Residual drift bias (DDM p.e.)')
    ax=gca;
    ax.FontSize = 8;
     
    if strcmp(gammafocus, 'transient')
%         if iplot == 1; ax.XLim = [-20 15]; else ax.XLim = [-25 15]; end
%         ax.XLim = [-25 15];
        ax.XLim = [-20 20];
        ax.YLim = [-1.7 1.5];
        ax.YTick = -2:2;
    else
        ax.XLim = [-2 2];
    end
%     if icond == 2; ax.YLim = [-1.5 1.1]; else ax.YLim = [-2 2]; end
    
    [ outdat ] = export_for_rmcorr(corrdat(:,:,1), corrdat(:,:,2), [respavg.behav_conds{icond} '_gamma'], 'driftbias', '/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/data_rmcorr/gammavsdriftbias');
    correlationdata.dat(:,:,:,icond) = corrdat;

end

correlationdata.dimord = 'subj_bins_gammadrift_conslib'

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prealpha');
    
    mkdir(outpath)
    %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    outfile = fullfile(outpath, sprintf('gammavsdcsubjavg_%s_%s_%dbins_zacrosstrials', respavg.sens.leg{isoi}, respavg.sdt_conds{istim,iresp}, nbins) );
    disp(outfile)
    
    print(outfile, '-dpng') %'-png',  '-pdf',
    f.Renderer = 'Painters';

        export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
%     print(outfile, '-dpdf', '-fillpage') %'-png',  '-pdf',
    print(outfile, '-depsc2') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% plot single subj to look at shape of response
f = figure;
for icond = 4
    for isub = 1:14
        subplot(4,4,isub); axis square
        scatter(corrdat(isub,:,1), corrdat(isub,:,2))
    end
end



%% 3 way ANOVA task type (liberal / conservative) X measure (pre-stim alpha / post-stim gamma) x bins?
% run with bins
% 1D vectors: responses (14x10 blocks), do in spss
anovadat = [];         
figure; iplot=0;
for icond = 1:2
%     for im = [1,3]
%     for im = [1,11]
    for im = [1,8]
        iplot=iplot+1;
        temp = squeeze(dat_bins(:,4,im,:,icond));
        anovadat = [anovadat temp];
        subplot(2,2,iplot)
        plot(mean(temp))
        %         
    end
end

%% bins 1 way ANOVA on liberal-conservative 
% run with bins
% 1D vectors: responses (14x10 blocks), do in spss
anovadat = [];         
figure; iplot=0;
% for icond = 1:2
    for im = [1,8]
        iplot=iplot+1;
        temp = squeeze(dat_bins(:,4,im,:,2)) - squeeze(dat_bins(:,4,im,:,1));
        anovadat = [anovadat temp];
        subplot(2,2,iplot)
        plot(mean(temp))
        %         
    end
% end




%% old, not used
%     coeffs = polyfit(corrdat(:,1), corrdat(:,2), 1);    
%     [~,~,ci] = ttest(yvals); % plot conf int bowtie
%     % %     ci(1,:) = ci(1,:) - polyval(coeffs, xvals);
%     % %     ci(2,:) = ci(2,:) + polyval(coeffs, xvals);
%     %
%     %     shadedErrorBar(xvals, polyval(coeffs, xvals), flipud(ci)  )  %
%     % %     figure; plot(1:11,ci )
%     
%     pat = patch([xvals(1) xvals(end) xvals(end) xvals(1)]', ...
%         [ci(1,1) ci(2,end) ci(1,end)  ci(2,1)]', 'k' );
%     pat.FaceColor = linecol(icond,:);
%     pat.EdgeColor = 'none'; %[0.5 0.5 0.5];
%     pat.FaceAlpha = 0.25;
% 
%     pl = plot(xvals, polyval(coeffs, xvals));
%     pl.Color = linecol(icond,:);
%     pl.LineWidth = 3;

%     sc = scatter( corrdat(:,1),  corrdat(:,2) );
%     sc.MarkerEdgeColor = 'none';
%     l = lsline;
%     for isub = 1:14
%         l(isub).Color = colors(isub,:);
%         l(isub).LineWidth = 0.5;
%     end
%     l(15).LineWidth = 2;
%     l(15).Color = 'k';    
    %         if im == 11
%             yyaxis right
%         end
%         s(ctr) = shadedErrorBar(xdat, mean(ydat), std(ydat)/sqrt(length(ydat)), [], 1);
%         
%         s(ctr).mainLine.Color = linecol(icond,:);
%         s(ctr).mainLine.LineStyle = linestyle{ctr};
%         s(ctr).mainLine.LineWidth = 2;
%         s(ctr).patch.FaceColor = linecol(icond,:);
%         
%         ax = gca;
%         ax.YColor = [0 0 0];
        
        corrdat(:,:,ctr) = temp(:,:,icond);
%         if im == 8
% %             ylabel('Gamma power (z-score)')
%             ylabel('Alpha power (z-score)')
%         elseif im == 11
%             ylabel('Drift bias (DDM estimate)')
%             if icond == 2
%                 ylim([2 3])
%             else
%                 ylim([-0.5 0.5])
%             end 
%         end

    %     legend([s.mainLine], {'Gamma power', 'Drift bias'}, 'Location', 'South'); legend boxoff
%     %     legend(sc, {'Gamma power', 'Drift bias'}, 'Location', 'Southwest'); legend boxoff
%     xlabel('8 - 12 Hz power (log)')
% %     xlim([-16.7 -14.9])
%     xlim([0 11])

    corrdat = reshape(corrdat, [],2);    
    scatter(corrdat(outliers,1), corrdat(outliers,2), 'x')
    
    corrdat(outliers,:) = NaN;
    mdl = fitlm( corrdat(:,1), corrdat(:,2) );
    corrdat = reshape(corrdat, 14,10,2);
    
    %     scatter(mean(corrdat(:,:,1)), mean(corrdat(:,:,2)))
    
    %     l(15).LineWidth = 2;
    %     l(15).Color= 'k';
    
%     % compute pval given r and df    
    rho_outliersremoved = sqrt(mdl.Rsquared.Ordinary);
%     nrepl = nbins * 14; % 140 for 10 bin
%     % convert correlation coefficient to t-statistic (for MCP correction): t^2 = DF*R^2 / (1-R^2)
%     t = rho*(sqrt(nrepl-2))./sqrt((1-rho.^2));    
%     
%     v = nrepl - 14 - 1 ;
%     %     tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution, where ?t? is the t-statistic and ?v? the degrees-of-freedom.
%     tail2P = 2*tcdf(-abs(t),v)
%     
%     title(sprintf('r(%d) = %g, p = %g', v, rho, tail2P))
    title(sprintf('r = %g (%d removed: r = %g)', rho_incoutliers, length(outliers), rho_outliersremoved))

%     % plot scatter gamma vs db after averaging subj
%     correlationdata.dat(:,:,:,icond) = corrdat;
%     corravg = squeeze(nanmean(corrdat));
%     corrkeep(:,:,icond) = corravg;
%     xsem = nanstd(corrdat(:,:,1)) ./ sqrt(14) /2;
%     ysem = nanstd(corrdat(:,:,2)) ./ sqrt(14) /2;
    
%     iplot = iplot+1;
%     subplot(2,2,iplot); hold on; box on; axis square;
    
    
%     
%     e = errorbar(corravg(:,1),corravg(:,2), -ysem, ysem, -xsem, xsem, '.');
%     e.Color = linecol(icond,:); % same color for error bars
%     %     e.Color = 'k'; % same color for error bars
%     
%     [r,p] = corr(corravg(:,1),corravg(:,2), 'type', corrtype );
%     [b,bint,res,rint,stats] = regress(corravg(:,2),[corravg(:,1) ones(size(corravg(:,1)))] )     % do regression
%     %     text(-0.6, 0.1, sprintf('R^2 = %1.2f\np = %1.4f', stats(1), stats(3)))
%     
%     % plot points, Color points by excitability level
%     for ibin=1:nbins
%         e = errorbar(corravg(ibin,1),corravg(ibin,2), -ysem(ibin), ysem(ibin), -xsem(ibin), xsem(ibin), '.');
%         e.Color = cmap(ibin,:);
%     end
%     colormap(cmap)
%     scatter(corravg(:,1),corravg(:,2), dotsize, 1:nbins, 'filled');
% 
%     
% %     scatter(corravg(:,1),corravg(:,2), 'filled', 'Markerfacecolor', linecol(icond,:),  ...
% %         'Markeredgecolor', [1 1 1], 'sizedata', dotsize, 'linewidth', 0.5) ;
%     
%     %     title(sprintf('%s r = %g, p = %g', respavg.behav_conds{icond}, r, p))
%     title(sprintf('%s R^2 = %1.2f\np = %1.4f', respavg.behav_conds{icond}, stats(1), stats(3)))
    
    %     xlabel('Gamma activity (Z-score)')
    %     ylabel('Drift bias (Z-score)')
%     xlabel('Gamma power (z-score)')
    xlabel('Gamma power (psc)')
%     xlabel('Alpha power (z-score)')
    ylabel('Drift bias (DDM estimate)')
    
%     xlim([-0.5 0.5]);
%     ylim([-1 1]);
    
%     if strcmp(corrtype, 'Pearson')
%         if p < 0.05
%             l = lsline;
%             l.Color = 'k';
%             l.LineWidth = 1;
%         end
%     end


%% export for R in column format

% dat = squeeze(dat_bins(:,4,[1,8],:,:)); % subj alphagamma bin cond
% subjcol = repmat(transpose(1:14), 1, size(dat,2), size(dat,3), size(dat,4)  );
% condcol = repmat(1:2, size(dat,1), 1, size(dat,3), size(dat,4) );

% alpha, gamma both conds
outdat = [];
dat = squeeze(dat_bins(:,4,1,:,:)); % subj alphagamma bin cond
subjcol = repmat(transpose(1:14), 1, size(dat,2), size(dat,3) );
condcol = repmat(1:2, size(dat,1), size(dat,2), 1 );
bincol = repmat(1:nbins, size(dat,1), 1, size(dat,3) );
outdat = [subjcol(:) condcol(:) bincol(:) dat(:)];
dat = squeeze(dat_bins(:,4,8,:,:)); % subj alphagamma bin cond
outdat(:,5) = dat(:);

exportpath = '/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/data_R';
mkdir(exportpath)
% write header and matrix to file
outfile = sprintf('critEEGalpha_gamma_cond_bin2.csv');
fid = fopen(fullfile(exportpath, outfile), 'w') ;
fprintf(fid, '%s,%s,%s,%s,%s\n', 'subj_idx', 'Cond_Cons1_Lib2', 'alphabin', 'alpha', 'gamma');
dlmwrite(fullfile(exportpath, outfile), outdat, '-append', 'precision', 10); % outmat2 !!!
fclose(fid);

fprintf('Wrote %s to %s\n', outfile, exportpath)








%% 2 line fitting 
if nbins == 10
    leftbins = 1:5; rightbins = 6:10;
elseif  nbins == 5
    leftbins = 1:3; rightbins = 3:5;
end

im = 8; % gamma
pol = [];
for isub=1:14
    for icond = 1:2
        temp = squeeze(dat_bins(:,4,im,:,icond)); % -  squeeze(mean(dat_bins(:,4,im,:,icond),4))
        % temp = squeeze(dat_bins(:,4,im,:,2)) - squeeze(dat_bins(:,4,im,:,1)); % -  squeeze(mean(dat_bins(:,4,im,:,icond),4))
%         figure; plot(temp')

        pol(isub,icond,1,:) = polyfit(1:length(leftbins), temp(isub,leftbins),1);
        pol(isub,icond,2,:) = polyfit(1:length(rightbins), temp(isub,rightbins),1);
        
        
    end
end
pol = pol(:,:,:,1); % only keep slope
pol = reshape(pol, 14, [])

% % 
% % left = mean(pol(:,1))
% % pl = permtest(pol(:,1))
% % right = mean(pol(:,2))
% % pr = permtest(pol(:,2))
% 
% [ outdat ] = export_for_rmcorr(squeeze(dat_bins(:,4,im,3:5,1)), squeeze(dat_bins(:,4,im,3:5,2)), 'cons_binnedgamma', 'lib_binnedgamma', '/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/data_rmcorr/gammabinnedbyalpha');

%% correlate dc with gamma , corr within subj
close all
for icond= 1:2
    corrdat = [];
    corrdat(:,1,:) = squeeze(dat_bins(:,4,8,:,icond)); % subj  bins gamma/dc
    corrdat(:,2,:) = squeeze(dat_bins(:,4,11,:,icond));
    
    
    r = [];
    f = figure;
    f.Position=[ 680         112        1019         986]
    for isub = 1:14
        subplot(4,4,isub); axis square; box on
        temp = [squeeze(corrdat(isub,1,:)) squeeze(corrdat(isub,2,:))];
%                 temp = (temp - nanmean(temp)) ./ nanstd(temp); % zscore
        
        %     temp = [squeeze(corrdat(isub,1,2:9)) squeeze(corrdat(isub,2,2:10))];
        %     temp = temp(~isnan(temp(:,2)) , :)
        scatter( temp(:,1), temp(:,2));
        
        r(isub,:) = corr( temp(:,1), temp(:,2), 'type', 'Spearman');
        title(r(isub,:))
        lsline
        
    end
    nanmean(r)
    permtest(r)
end
% corrdat = corrdat(~isnan(corrdat(:,2)),:)
%
% figure; scatter(corrdat(:,1), corrdat(:,2))
% title(corr(corrdat(:,1), corrdat(:,2)))

%% fit gaussian to single ssubj data

close all
sigma = nan(14,2);
for isub = 1:14
    for icond= 1:2 %1:2
        
        x = squeeze(dat_bins(isub,4,1,:,icond));
        y = squeeze(dat_bins(isub,4,8,:,icond)); % subj  bins gamma/dc
        y = y - min(y);
        %         figure; plot(subjdat)
        try
            %         f = fit( (-4.5:4.5)', subjdat, 'gauss1');
            f = fit( x, y, 'gauss1');
            sigma(isub,icond) = f.c1;
            figure; plot(f,x,y)
            title(f.c1)
        catch
            disp ('fit no')
            figure; scatter(x,y)
            
        end
        
        % [mu, sh] = normfit(  subjdat)
% sigmahat(isub, icond) = sh;
        
%         %     r = [];
%         %     f = figure;
%         %     f.Position=[ 680         112        1019         986]
%         subplot(4,4,isub); axis square; box on
%         temp = [squeeze(corrdat(isub,1,:)) squeeze(corrdat(isub,2,:))];
%         %                 temp = (temp - nanmean(temp)) ./ nanstd(temp); % zscore
%         
%         %     temp = [squeeze(corrdat(isub,1,2:9)) squeeze(corrdat(isub,2,2:10))];
%         %     temp = temp(~isnan(temp(:,2)) , :)
%         scatter( temp(:,1), temp(:,2));
%         
%         r(isub,:) = corr( temp(:,1), temp(:,2), 'type', 'Spearman');
%         title(r(isub,:))
%         lsline
        
    end
%     nanmean(r)
%     permtest(r)
end

 sigma(sigma(:,1)>20,1) = nan;
 sigma(sigma(:,2)>20,2) = nan;
 sigma
nanmean(sigma)








%% correlate alphapow and gammapow with crit across subjects
%%% behav.criterion(:,:,4) = behav.criterion(:,:,2) - behav.criterion(:,:,1); % Nomiss -  NOFA
close all

% critdat = squeeze(respavg.behavior.criterion([1:13, 15,16],4,1:2)); % Steph out
% critdat = (critdat - critdat(:,1)) ./ critdat(:,1) * 100;
% take ddm dc
critdat = squeeze(dat_bins(:, 4, 11, 1:5, 2));
critdat = nanmean(critdat,2);

% take ddm dc, take psc
critdat = squeeze(dat_bins(:, 4, 11, 1:5, 1:2));
critdat = squeeze(nanmean(critdat,2));
critdat = critdat(:,2) ./ critdat(:,1);

% pow = (avgalphapow - avgalphapow(:,1)) ./ avgalphapow(:,1) * 100; % prestim alpha!
pow = gammapow; % (gammapow - gammapow(:,1)) ./ gammapow(:,1) * 100; % poststim gamma!
% pow = log(avgalphapow); % (gammapow - gammapow(:,1)) ./ gammapow(:,1) * 100; % poststim gamma!

% incsubj = [1:7, 9:12, 14,15]; % 1 outlier out + the others
% incsubj = [1:4, 6, 7, 9:12, 14,15]; % 1 outlier out + the others
incsubj = [1:7, 9:15]; % 1 outlier out + the others
pow = pow(incsubj,:);
% pow = (pow(:,2) -  pow(:,1)) ./ mean(pow,2) * 100;
pow = (pow(:,2) -  pow(:,1)) ./ pow(:,1) * 100;
% pow = (pow - mean(pow,2)) ./ std(pow,0,2)
% pow = pow(:,2) - pow(:,1);
%
% in = critdat < 15;
% pow = pow(in);
% critdat = critdat(in);
%
% in = pow < 40;
% pow = pow(in);
% critdat = critdat(in);

% [r,p] = corr(critdat([1:7, 9:15],2) - critdat([1:7, 9:15],1), pow([1:7, 9:15],2) -  pow([1:7, 9:15],1), 'type', 'spearman')
[r,p] = corr(critdat, pow, 'type', 'pearson')

figure; hold on; box on; axis square
% scatter(critdat([1:7, 9:15],2) - critdat([1:7, 9:15],1), pow([1:7, 9:15],2) -  pow([1:7, 9:15],1), 400, ...
%     'MarkerEdgecolor', [1 1 1], 'MarkerFacecolor', [0 0 0], 'Linewidth', 3)
scatter(critdat , pow, 400, ...
    'MarkerEdgecolor', [1 1 1], 'MarkerFacecolor', [0 0 0], 'Linewidth', 3)
% text(behavdat, alphapowdat, respavg.SUBJ );

ax=gca
ax.FontSize=20;
title(sprintf('r = %g, p = %g', r,p))

xlabel('criterion change lib vs cons (%)')
% xlabel('dprime change lib vs cons (%)')
% ylabel('alpha power change lib vs cons (%)')
ylabel('gamma power change lib vs cons (%)')
lsline


%% put matrix in spss format: copy paste to spss
% kol 1: subject ID: 1 1 1 1 1 ...
% kol 2: binno: 1 2 3 4 5 6 ...
% kol 3: gamma values
% kol 4: dc values
close all
spssdat = [];
for isub = 1:14
    ctr=0;
    dat = [];
    for iz = 1:2%:2 % raw / zscored
        for icond = 1:2
            for im = [8 11]
                
                temp = squeeze(dat_bins(isub,4,im,:,icond)); % subj gamma/dc bins
                if iz == 2
                    temp = (temp - nanmean(temp)) ./ nanstd(temp); % zscore
                end
                
                ctr = ctr + 1;
                dat(:,ctr+2) = temp;
                dat(:,1) = isub;
                dat(:,2) = 1:nbins;
                %                 dat(:,2) = 1:10;
            end
        end
    end
    spssdat = [spssdat; dat];
end

% figure; scatter(spssdat(:,3), spssdat(:,4))
% figure; scatter(spssdat(:,5), spssdat(:,6))
%
% [r,p] = corr(spssdat(:,3), spssdat(:,4))

%% plot ratio lib vs cons: gain increase?
% close all
figure; hold on; axis square; box on
plotdat = keepdat; % + abs(min(min(keepdat)));
uturn  = (plotdat(:,2) ./ plotdat(:,1))';
% plot(uturn)
s = scatter(1:nbins, uturn, 'filled'  );
s.MarkerFaceColor = 'k';
s.MarkerEdgeColor = [1 1 1];
s.SizeData = 250

p = polyfit(transpose(1:nbins), uturn', 1);
y = polyval(p, 1:nbins);
pl = plot(1:nbins, y);
pl.Color = 'k';
pl.LineStyle = '--';

p = polyfit(transpose(1:nbins), uturn', 2);
y = polyval(p, 1:nbins/20:nbins);
pl = plot(1:nbins/20:nbins, y);
pl.Color = 'k';

xlabel('Neural Excitability')
ylabel('Gain ratio (liberal / conservative)')
ax=gca;

xlim([1 10])
ax.XTick = [1 5.5 10]
ax.YTick = 0.9:0.1:1.4

% xlim([1 5])
% ax.XTick = [1 3 5];
% ax.YTick = 0.9:0.1:1.4;
% order = sequential_regression(1:10, uturn) % do in python

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prealpha');
    
    mkdir(outpath)
    %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    outfile = fullfile(outpath, sprintf('gainratio_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim,iresp}) );
    disp(outfile)
    %     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    print(outfile, '-dpdf') %'-png',  '-pdf',
    print(outfile, '-dpng') %'-png',  '-pdf',
    print(outfile, '-depsc2') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% 2 way ANOVA task type (liberal / conservative) X measure (pre-stim alpha / post-stim gamma)?
% run without bins
% mat: all lib, all con x alpha / gamma: 28 x 2
anovadat = [];
ims = [1,8] % alpha gamma
for im = 1:length(ims)
    temp = squeeze(mean(dat_bins(:,4,im,:,1:2),4)); %avg over bins
    anovadat(:,im) = temp(:);
end

[p,tbl] = anova2(anovadat,14)


%% bin_fits: test cubic fit coeff high vs low crit. Higher gain for low crit?
%  bin_fits(isub,:,im,icond)
close all
SAV=1
linecol = cbrewer('qual', 'Set1',3);
im = 8;
% ylims = [-1.5e-3 0.5e-3; -0.005 0.02; -0.0125 0.05] % not in psc yet
ylims = [-1.5e-1 0.5e-1; -0.5 2; -1.25 5]
% figure; bar(mean(fitdat))
f = figure; ipl = 0;
f. Position= [2095         556 240*3         120*3]
% f. Position= [2095         556 240         120]
for coeff = 1:3 % 1 = quadratic fit, 2 = linear, 3 = intercept
    fitdat = squeeze(bin_fits(:,coeff,im,1:2))
    ipl = ipl+1;
    subplot(1,3,ipl);
    b = barweb(mean(fitdat), std(fitdat) / sqrt(14), 0.5)
    
    %     b.ax.YLim = ylims(coeff, :);
    b.bars(1).FaceColor = linecol(1,:);
    b.bars(2).FaceColor  = linecol(2,:);
    %     box on
    %     randtest1d(fitdat(:,1) - fitdat(:,2), zeros(16,1), 0, 10000)
    p= []; stats=[];
    %     [~,p(1),~,stats{1}]= ttest(fitdat(:,1), fitdat(:,2))
    %     [~,p(2),~,stats{2}]= ttest(fitdat(:,1))
    %     [~,p(3),~,stats{3}]= ttest(fitdat(:,2))
    [p(1)]= permtest(fitdat(:,1), fitdat(:,2))
    [p(2)]= permtest(fitdat(:,1))
    [p(3)]= permtest(fitdat(:,2))
    title(sprintf('coeff%d p = %g\n%g  %g', coeff, p ))
    ax=gca;
    ax.FontSize= 12;
    %      axis square;
    
    %     % correlate gain to prestim alpha, hi-lo crit
    %     alphadat = mean(dat_bins(:, 1, 1:nbins, 4),3); % only sig bins
    %     ipl = ipl+1;
    %     subplot(3,2,ipl); hold on
    %     scatter(alphadat, fitdat(:,1)-fitdat(:,2))
    %     axis square; box on
    %     [rho, pval] = corr(alphadat, fitdat(:,1)-fitdat(:,2), 'type', 'Pearson')
    %     title(sprintf('r = %g, p = %g', rho, pval ))
    %         ax=gca;
    %     ax.FontSize= 12;
    %
    %     lsline
    
end

%% correlate gain to gamma
im =8;
f = figure;
for coeff = 1:3
    fitdat = squeeze(bin_fits(:,coeff,im,1:2));
    dc = squeeze(dat_bins(:,4,11,:,1:2));
    dc = squeeze(mean(dc,2));
    [r,p] =corr(fitdat(:,2)-fitdat(:,1), dc(:,2)-dc(:,1) , 'type', 'Spearman' )
    subplot(1,3,coeff);
    scatter(fitdat(:,2)-fitdat(:,1), dc(:,2)-dc(:,1)  )
    lsline
    axis square; box on
    xlabel('coeff lib vs cons')
    ylabel('dc lib vs cons')
    title(sprintf('r = %g, p = %g', r, p))
end

%%
if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prealpha');
    
    mkdir(outpath)
    %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    outfile = fullfile(outpath, sprintf('Quadrfitcoeff_%s', respavg.sens.leg{isoi}) );
    disp(outfile)
    %     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    %     export_fig(outfile, '-eps', '-transparent') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    print(outfile, '-dpdf')
    cd(outpath)
end %ises


% %% bin_fits: test cubic fit coeff high vs low crit. Higher gain for low crit?
% % all coeff in 1 barweb plot
% im = 8;
% % figure; bar(mean(fitdat))
% f = figure; ipl = 0;
% f. Position = [680 472 1055 626]
% fitdat = squeeze(bin_fits(:,:,im,1:2))
% % ipl = ipl+1;
% % subplot(3,2,ipl);
% barweb(squeeze(mean(fitdat)), squeeze(std(fitdat)) / sqrt(16), 0.5)
% box on
% %     randtest1d(fitdat(:,1) - fitdat(:,2), zeros(16,1), 0, 10000)
% p= [];
% [~,p(1)]= ttest(fitdat(:,1), fitdat(:,2))
% [~,p(2)]= ttest(fitdat(:,1))
% [~,p(3)]= ttest(fitdat(:,2))
% title(sprintf('coeff%d p = %g\n%g  %g', coeff, p ))
% ax=gca;
% ax.FontSize= 12;
% axis square;
%

%%

%% plot single subj dat_bins
%
% %dimord subj measureno binno cond
% figure
% for isub = 1:16
%     subplot(4,4,isub)
%     dat =  squeeze(dat_bins(isub,1,:,1));
%     plot(1:nbins, dat)
% end






% %% Plotting dc for the bins, along with across bin correlations
% close all
% SAV = 1;
% nbins = size(bin_ranges,1);
% % figure; plot(squeeze(mean(dat_bins)))
% % close all
% % linecol = {'r' 'b' 'k'};
% linecol = cbrewer('qual', 'Set1',3);
% linecol(4,:) = [0 0 0];
% set(0, 'Defaultaxesfontsize', 12)
% f = figure; hold on; iplot =0;
% % f.Position = [ 1989         642         600         300];  %  567         314        1354         628
% f.Position = [ 2088         166         600         600];
% keepdat = []; % for ratio lib vs cons
% for im = [8 11] % [1,8,9] % [1:4, 6:9] 1:8% 1:length(meas_leg)  % 1:5 %2:5  % ssvep response %2:5
%     for icond = 2:-1:1 % 1:2%:2
%         iplot=iplot+1;
%         subplot(2,2,iplot); hold on; axis square; box on
%         xlim([0 nbins-1])
%         %         if iplot > 1
%         %             r = refline(0,0);
%         %             r.LineStyle = '--';
%         %             r.Color = 'k';
%         %         end
%         h = []; dat = []; clear s;
%         %             % within subj error bars
%         temp = squeeze(dat_bins(:,4,im,:,1:2)); % -  squeeze(mean(dat_bins(:,4,im,:,icond),4))
%         condavg = squeeze(nanmean(temp,3));
%         grandavg = nanmean(condavg);
%         dat = temp(:,:,icond) - condavg + grandavg;
%         
%         %             dat(:,icond) = squeeze(mean(dat_bins(:,4,im,:,icond))); %             actual data
%         %             sem = squeeze(std(dat_bins(:,4,im,1:nbins,icond))) / sqrt(size(dat_bins(:,4,im,:,icond),1));
%         %             dat(:,icond) = flipud(dat(:,icond));
%         %             s(icond) = shadedErrorBar(0:nbins-1, dat(1:nbins,icond), sem, {'Color', linecol(icond,:)}, 1);
%         
%         dat = fliplr(dat);
%         sem = nanstd(dat,0,1) / sqrt(size(dat,1));
%         s(icond) = shadedErrorBar(0:nbins-1, nanmean(dat(:,1:nbins)), sem, {'Color', linecol(icond,:)}, 0);
%         scatter(0:nbins-1, nanmean(dat(:,1:nbins)), 'filled', 'MarkerFaceColor', linecol(icond,:))
%         keepdat(:,icond) = nanmean(dat);
%         
%         ax=gca;
%         ax.XLim = [0 nbins-1];
%         ax.XTick = [0 (nbins-1)/2 nbins-1];
%         ax.XTickLabel = {'Low' 'Medium' 'High'};
%         if icond == 1
%             %             ylim([0 1.5])
%             %             ylim([0 1.25])
%             %             ax.YTick = [0:0.5:3];
%         else
%             %             ylim([2 2.75])
%             %             ax.YTick = [2:0.25:3];
%         end
%         
%         
%         legend([s.mainLine], respavg.behav_conds(icond)); legend boxoff
%         xlabel('Neural excitability')
%         %         ylabel('DDM drift bias')
%         ylabel(meas_leg{im})
%         
%     end
%     
%     out = randtest_corr(corrdatcond{2},corrdatcond{1},1,1000, 'Pearson')
% end
% 
% if SAV
%     %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
%     outpath = fullfile(respavg.PREOUT, 'prealpha');
%     
%     mkdir(outpath)
%     %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
%     outfile = fullfile(outpath, sprintf('gammavsdc_%s_%s_%dbins', respavg.sens.leg{isoi}, respavg.sdt_conds{istim,iresp}, nbins) );
%     disp(outfile)
%     %     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
%     print(outfile, '-dpdf', '-fillpage') %'-png',  '-pdf',
%     print(outfile, '-dpng') %'-png',  '-pdf',
%     print(outfile, '-depsc2') %'-png',  '-pdf',
%     %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
%     cd(outpath)
% end %ises

%% old
%     dat_bins(:,:,8,:,:) = (dat_bins(:,:,8,:,:) - meandat(:,:,8)) ./ sddat(:,:,8);

% dat_bins = dat_bins - mean(dat_bins(:,:,:,:,1),4); % only subtract
% dat_bins = dat_bins - max(dat_bins(:,:,:,:,2),[],4); % take max for lib to normalize


% dat_bins = dat_bins - min(dat_bins(:,:,:,:,1),[],4); % take min for cons to normalize

% dat_bins = dat_bins ./ max(dat_bins(:,:,:,:,2),[],4); % divide by max for lib to normalize, stuff flipped bc negative

% dat_bins(:,:,8,:,:) = dat_bins(:,:,8,:,:) - min(dat_bins(:,:,8,:,1),[],4); % only subtract for gamm
% dat_bins(:,:,11,:,:) = dat_bins(:,:,11,:,:) - min(dat_bins(:,:,11,:,1),[],4); % only subtract for gamm

% dat_bins = dat_bins - min(dat_bins(:,:,:,:,1),[],4); % take min for cons to normalize


% permtest(bin_fits(:,1,8,2)) %=    0.0307 for lib!
% permtest(bin_fits(:,1,8,1)) %=    0.4195 for cons
% 
% leftslope =  bin_fits(:,2,8,2) + (2*bin_fits(:,1,8,2).*squeeze(dat_bins(:, 4, 8, 1, 2)) ); % ?1?+?2?2XL
% permtest(leftslope)  %  0.0364
% rightslope = bin_fits(:,2,8,2) + (2*bin_fits(:,1,8,2).*squeeze(dat_bins(:, 4, 8, nbins, 2)) ); %  ?1?+?2?2XH
% permtest(leftslope)  %  0.0388
% turningpoint = -1 * ( bin_fits(:,2,8,2)  ./2 .* bin_fits(:,1,8,2)); % ???1/2?2.
% [h,p,ci] = ttest(turningpoint)
% 

% fooling with the log for normalization
%                         % convert to nV/m2 and take log
%                         prestimdat = log1p(prestimdat * 1e9);
%                         poststimdat = log1p(poststimdat * 1e9);
% %                         %                         % convert to nV/m2 and take log
%                         prestimdat = log1p(prestimdat );
%                         poststimdat = log1p(poststimdat );
%                         % convert to nV/m2 and take log
%                         prestimdat = log(1+prestimdat );
%                         poststimdat = log(1+poststimdat );
                        
                        % TODO reject outliers before taking avg for psc
                        % old:
                        %                             zdatpre = (log(prestimdat)  - nanmean(log(prestimdat) )) ./ nanstd(log(prestimdat )); % take log to make normal
                        %                             zdatpost = (log(poststimdat)  - nanmean(log(poststimdat) )) ./ nanstd(log(poststimdat )); % take log to make normal
                        %                         prestimdat = prestimdat(abs(zdatpre) < outlierzthr,:); % remove outliers, take abs to remove both sides
                        %                         poststimdat = poststimdat(abs(zdatpost) < outlierzthr,:); % remove outliers, take abs to remove both sides
                        %                         temp = temp(abs(zdatpost) < outlierzthr,1);
                        
%                         % take log and shift to positive range
%                         poststimdat = log(poststimdat);
%                         poststimdat(isinf(poststimdat))  = NaN;
%                         prestimdat = log(prestimdat);
%                         prestimdat(isinf(prestimdat)) = NaN;
%                         
%                         postdatminimum = min(poststimdat);
%                         poststimdat = poststimdat - postdatminimum;
%                         prestimdat = prestimdat - postdatminimum;
%                         %      