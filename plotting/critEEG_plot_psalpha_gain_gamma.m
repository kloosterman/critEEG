%% PLot prestim alpha TFR + topo, then sorted alpha power, then gain in gamma + bars of quadratic coeff

% respavg = critEEG_load_respavg()

%% plot binned by alpha: alpha itself and gamma, 10 bins
% respavg = critEEG_load_respavg()

meas_leg = {'Alpha pre' 'Abeta' 'SSVEP' 'P1' 'N1' 'CPP' 'SSVEP high' 'gamma' 'allhighfreq' 'sdt' 'ddm_dc'}

nbins = 10;
if nbins == 10
    bin_ranges = [0    15;   5    25;   15    35;    25    45;    35    55;    45    65;   55 75; 65 85; 75 95;  85 100];
elseif nbins == 5
    bin_ranges = [0 33;  16 50; 33 67; 50 83; 66 100 ];
elseif nbins == 3
    bin_ranges = [0 33; 33 66; 66 100];    
end

dat_bins = nan(length(respavg.SUBJ), 3, 11, nbins, 2); %dimord subj ses measureno binno cond
dat_bin_count  = nan(length(respavg.SUBJ), 3, 11, nbins, 2); %for histogram of alpha vals
dat_bin_edges  = nan(length(respavg.SUBJ), 3, 11, nbins, 2, 2); %for histogram of alpha vals

istim = 3; iresp = 3;

isoi = 1; %     'occipital'    'motor'    'allsens'    'occpar'    'frontal'    'occlatr'    'posterior'    'Pz' 'POz'
alpha_soi = 1;
alpha_var = [];
avgalphapow = []; gammapow = [];
for isub = 1:15 % 1 :nsub %   % 1:nsub  %  
    for ises = 4 %1:3
        fprintf('SUBJ %s %d\n', respavg.SUBJ{isub}, ises)
        for im = [1 8] %1:length(meas_leg)  % 1:5 %2:5  % ssvep response %2:5
            if im < 7; iband = 1; else iband = 2; end
            
            for icond = 1:2%:3
                if ises < 4
                    %                     pre_alpha_dat = respavg.pow_singletrial{isub,1}{1,ises,icond,istim,iresp}(isoi,:);
                    %                     postdat = respavg.pow_singletrial{isub,iband}{im,ises,icond,istim,iresp};
                else % concatenate all trials across sessions
                    pre_alpha_dat = cat(1, respavg.pow_singletrial{isub,1}{1,1:3,icond,istim,iresp});
                    postdat = cat(1, respavg.pow_singletrial{isub,iband}{im,1:3,icond,istim,iresp});
                end
                
                if im == 10 % SDT
                    dat = [pre_alpha_dat(:,alpha_soi) postdat(:,1:2)]; % stim, resp cols
                else
                    dat = [pre_alpha_dat(:,alpha_soi) postdat(:,isoi)];
                end
                
                if isempty(dat)
                    fprintf('dat not found\n')
                    continue;
                end
                dat = dat(~isnan(dat(:,1)),:); % remove nans 
                
                ntrials(isub,icond) = length(dat(:,1));
                if im == 1 % for corr with crit
%                     avgalphapow(isub, icond) = mean(dat(:,1));
                    avgalphapow(isub, icond) = mean(dat(:,1));
                end
                
                if im == 10 %SDT
                    dat(:,1) = log(dat(:,1)); % take log only of alpha pow
                elseif im == 1 % alpha
                    dat(:,1) = log(dat(:,1)); % take log only of alpha pow: col1 used for binning
                    dat(:,2) = log(dat(:,2)); % take log only of tobeplotted alpha pow in col 2
                else
                    %                     dat = log(dat); % take log
                    %                     dat = log(dat(:,2)); % take log
                    dat(:,1) = log(dat(:,1)); % take log only of alpha pow: col1 used for binning
                    dat(:,2) = log(dat(:,2)); % take log only of tobeplotted alpha pow in col 2

                    %                     dat(:,2) = log(dat(:,2)); % take log only of gamma pow
                end
                if im == 8
                    gammapow(isub, icond) = mean(dat(:,2));
%                     dat(:,2) = (dat(:,2) - gammapow(isub, 1)) ./ gammapow(isub, 1) * 100;
                end

%                 %figure; scatter(log(dat(:,1)), dat(:,2))
%                 if im == 8 % fit single trials
%                     pf(isub,icond,:) = polyfit(dat(:,1), dat(:,2), 2);
%                     % y = polyval(p, 1:nbins); 
%                 end
                    
                alpha_var(isub, icond) = var(dat(:,1));
                
                % control bin_ranges, e.g. to get 30 % of the observations in each bin, cf Rajagovindan etal
                data_range = linspace(min(dat(:,1)), max(dat(:,1)), 100); 
                for ib = 1:size(bin_ranges,1)
%                     bin_edges = prctile(dat(:,1), bin_ranges(ib,:)); % use  percentiles itself
                    bin_edges = prctile(data_range, bin_ranges(ib,:)); % use data range, not percentiles itself
                    bin_ind = discretize(dat(:,1), bin_edges, 1);
                    
                    if im == 10 % SDT
                        ydat = dat(bin_ind == 1, 2:3);
                        hitrate = length(find(ydat(:,1)==1 & ydat(:,2)==1)) / length(find(ydat(:,1)==1)) ;
                        farate = length(find(ydat(:,1)==2 & ydat(:,2)==1)) / length(find(ydat(:,1)==1));
                        if hitrate == 0; hitrate = 0.0001; end
                        if hitrate == 1; hitrate = 0.9999; end
                        if farate == 0; farate = 0.0001; end
                        if farate == 1; farate = 0.9999; end
                        % criterion:
                        dat_bins(isub, ises, im, ib, icond) =  -0.5 * (norminv(hitrate) + norminv(farate));
                        % dprime:
%                         dat_bins(isub, ises, im, ib, icond) =  norminv(hitrate) - norminv(farate);
                    else
                        ydat = dat(bin_ind == 1, 2);
                        dat_bins(isub, ises, im, ib, icond) = mean(ydat);
                        dat_bin_count(isub, ises, im, ib, icond) = length(ydat);
                        dat_bin_edges(isub, ises, im, ib, icond,:) = bin_edges;
                    end
                end
                
            end
        end
    end
end

% % add ddm dc 10 bins
if nbins == 10
    dat_bins(:, 4, 11, 1:10, 1:2) = respavg.ddmpars10bins.dc([1:13, 15:end],:,: ); % remove 8 and 14
elseif nbins == 5
    % add ddm dc 5
    dat_bins(:, 4, 11, 1:5, 1:2) = respavg.ddmpars5bins.dc([1:13, 15:end],:,: ); % remove 14
elseif nbins == 3
    dat_bins(:, 4, 11, 1:3, 1:2) = respavg.ddmpars3bins.dc([1:13, 15:end],:,: ); % remove 14
end

dat_bins = dat_bins - min(dat_bins(:,:,:,:,1),[],4); % take min for cons to normalize

% polyfit 2nd order on single subjects, get quadratic fit (gain)
for isub = 1:15
    for im = 8% [1,8,9]
        for icond = 1:2
            bin_fits(isub,:,im,icond) = polyfit(transpose(1:nbins), squeeze(dat_bins(isub, ises, im, :, icond)), 2);
        end
    end
end
bin_fits = bin_fits([1:7, 9:15],:,:,:,:);

% remove subj 8, very low gamma
dat_bins = dat_bins([1:7, 9:15],:,:,:,:);


%% Plotting bins

SAV = 1;
nbins = size(bin_ranges,1);
% figure; plot(squeeze(mean(dat_bins)))
% close all
% linecol = {'r' 'b' 'k'};
linecol = cbrewer('qual', 'Set1',3);
linecol(4,:) = [0 0 0];
set(0, 'Defaultaxesfontsize', 12)
f = figure; hold on; iplot =0;
% f.Position = [ 1989         642         600         300];  %  567         314        1354         628
f.Position = [ 567         314        1354         628];
keepdat = []; % for ratio lib vs cons
for im =  [1,8] % [1,8, 11] % [1,8,9] % [1:4, 6:9] 1:8% 1:length(meas_leg)  % 1:5 %2:5  % ssvep response %2:5
    iplot=iplot+1;
    subplot(1,3,iplot); hold on; axis square; box on
    xlim([0 nbins-1])
    if iplot > 1
        r = refline(0,0);
        r.LineStyle = '--';
        r.Color = 'k';
    end
    h = []; dat = []; clear s;
    for icond = 1:2%:2
        if icond<3 || icond == 4
%             % within subj error bars
            temp = squeeze(dat_bins(:,4,im,:,1:2)); % -  squeeze(mean(dat_bins(:,4,im,:,icond),4))
            condavg = squeeze(nanmean(temp,3));
            grandavg = nanmean(condavg);
            dat = temp(:,:,icond) - condavg + grandavg;
            
%             dat(:,icond) = squeeze(mean(dat_bins(:,4,im,:,icond))); %             actual data
%             sem = squeeze(std(dat_bins(:,4,im,1:nbins,icond))) / sqrt(size(dat_bins(:,4,im,:,icond),1));
%             dat(:,icond) = flipud(dat(:,icond));
%             s(icond) = shadedErrorBar(0:nbins-1, dat(1:nbins,icond), sem, {'Color', linecol(icond,:)}, 1);

            dat = fliplr(dat);
            sem = nanstd(dat,0,1) / sqrt(size(dat,1));
            s(icond) = shadedErrorBar(0:nbins-1, nanmean(dat(:,1:nbins)), sem, {'Color', linecol(icond,:)}, 1);
            scatter(0:nbins-1, mean(dat(:,1:nbins)), 'filled')
            keepdat(:,icond) = mean(dat);
        else
            dat(:,icond) = dat(:,2) ./ dat(:,1);
            scatter(1:nbins, dat(:,icond))
        end
    end
   
    ax=gca;
    ax.XLim = [0 nbins-1];
    ax.XTick = [0 (nbins-1)/2 nbins-1];
    ax.XTickLabel = {'Low' 'Medium' 'High'};
    
    legend([s.mainLine], respavg.behav_conds(1:2)); legend boxoff
    xlabel('Neural excitability')
end

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prealpha');
    
    mkdir(outpath)
    %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    outfile = fullfile(outpath, sprintf('alphagroupfits_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim,iresp}) );
    disp(outfile)
%     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    print(outfile, '-dpdf') %'-png',  '-pdf',
    print(outfile, '-dpng') %'-png',  '-pdf',
    print(outfile, '-depsc2') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% Average zscored dc and gamma over subj, corr bins
close all
SAV=1;
f = figure;
f.Position = [2080         438 300 150];
% f.Position = [2080         438 1000 375];
corrkeep = [];
% corrtype = 'Spearman';
corrtype = 'Pearson';
iplot=0;
correlationdata = [];
for icond = 2:-1:1 % [2,1,4] %
    iplot = iplot+1;
    ctr=0;
    corrdat = [];
    for im = [8 11]
        ctr=ctr+1;
        if icond < 3
            temp = squeeze(dat_bins(:,4,im,:,icond)); % subj gamma/dc bins
            temp = (temp - nanmean(temp,2)) ./ nanstd(temp,0,2); % zscore
        else
%             %psc wrt cons
%             temp = (squeeze(dat_bins(:,4,im,:,2)) - squeeze(dat_bins(:,4,im,:,1))) ./ ...
%                 squeeze(dat_bins(:,4,im,:,2)) * 100; % subj gamma/dc bins
            
            temp = squeeze(dat_bins(:,4,im,:,2)) - squeeze(dat_bins(:,4,im,:,1)); % subj gamma/dc bins
%             temp = (temp - nanmean(temp,2)) ./ nanstd(temp,0,2); % zscore

%             % zscore conds separately
%             temp2 = squeeze(dat_bins(:,4,im,:,2)); % subj gamma/dc bins
%             temp2 = (temp2 - nanmean(temp2,2)) ./ nanstd(temp2,0,2); % zscore
%             temp1 = squeeze(dat_bins(:,4,im,:,1));
%             temp1 = (temp1 - nanmean(temp1,2)) ./ nanstd(temp1,0,2); % zscore
%             temp = temp2-temp1;
        end
        
        
        corrdat(:,:,ctr) = fliplr(temp);
        %     corrdat = [corrdat temp(:)]
    end
    correlationdata.dat(:,:,:,icond) = corrdat;
    corravg = squeeze(nanmean(corrdat));
    corrkeep(:,:,icond) = corravg;
    xsem = nanstd(corrdat(:,:,1)) ./ sqrt(14) /2;
    ysem = nanstd(corrdat(:,:,2)) ./ sqrt(14) /2;
    
    s = subplot(1,2,iplot); hold on
    box on; axis square;

    e = errorbar(corravg(:,1),corravg(:,2), -ysem, ysem, -xsem, xsem, '.');
    e.Color = linecol(icond,:); % same color for error bars
%     e.Color = 'k'; % same color for error bars
       
    [r,p] = corr(corravg(:,1),corravg(:,2), 'type', corrtype );
    % do regression
    [b,bint,res,rint,stats] = regress(corravg(:,2),[corravg(:,1) ones(size(corravg(:,1)))] )
    %     text(-0.6, 0.1, sprintf('R^2 = %1.2f\np = %1.4f', stats(1), stats(3)))
    
    scatter(corravg(:,1),corravg(:,2), 'filled', 'Markerfacecolor', linecol(icond,:),  ...
        'Markeredgecolor', [1 1 1], 'sizedata', 50, 'linewidth', 0.5) ;
    
%     title(sprintf('%s r = %g, p = %g', respavg.behav_conds{icond}, r, p))
    title(sprintf('%s R^2 = %1.2f\np = %1.4f', respavg.behav_conds{icond}, stats(1), stats(3)))
    
    xlabel('Gamma activity (Z-score)')
    ylabel('Drift bias (Z-score)')
%     xlabel('Gamma raw')
%     ylabel('Drift criterion raw')

%     s.XTick = -1:0.5:1;
%     s.YTick = -1:0.5:1;
%     s.XLim = [-1 0.75];
%     s.YLim = [-1 0.75];

% cond 4:
%     s.XLim = [-0.1 0.1];
%     s.YLim = [1 2.5];

    %     s.XLim = [-0.75 0.75];
%     s.YLim = [-0.75 0.75];
%     if iplot == 2
%         s.XLim = [-0.75 0.5];
%         s.YLim = [-0.75 0.75];        
%     elseif iplot == 1
%         s.XLim = [-0.825 0.5];
%         s.YLim = [-0.75 0.5];
%     end
    if strcmp(corrtype, 'Pearson')
        if p < 0.05
            l = lsline;
            l.Color = 'k';
            l.LineWidth = 2;
        end
    end
    ax=gca;
    ax.FontSize = 8;
%     if iplot==2
%         c = colorbar;
%         c.Ticks = [1 5.5 10];
%         c.TickLabels = {'Low' 'Medium' 'High'};
%         c.Label.String = 'Neural excitability';
%         %         c.Location = 'East'
%         %         c.Position(4) = c.Position(4) * 0.5
%         c.Position = [ 0.91    0.3281    0.0200    0.3760] ;
%         c.Label.FontSize=12;
%     end
end
%     out = randtest_corr( corrkeep(:,:,2), corrkeep(:,:,1) ,1, 10000,  corrtype)

correlationdata.dimord = 'subj_bins_gammadrift_conslib'

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prealpha');
    
    mkdir(outpath)
    %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    outfile = fullfile(outpath, sprintf('gammavsdcsubjavg_%s_%s_%dbins', respavg.sens.leg{isoi}, respavg.sdt_conds{istim,iresp}, nbins) );
    disp(outfile)
%     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    print(outfile, '-dpdf', '-fillpage') %'-png',  '-pdf',
    print(outfile, '-dpng') %'-png',  '-pdf',
    print(outfile, '-depsc2') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises
