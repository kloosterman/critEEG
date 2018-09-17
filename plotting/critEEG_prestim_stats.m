% load in basespec and plot

% respavg = critEEG_load_respavg()

%%
% nsub = length(respavg.SUBJ);
% sub_ind = [1:15];
nsub = 14;
sub_ind = [1:7, 9:15];

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
% cfg.statistic        = 'correlationT';
cfg.correctm         = 'cluster';
% cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05; % for evoked data
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;
% % cfg.neighbours       = []; %in case no channel data present
% if ismac
%     cfg0_neighb.template  = 'elec1010_neighb.mat';
%     cfg0_neighb.elecfile  = 'standard_1020.elc';
% else
%     cfg0_neighb.template  = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20161220/template/neighbours/elec1010_neighb.mat';
%     cfg0_neighb.elecfile  = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20161220/template/electrode/standard_1020.elc';
% end
% cfg0_neighb.feedback = 'no';
% cfg0_neighb.method = 'template';
% cfg.neighbours       = ft_prepare_neighbours(cfg0_neighb);
cfg.neighbours       = []

design = zeros(2,2*nsub);
for i = 1:nsub
    design(1,i) = i;
end
for i = 1:nsub
    design(1,nsub+i) = i;
end
design(2,1:nsub)        = 1;
design(2,nsub+1:2*nsub) = 2;
cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

cfg.minnbchan = 0;

%%% correlation design 
% design = respavg.behavior.criterion(:,4,4)';
% cfg.design   = design;
% cfg.uvar     = [];
% cfg.ivar     = 1;
% cfg.type = 'Spearman';

cfg.latency = [-0.8 -0.2]; % -0.8 to 0 : p=0.02

% cfg.latency = [-0.8 -0.25];

% % % just occ
% cfg.channel = respavg.sens.ind{1};
% cfg.avgoverchan = 'yes';

% just occ, spectrum
cfg.channel = respavg.sens.ind{1};
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
% 
% %just for alpha, space by time
% cfg.frequency = [8 13];
% cfg.avgoverfreq = 'yes';

% % just for alpha occipital, time
% cfg.frequency = [9 13];
% cfg.avgoverfreq = 'yes';
% % cfg.channel = respavg.sens.ind{1};  % Pz, PO3, O1, Oz, O2, PO4, P1, P2, POz)
% cfg.channel = respavg.sens.ind{1};
% cfg.avgoverchan = 'yes';
% cfg.latency = [-0.6 0];
% cfg.avgovertime = 'yes';


ises=4;
isoi=1;
poolstat = []; stat = [];
for icond = 4 %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
    for istim = 3%:2%:3
        for iresp = 3%:2
            for iband = 1%:2
                for itrig = 1
                    
                    freq = [];
                    freq.dimord = 'subj_chan_freq_time';
                    freq.label = respavg.label;
                    freq.freq = respavg.freq{iband};
                    freq.time = respavg.time{itrig};
                    
                    freq.powspctrm  = squeeze(respavg.pow(sub_ind,:, 1:length(respavg.freq{iband}), 1:length(freq.time),   iband, itrig,     ises,icond,istim,iresp ));
%                     freq.powspctrm = log(squeeze(respavg.pow(sub_ind,:, 1:length(respavg.freq{iband}), 1:length(freq.time),   iband, itrig,     ises,2,istim,iresp ))) - ...
%                         log(squeeze(respavg.pow(:,:, 1:length(respavg.freq{iband}), 1:length(freq.time),   iband, itrig,     ises,1,istim,iresp )));
%                     dum = squeeze(mean(respavg.pow(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
%                         iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                    
                    %                 freq.powspctrm  = squeeze(mean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                    %                     iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                    
                    freqzero = freq; %create zero freq to test against
                    freqzero.powspctrm = zeros(size(freq.powspctrm));
                    poolstat{ iband, icond, istim, iresp } = ft_freqstatistics(cfg, freq, freqzero);

                    %                     poolstat{ iband, icond, istim, iresp } = ft_freqstatistics(cfg, freq);
                    
                    stat{iband} = poolstat{ iband, icond, istim, iresp };
                    
                    tind = respavg.time{itrig} > cfg.latency(1) & respavg.time{itrig} < cfg.latency(2);
                    stat{iband}.powspctrm = mean(freq.powspctrm(:,:,:,tind),4);
%                     stat.powspctrm = squeeze(mean(stat.powspctrm,2))
                end
            end
        end
    end
end
stat{iband}.negclusters(1).prob

%% plot significant cluster in freq (space)
iband = 1;

% chans = any(stat.mask,2);
chans = respavg.sens.ind{1};
% dum = squeeze(trapz(stat.powspctrm(:,chans,:),2));
dum = squeeze(mean(stat{iband}.powspctrm(:,chans,:),2));
% dum = dum([1:7, 9:15],:)
% close all
f = figure; hold on
f.Position = [680 873 345 225];
h =shadedErrorBar(stat{iband}.freq, mean(dum ), std(dum)/sqrt(size(dum,1)), 'k', 1 );
h.mainLine.LineWidth = 2;
r= refline(0,0);
r.Color = 'k';
r.LineStyle = '--';
% plos sig bar
sigfreqs = any(stat{iband}.mask,1);
if any(sigfreqs)
    barind = find(diff(sigfreqs));
    barind(1) = barind(1)+1
    if mod(length(barind),2), barind = [barind length(any(mask,1))]; end
    % ylim([-2e-8 0.5e-8])
    % ylim([-15 10])
    ax=gca;
    ax.FontSize = 12;
    % ypos = -12.5;% diff(ax.YLim)*-0.1;
    ypos = diff(ax.YLim)*-0.1;
    pl = plot( stat{iband}.freq(barind), [ypos ypos] );
    pl.LineWidth = 10;
    pl.Color = 'k';
    xlabel('Frequency (Hz)')
    ylabel('Lib - Cons power')
    xlim([0 35])
end

% topo's for cons, lib en lib-cons
f = figure; hold on
f.Position = [        1026         873         345         225];
cfg = [];
cfg.layout = 'elec1010.lay';
cfg.comment = 'no';
cfg.marker = 'on';
cfg.shading = 'flat';
cfg.style = 'straight'; %both  straight
cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
cfg.markersize = 4;
cfg.highlightchannel = chans; %respavg.sens.ind{1};
cfg.highlightsymbol = '.';
cfg.highlightsize = 20;
cfg.highlightcolor = 'w';
load colormap_jetlightgray.mat
cfg.colormap = cmap;
cfg.highlight = 'on';

%topo
% cfg.xlim = times(itop,:);
% cfg.ylim = freqs(itop,:);
clear dum
freq=[];
freq.label = respavg.label;
freq.dimord = 'chan_freq';
freq.freq = mean(stat{iband}.freq(sigfreqs));
dum = squeeze(mean(stat{iband}.powspctrm(:,:,sigfreqs),3));
% dum = squeeze(trapz(stat{iband}.powspctrm(:,:,freqs),3));

freq.powspctrm = mean(dum)';
% cfg.zlim = [-10 10];
                    cfg.zlim = 'maxabs';
%                     cfg.zlim = scales(itop,:);
ft_topoplotTFR(cfg, freq);
colorbar
ax=gca;
ax.FontSize = 12;


%% correlate alpha with  crit dprime db drift
close all
SAV = 1;
ises=4; isoi=1; icond = 1:2; istim = 3; iresp = 3;  itrig = 1;

incsubj = [1:7, 9:15]; % 14 already dropped, only drop 8 here

freq = [];
freq.dimord = 'subj_chan_freq_time';
freq.label = respavg.label;
iband = 2;
freq.freq = respavg.freq{iband};
freq.time = respavg.time{itrig};
cfg=[];
% cfg.freq = [8 12]; cfg.avgoverfreq = 'yes';
cfg.freq = [59 100]; cfg.avgoverfreq = 'yes';
% cfg.latency = [-0.8 -0.2]; cfg.avgovertime = 'yes';
cfg.latency = [0.2 0.6]; cfg.avgovertime = 'yes';
cfg.channel = respavg.sens.ind{isoi}; cfg.avgoverchan = 'yes';

alphapowdat=[];
for icond = 1:2
    freq.powspctrm  = squeeze(respavg.pow(incsubj,:, 1:length(respavg.freq{iband}), 1:length(freq.time),   iband, itrig,     ises,icond,istim,iresp ));
%     freq.powspctrm  = squeeze(respavg.dat(incsubj,:, 1:length(respavg.freq{iband}), 1:length(freq.time),   iband, itrig,     ises,icond,istim,iresp ));
    freqout = ft_selectdata(cfg, freq);
    alphapowdat(:,icond) = freqout.powspctrm;
end

% %%
% incsubj = [1:7, 9:15]; % 14 already dropped, only drop 8 here
% dum2 = squeeze(mean(respavg.pow(incsubj,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
%     iband, itrig, ises,1:2, istim, iresp), 2)); %average over sens
% dum2 = mean(dum2(:,sigfreqs,:,:),2);
% alphapowdat = squeeze(mean(dum2(:,:,tind,:),3));

% alphapowdat = log(alphapowdat);
alphapowdat(:,3) = alphapowdat(:,2) - alphapowdat(:,1);

incsubj = [1:7, 9:13, 15:16]; % drop 8 and 14

% load in DDM csv
PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/critEEG_analysis/ddm/';
PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/plots/';
% csvin = 'results_0_-016s.csv';  % 250 ms cutoff, counted from target onset
csvin = 'results_0_200ms_cutoff.csv';
ddmtemp = csvread(fullfile(PREIN, csvin), 1, 0);
ddmtemp = ddmtemp(:,[2:7, 9:11]);

behavdat =[]; % crit dprime db drift
behavdat(:,1:2,1) = squeeze(respavg.behavior.criterion(incsubj,4,1:2));
behavdat(:,1:2,2) = squeeze(respavg.behavior.dprime(incsubj,4,1:2));
behavdat(:,1:2,3) = [ddmtemp(incsubj,8) ddmtemp(incsubj,9)];
behavdat(:,1:2,4) = [ddmtemp(incsubj,3) ddmtemp(incsubj,4)];
behavdat(:,3,:) = behavdat(:,2,:) - behavdat(:,1,:) ; % lib - cons

measleg = { 'SDT criterion' 'SDT dprime' 'DDM drift bias' 'DDM drift rate'};

linecol(3,:) = [0 0 0];
f = figure; iplot=0;
f.Position =[   680   278   900 600   ];
for icond = 1:3
    for im = 1:4 % crit dprime db drift
        iplot=iplot+1;
        subplot(3,4,iplot);     hold on; box on; axis square
        
        scatter(alphapowdat(:,icond), behavdat(:,icond,im), ...
            'MarkerEdgecolor', [1 1 1], 'MarkerFacecolor', linecol(icond,:), 'Linewidth', 1, 'Sizedata', 100)
        [r,p] = corr(alphapowdat(:,icond), behavdat(:,icond,im), 'type', 'Pearson'); % Pearson  Spearman
        title(sprintf('r = %1.2f, p = %1.2f', r,p))
%         if  im == 1
            ylabel(measleg{im});
%         end
        if icond == 3
            xlabel('8-12 Hz power (log)')
%             xlabel('59 - 100 Hz power (log)')
        end
        lsline
        
    end
end

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prealpha');
    
    mkdir(outpath)
    %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    outfile = fullfile(outpath, sprintf('acrosssubcorr%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim,iresp} ));
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





%%



figure; 
for ises = 1:3
    subplot(1,3,ises)
    hold on; box on; axis square
    temp = [behavdat(:,ises), alphapowdat(:,ises)];
    [r,p] = corr(temp(~isnan(temp(:,1)),1), temp(~isnan(temp(:,1)),2), 'type', 'spearman');
    scatter(temp(~isnan(temp(:,1)),1), temp(~isnan(temp(:,1)),2), 400, ...
        'MarkerEdgecolor', [1 1 1], 'MarkerFacecolor', [0 0 1/ises], 'Linewidth', 3)
    % text(behavdat, alphapowdat, respavg.SUBJ );
    title(sprintf('r = %g, p = %g', r,p))

xlabel('criterion change lib vs cons (%)')
% xlabel('dprime change lib vs cons (%)')
ylabel('alpha power change lib vs cons (%)')
lsline

end


%% make 2 bins obv crit shift
lowgroup = behavdat < median(behavdat);
highgroup = behavdat > median(behavdat);
permtest(alphapowdat(lowgroup), alphapowdat(highgroup))
figure; bar([mean(alphapowdat(lowgroup)) mean(alphapowdat(highgroup))])


%% corr dc and crit

% incsubj = [1:13,15,16];
incsubj = [1:16];
behavdat = (respavg.behavior.criterion(incsubj,4,2) - respavg.behavior.criterion(incsubj,4,1)); % ./ ...
%    respavg.behavior.criterion(incsubj,4,1) * 100; %(isub, ises, icond) 
dcdat = (ddmdat(incsubj,2,4) - ddmdat(incsubj,1,4));
[r,p] = corr(behavdat, dcdat)
figure; hold on; box on; axis square
scatter(behavdat, dcdat, 400, ...
    'MarkerEdgecolor', [1 1 1], 'MarkerFacecolor', [0 0 0], 'Linewidth', 3)
% text(behavdat, dcdat, respavg.SUBJ );

ax=gca
ax.FontSize=20;
title(sprintf('r = %g, p = %g', r,p))

xlabel('criterion lib - cons')
% xlabel('dprime change lib vs cons (%)')
ylabel('DDM dc lib - cons')
lsline

SAV=1;
if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'ddm_vs_dc');
    
    mkdir(fullfile(outpath))
    %                 if showstats
    %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
    %                 else
    outfile = fullfile(outpath, 'ddm_vs_dc');
    %                 end
    disp(outfile)
    export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% make 3 bins obv crit shift
% alldat = [behavdat alphapowdat]
[~,ind] = sort(behavdat )

lowgroup = alphapowdat(ind(1:5));
mediumgroup = alphapowdat(ind(6:10));
highgroup = alphapowdat(ind(11:15));

% permtest(, alphapowdat(highgroup))
figure; bar([mean(lowgroup) mean(mediumgroup) mean(highgroup)])
ax=gca
ax.FontSize=20;


%% Plot TFR for occipital pooling

close all
SAV = 1;

% addpath(genpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/plotting'));
showstats = 1;
load colormap_jetlightgray.mat
% cmap = flipud(cmap) % !!!! flip colorbar to get liberal - cons
ncols = 1;
ises = 4;
for istim = 3%:2
    for iresp = 3%:2        %1:2;
        for isoi = 1 %1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
            
            f = figure;    iplot=0; hold on
            f.Position = [2029         314         360         500]; %
            set(0,'DefaultAxesFontsize',30)
            pl = [];
            for iband = 1 %2:-1:1 % 1:2 %
                for icond = 4 %[1:2, 4] %1:4 %1:4%  % 1:2
                    for itrig = 1%:2
                        
                        ZLIMS = [-30 30; -10 10];
                        if icond == 4; ZLIMS = ZLIMS/2; end
%                         ZLIM = ZLIMS(iband,:);
                        
                        iplot = iplot+1;
                        pl = subplot(2,ncols,iplot); hold on
                        pl.Units = 'pixels';
                        if iplot == 1; yl = ylabel('Frequency (Hz)'); end
                        colormap(cmap)
                        
%                         dum = squeeze(mean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
%                             iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        dum = squeeze(mean(respavg.pow(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:17, ... ipv length(respavg.time{itrig}), ...
                            iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        dumsubj = dum; % keep singlesub
                        dum = squeeze(mean(dumsubj));
%                         ZLIM = [min(min(dum)) -min(min(dum))];
%                         ZLIM = [-1.5e-08 1.5e-08];
%                         ZLIM = [-1e-08 1e-08]; % raw
                        ZLIM = [-0.15 0.15];

                        if showstats
%                             mask = double(squeeze( poolstat{ iband, icond, istim, iresp }.mask));
                            mask = squeeze(double(stat.mask));
                            mask(mask==0) = 0.25;
                            %                             mask = double(squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.prob < 0.15));
                            %                             mask(mask==0) = 0.25;
                            
                            ft_plot_matrix(respavg.time{itrig}(1:17), respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no', ...
                                'highlight', mask, 'highlightstyle', 'opacity'); % opacity
                        else
                            ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no' ); % opacity
%                             ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum, 'box', 'no' ); % opacity
                        end
                        
%                         XLIM = [respavg.time{itrig}(1) respavg.time{itrig}(end)];
                        XLIM = [-0.8 0];
                        YLIM = [respavg.freq{iband}(1) respavg.freq{iband}(end)];
                        xlim(XLIM)
                        ylim(YLIM)
                        
%                         binwidth = 0.01; %on screen
%                         %                         pl.Position(3) = pl.Position(3) * diff(XLIM);
%                         pl.Position(3) = length(respavg.time{itrig}) * binwidth;
%                         
                        ax = gca;
                        ax.FontSize = 12;
                        ax.Box = 'on';
                        ax.XTick = -1:0.25:2;
                        if iband==1; ax.YTick = [0:5:200]; else ax.YTick = [0:10:200]; end
                        
                        if iplot == 1
                            title( sprintf('%s, %s, psc Lo: %g, Hi: %g\n%s', respavg.sdt_conds{istim, iresp}, respavg.sens.leg{isoi}, ZLIMS(:,2), respavg.behav_conds{icond} ))
%                             pl.XTickLabel = [];
                        else
%                             pl.Position(2) = pl.Position(2)*1.95;
%                             pl.Position(4) = pl.Position(4)*1.5;
%                             xlabel(sprintf('Time from %s (s)', respavg.trigger_leg{itrig} ));
                            xlabel('Time from sequence onset (s)');

%                             title( sprintf('%s, [%g-%g%%]', respavg.behav_conds{icond}, ZLIM*100))
                        end
                        
                        yaxis = [respavg.freq{iband}(1), respavg.freq{iband}(end)];
                        plot([0,0],yaxis,'k',[0,0],yaxis,'k');
                        if strcmp(respavg.trigger_leg{itrig}, 'stim')
                            stimonset = 0.16;
                            plot([stimonset,stimonset],yaxis,'--k',[stimonset,stimonset],yaxis,'--k');
                            plot([1,1],yaxis,'--w',[1,1],yaxis,'--w', 'Linewidth', 2);

                            RT = mean(respavg.behavior.RT(:, ises, istim, icond)); %RT
                            if icond<4
                                plot(RT, respavg.freq{iband}(end-1), 'vk', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
                            end
                        end
                        h=colorbar;
%                         h.Position(3) = h.Position(3) * 0.75;
%                         h.Position(4) = h.Position(4) * 0.5;
% %                         h.Position(1) = 0.95
% %                         h.Position(2) = 0.5
%                         h.Ticks = [h.Limits(1) 0 h.Limits(2)];
%                         % Specify pos as a four-element vector of the form [x y w h] in data units.
                        if iplot == 1
                            boxlocs = [-0.35 8 0.35 5; ...
                                       0.2 42 0.2 16];
                        else
                            boxlocs = [0.25 11 0.2 11; ...
                                       0.25 24 0.2 3];
                        end
                        for ib = 1%:2
%                             if ib < 3; subplot(1,2,1); else subplot(1,2,2); end
%                             r=rectangle('Position', boxlocs(ib,:));
%                             r.LineStyle = '--';
                        end

                    end %ises
                end
            end
            set(gcf, 'Color', 'white')
            
            if SAV
                %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
                outpath = fullfile(respavg.PREOUT, 'poolings');
                
                mkdir(fullfile(outpath, 'stats'))
                %                 if showstats
                %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
                %                 else
                outfile = fullfile(outpath, sprintf( 'TFRrawpow_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
                %                 end
                disp(outfile)
                export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
                export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
                cd(outpath)
            end %ises
        end
    end
end %ises



%% plot multiplot for significant clusters

close all

istim=3; iresp=3; 
for icond = 4% 1:4 %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
    for iband = 1%:2%:2
        for itrig = 1 %1:2
            freq=[];
% %             freq2stats.mask = poolstat{isoi, iband, itrig, icond, istim, iresp}.mask;
            freq.mask = stat.negclusterslabelmat == 1;
% %             freq.mask = poolstat{isoi, iband, itrig, icond, istim, iresp}.negclusterslabelmat == 2;
%             if ~any(freq.mask)
%                 continue;
%             end
%             freq.time = poolstat{iband, icond, istim, iresp}.time;
            freq.label = respavg.label;
            freq.dimord = 'chan_freq_time';
            freq.time = respavg.time{itrig};
            freq.freq = respavg.freq{iband};
%             freq.powspctrm  = squeeze(respavg.pow(:,:, 1:length(respavg.freq{iband}), 1:23, ...
%                 iband, itrig, 4,icond, istim, iresp));
            freq.powspctrm  = squeeze(mean(respavg.pow(:,:, 1:length(respavg.freq{iband}), :, ...
                iband, itrig, 4,icond, istim, iresp),1));
            
            cfg = [];
            cfg.xlim = [-0.8 0];
            cfg.shading = 'flat';
%             cfg.layout = 'biosemi64incI1I2.lay';  % biosemi64  elec1010
            % cfg.layout = 'biosemi64.lay';  % biosemi64  elec1010
            cfg.layout = 'elec1010.lay';  % biosemi64  elec1010
            cfg.layout = ft_prepare_layout(cfg);
            cfg.layout.width(:) = 0.075;
            cfg.layout.height(:) = 0.075;
            
%             cfg.zlim = [-0.1 0.1];
            cfg.zlim = 'maxabs';
            % cfg.zlim = [-1 1];
            load( 'colormap_jetlightgray.mat')
            cfg.colormap = cmap;
            cfg.maskparameter = 'mask';
            cfg.hotkeys = 'yes';
            cfg.fontsize = 18;
            cfg.colorbar = 'yes';
            
            f = figure;
            f.Position = [ 680   678   1200   1200];
            
            ft_multiplotTFR(cfg, freq);
            
            title( sprintf('%s, %s, %s, %s',  respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}, respavg.trigger_leg{itrig} ))
        end
    end
end


%% Plot spectra for occipital pooling

TIM = [-0.8 -0.2];
tind = respavg.time{itrig} > TIM(1) & respavg.time{itrig} < TIM(2)

% close all
SAV = 1;

% addpath(genpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/plotting'));
showstats = 0;
% load colormap_jetlightgray.mat
% cmap = flipud(cmap) % !!!! flip colorbar to get liberal - cons
ncols = 1;
ises = 4;
for istim = 3%:2
    for iresp = 3%:2        %1:2;
        for isoi = 1 %1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
            
            f = figure;    iplot=0; hold on
            f.Position = [ 2092         504         894         314]; %
            set(0,'DefaultAxesFontsize',12)
            pl = [];
            for iband = 1 %2:-1:1 % 1:2 %
                for icond = [1:2, 4] %1:4 %1:4%  % 1:2
                    for itrig = 1%:2
                        
                        ZLIMS = [-30 30; -10 10];
                        if icond == 4; ZLIMS = ZLIMS/2; end
%                         ZLIM = ZLIMS(iband,:);
                        
if icond==2; iplot=0; end
                        iplot = iplot+1;
                        pl = subplot(1,2,iplot); hold on
                        
                        xlabel('Frequency (Hz)');

                        colormap(cmap)
                        
%                         dum = squeeze(mean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
%                             iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        dum = squeeze(mean(respavg.pow(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                            iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        dum = mean(dum(:,:,tind),3);
%                         dum = log(dum);
                        
                        %                       close all
                        %                         f = figure; hold on
                        %                         f.Position =[        2001         172        1202         896];
                        if iplot <2
                            pl = plot(respavg.freq{1}, mean(dum))
                            pl.LineWidth = 2;
                            legend(respavg.behav_conds(1:2)); legend boxoff
                            ylabel(sprintf('EEG scalp current density (V/m^2)'))
                        else
                            sh =shadedErrorBar(respavg.freq{1}, mean(dum), std(dum)/sqrt(nsub), 'k', 0)
                            sh.mainLine.LineWidth=2;
                            r = refline(0,0);
                            r.Color = 'k';
                            r.LineStyle = '--';
                            ylabel('Power modulation (%)')
%                             legend(sh.mainLine, 'Liberal vs. conservative'); legend boxoff
                            title('Liberal vs. conservative')
                            % plos sig bar
                            sigfreqs = any(stat.mask,1);
                            barind = find(diff(sigfreqs));
                            barind(1) = barind(1)+1;
                            if mod(length(barind),2), barind = [barind length(any(mask,1))]; end
                            ax=gca;
                            %                             ypos = -13;% diff(ax.YLim)*-0.1;
                            ypos = diff(ax.YLim)*-0.1;

                            pl = plot( stat.freq(barind), [ypos ypos] );
                            pl.LineWidth = 10;
                            pl.Color = 'k';
                            text(14, ypos, 'p < 0.05, corrected')

                        end
                        xlim([0 35])
                        ax=gca;
                        ax.XTick = 0:5:35;
%                         ax.YLim = [-0.3e-7 0.11e-7];
                        if iplot == 2, ax.YTick = -0.3e-7:0.1e-7:1e-7; end
                        % plot single subj
%                         for isub = 1:15 
%                             pl = plot(respavg.freq{iband},   dum(isub,:)  )
%                             pl.LineWidth = 2;
%                              pl.LineStyle = '-';
%                             if isub > 7
%                                 pl.LineStyle = '--';
%                             end
%                         end
%                         legend(respavg.SUBJ(1:15))

                    end %ises
                end
            end
            set(gcf, 'Color', 'white')
            
            if SAV
                %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
                outpath = fullfile(respavg.PREOUT, 'spectra');
                
                mkdir(fullfile(outpath))
                %                 if showstats
                %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
                %                 else
                outfile = fullfile(outpath, sprintf( 'spectrarawpow_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
                %                 end
                disp(outfile)
                export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
                print(outfile, '-dpdf') %'-png',  '-pdf',
                cd(outpath)
            end %ises
        end
    end
end %ises

%% plotting: critEEG TFR NoFA vs NoMiss: motor and visual cortex + stats
close all
SAV = 1;


showstats = 0;
%
% cmap = cbrewer('div', 'RdBu',256);
% cmap = flipud(cmap);
% load('colormap170613.mat');
% load('colormap_jetgray', 'cmap');
% cmap = jet(256);
load colormap_jetlightgray.mat

ncols = 1;
ises = 4;
for istim = 3 %1:2
    for iresp = 3 %1:2        %1:2;
        for isoi = 1 %1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
            
            f = figure;    iplot=0; hold on
            f.Position = [2029         314         360         500]; %
            set(0,'DefaultAxesFontsize',30)
            pl = [];
            for iband = 2:-1:1 % 1:2 %
                for icond = 4 %[1:2, 4] %1:4 %1:4%  % 1:2
                    for itrig = 1%:2
                        
                        ZLIMS = [-30 30; -10 10];
                        if icond == 4; ZLIMS = ZLIMS/2; end
                        ZLIM = ZLIMS(iband,:);

                        iplot = iplot+1;
                        pl = subplot(2,ncols,iplot); hold on
                        pl.Units = 'pixels';
                        if iplot == 1; yl = ylabel('Frequency (Hz)'); end
                        colormap(cmap)
                        
%                         dum = squeeze(mean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
%                             iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        dum = squeeze(mean(respavg.pow(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                            iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        dumsubj = dum; % keep singlesub
                        dum = squeeze(mean(dumsubj));
                        
                        %                         ZLIMS = [-30 30; -10 10];
                        %                         if icond == 4; ZLIMS = ZLIMS/2; end
                        %                         ZLIM = ZLIMS(iband,:);
%                         ZLIM = [-max(max(dum)) max(max(dum))];
%                         ZLIM = [min(min(dum)) -min(min(dum))];
%                         ZLIM = [-1e-8 1e-8];
                        ZLIM = [-0.15 0.15];
                        
                        if showstats
                            mask = double(squeeze( poolstat{isoi, iband, itrig, icond, istim, iresp}.mask));
                            mask(mask==0) = 0.25;
                            %                             mask = double(squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.prob < 0.15));
                            %                             mask(mask==0) = 0.25;
                            
                            ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no', ...
                                'highlight', mask, 'highlightstyle', 'opacity'); % opacity
                        else
                            ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no' ); % opacity
%                             ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum, 'box', 'no' ); % opacity
                        end
                        
%                         XLIM = [respavg.time{itrig}(1) respavg.time{itrig}(end)];
                        XLIM = [-0.8 0.3];
                        YLIM = [respavg.freq{iband}(1) respavg.freq{iband}(end)];
                        xlim(XLIM)
                        ylim(YLIM)
                        
%                         binwidth = 0.01; %on screen
%                         %                         pl.Position(3) = pl.Position(3) * diff(XLIM);
%                         pl.Position(3) = length(respavg.time{itrig}) * binwidth;
%                         
                        ax = gca;
                        ax.FontSize = 12;
                        ax.Box = 'on';
                        ax.XTick = -1:0.25:2;
                        if iband==1; ax.YTick = [0:5:200]; else ax.YTick = [0:10:200]; end
                        
                        if iplot == 1
                            title( sprintf('%s, %s, psc Lo: %g, Hi: %g\n%s', respavg.sdt_conds{istim, iresp}, respavg.sens.leg{isoi}, ZLIMS(:,2), respavg.behav_conds{icond} ))
%                             pl.XTickLabel = [];
                            xlabel('Time from stimulus train onset (s)');
                        else
                            pl.Position(2) = pl.Position(2)*1.95;
%                             pl.Position(4) = pl.Position(4)*1.5;
%                             xlabel(sprintf('Time from %s (s)', respavg.trigger_leg{itrig} ));
                            xlabel('Time from stimulus train onset (s)');

%                             title( sprintf('%s, [%g-%g%%]', respavg.behav_conds{icond}, ZLIM*100))
                        end
                        
                        yaxis = [respavg.freq{iband}(1), respavg.freq{iband}(end)];
                        plot([0,0],yaxis,'k',[0,0],yaxis,'k');
                        if strcmp(respavg.trigger_leg{itrig}, 'stim')
                            stimonset = 0.16;
                            plot([stimonset,stimonset],yaxis,'--k',[stimonset,stimonset],yaxis,'--k');
                            plot([1,1],yaxis,'--w',[1,1],yaxis,'--w', 'Linewidth', 2);

                            RT = mean(respavg.behavior.RT(:, ises, istim, icond)); %RT
                            if icond<4
                                plot(RT, respavg.freq{iband}(end-1), 'vk', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
                            end
                        end
                        h=colorbar;
                        h.Position(3) = h.Position(3) * 0.75;
                        h.Position(4) = h.Position(4) * 0.5;
                        h.Ticks = [h.Limits(1) 0 h.Limits(2)];
                        % Specify pos as a four-element vector of the form [x y w h] in data units.
                        if iplot == 1
                            boxlocs = [-0.4 8 0.35 5; ...
                                       0.2 42 0.2 16];
                        else
                            boxlocs = [0.25 11 0.2 11; ...
                                       0.25 24 0.2 3];
                        end
                        for ib = 1%:2
%                             if ib < 3; subplot(1,2,1); else subplot(1,2,2); end
                            r=rectangle('Position', boxlocs(ib,:));
                            r.LineStyle = '--';
                        end

                    end %ises
                end
            end
            set(gcf, 'Color', 'white')
            
            if SAV
                %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
                outpath = fullfile(respavg.PREOUT, 'prestim');
                
                mkdir(fullfile(outpath, 'stats'))
                %                 if showstats
                %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
                %                 else
                outfile = fullfile(outpath, sprintf( 'TFRstats_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
                %                 end
                disp(outfile)
                export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
                export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
                cd(outpath)
            end %ises
        end
    end
end %ises

%% plot topo's fpr tf roi
SAV=1;

close all
cfg = [];
cfg.layout = 'elec1010.lay';
% cfg.comment = 'no';
cfg.marker = 'on';
cfg.shading = 'flat';
cfg.style = 'straight'; %both  straight
cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
cfg.markersize = 4;
cfg.highlightchannel = respavg.sens.ind{1};
cfg.highlightsymbol = '.';
cfg.highlightcolor = 'w';
cfg.highlightsize = 30;
cfg.colormap = cmap;
cfg.highlight = 'on';
%                             boxlocs = [-0.35 8 0.35 5; ...
%                                        0.2 42 0.2 16];

% times = [-0.4 0];
times = [-0.8 -0.2];
freqs = [8 12];
topotitles = {'Alpha'};
% scales =  [-0.05 0.05]; %[-1e-8 1e-8];
% scales =  [-0.1 0.1]; %[-1e-8 1e-8];
scales =  [-2e-8 2e-8];

ncols = 2;
ises = 4;
iband = 1;
for istim = 3
    for iresp = 3        %1:2;
        
        f = figure;    iplot=0; hold on
        f.Position = [ 2063         121        600         600]; % 
        set(0,'DefaultAxesFontsize',30)
        pl = [];
        for icond = 3:4 %[1:2, 4] %1:4 %1:4%  % 1:2
            if icond==4, scales =  [-2e-8 2e-8]; end

            for itrig = 1%:2
                for itop = 1%:4
                    
                    ZLIMS = [-30 30; -10 10];
                    if icond == 4; ZLIMS = ZLIMS/2; end
                    ZLIM = ZLIMS(iband,:);
                    
                    iplot = iplot+1;
                    pl = subplot(2,ncols,iplot); hold on
                    %                     pl.Units = 'pixels';
                    if iplot == 1; yl = ylabel('Frequency (Hz)'); end
                    colormap(cmap)
                    
                    %topo
                    cfg.xlim = times(itop,:);
                    cfg.ylim = freqs(itop,:);
                    freq=[];
                    freq.label = respavg.label;
                    freq.dimord = 'chan_freq_time';
                    freq.time = respavg.time{itrig};
                    freq.freq = respavg.freq{iband};
%                     dum = squeeze(mean(respavg.pow(:, :, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
%                         iband, itrig, 4,icond, istim, iresp), 1)); %average over sens
                    dum = squeeze(respavg.pow(:, :, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                        iband, itrig, 4,icond, istim, iresp)); %average over sens
                    % dum = zscore(dum); % over subj

                    freq.powspctrm = squeeze(mean(dum));
                    if iplot < 4
                        cfg.zlim =  [-5e-7 5e-7];
                    else
                        cfg.zlim = [-15 15];
                    end
                    cfg.zlim = 'maxabs';
                    cfg.comment = 'no';
%                     cfg.zlim = scales(itop,:);
                    ft_topoplotTFR(cfg, freq);
                    c = colorbar;
%                     c.Ticks = [scales(itop,1) 0 scales(itop,2)];
                    c.Position(4) = c.Position(4)*0.5;
                    c.Position(1) = c.Position(1)*1.1;
                    c.Label.String = 'Modulation (%)';
                    ax=gca;
                    ax.FontSize = 12;
%                     title(topotitles{itop})
                    title( respavg.behav_conds{icond} )
                end
            end
        end
    end
end

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prestim');
    
    mkdir(fullfile(outpath, 'stats'))
    %                 if showstats
    %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
    %                 else
    outfile = fullfile(outpath, sprintf( 'topo_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    %                 end
    disp(outfile)
    export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% plot spectra all SDT lib-cons
close all

for icond = [1:2, 4] %1:4 %1:4%  % 1:2
    f = figure;    iplot=0; hold on
    f.Position = [ 2063         121        800         600]; %
    set(0,'DefaultAxesFontsize',10)
    pl = [];
    ncols=1;
    isoi=1;
    iplot = iplot+1;
    pl = subplot(2,ncols,iplot); hold on
    for iresp = 1:2        %1:2;
        for istim = 1:2

            for itrig = 1%:2
                
                ZLIMS = [-30 30; -10 10];
                if icond == 4; ZLIMS = ZLIMS/2; end
                ZLIM = ZLIMS(iband,:);
                
                %                     pl.Units = 'pixels';
                if iplot == 1; yl = xlabel('Frequency (Hz)'); end
                colormap(cmap)
                
                freq=[];
                %                     freq.label = respavg.label;
                %                     freq.dimord = 'chan_freq_time';
                %                     freq.time = respavg.time{itrig};
                %                     freq.freq = respavg.freq{iband};
                dum = squeeze(mean(respavg.pow(:, :, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                    iband, itrig, 4,icond, istim, iresp), 1)); %average over sens
                dum = squeeze(mean(dum(respavg.sens.ind{isoi},:,:))); % avg over chan
                tind = respavg.time{itrig} >= -0.4 & respavg.time{itrig} <= 0;
                dum = mean(dum(:,tind),2);
                
                pl=plot(respavg.freq{iband}, dum);
                pl.LineWidth = 2;
                legnames = respavg.sdt_conds(1:2,1:2);
                legend(legnames(:))

                
                
%                 cfg.zlim = scales(itop,:);
%                 ft_topoplotTFR(cfg, freq);
%                 c = colorbar;
%                 c.Ticks = [scales(itop,1) 0 scales(itop,2)];
%                 c.Position(4) = c.Position(4)*0.5;
%                 c.Position(1) = c.Position(1)*1.1;
%                 c.Label.String = 'Modulation (%)';
%                 ax=gca;
%                 ax.FontSize = 12;
%                 title(topotitles{itop})
            end
        end
    end
end
ref = refline(0,0);
ref.Color = 'k';
ref.LineWidth = 1;
ref.LineStyle = '--';
title('Liberal - conservative')

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'prestim');
    
    mkdir(fullfile(outpath, 'stats'))
    %                 if showstats
    %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
    %                 else
    outfile = fullfile(outpath, sprintf( 'spectra_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    %                 end
    disp(outfile)
    export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

