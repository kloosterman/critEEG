% load in basespec and plot

respavg = critEEG_load_respavg()

%%
nsub = length(respavg.SUBJ);

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.025; % for evoked data
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;
cfg0_neighb.template  = 'elec1010_neighb.mat';
cfg0_neighb.elecfile  = 'standard_1020.elc';
cfg0_neighb.feedback = 'no';
cfg0_neighb.method = 'template';
cfg.neighbours       = ft_prepare_neighbours(cfg0_neighb);

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

cfg.latency = [-0.8 -0.2]; % -0.8 to 0 : p=0.02

% just occ, spectrum
cfg.channel = respavg.sens.ind{1};
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';

ises=4;
isoi=1;
poolstat = [];
for icond = 4 %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
    for istim = 3%:2%:3
        for iresp = 3%:2
            for iband = 1
                for itrig = 1
                    
                    freq = [];
                    freq.dimord = 'subj_chan_freq_time';
                    freq.label = respavg.label;
                    freq.freq = respavg.freq{iband};
                    freq.time = respavg.time{itrig};
                    
                    freq.powspctrm  = squeeze(respavg.pow(:,:, 1:length(respavg.freq{iband}), 1:length(freq.time),   iband, itrig,     ises,icond,istim,iresp ));
                    
                    freqzero = freq; %create zero freq to test against
                    freqzero.powspctrm = zeros(size(freq.powspctrm));
                    poolstat{ iband, icond, istim, iresp } = ft_freqstatistics(cfg, freq, freqzero);
                    
                    stat = poolstat{ iband, icond, istim, iresp };
                    
                    tind = respavg.time{itrig} > cfg.latency(1) & respavg.time{itrig} < cfg.latency(2);
                    stat.powspctrm = mean(freq.powspctrm(:,:,:,tind),4);
                end
            end
        end
    end
end
stat.negclusters(1).prob

%% plot significant cluster in freq (space)

% chans = any(stat.mask,2);
chans = respavg.sens.ind{1};
% dum = squeeze(trapz(stat.powspctrm(:,chans,:),2));
dum = squeeze(mean(stat.powspctrm(:,chans,:),2));
% dum = dum([1:7, 9:15],:)
close all
f = figure; hold on
f.Position = [680 873 345 225];
h =shadedErrorBar(stat.freq, mean(dum ), std(dum)/sqrt(size(dum,1)), 'k', 1 );
h.mainLine.LineWidth = 2;
r= refline(0,0);
r.Color = 'k';
r.LineStyle = '--';
% plos sig bar
sigfreqs = any(stat.mask,1);
barind = find(diff(sigfreqs));
barind(1) = barind(1)+1
if mod(length(barind),2), barind = [barind length(any(mask,1))]; end
% ylim([-2e-8 0.5e-8])
% ylim([-15 10])
ax=gca;
ax.FontSize = 12;
% ypos = -12.5;% diff(ax.YLim)*-0.1;
ypos = diff(ax.YLim)*-0.1;
pl = plot( stat.freq(barind), [ypos ypos] );
pl.LineWidth = 10;
pl.Color = 'k';
xlabel('Frequency (Hz)')
ylabel('Power change from conservative (%)')
xlim([0 35])

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
freq=[];
freq.label = respavg.label;
freq.dimord = 'chan_freq';
freq.freq = mean(stat.freq(sigfreqs));
dum = squeeze(mean(stat.powspctrm(:,:,sigfreqs),3));
% dum = squeeze(trapz(stat.powspctrm(:,:,freqs),3));

freq.powspctrm = mean(dum)';
% cfg.zlim = [-10 10];
                    cfg.zlim = 'maxabs';
%                     cfg.zlim = scales(itop,:);
ft_topoplotTFR(cfg, freq);
colorbar
ax=gca;
ax.FontSize = 12;

%% correlate alpha with crit
ises=1:3;
isoi=1;
poolstat = [];
icond = 4; %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
istim = 3;%:2%:3
iresp = 3;%:2
iband = 1;
itrig = 1;
dum2 = squeeze(mean(respavg.pow(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
    iband, itrig, ises,icond, istim, iresp), 2)); %average over sens
dum2 = mean(dum2(:,sigfreqs,:,:),2);
alphapowdat = squeeze(mean(dum2(:,:,tind,:),3));
behavdat = respavg.behavior.criterion([1:13,15,16],1:3,2) - respavg.behavior.criterion([1:13,15,16],1:3,1);

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

%% Plot TFR for occipital pooling

close all
SAV = 1;

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
                        
                        dum = squeeze(mean(respavg.pow(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:17, ... ipv length(respavg.time{itrig}), ...
                            iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        dumsubj = dum; % keep singlesub
                        dum = squeeze(mean(dumsubj));
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

%% Plot spectra for occipital pooling

TIM = [-0.8 -0.2];
tind = respavg.time{itrig} > TIM(1) & respavg.time{itrig} < TIM(2)

close all
SAV = 1;

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
                        ax.YLim = [-0.3e-7 0.11e-7];
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

