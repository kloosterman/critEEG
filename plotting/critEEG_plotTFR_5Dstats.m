
% respavg = critEEG_load_respavg()

%% freqstatistics on TFR:
% motor and occ NoMiss vs NoFA
nsub = length(respavg.SUBJ);

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
% cfg.correctm         = 'no';
% cfg.clusteralpha     = 0.0001; % for evoked data
cfg.clusteralpha     = 0.05; % for evoked data
% cfg.clusteralpha     = 0.1;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% cfg.neighbours       = []; %in case no channel data present
if ismac
    cfg0_neighb.template  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20161220/template/neighbours/elec1010_neighb.mat';
    cfg0_neighb.elecfile  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20161220/template/electrode/standard_1020.elc';
else
    cfg0_neighb.template  = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20161220/template/neighbours/elec1010_neighb.mat';
    cfg0_neighb.elecfile  = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20161220/template/electrode/standard_1020.elc';
end
cfg0_neighb.feedback = 'yes'
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

cfg.minnbchan = 1;

poolstat=cell(2,2,2,4,2,2);
% for isoi = 1%:2 % 2%:2 % occ and motor
    for ises = 4 % 4 = collapsed over session
        for icond = 1:4 %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
            for istim = 1%:2%:3
                for iresp = 1%:2
                    for iband = 1:2
                        for itrig = 1:2
                            
                            freq2stats = [];
                            freq2stats.dimord = 'subj_chan_freq_time';
                            freq2stats.label =respavg.label;
                            freq2stats.time = respavg.time{itrig};
                            freq2stats.freq = respavg.freq{iband};
                            
                            freq2stats.powspctrm  = squeeze(respavg.dat(:,:, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                                iband, itrig, 4,icond, istim, iresp)); 
                            
                            freq2statszero = freq2stats; %create zero freq to test against
                            freq2statszero.powspctrm = zeros(size(freq2stats.powspctrm));
                            poolstat{isoi, iband, itrig, icond, istim, iresp} = ft_freqstatistics(cfg, freq2stats, freq2statszero);
                        end
                    end
                end
            end
        end
    end
% end

%% plot multiplot for significant clusters

close all

istim=1; iresp=1; 
for icond = 4% 1:4 %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
    for iband = 1:2%:2
        for itrig = 2 %1:2
%             freq2stats.mask = poolstat{isoi, iband, itrig, icond, istim, iresp}.mask;
%             freq2stats.mask = poolstat{isoi, iband, itrig, icond, istim, iresp}.posclusterslabelmat == 1;
            freq2stats.mask = poolstat{isoi, iband, itrig, icond, istim, iresp}.negclusterslabelmat == 2;
            if ~any(freq2stats.mask)
                continue;
            end
            freq2stats.time = respavg.time{itrig};
            freq2stats.freq = respavg.freq{iband};
            freq2stats.powspctrm  = squeeze(respavg.dat(:,:, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                iband, itrig, 4,icond, istim, iresp));
            
            cfg = [];
%             cfg.layout = 'biosemi64incI1I2.lay';  % biosemi64  elec1010
            % cfg.layout = 'biosemi64.lay';  % biosemi64  elec1010
            cfg.layout = 'elec1010.lay';  % biosemi64  elec1010
            cfg.layout = ft_prepare_layout(cfg);
            cfg.layout.width(:) = 0.075;
            cfg.layout.height(:) = 0.075;
            
            cfg.zlim = [-0.1 0.1];
            % cfg.zlim = [-1 1];
            load( 'colormap_jetlightgray.mat')
            cfg.colormap = cmap;
            cfg.maskparameter = 'mask';
            cfg.hotkeys = 'yes';
            cfg.fontsize = 18;
            
            f = figure;
            f.Position = [ 680   678   1200   1200];
            
            ft_multiplotTFR(cfg, freq2stats);
            
            title( sprintf('%s, %s, %s, %s',  respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}, respavg.trigger_leg{itrig} ))
        end
    end
end
%% old
% %% plotting: critEEG TFR NoFA vs NoMiss: motor and visual cortex + stats
% close all
% SAV = 1;
% 
% % set(0,'DefaultAxesFontsize',10)
% addpath(genpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/plotting'));
% % addpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/cbrewer')
% 
% showstats = 1;
% % istim = 2;
% % iresp = 3;
% %
% % cmap = cbrewer('div', 'RdBu',256);
% % cmap = flipud(cmap);
% 
% % load('colormap170613.mat');
% load('colormap_jetgray', 'cmap');
% 
% % cmap = jet(256);
% % for topo:
% cfg = [];
% cfg.layout = 'elec1010.lay';
% cfg.comment = 'no';
% cfg.marker = 'on';
% cfg.shading = 'flat';
% cfg.style = 'straight'; %both  straight
% cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
% cfg.markersize = 4;
% cfg.highlight = 'on';
% cfg.highlightsymbol = '.';
% cfg.highlightsize = 20;
% cfg.highlightchannel = respavg.sens.ind{1};
% cfg.colormap = cmap;
% 
% ncols = 4;
% 
% ises = 4;
% for istim = 1%:2
%     for iresp = 1%:2        %1:2;
%         for isoi = 1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
%             
%             figure;    iplot=0; hold on
% %             set(gcf, 'Position', [0 -200 1000 500]) % 2,3 subplot
%             set(gcf, 'Position', [0 -200 1000/3*ncols 500]) %
%             set(0,'DefaultAxesFontsize',30)
%             pl = [];
%             for iband = 2:-1:1 % 1:2 %
%                 for icond = [1:2, 4] %1:4 %1:4%  % 1:2
%                     for itrig = 1%:2
%                         
%                         ZLIMS = [-0.3 0.3; -0.1 0.1];
% %                         ZLIMS = [-0.4 0.4; -0.05 0.05];
% %                         ZLIMS = [-0.5 0.5; -0.2 0.2];
%                         if icond == 4; ZLIMS = ZLIMS/2; end
%                         ZLIM = ZLIMS(iband,:);
%                         
%                         iplot = iplot+1;
%                         pl = subplot(2,ncols,iplot); hold on
%                         pl.Units = 'pixels';
%                         if iplot == 1 || iplot == 5 %~mod(iplot,2)
% %                             pl.YTickLabel = [];
%                             ylabel('Frequency (Hz)');
%                         else
%                         end
%                         colormap(cmap)
%                         
%                         dum = squeeze(mean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
%                             iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
%                         dumsubj = dum; % keep singlesub
%                         dum = squeeze(mean(dumsubj));
%                         
%                         if showstats
%                             mask = double(squeeze( poolstat{isoi, iband, itrig, icond, istim, iresp}.mask));
%                             mask(mask==0) = 0.25;
%                             %                             mask = double(squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.prob < 0.15));
%                             %                             mask(mask==0) = 0.25;
%                             
%                             ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no', ...
%                                 'highlight', mask, 'highlightstyle', 'opacity'); % opacity
%                         else
%                             ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no' ); % opacity
%                         end
%                         
%                         XLIM = [respavg.time{itrig}(1) respavg.time{itrig}(end)];
%                         YLIM = [respavg.freq{iband}(1) respavg.freq{iband}(end)];
%                         xlim(XLIM)
%                         ylim(YLIM)
%                         
% %                         binwidth = 0.01; %on screen
% %                         %                         pl.Position(3) = pl.Position(3) * diff(XLIM);
% %                         pl.Position(3) = length(respavg.time{itrig}) * binwidth;
% %                         
%                         ax = gca;
%                         ax.FontSize = 12;
%                         ax.Box = 'on';
%                         ax.XTick = -1:0.5:2;
%                         if iband==1
%                             ax.YTick = [0:5:200];
% %                             ax.YTick = [0:10:200];
%                         else
%                             ax.YTick = [0:10:200];
% %                             ax.YTick = [0:20:200];
%                         end
%                         
%                         if iplot < 4
%                             title( sprintf('%s, %s, psc Lo: %g, Hi: %g\n%s', respavg.sdt_conds{istim, iresp}, respavg.sens.leg{isoi}, ZLIMS(:,2)*100, respavg.behav_conds{icond} ))
%                             pl.XTickLabel = [];
%                         else
%                             pl.Position(2) = pl.Position(2)*1.95;
% %                             pl.Position(4) = pl.Position(4)*1.5;
%                             xlabel(sprintf('Time from %s (s)', respavg.trigger_leg{itrig} ));
% 
% %                             title( sprintf('%s, [%g-%g%%]', respavg.behav_conds{icond}, ZLIM*100))
%                         end
%                         
%                         %
%                         %                         title({sprintf('%s %s', respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}) ;
%                         %                             sprintf('%s, scale [%g-%g]', respavg.sens.leg{isoi}, ZLIM )}); % 'FontWeight','bold'
%                         
%                         yaxis = [respavg.freq{iband}(1), respavg.freq{iband}(end)];
%                         plot([0,0],yaxis,'k',[0,0],yaxis,'k');
%                         if strcmp(respavg.trigger_leg{itrig}, 'stim')
%                             stimonset = 0.16;
%                             plot([stimonset,stimonset],yaxis,'--k',[stimonset,stimonset],yaxis,'--k');
%                             plot([1,1],yaxis,'--w',[1,1],yaxis,'--w', 'Linewidth', 2);
% 
%                             RT = mean(respavg.behavior.RT(:, ises, istim, icond)); %RT
%                             if icond<4
%                                 plot(RT, respavg.freq{iband}(end-1), 'vk', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
%                             end
%                         end
%                         % h=colorbar;
%                       
%                         %topo
%                         if iplot == 3 && showstats
%                             cond2plot = 4;
%                             band2plot = 1;
%                             iplot = iplot+1;
%                             pl = subplot(2,4,iplot); hold on
%                             freq=[];
%                             freq.label = respavg.label;
%                             freq.dimord = 'chan';
%                             dum = squeeze(mean(respavg.dat(:, :, 1:length(respavg.freq{band2plot}), 1:length(respavg.time{itrig}), ...
%                                 band2plot, itrig, 4,cond2plot, istim, iresp), 1)); %average over sens
%                             mask = squeeze( poolstat{isoi, band2plot, itrig, cond2plot, istim, iresp}.mask);
%                             % inspect labelmat to see indexing:
%                             %                         figure; imagesc(poolstat{isoi, iband, itrig,
%                             %                         icond, istim, iresp}.time, poolstat{isoi, iband,
%                             %                         itrig, icond, istim, iresp}.freq,
%                             %                         squeeze(poolstat{isoi, iband, itrig, icond,
%                             %                         istim, iresp}.posclusterslabelmat)); colorbar
%                             %                         mask = squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.posclusterslabelmat == 2);
%                             freq.avg = squeeze(mean(dum(:,mask),2));
%                             cfg.zlim = [-0.15 0.15];
%                             ft_topoplotTFR(cfg, freq);
% 
%                         end
%                         
%                         % spectrum
%                         if iplot == 7 && showstats
%                             iplot = iplot+1;
%                             pl = subplot(2,4,iplot); hold on
%                             
%                             mask = double(squeeze( poolstat{isoi, iband, itrig, icond, istim, iresp}.mask));
%                             dumsubj(:,~mask) = NaN;
%                             dum = squeeze(nanmean(dumsubj,3)) * 100; % avg over time
%                             h = shadedErrorBar(respavg.freq{iband}, mean(dum), std(dum)/sqrt(nsub), 'k' );
% %                             h.mainLine.Color = cmap(1,:); % fix color? 
% %                             h.patch.FaceColor = [0 0 1]; % fix color? 
% %                             h.patch.FaceAlpha = 0.25;
% %                             % pl.ColorOrder = cmap(end,:); % fix color?
% %                             xlim([ respavg.freq{iband}(1) respavg.freq{iband}(end)])
% %                             pl.XAxisLocation = 'top'; %freq axis
% %                             camroll(90);  %camorbit(0,180)
%                             pl.Position(2) = pl.Position(2)*1.95;
%                             ax = gca;
%                             ax.FontSize = 12;
%                             ax.Box = 'on';
%                             ax.XTick = [0:5:35];
% %                             ax.YDir = 'reverse';
%                             xlabel('Frequency (Hz)')
%                             ylabel( 'Modulation (%)' )
%                             xlim([3 25])
% 
%                         end
% 
%                     end %ises
%                 end
%             end
%             set(gcf, 'Color', 'white')
%             
%             if SAV
%                 %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
%                 outpath = fullfile(respavg.PREOUT, 'poolings');
%                 
%                 mkdir(fullfile(outpath, 'stats'))
%                 %                 if showstats
%                 %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
%                 %                 else
%                 outfile = fullfile(outpath, sprintf( 'TFRstats_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
%                 %                 end
%                 disp(outfile)
%                 export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
%                 export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
%                 cd(outpath)
%             end %ises
%         end
%     end
% end %ises
% 
% %% plot topo's based on stat.mask field
% % close all
% load('colormap_jetgray', 'cmap');
% SAV=1;
% ises = 4;
% for istim = 1%:2
%     for iresp = 1%:2        %1:2;
%         for isoi = 1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
%             
%             figure;    iplot=0; hold on
%             set(gcf, 'Position', [0 -200 1000 500])
%             set(0,'DefaultAxesFontsize',30)
%             
%             for iband = 2% 2:-1:1 % 1:2 %
%                 for icond = 3 %[1,3, 4]  %[1:2, 4] %1:4 %1:4%  % 1:2
%                     for itrig = 1%:2
%                         
%                         ZLIMS = [-0.4 0.4; -0.05 0.05];
% %                         if icond == 4; ZLIMS = ZLIMS/4; end
%                         ZLIM = ZLIMS(iband,:);
%                         
%                         iplot = iplot+1;
%                         pl = subplot(1,1,iplot); hold on
% %                         if iplot == 1 || iplot == 4 %~mod(iplot,2)
% %                             ylabel('Frequency (Hz)');
% %                         else
% %                         end
% %                         colormap(cmap)
%                         
%                         freq=[];
%                         freq.label = respavg.label;
%                         freq.dimord = 'chan';
%                         dum = squeeze(mean(respavg.dat(:, :, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
%                             iband, itrig, 4,icond, istim, iresp), 1)); %average over sens
%                         mask = squeeze( poolstat{isoi, iband, itrig, icond, istim, iresp}.mask);
%                         % inspect labelmat to see indexing:
% %                         figure; imagesc(poolstat{isoi, iband, itrig,
% %                         icond, istim, iresp}.time, poolstat{isoi, iband,
% %                         itrig, icond, istim, iresp}.freq,
% %                         squeeze(poolstat{isoi, iband, itrig, icond,
% %                         istim, iresp}.posclusterslabelmat)); colorbar
% %                         mask = squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.posclusterslabelmat == 2);
%                         freq.avg = squeeze(mean(dum(:,mask),2));
%                       
%                         cfg = [];
%                         cfg.layout = 'elec1010.lay';
%                         cfg.comment = 'no';
%                         cfg.marker = 'on';
%                         cfg.shading = 'flat';
%                         cfg.style = 'straight'; %both  straight
%                         cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
%                         cfg.markersize = 6;
%                         cfg.highlight = 'on';
%                         cfg.highlightsymbol = '.';
%                         cfg.highlightsize = 30;
%                         cfg.highlightchannel = respavg.sens.ind{1};
%                         cfg.colormap = cmap;
%                         cfg.zlim = ZLIM;
%                         ft_topoplotTFR(cfg, freq);
% 
%                         colorbar
% 
%                     end %ises
%                 end
%                 
%                 set(gcf, 'Color', 'white')
%             end
%             
%             if SAV
%                 %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
%                 outpath = fullfile(respavg.PREOUT, 'poolings');
%                 
%                 mkdir(fullfile(outpath, 'stats'))
%                 %                 if showstats
%                 %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
%                 %                 else
%                 outfile = fullfile(outpath, sprintf( 'topo_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
%                 %                 end
%                 disp(outfile)
%                 export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
%                 export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
%                 cd(outpath)
%             end %ises
%         end
%     end
% end %ises
% 
% %% Plot spectrum of high-low crit power
% cfg = [];
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'depsamplesT';
% cfg.correctm         = 'cluster';
% % cfg.correctm         = 'no';
% % cfg.clusteralpha     = 0.0001; % for evoked data
% cfg.clusteralpha     = 0.05; % for evoked data
% % cfg.clusteralpha     = 0.1;
% cfg.clusterstatistic = 'maxsum';
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.025;
% cfg.numrandomization = 1000;
% cfg.neighbours       = []; %in case no channel data present
% 
% design = zeros(2,2*nsub);
% for i = 1:nsub
%     design(1,i) = i;
% end
% for i = 1:nsub
%     design(1,nsub+i) = i;
% end
% design(2,1:nsub)        = 1;
% design(2,nsub+1:2*nsub) = 2;
% 
% cfg.design   = design;
% cfg.uvar     = 1;
% cfg.ivar     = 2;
% 
% close all
% SAV=1;
% ises = 4;
% YLIMS = [5; 2.5];
% for istim = 1%:2
%     for iresp = 1%:2        %1:2;
%         for isoi = 1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
%             
%             figure;    iplot=0; hold on
%             set(gcf, 'Position', [500 500 800 300])
%             set(0,'DefaultAxesFontsize',12)
%             stat = {};
%             for iband = 1:2% 2:-1:1 % 1:2 %
%                 iplot = iplot + 1;
%                 subplot(1,2,iplot)
%                 for icond = 4 %[1,3, 4]  %[1:2, 4] %1:4 %1:4%  % 1:2
%                     for itrig = 1%:2
%                         dum = squeeze(mean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
%                             iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
%                         % take 0.16 - 0.4 s (stim response)
%                         TIM = [0.15 0.3];
%                         tinds = respavg.time{itrig} >= TIM(1) & respavg.time{itrig} <= TIM(2);
%                         dum = mean(dum(:,:,tinds),3);
%                         dumsem = std(dum * 100) / sqrt(nsub);
% %                         dum = squeeze(mean(dum));
%                         %                         % use mask:
%                         % %                         mask = squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.posclusterslabelmat == 1);
%                         %                         mask = squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.mask);
%                         %                         dum(~mask) = NaN;
%                         
%                         freq = [];
%                         freq.label = {'occ'};
% %                         freq.freq = [respavg.freq{1} respavg.freq{2}];
%                         freq.freq = respavg.freq{iband};
%                         freq.dimord = 'chan_subj_freq';
%                         freq.powspctrm = shiftdim(dum, -1);
%                         
%                         freqzero = freq; %create zero freq to test against
%                         freqzero.powspctrm = zeros(size(freq.powspctrm));
%                         stat{iband} = ft_freqstatistics(cfg, freq, freqzero);
%  
% 
%                         
%                         
% %                         plot(respavg.freq{iband}, nanmean(dum,2) * 100);
%                         shadedErrorBar(respavg.freq{iband}, mean(dum) * 100, dumsem);
%                         xlim([ respavg.freq{iband}(1) respavg.freq{iband}(end)])
%                         ylim([-YLIMS(iband) YLIMS(iband)])
%                         line(get(gca, 'XLIM'), [0 0])
%                         title(TIM)
%                         
%                         %                         shadedErrorBar
%                         
%                     end %ises
%                 end
%             end
%         end %ises
%     end
% end %ises




