
% respavg = critEEG_load_respavg()

%% freqstatistics on TFR:
% motor and occ NoMiss vs NoFA
% nsub = length(respavg.SUBJ);
nsub = 14; % subj 8 removed
sub_ind = [1:7, 9:15];

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
% cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05; 
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

% cfg.neighbours       = []; %in case no channel data present

% prepare_neighbours determines what sensors may form clusters
cfg0_neighb = [];
cfg0_neighb.method    = 'template';
if ismac
    cfg0_neighb.template  = 'elec1010_neighb.mat';
    cfg0_neighb.elecfile  = 'standard_1020.elc';
else
    cfg0_neighb.template  = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20161220/template/neighbours/elec1010_neighb.mat';
    cfg0_neighb.elecfile  = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20161220/template/electrode/standard_1020.elc';
end
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

toi = [-0.2 0.9];
cfg.latency = toi;
tind = respavg.time{1} >= cfg.latency(1) & respavg.time{1} <= cfg.latency(2);

poolstat=cell(2,2,2,4,2,2);
freqexp = {}; % for source

for isoi = 1%:2 % 2%:2 % occ and motor
    cfg.channel = respavg.sens.ind{isoi};
    cfg.avgoverchan = 'yes';
    for ises = 4 % 4 = collapsed over session
        for icond = 1:4 %1:4 %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
            for istim = 3%:3%1:2%:3
                for iresp = 3%:3% 1:2
                    for iband = 1:2
                        for itrig = 1%:2 %1%:2
                            
                            freq2stats = [];
                            freq2stats.dimord = 'subj_chan_freq_time';
                            freq2stats.label = respavg.label;
                            freq2stats.time = respavg.time{itrig};
                            freq2stats.freq = respavg.freq{iband};
                            
%                             freq2stats.powspctrm  = squeeze(respavg.dat(:,:, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
%                                 iband, itrig, 4,icond, istim, iresp)); %average over sens
                            freq2stats.powspctrm  = squeeze(respavg.dat(sub_ind,:, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                                iband, itrig, 4,icond, istim, iresp)); %average over sens
                            

                            freqexp{iband} = freq2stats; % for source
                            freq2statszero = freq2stats; %create zero freq to test against
                            freq2statszero.powspctrm = zeros(size(freq2stats.powspctrm));
                            poolstat{isoi, iband, itrig, icond, istim, iresp} = ft_freqstatistics(cfg, freq2stats, freq2statszero);
                        end
                    end
                end
            end
        end
    end
end

cfg2=[];
cfg2.latency = [-0.2 0.9];
for iband=1:2
    freq = ft_selectdata(cfg2, freqexp{iband});
    freqexp{iband} = freq;
end

%% plotting: critEEG TFR NoFA vs NoMiss: motor and visual cortex + stats
close all
SAV = 1;

% set(0,'DefaultAxesFontsize',10)
% addpath(genpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/plotting'));
% addpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/cbrewer')

showstats = 1;
% istim = 2;
% iresp = 3;
%
% cmap = cbrewer('div', 'RdBu',256);
% cmap = flipud(cmap);

% load('colormap170613.mat');
% load('colormap_jetgray', 'cmap');
load colormap_jetlightgray.mat
% cmap = jet(256);

ncols = 3;

ises = 4;
for istim = 3%:3 %1:2
    for iresp = 3%:3%:3 %1:2        %1:2;
        for isoi = 1 %1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
            
            for icond = 3% 1:4%:4 %[1:2, 4] %1:4 %1:4%  % 1:2
                f = figure;     hold on
                %             set(gcf, 'Position', [0 -200 1000 500]) % 2,3 subplot
%                 f.Position = [2029         314         1080        500]; %
                f.Position = [2029         314         600        350]; %
                set(0,'DefaultAxesFontsize',8)
                pl = [];
                for iband = 2:-1:1 % 1:2 %
                    for itrig = 1%:2
                        
                        %                         ZLIMS = [-30 30; -10 10];
                        ZLIMS = [-20 20; -7 7];
                        if icond == 4
                            ZLIMS = [-10 10; -10 10];
                        end
                        ZLIM = ZLIMS(iband,:);
                        
                        if iband == 1; iplot = 4; else iplot = 1; end
                        pl = subplot(2,ncols,iplot); hold on
%                         pl.Units = 'pixels';
                    if iplot == 1; yl = ylabel('Frequency (Hz)'); yl.Units = 'Normalized'; yl.Position(2) = 0;  end
                        colormap(cmap)
                        
                        dum = squeeze(mean(respavg.dat(sub_ind,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), tind, ...
                            iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        %                         dum = squeeze(mean(respavg.pow(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                        %                             iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        dumsubj = dum; % keep singlesub
                        dum = squeeze(mean(dumsubj));
                        
                        if showstats
                            mask = double(squeeze( poolstat{isoi, iband, itrig, icond, istim, iresp}.mask));
                            mask(mask==0) = 0.25;
                            %                             mask = double(squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.prob < 0.15));
                            %                             mask(mask==0) = 0.25;
                            
                            ft_plot_matrix(respavg.time{itrig}(tind), respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no', ...
                                'highlight', mask, 'highlightstyle', 'opacity'); % opacity
                        else
                            ft_plot_matrix(respavg.time{itrig}(tind), respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no' ); % opacity
                            %                             ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum); % opacity
                        end
                        if iplot == 1; pl.Position(2) = pl.Position(2)-0.05;
                        else pl.Position(2) = pl.Position(2)+0.05;
                        end
                        pl.Position(3) = pl.Position(3)*1.2;
                        pl.Position(1) = pl.Position(1) - 0.05;
                        %                         XLIM = [respavg.time{itrig}(1) respavg.time{itrig}(end)];
                        xlim(toi);
                        %                         xlim([respavg.time{itrig}(1) respavg.time{itrig}(end)])
                        YLIM = [respavg.freq{iband}(1) respavg.freq{iband}(end)];
                        
                        ylim(YLIM)
                        
                        %                         binwidth = 0.01; %on screen
                        %                         %                         pl.Position(3) = pl.Position(3) * diff(XLIM);
                        %                         pl.Position(3) = length(respavg.time{itrig}) * binwidth;
                        %
                        ax = gca;
                        ax.FontSize = 8;
                        ax.Box = 'on';
                        ax.XTick = -1:0.25:2;
                        if iband==1; ax.YTick = [0:10:200]; else ax.YTick = [0:20:200]; end
                        
                        if iplot == 1
%                             title( sprintf('%s, %s, psc Lo: %g, Hi: %g\n%s', respavg.sdt_conds{istim, iresp}, respavg.sens.leg{isoi}, ZLIMS(:,2), respavg.behav_conds{icond} ))
                            pl.XTickLabel = [];
                        else
%                             pl.Position(2) = pl.Position(2)*1.95;
                            %                             pl.Position(4) = pl.Position(4)*1.5;
                            %                             xlabel(sprintf('Time from %s (s)', respavg.trigger_leg{itrig} ));
                            if itrig == 1
                                xlabel('Time from trial onset (s)');
                            else
                                xlabel('Time from response (s)');
                            end
                            
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
                        c=colorbar;

                        c.Label.String = 'Modulation (%)';
                        c.Position(4) = c.Position(4)*0.5;
                        c.Position(3) = c.Position(3)*0.5;
                        c.Position(1) = c.Position(1)+0.075;
                        c.Position(2) = c.Position(2)+0.075;
                        c.Box = 'off';

                        % Specify pos as a four-element vector of the form [x y w h] in data units.
                        if iplot == 1
                            boxlocs = [0.2 59 0.4 40; ...
                                0.2 42 0.4 16];
                        else
                            boxlocs = [0.2 11 0.4 11; ...
                                0.2 23 0.4 4];
                        end
                        for ib = 1:2
                            %                             if ib < 3; subplot(1,2,1); else subplot(1,2,2); end
                            r=rectangle('Position', boxlocs(ib,:));
                            r.LineStyle = '--';
                        end
                        
                    end %ises
                end
                set(gcf, 'Color', 'white')
                
            end
        end
    end
end %ises


% plot topo's for beta suppr, ssvep, ssvep harmonic and gamma
load colormap_jetlightgray.mat

% close all
cfg = [];
cfg.layout = 'elec1010.lay';
cfg.comment = 'no';
cfg.marker = 'on';
cfg.shading = 'flat';
cfg.style = 'straight'; %both  straight
cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
cfg.markersize = 3;
cfg.highlightchannel = respavg.sens.ind{1};
cfg.highlightsymbol = '.';
cfg.highlightsize = 20;
cfg.highlightcolor = [1 1 1];
cfg.colormap = cmap;

% boxlocs = [0.2 59 0.2 40; ...
%     0.2 42 0.2 16];
% boxlocs = [0.25 11 0.2 11; ...
%     0.25 24 0.2 3];
% times = [0.2 0.4; 0.2 0.4; 0.2 0.4; 0.2 0.4];
% times = [0.2 0.6; 0.2 0.6; 0.2 0.6; 0.2 0.6];
times = [0.2 0.5; 0.2 0.5; 0.2 0.5; 0.2 0.5];
freqs = [42 58; 59 100; 23 27; 11 22];
topotitles = {'SSVEP harmonic', 'Gamma', 'SSVEP', 'Beta'};
topotitlesfreq = {'42-58 Hz', '59-100 Hz', '23-27 Hz', '11-22 Hz'};
scales = [-12 12; -5 5; -25 25;  -10 10];

% ncols = 2;
ises = 4;
for istim = 3
    for iresp = 3        %1:2;
        
%         f = figure;    iplot=0; hold on
%         f.Position = [ 2063         121        1063         849]; % 
%         set(0,'DefaultAxesFontsize',30)
        pl = [];
        for icond = 3 %[1:2, 4] %1:4 %1:4%  % 1:2
            for itrig = 1%:2
                for itop = 1:4
                    if itop < 3
                        iband = 2; 
                        cfg.highlight = 'on';
                    else
                        iband = 1; 
                        cfg.highlight = 'off';
                    end
                    
                    ZLIMS = [-30 30; -10 10];
                    if icond == 4; ZLIMS = ZLIMS/2; end
                    ZLIM = ZLIMS(iband,:);
                    
                    %                     iplot = iplot+1;
                    if itop == 1; iplot = 2;
                    elseif itop == 2; iplot = 3;
                    elseif itop == 3; iplot = 5;
                    elseif itop == 4; iplot = 6;
                    end
                    pl = subplot(2,ncols,iplot); hold on
                    if iplot < 4
                        pl.Position(2) =  pl.Position(2) - 0.02;
                    else
                        pl.Position(2) =  pl.Position(2) + 0.02;
                    end
                    %                     pl.Units = 'pixels';
                    colormap(cmap)
                    
                    %topo
                    cfg.xlim = times(itop,:);
                    cfg.ylim = freqs(itop,:);
                    freq=[];
                    freq.label = respavg.label;
                    freq.dimord = 'chan_freq_time';
                    freq.time = respavg.time{itrig};
                    freq.freq = respavg.freq{iband};
                    dum = squeeze(mean(respavg.dat(sub_ind, :, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                        iband, itrig, 4,icond, istim, iresp), 1)); %average over sens
                    %                             freq.avg = squeeze(mean(dum(:,mask),2));
                    freq.powspctrm = dum;
%                     cfg.zlim = [-0.15 0.15];
%                     cfg.zlim = 'maxabs';
                    cfg.zlim = scales(itop,:);
                    ft_topoplotTFR(cfg, freq);
                    c = colorbar;
                    c.Ticks = [scales(itop,1) 0 scales(itop,2)];
                    c.Position(4) = c.Position(4)*0.5;
                    c.Position(3) = c.Position(3)*0.5;
                    c.Position(1) = c.Position(1)+0.065;
                    c.Position(2) = c.Position(2)+0.03;
                    c.Label.String = 'Modulation (%)';
                    c.Box = 'off';
                    ax=gca;
                    ax.FontSize = 8;
                    title(sprintf('%s\n%s', topotitles{itop}, topotitlesfreq{itop}))
                end
            end
        end
    end
end

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'topos');
    
    mkdir(fullfile(outpath, 'stats'))
    %                 if showstats
    %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
    %                 else
    outfile = fullfile(outpath, sprintf( 'Topo_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    %                 end
    disp(outfile)
%     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
        f.Renderer = 'Painters';
        print(outfile, '-depsc2') %'-png',  '-pdf',
        %                     print(outfile, '-depsc2') %'-png',  '-pdf',
        print(outfile, '-dpdf') %'-png',  '-pdf',

    cd(outpath)
end %ises

%% plot bars and time courses of gamma for lib and cons
% gamma
SAV = 1;
close all
conds = 1:2;  itrig = 1;
cfg = [];
cfg.channel = respavg.sens.ind{1}; cfg.avgoverchan = 'yes';

% whattoplot = 'bars';
whattoplot = 'timecourses';
stats = 1;

foi = [42 58; 59 100; 23 27;  11 22; ];
% foi = [42 58; 59 100; 23 27;  3 6; ]; % ERP marginally sig
cfg.avgoverfreq = 'yes';
switch whattoplot
    case 'bars'
        toi = [ 0.2 0.6; 0.2 0.6;  0.2 0.6; 0.25 0.6; ];
        % toi = [ 0.25 0.65; 0.25 0.65; 0.25 0.65; 0.25 0.65; ];
        cfg.latency = toi(1,:);
        cfg.avgovertime = 'yes';
    case 'timecourses'
        cfg.latency = [-0.4 0.9];
        cfg.avgovertime = 'no';
end

linecol = cbrewer('qual', 'Set1', 3);
linecol(4,:) = [0 0 0];
titleg = {'SSVEP harmonic', 'Gamma', 'SSVEP', 'Beta'};

f = figure; 
iplot= 0;
linestyle = {'-', '--'};
stat = {};
for ifoi = 1:4
    cfg.frequency = foi(ifoi,:);
    if min(cfg.frequency) > 35; iband = 2; else iband = 1; end;
    
    for istim = 3
        for iresp = 3
%             iplot = iplot+1; sp = subplot(2,2, iplot);
            iplot = iplot+1; sp = subplot(1,4, iplot);
            
            clear freqdat;
            for icond = 1:2
                freq = [];
                freq.label = respavg.label;
                freq.time = respavg.time{itrig};
                freq.freq = respavg.freq{iband};
                freq.dimord = 'subj_chan_freq_time';
                freq.powspctrm = squeeze(respavg.dat(sub_ind,:, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                    iband, itrig, 4, icond, istim, iresp));
                freqsel = ft_selectdata(cfg, freq);
                freqdat(:,:,icond) = squeeze(freqsel.powspctrm);
            end

            %         iplot = iplot+1; sp = subplot(3,3, iplot);
            switch cfg.avgovertime
                case 'yes'
                    f.Position = [   204         346        500         750 ];
                    axis square
                    b = barweb( squeeze(mean(freqdat)), squeeze(std(freqdat)/sqrt(14)), 0.5,[],[],[],[], linecol(1:2,:) );
                                    pval = permtest(freqdat(:,2) - freqdat(:,1));
%                     [~, pval] = ttest(freqdat(:,2), freqdat(:,1));
                    title(sprintf('%s, p = %g', titleg{ifoi}, pval));
                    ylabel(sprintf('%d-%d Hz power (psc)', cfg.frequency))
                otherwise
%                     f.Position = [   204         346        800         600 ];
                    f.Position = [    2048         600        1250         250 ];
                    hold on; box on
                    for icond = 1:2
                        s = shadedErrorBar(freqsel.time, squeeze(mean(freqdat(:,:,icond))), squeeze(std(freqdat(:,:,icond))/sqrt(14)), [], 1  );
                        s.mainLine.Color =  linecol(icond,:); s.patch.FaceColor =  linecol(icond,:);
                        %                     s.mainLine.LineStyle = linestyle{iresp}
                    end
                    if stats
                        freq.powspctrm = squeeze(respavg.dat(sub_ind,:, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                            iband, itrig, 4, 4, istim, iresp)); % take icond 4: difference
                        %                     freqsel = ft_selectdata(cfg, freq);
                        cfg2 = cfg;
                        cfg2.method           = 'montecarlo';
                        cfg2.statistic        = 'depsamplesT';
                        cfg2.correctm         = 'cluster';
                        % cfg2.correctm         = 'no';
                        cfg2.clusteralpha     = 0.05;
                        cfg2.clusterstatistic = 'maxsum';
                        cfg2.tail             = 0;
                        cfg2.clustertail      = 0;
                        cfg2.alpha            = 0.05;
                        cfg2.numrandomization = 1000;
                        cfg2.neighbours       = []; %in case no channel data present
                        
                        design = zeros(2,2*nsub);
                        for i = 1:nsub
                            design(1,i) = i;
                        end
                        for i = 1:nsub
                            design(1,nsub+i) = i;
                        end
                        design(2,1:nsub)        = 1;
                        design(2,nsub+1:2*nsub) = 2;
                        
                        cfg2.design   = design;
                        cfg2.uvar     = 1;
                        cfg2.ivar     = 2;                        
                        
                        freq0 = freq;
                        freq0.powspctrm = zeros(size(freq.powspctrm));
                        stat{ifoi} = ft_freqstatistics(cfg2, freq, freq0);
                        plot_sig_bar(squeeze(stat{ifoi}.mask)', stat{ifoi}.time, -8, 8, 'k');

                    end
                    
                    %                     title(sprintf('%s', respavg.sdt_conds{istim, iresp}))
                    title(sprintf('%s', titleg{ifoi}));

                    %                 ref = refline(0,0); ref.LineStyle = '--'; ref.Color = 'k';
                    sp.XLim = cfg.latency;
                    sp.XTick = [-0.25:0.25:1];
                    line([0 0], sp.YLim, 'Color', 'k')
%                     line(sp.XLim, [0 0], 'Color', 'k')
                    line([0.16 0.16], sp.YLim, 'Color', 'k', 'Linestyle', '--')
                    ylabel(sprintf('%d-%d Hz power (psc)', cfg.frequency))
                    xlabel('Time from trial start (s)')
                    if iplot == 1; text(0.2, -6, 'P < 0.05, corrected'); ylim([-10 15])
                    elseif iplot == 2; ylim([-10 10]);
                    elseif iplot == 3; ylim([-20 30]);
                    elseif iplot == 4; ylim([-20 30]); 
                    end
            end
%             sp.YLim = [-5 12];
            sp.FontSize = 14;
        end
    end
end
if SAV
    outpath = fullfile(respavg.PREOUT, 'poolings');
    outfile = fullfile(outpath, sprintf( 'power4foi_avgtime%s_%s_%s', cfg.avgovertime, respavg.sens.leg{isoi}, respavg.freqband{iband}));
    disp(outfile)
    %     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    print(outfile, '-dpng') %'-png',  '-pdf',
    %     print(outfile, '-dpdf') %'-png',  '-pdf',
    f.Renderer = 'Painters';
    print(outfile, '-depsc2') %'-png',  '-pdf',
    cd(outpath)
end 



%% plotting: TFR's for all SDT conditions to show that SSVEP mod is always there
close all
SAV = 1;
showstats = 0;
load colormap_jetlightgray.mat
ncols = 2;

ises = 4;
f = figure;    iplot=0; hold on
%             set(gcf, 'Position', [0 -200 1000 500]) % 2,3 subplot
f.Position = [ 2029         419         600         395]; %
set(0,'DefaultAxesFontsize',30)
pl = [];
for istim = 1:2
    for iresp = 1:2        %1:2;
        for isoi = 1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
            
            for iband = 1 % 2:-1:1 % 1:2 %
                for icond = 3 %[1:2, 4] %1:4 %1:4%  % 1:2
                    for itrig = 1%:2
                        
                        ZLIMS = [-30 30; -10 10];
                        if icond == 4; ZLIMS = ZLIMS/2; end
                        ZLIM = ZLIMS(iband,:);
                        
                        iplot = iplot+1;
                        pl = subplot(2,ncols,iplot); hold on
%                         pl.Units = 'pixels';
                        if iplot == 1; yl = ylabel('Frequency (Hz)'); end
                        colormap(cmap)
                        
                        dum = squeeze(mean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                            iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        dumsubj = dum; % keep singlesub
                        dum = squeeze(mean(dumsubj));
                        
                        if showstats
                            mask = double(squeeze( poolstat{isoi, iband, itrig, icond, istim, iresp}.mask));
                            mask(mask==0) = 0.25;
                            %                             mask = double(squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.prob < 0.15));
                            %                             mask(mask==0) = 0.25;
                            
                            ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no', ...
                                'highlight', mask, 'highlightstyle', 'opacity'); % opacity
                        else
                            ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no' ); % opacity
                        end
                        
%                         XLIM = [respavg.time{itrig}(1) respavg.time{itrig}(end)];
                        XLIM = [-0.2 0.9];
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
                        if iband==1; ax.YTick = [0:10:200]; else ax.YTick = [0:20:200]; end
                        
%                         if iplot == 1
%                             title( sprintf('%s, %s, psc Lo: %g, Hi: %g\n%s', respavg.sdt_conds{istim, iresp}, respavg.sens.leg{isoi}, ZLIMS(:,2), respavg.behav_conds{icond} ))
                            pl.XTickLabel = [];
%                         else
%                             pl.Position(2) = pl.Position(2)*1.95;
%                             pl.Position(4) = pl.Position(4)*1.5;
%                             xlabel(sprintf('Time from %s (s)', respavg.trigger_leg{itrig} ));
                            xlabel('Time from trial onset (s)');

                            title( sprintf('%s', respavg.sdt_conds{istim, iresp} ))
%                         end
                        
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
%                         h=colorbar;
%                         h.Position(3) = h.Position(3) * 0.75;
%                         h.Position(4) = h.Position(4) * 0.5;
% %                         h.Position(1) = 0.95
% %                         h.Position(2) = 0.5
%                         h.Ticks = [h.Limits(1) 0 h.Limits(2)];
%                         % Specify pos as a four-element vector of the form [x y w h] in data units.
%                         if iplot == 1
%                             boxlocs = [0.2 59 0.2 40; ...
%                                        0.2 42 0.2 16];
%                         else
%                             boxlocs = [0.25 11 0.2 11; ...
%                                        0.25 24 0.2 3];
%                         end
%                         for ib = 1:2
% %                             if ib < 3; subplot(1,2,1); else subplot(1,2,2); end
%                             r=rectangle('Position', boxlocs(ib,:));
%                             r.LineStyle = '--';
%                         end

                    end %ises
                end
            end
            set(gcf, 'Color', 'white')
        end
    end %ises
end

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'poolings');
    
    mkdir(fullfile(outpath, 'stats'))
    %                 if showstats
    %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
    %                 else
    outfile = fullfile(outpath, sprintf( 'TFRstatsALLSDT_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    %                 end
    disp(outfile)
    export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% plotting: TFR's for all SDT conditions to show that SSVEP mod is always there SINGLE SUBJECTS
close all
SAV = 1;
showstats = 0;
load colormap_jetlightgray.mat
ncols = 3;

ises = 4;
f = figure;    iplot=0; hold on
%             set(gcf, 'Position', [0 -200 1000 500]) % 2,3 subplot
f.Position = [ 2029         419         900         1052]; %
set(0,'DefaultAxesFontsize',30)
pl = [];
for istim = 3
    for iresp = 3        %1:2;
        for isoi = 1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
            
            for iband = 1 % 2:-1:1 % 1:2 %
                for icond = 3 %[1:2, 4] %1:4 %1:4%  % 1:2
                    for itrig = 1%:2
                        
                        dum = squeeze(mean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                            iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        
                        for isub = 1:15
                            
                            ZLIMS = [-30 30; -10 10];
                            if icond == 4; ZLIMS = ZLIMS/2; end
                            ZLIM = ZLIMS(iband,:);
                            
                            iplot = iplot+1;
                            pl = subplot(5,ncols,iplot); hold on
                            %                         pl.Units = 'pixels';
                            if iplot == 1; yl = ylabel('Frequency (Hz)'); end
                            colormap(cmap)
                            %                         dumsubj = dum; % keep singlesub
%                             dum = squeeze(mean(dumsubj));
                            
                            if showstats
                                mask = double(squeeze( poolstat{isoi, iband, itrig, icond, istim, iresp}.mask));
                                mask(mask==0) = 0.25;
                                %                             mask = double(squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.prob < 0.15));
                                %                             mask(mask==0) = 0.25;
                                
                                ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, dum, 'clim', ZLIM, 'box', 'no', ...
                                    'highlight', mask, 'highlightstyle', 'opacity'); % opacity
                            else
                                ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, squeeze(dum(isub,:,:)), 'clim', ZLIM, 'box', 'no' ); % opacity
                            end
                            
                            %                         XLIM = [respavg.time{itrig}(1) respavg.time{itrig}(end)];
                            XLIM = [-0.2 0.9];
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
                            
                            %                         if iplot == 1
                            %                             title( sprintf('%s, %s, psc Lo: %g, Hi: %g\n%s', respavg.sdt_conds{istim, iresp}, respavg.sens.leg{isoi}, ZLIMS(:,2), respavg.behav_conds{icond} ))
                            pl.XTickLabel = [];
                            %                         else
                            %                             pl.Position(2) = pl.Position(2)*1.95;
                            %                             pl.Position(4) = pl.Position(4)*1.5;
                            %                             xlabel(sprintf('Time from %s (s)', respavg.trigger_leg{itrig} ));
                            xlabel('Time from trial onset (s)');
                            
%                             title( sprintf('%s', respavg.sdt_conds{istim, iresp} ))
                            title( sprintf('Subj %d', isub ))
                            %                         end
                            
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
                                                    h.Position(1) = 0.95
                                                    h.Position(2) = 0.5
                            %                         h.Ticks = [h.Limits(1) 0 h.Limits(2)];
                            %                         % Specify pos as a four-element vector of the form [x y w h] in data units.
                            %                         if iplot == 1
                            %                             boxlocs = [0.2 59 0.2 40; ...
                            %                                        0.2 42 0.2 16];
                            %                         else
                            %                             boxlocs = [0.25 11 0.2 11; ...
                            %                                        0.25 24 0.2 3];
                            %                         end
                            %                         for ib = 1:2
                            % %                             if ib < 3; subplot(1,2,1); else subplot(1,2,2); end
                            %                             r=rectangle('Position', boxlocs(ib,:));
                            %                             r.LineStyle = '--';
                            %                         end
                            
                        end %ises
                    end %ises
                end
            end
            set(gcf, 'Color', 'white')
        end
    end %ises
end

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(respavg.PREOUT, 'poolings');
    
    mkdir(fullfile(outpath, 'stats'))
    %                 if showstats
    %                     %                     outfile = fullfile(outpath, 'stats',   sprintf( 'TFRstats_%s_%s_%s', respavg.sens.leg{isoi}, respavg.sdt_conds{istim, iresp}));
    %                 else
    outfile = fullfile(outpath, sprintf( 'TFRstatssinglesub_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    %                 end
    disp(outfile)
    export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises




%% test in topo where gamma is significant

iband = 2;
itrig = 1;
    
freq=[];
freq.label = respavg.label;
freq.dimord = 'subj_chan_freq_time';
freq.time = respavg.time{itrig};
freq.freq = respavg.freq{iband};
dum = squeeze(respavg.dat(:, :, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
    iband, itrig, 4,icond, istim, iresp)); %average over sens
%                             freq.avg = squeeze(mean(dum(:,mask),2));
freq.powspctrm = dum;
          
itop = 1;
cfg = [];
cfg.latency = times(itop,:);
cfg.frequency = freqs(itop,:);
cfg.avgovertime = 'yes';
cfg.avgoverfreq = 'yes';

freq = ft_selectdata(cfg, freq);
freq.dimord = 'subj_chan'

nsub = length(respavg.SUBJ);

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; % for evoked data
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% prepare_neighbours determines what sensors may form clusters
cfg0_neighb = [];
cfg0_neighb.method    = 'template';
cfg0_neighb.template  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20161220/template/neighbours/elec1010_neighb.mat';
cfg0_neighb.elecfile  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20161220/template/electrode/standard_1020.elc';
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

freqzero=freq;
freqzero.powspctrm = zeros(size(freq.powspctrm));
stat = ft_freqstatistics(cfg, freq, freqzero);



%%                         
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


%% plot topo's based on stat.mask field
% close all
load('colormap_jetgray', 'cmap');
SAV=1;
ises = 4;
for istim = 1%:2
    for iresp = 1%:2        %1:2;
        for isoi = 1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
            
            figure;    iplot=0; hold on
            set(gcf, 'Position', [0 -200 1000 500])
            set(0,'DefaultAxesFontsize',30)
            
            for iband = 2% 2:-1:1 % 1:2 %
                for icond = 3 %[1,3, 4]  %[1:2, 4] %1:4 %1:4%  % 1:2
                    for itrig = 1%:2
                        
                        ZLIMS = [-0.4 0.4; -0.05 0.05];
%                         if icond == 4; ZLIMS = ZLIMS/4; end
                        ZLIM = ZLIMS(iband,:);
                        
                        iplot = iplot+1;
                        pl = subplot(1,1,iplot); hold on
%                         if iplot == 1 || iplot == 4 %~mod(iplot,2)
%                             ylabel('Frequency (Hz)');
%                         else
%                         end
%                         colormap(cmap)
                        
                        freq=[];
                        freq.label = respavg.label;
                        freq.dimord = 'chan';
                        dum = squeeze(mean(respavg.dat(:, :, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                            iband, itrig, 4,icond, istim, iresp), 1)); %average over sens
                        mask = squeeze( poolstat{isoi, iband, itrig, icond, istim, iresp}.mask);
                        % inspect labelmat to see indexing:
%                         figure; imagesc(poolstat{isoi, iband, itrig,
%                         icond, istim, iresp}.time, poolstat{isoi, iband,
%                         itrig, icond, istim, iresp}.freq,
%                         squeeze(poolstat{isoi, iband, itrig, icond,
%                         istim, iresp}.posclusterslabelmat)); colorbar
%                         mask = squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.posclusterslabelmat == 2);
                        freq.avg = squeeze(mean(dum(:,mask),2));
                      
                        cfg = [];
                        cfg.layout = 'elec1010.lay';
                        cfg.comment = 'no';
                        cfg.marker = 'on';
                        cfg.shading = 'flat';
                        cfg.style = 'straight'; %both  straight
                        cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
                        cfg.markersize = 6;
                        cfg.highlight = 'on';
                        cfg.highlightsymbol = '.';
                        cfg.highlightsize = 30;
                        cfg.highlightchannel = respavg.sens.ind{1};
                        cfg.colormap = cmap;
                        cfg.zlim = ZLIM;
                        ft_topoplotTFR(cfg, freq);

                        colorbar

                    end %ises
                end
                
                set(gcf, 'Color', 'white')
            end
            
            if SAV
                %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
                outpath = fullfile(respavg.PREOUT, 'poolings');
                
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
        end
    end
end %ises

%% Plot spectrum of high-low crit power
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
cfg.neighbours       = []; %in case no channel data present

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

close all
SAV=1;
ises = 4;
YLIMS = [5; 2.5];
for istim = 1%:2
    for iresp = 1%:2        %1:2;
        for isoi = 1%:2%:2%:4%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
            
            figure;    iplot=0; hold on
            set(gcf, 'Position', [500 500 800 300])
            set(0,'DefaultAxesFontsize',12)
            stat = {};
            for iband = 1:2% 2:-1:1 % 1:2 %
                iplot = iplot + 1;
                subplot(1,2,iplot)
                for icond = 4 %[1,3, 4]  %[1:2, 4] %1:4 %1:4%  % 1:2
                    for itrig = 1%:2
                        dum = squeeze(mean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{iband}), 1:length(respavg.time{itrig}), ...
                            iband, itrig, 4,icond, istim, iresp), 2)); %average over sens
                        % take 0.16 - 0.4 s (stim response)
                        TIM = [0.15 0.3];
                        tinds = respavg.time{itrig} >= TIM(1) & respavg.time{itrig} <= TIM(2);
                        dum = mean(dum(:,:,tinds),3);
                        dumsem = std(dum * 100) / sqrt(nsub);
%                         dum = squeeze(mean(dum));
                        %                         % use mask:
                        % %                         mask = squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.posclusterslabelmat == 1);
                        %                         mask = squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.mask);
                        %                         dum(~mask) = NaN;
                        
                        freq = [];
                        freq.label = {'occ'};
%                         freq.freq = [respavg.freq{1} respavg.freq{2}];
                        freq.freq = respavg.freq{iband};
                        freq.dimord = 'chan_subj_freq';
                        freq.powspctrm = shiftdim(dum, -1);
                        
                        freqzero = freq; %create zero freq to test against
                        freqzero.powspctrm = zeros(size(freq.powspctrm));
                        stat{iband} = ft_freqstatistics(cfg, freq, freqzero);
 

                        
                        
%                         plot(respavg.freq{iband}, nanmean(dum,2) * 100);
                        shadedErrorBar(respavg.freq{iband}, mean(dum) * 100, dumsem);
                        xlim([ respavg.freq{iband}(1) respavg.freq{iband}(end)])
                        ylim([-YLIMS(iband) YLIMS(iband)])
                        line(get(gca, 'XLIM'), [0 0])
                        title(TIM)
                        
                        %                         shadedErrorBar
                        
                    end %ises
                end
            end
        end %ises
    end
end %ises

