
%
% behav=respavg.behavior;
behav=mseavg.behavior;
SUBJ = mseavg.SUBJ;
behav.condleg = {'Conservative', 'Liberal', '', 'Cons-lib'};

% plot criterion
condcolors = cbrewer('qual', 'Set1',3);
condcolors = [condcolors(2,:); condcolors(1,:); condcolors(3,:)];

SAV = 1;
close all
set(0, 'DefaultAxesFontsize', 6)
figh = figure; set(gcf, 'Position', [ 2055   320 600 150]) % 240 px wide
% figh = figure; set(gcf, 'Position', [ 2055   320 240 240]) % 240 px wide
iplot = 0;
iconds = 2:-1:1;
ncols = 3;
for im =1:length(behav.behav_measures)
    iplot = iplot + 1;
    pl = subplot(1,ncols,iplot);
    h = barweb( squeeze(nanmean(behav.(behav.behav_measures{im})(:,4,iconds))), squeeze(nanstd(behav.(behav.behav_measures{im})(:,4,iconds))/sqrt(length(SUBJ))), ...
        0.5, [], sprintf('%s: p = %g p = %g p = %g', behav.behav_measures{im}, behav.stat.(behav.behav_measures{im})));
    ylabel(behav.behav_measures{im}, 'Fontsize', 7)
    %     axis square
    %     if im==2; legend(behav.condleg(iconds), 'Location', 'NorthEast'); legend boxoff; end
    shading flat
    h.bars(1).FaceColor = condcolors(1,:);
    h.bars(2).FaceColor = condcolors(2,:);
    h.bars(1).EdgeAlpha = 0;
    h.bars(2).EdgeAlpha = 0;
    
    ax=gca;
    ax.Position(3) = ax.Position(3) * 0.6;
    %     ax.YTick = [-1:0.25:1];
    
    % frequentist stats
    [h,pconsv0] = ttest(behav.criterion(:, 4, 1));  % (isub, ises, icond)
    [h,plibvs0] = ttest(behav.criterion(:, 4, 2)) ; % (isub, ises, icond)
    [h,pconsvslib] = ttest(behav.criterion(:, 4, 1), behav.criterion(:, 4, 2)) ; % (isub, ises, icond)
end
if SAV
    disp(behav.PREOUT)
    %     export_fig(fullfile(behav.PREOUT, 'Fig_behavior'), '-pdf')
    print(fullfile(behav.PREOUT, 'Fig_behavior_c_dprime'), '-dpdf') % , '-fillpage'
end

% print sourcedat for Figure 2
sourcedat = [squeeze(behav.hitrate(:,4,1:2)) squeeze(behav.farate(:,4,1:2)) squeeze(behav.dprime(:,4,1:2)) squeeze(behav.criterion(:,4,1:2)) ]


%%

% plot H and FA separately
f = figure; iplot=0;
f.Position = [  680   955   800   150]
iplot = iplot + 1;  subplot(1,ncols,iplot)
bardat = [ squeeze(nanmean(behav.hitrate(:,4,iconds))) squeeze(nanmean(behav.farate(:,4,iconds))) ...
    squeeze(nanmean(behav.hitrate(:,4,iconds))) - squeeze(nanmean(behav.farate(:,4,iconds))) ];
semdat = [ squeeze(nanstd(behav.hitrate(:,4,iconds))/sqrt(length(SUBJ))) squeeze(nanstd(behav.farate(:,4,iconds))/sqrt(length(SUBJ))) ...
    squeeze(nanstd(behav.hitrate(:,4,iconds))/sqrt(length(SUBJ))) - squeeze(nanstd(behav.farate(:,4,iconds))/sqrt(length(SUBJ))) ];

% barweb( bardat, semdat, 0.5, behav.condleg )
h = barweb( bardat', semdat', 0.75, {'Hits' 'False alarms' 'Hits-False alarms'});
% legend({'Hitrate' 'FArate' 'Hit-FArate' }, 'Location', 'NorthWest'); legend boxoff;
% legend(behav.condleg(iconds), 'Location', 'NorthEast'); legend boxoff;
h.bars(1).FaceColor = condcolors(1,:);
h.bars(2).FaceColor = condcolors(2,:);
h.bars(1).EdgeAlpha = 0;
h.bars(2).EdgeAlpha = 0;

ax=gca;
ax.YTick = [0:0.25:1];
ax.Position(3) = ax.Position(3) * 0.6;

% p=[];
% difference = squeeze(behav.hitrate(:,4,1)) - squeeze(behav.hitrate(:,4,2)); % difference between H and FA
% p(1)= randtest1d( difference, zeros(size(difference(:,1))), 0, 1000);
%
% difference = squeeze(behav.farate(:,4,1)) - squeeze(behav.farate(:,4,2)); % difference between H and FA
% p(2)= randtest1d( difference, zeros(size(difference(:,1))), 0, 1000);
%
% difference = squeeze(behav.hitrate(:,4,iconds)) - squeeze(behav.farate(:,4,iconds)); % difference between H and FA
% p(3)= randtest1d( difference(:,1) - difference(:,2), zeros(size(difference(:,1))), 0, 1000);

[p_hit] = permtest( squeeze(behav.hitrate(:,4,1)) , squeeze(behav.hitrate(:,4,2)))

difference = squeeze(behav.farate(:,4,1)) - squeeze(behav.farate(:,4,2)); % difference between H and FA
p(2)= permtest( difference, zeros(size(difference(:,1))));

difference = squeeze(behav.hitrate(:,4,iconds)) - squeeze(behav.farate(:,4,iconds)); % difference between H and FA
p(3)= permtest( difference(:,1) , difference(:,2));

title(sprintf('pvals: %g %g %g', p ))
% axis square
ylabel('Proportion of reports'); ylim([0 1])
% sigstar([1,2])  % TODO add stars
shading flat
legend(behav.condleg(iconds), 'Location', 'NorthEast'); legend boxoff;

% % plot RT distributions
% iplot = iplot + 1;  subplot(1,ncols,iplot); hold on
% stimleg = {' Hit', ' False alarm'};
% stimlines = {'-', '--'}; % fig hom
% % condcolors = ['b', 'r']; % {'NoMiss', 'NoFA'}
% clear leg
% for istim=1:2
%     for icond = 1:2
% %         handle(istim, icond) = shadedErrorBar(behav.RT_edges, squeeze(mean(behav.RThist(:,4,istim,icond,:))), ...
% %             squeeze(std(behav.RThist(:,4,istim,icond,:))) /sqrt(length(SUBJ)), [stimlines{istim} condcolors(icond,:)], 1 );
%         handle(istim, icond) = shadedErrorBar(behav.RT_edges-0.16, squeeze(mean(behav.RThist(:,4,istim,icond,:))), ...
%             squeeze(std(behav.RThist(:,4,istim,icond,:))) /sqrt(length(SUBJ)), {'Linestyle', stimlines{istim}, 'Color', condcolors(icond,:)}, 0 );
%         leg{istim,icond} = [behav.condleg{icond} stimleg{istim}];
%     end
% end
% ax=gca;
% ax.YTick = [0:0.05:0.15];
% ax.XTick = [0:0.25:1];
% % xlabel('Time from train onset (s)', 'Fontsize', 10);
% xlabel('Time from (non-)target (s)', 'Fontsize', 10);
% ylabel('Proportion of reaction times', 'Fontsize', 10);
% xlim([behav.RT_edges(1) behav.RT_edges(end)] )
% legend([handle.mainLine], leg(:), 'Location', 'NorthEast'); legend boxoff
% xlim([0 1])
% ylim([-0.01 0.125])

% plot bars RT H and FA cons and lib
% plot H and FA separately
iplot = iplot + 1;
subplot(1,ncols,iplot)
bardat = [ squeeze(mean(behav.RT(:,4,1,iconds))) squeeze(mean(behav.RT(:,4,2,iconds))) ];
semdat = [ squeeze(std(behav.RT(:,4,1,iconds))/sqrt(length(SUBJ))) squeeze(std(behav.RT(:,4,1,iconds))/sqrt(length(SUBJ))) ];

% barweb( bardat, semdat, 0.5, behav.condleg )
h = barweb( bardat', semdat', 0.75, {'Hits' 'False alarms'})
% legend({'Hitrate' 'FArate' 'Hit-FArate' }, 'Location', 'NorthWest'); legend boxoff;
h.bars(1).FaceColor = condcolors(1,:);
h.bars(2).FaceColor = condcolors(2,:);
h.bars(1).EdgeAlpha = 0;
h.bars(2).EdgeAlpha = 0;

ax=gca;
ax.YTick = [0:0.05:1];
ax.Position(3) = ax.Position(3) * 0.6;

p=[];
p(1) = permtest(behav.RT(:,4,1,1), behav.RT(:,4,1,2)); % Hit lib vs cons
p(2) = permtest(behav.RT(:,4,2,1), behav.RT(:,4,2,2)); % Hit lib vs cons

title(sprintf('pvals: %g %g', p ))
% axis square
ylabel('Reaction time (s)'); ylim([0 1])
% sigstar([1,2])  % TODO add stars
shading flat

ylim([0.3 0.5])



if SAV
    disp(behav.PREOUT)
    %     export_fig(fullfile(behav.PREOUT, 'Fig_behavior'), '-pdf')
    print(fullfile(behav.PREOUT, 'Fig_behavior'), '-dpdf') % , '-fillpage'
end
%% Plotting overview
set(0, 'DefaultAxesFontsize', 8)
figure; set(gcf, 'Position', [0 -200 210*4 375*3])
iconds = 1:2;
for im = 1:2%length(behav.behav_measures)
    subplot(3,3,im)
    barweb( squeeze(nanmean(behav.(behav.behav_measures{im})(:,4,iconds))), squeeze(nanstd(behav.(behav.behav_measures{im})(:,4,iconds))/sqrt(length(SUBJ))), 0.5, [], sprintf('%s: HivsLocrit p = %g', behav.behav_measures{im}, behav.stat.(behav.behav_measures{im})))
    ylabel(behav.behav_measures{im})
    axis square
    if im==2; legend(behav.condleg, 'Location', 'East'); legend boxoff; end
end

dprime = squeeze(behav.(behav.behav_measures{1})(:,4,iconds) )
% test overall crit against 0
% crit = squeeze(mean(behav.(behav.behav_measures{2})(:,4,iconds),3) )
crit = squeeze(behav.(behav.behav_measures{2})(:,4,iconds) )

%%
% plot RT distributions
stimleg = {'Fig', 'Hom'};
stimlines = {'-', '--'}; % fig hom
condcolors = ['b', 'r']; % {'NoMiss', 'NoFA'}
subplot(3,3,4); hold on
for istim=1:2
    for icond = 1:2
        handle(istim, icond) = shadedErrorBar(behav.RT_edges, squeeze(mean(behav.RThist(:,4,istim,icond,:))), squeeze(std(behav.RThist(:,4,istim,icond,:))) /sqrt(length(SUBJ)), [stimlines{istim} condcolors(icond)], 1 );
        leg{istim,icond} = [behav.condleg{icond} stimleg{istim}];
    end
end
title('RT distributions, Fig', 'Fontsize', 14)
xlabel('Time from train onset (s)'); ylabel('Fraction of reaction times');
xlim([behav.RT_edges(1) behav.RT_edges(end)] )
legend([handle.mainLine], leg(:)); legend boxoff
axis square

% plot mean RT
for istim = 1:2
    subplot(3,3,4+istim)
    barweb(squeeze( nanmean(behav.RT(:,:,istim,1:2))), squeeze(nanstd(behav.RT(:,:,istim,1:2))/sqrt(length(SUBJ))), 0.5, behav.sesleg,  sprintf('%s RT: HivsLocrit p = %g', stimleg{istim}, behav.stat.RT(istim)))
    %     legend(behav.condleg, 'Location', 'NorthEast'); legend boxoff; axis square
    %     title(sprintf(' %s reaction times', stimleg{is}))
    ylabel('Reaction time (s)')
    axis square
    ylim([0 0.9])
end

% figure; set(gcf, 'Position', [0 -200 210*4 375*3])
subplot(3,3,7)
barweb( squeeze(nanmean(behav.hitrate(:,:,1:2))), squeeze(nanstd(behav.hitrate(:,:,1:2))/sqrt(length(SUBJ))), 0.5, behav.sesleg )
ylabel('Proportion of responses'); ylim([0 1])
title('Hitrate')
axis square

subplot(3,3,8)
barweb( squeeze(nanmean(behav.farate(:,:,1:2))), squeeze(nanstd(behav.farate(:,:,1:2))/sqrt(length(SUBJ))), 0.5, behav.sesleg )
ylabel('Proportion of responses'); ylim([0 1])
title('False alarm rate')
legend(behav.condleg, 'Location', 'NorthEast'); legend boxoff
axis square

subplot(3,3,9)
barweb( [ squeeze(nanmean(behav.hitrate(:,4,1:2))) squeeze(nanmean(behav.farate(:,4,1:2))) ...
    squeeze(nanmean(behav.hitrate(:,4,1:2))) - squeeze(nanmean(behav.farate(:,4,1:2))) ], ...
    [ squeeze(nanstd(behav.hitrate(:,4,1:2))/sqrt(length(SUBJ))) squeeze(nanstd(behav.farate(:,4,1:2))/sqrt(length(SUBJ))) ...
    squeeze(nanstd(behav.hitrate(:,4,1:2))/sqrt(length(SUBJ))) - squeeze(nanstd(behav.farate(:,4,1:2))/sqrt(length(SUBJ))) ] ...
    , 0.5, behav.condleg )
ylabel('Proportion of responses'); ylim([0 1])
legend({'Hitrate' 'FArate' 'Hit-FArate' }, 'Location', 'NorthWest'); legend boxoff;
difference = squeeze(behav.hitrate(:,4,1:2)) - squeeze(behav.farate(:,4,1:2)); % difference between H and FA
title(sprintf('Hit-FArate HivsLocrit: p = %g', randtest1d( difference(:,1) - difference(:,2), zeros(size(difference(:,1))), 0, 1000)))
axis square

set(gcf, 'Color', 'white')
%
disp(behav.PREOUT)
export_fig(fullfile(behav.PREOUT, 'behavior'), '-pdf')

%% Correlate behavioral measures with each other
% c vs d, c vs RT, d vs RT high crit
% c vs d, c vs RT, d vs RT low crit
% x_dv = {'criterion', 'criterion', 'dprime'};
% y_dv = {'dprime', 'RT', 'RT'};
x_dv = {'criterion'};
% y_dv = {'dprime'}; % hitminusfarate
y_dv = {'hitminusfarate'}; % hitminusfarate

behav.condleg = {'Conservative'  'Liberal'  'Condition average'  'Lib ? Cons'};

% close all
figure; iplot=0;
% set(gcf, 'Position', [0 -200 210*4 375*3])
set(gcf, 'Position', [0 -200 375*3 210*4 ])
for idv = 1:length(x_dv)
    for icond = 4 %1:4
        iplot=iplot+1;
        subplot(3,4,iplot);hold on; axis square; box on
        %     scatter(squeeze(behav.(x_dv{idv})(:,4,icond) , squeeze(behav.dprime(:,4,icond)) )   % isub, ises, icond
        datx = squeeze(behav.(x_dv{idv})(:,4,icond));
        if idv == 1
            daty = squeeze(behav.(y_dv{idv})(:,4,icond));
        else
            daty = squeeze(behav.(y_dv{idv})(:,4,3,icond)); % dimord subj_ses_stim_cond  take all stim for RT
        end
        
        scatter(datx, daty, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 200 )   % isub, ises, icond
        [rho, pval] = corr(datx, daty, 'type', 'Pearson' );
        %         [pval, rho] = permuteCorr(datx, daty, 10000, 0);
        
        title(sprintf('%s, r = %g, p = %g', behav.condleg{icond}, rho, pval))
        if pval < 0.05
            h = lsline;
            h.Color = 'k';
            h.LineWidth = 1.5;
        elseif pval < 0.01
            h = lsline;
            h.Color = 'k';
            h.LineWidth = 3;
        end
        xlabel(x_dv{idv})
        ylabel(y_dv{idv})
    end
end

export_fig(fullfile(behav.PREOUT, 'behavior_correlations'), '-pdf')

cd(behav.PREOUT)

%% plot RT distributions on their own
close all
f = figure;
f.Position = [ 680   952   250   150];
stimleg = {'Fig', 'Hom'};
stimlines = {'-', '--'}; % fig hom
condcolors = ['r', 'b']; % {'NoMiss', 'NoFA'}
condcolors = cbrewer('qual', 'Set1',3);

% subplot(3,3,4);
hold on
colcode = 0;
for istim=1:2
    for icond = 1:2
        colcode = colcode + 0.25;
        handle(istim, icond) = shadedErrorBar(behav.RT_edges, squeeze(mean(behav.RThist(:,4,istim,icond,:))), ...
            squeeze(std(behav.RThist(:,4,istim,icond,:))) /sqrt(length(SUBJ)), [], 0 );
        handle(istim, icond).mainLine.Color =  condcolors(icond,:);
        handle(istim, icond).patch.FaceColor =  [colcode colcode colcode];
        if istim == 2
            handle(istim, icond).mainLine.LineStyle = '--';
        end
        leg{istim,icond} = [behav.condleg{icond} ' ' stimleg{istim}];
    end
end
ax=gca;
ax.FontSize = 8;

xlabel('Time from target onset (s)');
ylabel('Fraction of reaction times');
xlim([behav.RT_edges(1) behav.RT_edges(end)] )
legend([handle.mainLine], leg(:), 'Location', 'Eastoutside'); legend boxoff
axis square

export_fig(fullfile(behav.PREOUT, 'RT_distribution'), '-pdf')
% print(fullfile(behav.PREOUT, 'RT_distribution'), '-depsc2')
print(fullfile(behav.PREOUT, 'RT_distribution'), '-dpdf')
cd(behav.PREOUT)

%% c_conservative - c_liberal > d?_conservative - d?_liberal (waarbij c = -.5* (Z(HR) + Z(FAR))  en d? = (Z(HR) - Z(FAR)))?

% critdiff = behav.criterion(:,4,2) - behav.criterion(:,4,1)
% dprimediff = behav.dprime(:,4,2) - behav.dprime(:,4,1)

% critdiff = zscore(behav.criterion(:,4,2)) - zscore(behav.criterion(:,4,1))
% dprimediff = zscore(behav.dprime(:,4,2)) - zscore(behav.dprime(:,4,1))
%
% critdiff = zscore(behav.criterion(:,4,2) - behav.criterion(:,4,1))
% dprimediff = zscore(behav.dprime(:,4,2) - behav.dprime(:,4,1))

critdiff = (behav.criterion(:,4,2) - behav.criterion(:,4,1)) ./ behav.criterion(:,4,1) * 100;
dprimediff = (behav.dprime(:,4,2) - behav.dprime(:,4,1)) ./ behav.dprime(:,4,1) * 100;

permtest(critdiff, dprimediff)
figure; bar([mean(critdiff) mean(dprimediff)])
[mean(critdiff) mean(dprimediff)]
[std(critdiff) std(dprimediff)]
%%
RTcons = behav.RT(:,4,3,1);
RTlib = behav.RT(:,4,3,2);
figure; bar([mean(RTcons) mean(RTlib)])
permtest(RTcons, RTlib)


%% Entropy paper: plot c and d' distributions
set(0, 'DefaultAxesFontsize', 10)

close all
SAV = 1;

condcolors = {'r', 'b', 'k', 'k'}; % {'NoMiss', 'NoFA'}
% condcolors = cbrewer('qual', 'Set1',3);
nrow = 3;
behav.condleg = {'Conservative'  'Liberal'  'Condition average'  'Liberal ? Conservative'};
nsub = length(mseavg.SUBJ);

f = figure;
f.Position = [          2081         463        1000         750];

iplot=0;
% behavmeas = {'criterion', 'dprime', 'hitminusfarate'};
% plotname = {'Bias', 'Sensitivity', 'hitminusfarate'};

behavmeas = {'criterion', 'dprime'}; 
% behavmeas = {'criterion', 'hitminusfarate'}; 
% behavmeas = {'hitplusfarate', 'hitminusfarate'};
% behavmeas = {'prop_yes', 'hitminusfarate'};
% behavmeas = {'Bprime', 'Aprime'};

plotname = {'Bias', 'Sensitivity'};

for ib = 1:2 %[1,2,4] %1:2 
    %histogram lib and cons
    iplot=iplot+1;    subplot(nrow,3,iplot); hold on
    temp = behav.(behavmeas{ib})(:,4,[1,2]);
    leftedge = min(temp(:));
    rightedge = max(temp(:));
    stepsize = (max(temp(:)) - min(temp(:))) / 10;
    edges = leftedge:stepsize:rightedge;
    for icond = [1,2]      
        h = histogram(behav.(behavmeas{ib})(:,4,icond), edges);
        h.FaceColor = condcolors{icond};
%         h.FaceColor = condcolors(icond,:);
    end
    title(sprintf('%s',  plotname{ib}))
    xlabel(behavmeas{ib})
    ylabel('N participants')
    if iplot == 1
        legend(behav.condleg(1:2), 'Location', 'Northwest'); legend boxoff
    end
    ax=gca;
    l = line([0,0], [0 ax.YLim(2)]);
    l.LineStyle = '--';
    l.Color = 'k';
    
    % scatter lib vs cons
    iplot=iplot+1; subplot(nrow,3,iplot);hold on; axis square; box on
    scatter(behav.(behavmeas{ib})(:,4,1), behav.(behavmeas{ib})(:,4,2), 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 200 )   % isub, ises, icond
    corrtype = 'Pearson'; 
%     if ib == 1; corrtype = 'Spearman';end
    [r,p] = corr(behav.(behavmeas{ib})(:,4,1), behav.(behavmeas{ib})(:,4,2), 'type', corrtype);
    if p <0.05; lsline; end
%     title(sprintf('%s\nr = %1.2f, p = %g',  behavmeas{ib}, r, p))
    title(sprintf('%ss r = %1.2f, p = %g', corrtype, r, p))
    xlabel(sprintf('%s conservative', behavmeas{ib}))
    ylabel(sprintf('%s liberal', behavmeas{ib}))
     
    % slopes cons to lib
    iplot=iplot+1; subplot(nrow,3,iplot);hold on; axis square; box on
    ax=gca;
    ax.ColorOrder = parula(nsub); %cbrewer('qual', 'Set1',3);
%     ax.ColorOrder = cbrewer('qual', 'Set2', nsub);
%     pl = plot([ones(nsub,1) ones(nsub,1)+1]', [behav.(behavmeas{ib})(:,4,1) behav.(behavmeas{ib})(:,4,2)]', 'Linewidth', 2);
    pl = plot([ones(nsub,1) ones(nsub,1)+1]', [behav.(behavmeas{ib})(:,4,1) behav.(behavmeas{ib})(:,4,2)]', 'Color', 'k', 'Linewidth', 1.5  );
    xlim([0.85 2.15])
    ax.XTick = [1 2];
    ax.XTickLabel = behav.condleg(1:2);
    ylabel(behavmeas{ib})
    title('Single subject slopes')
end

% correlate cond avg c and d' and lib -cons
for icond = 3:4
    iplot=iplot+1; subplot(nrow,3,iplot);hold on; axis square; box on
    scatter(behav.(behavmeas{1})(:,4,icond), behav.(behavmeas{2})(:,4,icond), 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 200 )   % isub, ises, icond
    [r,p] = corr(behav.(behavmeas{1})(:,4,icond), behav.(behavmeas{2})(:,4,icond));
    if p <0.05; lsline; end
    %     title(sprintf('%s\nr = %1.2f, p = %g',  behavmeas{ib}, r, p))
    title(sprintf('%s\nr = %1.2f, p = %g', behav.condleg{icond}, r, p))
    xlabel(sprintf('%s', behavmeas{1}))
    ylabel(sprintf('%s', behavmeas{2}))
end

if SAV
%     export_fig(fullfile(behav.PREOUT, 'behavior_quenchingpaper'), '-pdf')
    export_fig(fullfile(behav.PREOUT, 'behavior_quenchingpaper'), '-png')
    f.Renderer = 'Painters';
    print('-depsc2', 'behavior_quenchingpaper')
%     saveas(gcf, 'behavior_quenchingpaper', 'epsc' )
    cd(behav.PREOUT)
end

%% export data for rmcorr
xdat = behav.criterion(:, 1:3, 1); % cons
ydat = behav.criterion(:, 1:3, 2); % lib
export_for_rmcorr(xdat, ydat, 'critcons', 'critlib', '/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/data_rmcorr');

xdat = behav.dprime(:, 1:3, 1); % cons
ydat = behav.dprime(:, 1:3, 2); % lib
export_for_rmcorr(xdat, ydat, 'dprimecons', 'dprimelib', '/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/data_rmcorr');

xdat = behav.criterion(:, 1:3, 3); % cons
ydat = behav.dprime(:, 1:3, 3); % lib
export_for_rmcorr(xdat, ydat, 'critavg', 'dprimeavg', '/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/data_rmcorr');

xdat = behav.criterion(:, 1:3, 4); % cons
ydat = behav.dprime(:, 1:3, 4); % lib
export_for_rmcorr(xdat, ydat, 'critdiff', 'dprimediff', '/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/data_rmcorr');

xdat = behav.criterion(:, 1:3, 1); 
ydat = behav.dprime(:, 1:3, 1); 
export_for_rmcorr(xdat, ydat, 'critcons', 'dprimecons', '/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/data_rmcorr');

xdat = behav.criterion(:, 1:3, 2); 
ydat = behav.dprime(:, 1:3, 2); 
export_for_rmcorr(xdat, ydat, 'critlib', 'dprimelib', '/Users/kloosterman/Dropbox/PROJECTS/CriterionEEG/data_rmcorr');

