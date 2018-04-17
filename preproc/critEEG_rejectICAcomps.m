function critEEG_rejectICAcomps(SUBJ)
% select components to be rejected COSTRAP!

basepath = '/path/'; % on the cluster

for isub = 1:length(SUBJ)
    for ises = 1:3
        PREIN = fullfile(basepath, 'preproc', SUBJ{isub}, sprintf('ses%d', ises));  
        PREOUT = fullfile(basepath, 'preproc', SUBJ{isub}, sprintf('ses%d', ises));  
        if ~exist(PREIN, 'dir')
            continue
        end
        cd(PREIN)
%         inputfile = sprintf('%s_*ses%d_*comp.mat', SUBJ{isub}, ises);
        inputfile = sprintf('%s_*ses%d_*costrap_comp.mat', SUBJ{isub}, ises);
        
        inputfile = dir(inputfile);
        fprintf('Loading %s from...\n %s\n', inputfile.name, PREIN)
        load(inputfile.name)
        
        %     visualize components
        cfg = [];
        
        cfg.component = [1:50];       % specify the component(s) that should be plotted
%         cfg.layout = 'biosemi64.lay'; % specify the layout file that should be used for plotting
        cfg.layout = 'elec1010.lay'; % specify the layout file that should be used for plotting
        cfg.comment   = 'no';
        figure('units','normalized','outerposition', [0.9995 0.0367 1 0.8775] )
        ft_topoplotIC(cfg, comp)
        %
        cfg = [];
%         cfg.layout = 'biosemi64_incI1I2.lay'; % specify the layout file that should be used for plotting
%         cfg.layout = 'biosemi64.lay'; % specify the layout file that should be used for plotting
        cfg.layout = 'elec1010.lay'; % specify the layout file that should be used for plotting
        cfg.viewmode = 'vertical'; % component
        cfg.channel = 1:50;       % specify the component(s) that should be plotted
        % figure('Position', [69 58 774 1045])
        ft_databrowser(cfg, comp)
        f = gcf;
        f.Position = [69 58 774 1045];
        
        manual = 1;
        if manual
            opt=[];
            opt.WindowStyle = 'Normal';
            x = inputdlg('Enter space-separated components to be rejected:',...
                'ICA rejection', [1 50], {''}, opt);
            comps2reject = str2num(x{:});
            save comps2reject comps2reject
        else
            load comps2reject
        end
        close all
        
        cfg = [];
        cfg.component = comps2reject; % to be removed component(s)
        data_dirty = data;
        data = ft_rejectcomponent(cfg, comp, data);
        
        outputfile = sprintf('%s_ses%d_data_costrap_clean',  SUBJ{isub}, ises); %2 also remove muscle comps: focal topography
        
        fprintf('Saving %s to...\n %s\n', outputfile, PREOUT)
        save(outputfile, 'data');
        
        % % % % visualize result
        close all
        trialind = 1;
        figure; subplot(2,1,1)
        plot(data.time{trialind}, data.trial{trialind})
        subplot(2,1,2)
        plot(data_dirty.time{trialind}, data_dirty.trial{trialind})
        clear data_dirty data comp
    end
end


