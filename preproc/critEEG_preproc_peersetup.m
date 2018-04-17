% MIBexp_preproc_peersetup
% run from runMIBmeg_analysis

cfg1 = {};
cfg2 = {};
outputfile = {};
overwrite = runcfg.overwrite;

load('all_events_unbalanced.mat')
eventall = event; % JJ's file

%make cells for each subject, to analyze in parallel
ctr = 0;
for ibatch = 1:length(runcfg.batchlists)
    eval(runcfg.batchlists{ibatch}); %load in batchlist file, batch, PREOUT and PREIN come out
    for irun=1:length(batch)
        if isempty(batch(irun).dataset) % probably commented out
            continue
        end
        PREIN = fullfile(basepath, 'projectdata', 'critEEG', 'data', PRE);
        session_name = ['ses' int2str(batch(irun).sessionno)];
        PREOUT = fullfile(basepath, 'projectdata', 'critEEG', 'preproc', batch(irun).subj, session_name, filesep);
        
        outfile = sprintf('%s%s_ses%d_%s_%s_run%d_%s_data', PREOUT, batch(irun).subj, batch(irun).sessionno, batch(irun).sessiondate, ...
            batch(irun).type, batch(irun).exp, 'stim'); % runcfg.trigger{itrg} % ALWAYS stim, see trialfun
        if ~exist([outfile '.mat'], 'file') || overwrite % if the matfile does not yet exist, then add to the joblist
            ctr = ctr + 1;
            outputfile{ctr} = outfile;
            cfg1{ctr}.runcfg = runcfg;  %analysis specifics
            cfg1{ctr}.runcfg.PRE = PRE;  %analysis specifics
            cfg1{ctr}.runcfg.batch = batch(irun);  %analysis specifics
            
            cfg1{ctr}.dataset= fullfile(PREIN, batch(irun).dataset);
            cfg1{ctr}.channel = {'EEG' 'EXG*'}; % 1-2 horizontal, 3-4 vertical, 5-6 earlobes
            cfg1{ctr}.fsample = 256; %
            
            cfg1{ctr}.reref = 'yes';
            cfg1{ctr}.refchannel = { 'EXG5', 'EXG6' }; % earlobes

%             Biosemi: Because of BioSemi's unique setup with active electrodes,
%             battery power supply, and optic fiber signal link, 50 Hz interference 
%             is eliminated from the recordings. Therefore, reverting to signal 
%             corrupting patches like 50 Hz notch filters is never needed.
%             cfg1{ctr}.dftfilter = 'yes';
%             cfg1{ctr}.dftfreq = [49:0.1:51, 100];  % line noise removal using discrete fourier transform
%             cfg1{ctr}.padding   = 10;
            cfg1{ctr}.padding   = 0;
            
            cfg1{ctr}.continuous = 'yes';
%             cfg1{ctr}.demean = 'yes'; 
            cfg1{ctr}.detrend = 'yes'; % not good for timedomain analyses, good for gamma
            
            cfg1{ctr}.trialfun = 'sortTrials_critEEG';
            cfg1{ctr}.trialdef.begtim = -1;  % before stim onset
            cfg1{ctr}.trialdef.endtim = 1.25; % after report or stim if no report
            
            % check if event is present in JJ's file
            [~, runname] = fileparts(cfg1{ctr}.dataset);
            try
                cfg1{ctr}.event = eventall.(batch(irun).subj).(runname);
            catch
                fprintf('%s event not found\n', runname)
                ctr = ctr - 1;
                continue
            end

            cfg1{ctr}.artfrej = runcfg.preproc.artf_rejection; %do artifact rejection? yes
            cfg1{ctr}.artf_feedback = runcfg.preproc.artf_feedback; %feedback for inspection automatic artf detection
            cfg1{ctr}.loadartf = runcfg.preproc.loadartf;%load from file?
            %automatic artf detection: cfg specified in resp scripts
            cfg1{ctr}.musclethr = batch(irun).musclethr;
            cfg1{ctr}.eogverthr = batch(irun).eogverthr;
            cfg1{ctr}.eoghorthr = batch(irun).eoghorthr;
            
            % visual artifact rejection parameters
            cfg2{ctr}.method   = 'summary'; % channel trial summary
            cfg2{ctr}.channel = 'MEG';
            cfg2{ctr}.alim     = 1e-10;
            cfg2{ctr}.megscale = 1;
            cfg2{ctr}.eogscale = 5e-8;
            cfg2{ctr}.layout = 'biosemi64incI1I2.lay';
            
        else
            fprintf('%s exists! Skip it\n', outfile)
        end   %%% if ~exist([outfile '.mat'], 'file')
    end
    clear batch
end
% cfg1 = cfg1(1), cfg2 = cfg2(1), outputfile = outputfile(1) % for testing

fprintf('Running critEEG_preproc for %d cfgs\n', length(cfg1))

switch runcfg.preproc.parallel
    case 'local'
        cellfun(@critEEG_preproc, cfg1, cfg2, outputfile);
    case 'peer'
        peercellfun(@critEEG_preproc, cfg1, cfg2, outputfile);
    case 'torque'
        timreq = 15; %in minutes per run
        memreq = 2000; % in MB  memreq*1024^3  3 is added anyway

%         timreq = 2; %in minutes per run
        setenv('TORQUEHOME', 'yes')
        mkdir('~/qsub'); cd('~/qsub');
        switch runcfg.preproc.compile
            case 'no'
                nnodes = 30; % how many licenses available?
                stack = round(length(cfg1)/nnodes); % only used when not compiling 
%                 qsubcellfun(@MIBexp_preproc, cfg1, cfg2, cfg3, outputfile, 'memreq', memreq*1024^3, 'timreq', timreq*60, ...
%                     'stack', stack, 'StopOnError', true, 'backend', runcfg.preproc.parallel);

                qsubcellfun(@critEEG_preproc, cfg1, cfg2, outputfile, 'memreq', memreq, 'timreq', timreq*60, ...
                    'stack', stack, 'StopOnError', true, 'backend', runcfg.preproc.parallel, 'options', '-l nodes=1:ppn=1');
                
            case 'yes'
                compiledfun = qsubcompile({@critEEG_preproc @sortTrials_critEEG}, 'toolbox', {'signal', 'stats'});
                qsubcellfun(compiledfun, cfg1, cfg2, outputfile, 'memreq', memreq, 'timreq', timreq*60, ...
                    'stack', 1, 'StopOnError', true, 'backend', runcfg.preproc.parallel, 'options', '-l nodes=1:ppn=1');
        end

        
    case 'parfor'
        parfor ibatch = 1:length(cfg1)
            critEEG_preproc(cfg1{ibatch}, cfg2{ibatch}, outputfile{ibatch})
        end
    otherwise
        error('Unknown backend, aborting . . .\n')
end

