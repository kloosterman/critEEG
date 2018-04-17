function critEEG_rejectEOGartf(SUBJ)
% reject trials with blinks during stim

% parallel or not
compile = 'no';

parallel = 'local';  % torque
basepath = '/path/'; % on the cluster

cfg = {};
ctr=0;
for isub = 1:length(SUBJ)
    for ises= 1:3
        
        PREIN = fullfile(basepath, 'preproc', SUBJ{isub}, sprintf('ses%d', ises));
        if ~exist(PREIN, 'dir')
            continue
        end

        ctr = ctr + 1;
        cfg{ctr} = [];
        cfg{ctr}.PREIN = PREIN;
    end
end

fprintf('Running rejectEOGandCSD for %d cfgs\n', length(cfg))

switch parallel
    case 'local'
        cellfun(@rejectEOGandCSD, cfg(:));
    case 'peer'
        peercellfun(@rejectEOGandCSD, cfg(:));
    case {'torque' 'qsublocal'}
        timreq = 10*60; %in sec per run
        memreq = 8000; % in MB  memreq*1024^3  3 is added anyway

%         timreq = 2; %in minutes per run
        setenv('TORQUEHOME', 'yes')
        mkdir('~/qsub'); cd('~/qsub');
        switch compile
            case 'no'
                nnodes = 10000; % how many licenses available?
                stack = round(length(cfg(:))/nnodes); % only used when not compiling 

                qsubcellfun(@rejectEOGandCSD, cfg(:), 'memreq', memreq, 'timreq', timreq*60, ...
                    'stack', stack, 'StopOnError', true, 'backend', parallel, 'options', '-l nodes=1:ppn=1');
                
            case 'yes'
                compiledfun = qsubcompile(@rejectEOGandCSD, 'toolbox', {'signal', 'stats'});
                qsubcellfun(compiledfun, cfg(:), 'memreq', memreq, 'timreq', timreq*60, ...
                    'stack', 1, 'StopOnError', false, 'backend', parallel, 'options', '-l nodes=1:ppn=1');
        end

        
    case 'parfor'
%         parfor ibatch = 1:length(cfg(:))
%             rejectEOGandCSD(cfg1{ibatch}, cfg2{ibatch}, outputfile{ibatch})
%         end
    otherwise
        error('Unknown backend, aborting . . .\n')
end


