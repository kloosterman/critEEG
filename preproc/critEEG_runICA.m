function critEEG_runICA(SUBJ)

% parallel or not
if ismac
    parallel = 'local';
else
    parallel = 'torque';
end
compile = 'no';

basepath = '/path/'; % on the cluster

cfg = {};
ctr=0;
for isub = 1:length(SUBJ)
    for ises=1:3
        PREIN = fullfile(basepath, 'preproc', SUBJ{isub}, sprintf('ses%d', (ises)));
        PREOUT = fullfile(basepath, 'preproc', SUBJ{isub}, sprintf('ses%d', (ises)));
        if exist(PREIN, 'dir')
            ctr = ctr + 1;
            cfg{ctr} = [];
            cfg{ctr}.PREIN = PREIN;
            cfg{ctr}.PREOUT = PREOUT;
        else
            fprintf('%s dir not found\n', PREIN)
            continue
        end
    end
end

fprintf('Running doICA for %d cfgs\n', length(cfg))

switch parallel
    case 'local'
        cellfun(@doICA, cfg(:));
    case 'peer'
        peercellfun(@doICA, cfg(:));
    case {'torque' 'qsublocal'}
        timreq = 30; %in minutes per run
        memreq = 8000; % in MB  memreq*1024^3  3 is added anyway

%         timreq = 2; %in minutes per run
        setenv('TORQUEHOME', 'yes')
        mkdir('~/qsub'); cd('~/qsub');
        switch compile
            case 'no'
                nnodes = 30; % how many licenses available?
                stack = round(length(cfg(:))/nnodes); % only used when not compiling 
%                 qsubcellfun(@MIBexp_preproc, cfg1, cfg2, cfg3, outputfile, 'memreq', memreq*1024^3, 'timreq', timreq*60, ...
%                     'stack', stack, 'StopOnError', true, 'backend', runcfg.preproc.parallel);

                qsubcellfun(@doICA, cfg(:), 'memreq', memreq, 'timreq', timreq*60, ...
                    'stack', stack, 'StopOnError', true, 'backend', parallel, 'options', '-l nodes=1:ppn=1');
                
            case 'yes'
                compiledfun = qsubcompile({@doICA, @runica}, 'toolbox', {'signal', 'stats'});
                qsubcellfun(compiledfun, cfg(:), 'memreq', memreq, 'timreq', timreq*60, ...
                    'stack', 1, 'StopOnError', false, 'backend', parallel, 'options', '-l nodes=1:ppn=1');
        end

        
    case 'parfor'
%         parfor ibatch = 1:length(cfg(:))
%             doICA(cfg1{ibatch}, cfg2{ibatch}, outputfile{ibatch})
%         end
    otherwise
        error('Unknown backend, aborting . . .\n')
end




%%%%%%%%%%%%%%%%%%%%%%%%%%% doICA saved in separate function

