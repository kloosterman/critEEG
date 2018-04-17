function doICA(cfg)
% concatenate all runs per session
% remove EOG artifacts with ICA
% remove muscle artifacts per run

PREIN = cfg.PREIN;
PREOUT = cfg.PREOUT;
mkdir(PREOUT)

cd(PREIN)
runlist = dir('*data.mat');
data = cell(size(runlist,1),1);
for irun = 1:length(runlist)
    temp = load(runlist(irun).name);
    data{irun} = temp.data;   
    clear temp
end

if length(data) > 1
    data = ft_appenddata([], data{:});
else
    data = data{:};
end
data.fsample = 256;
data = rmfield(data, 'sampleinfo');

%% Detrending! Done here to avoid "same sample different value" ft issue
cfgtmp = [];
cfgtmp.detrend = 'yes';
data = ft_preprocessing(cfgtmp, data);

cfg.rEOG.channel    = {'EOGH','Fp1','Pz'};
cfg.rEOG.refchannel = {'Pz'};
cfg.rEOG.reref      = 'yes';
tmp.EOG = ft_preprocessing(cfg.rEOG, data);

cfg.tmp.channel    = { 'EOGH','Fp1' }; % get rid of Pz
tmp.EOG = ft_preprocessing(cfg.tmp,tmp.EOG);
cfg = rmfield(cfg,'tmp');

% compute radial EOG and add it to data
for iEOG = 1:length(tmp.EOG.trial)
    data.trial{iEOG}(end+1,:) = mean(tmp.EOG.trial{iEOG},1);
end
data.label(end+1) = {'rEOG'};
clear tmp

% perform COSTRAP-TSP-detection
data = nak_mwb_costrap_detection2(data,{'rEOG'},{'Fp1', 'Fp2'}, 2);

%% ICA
cfg = [];
cfg.channel = {'EEG' 'TSP*' };
cfg.method = 'runica';

%  cfg.trials = 1; % testing

comp = ft_componentanalysis(cfg, data);

[~, compfilename]= fileparts(runlist(1).name);
% outputfile = sprintf('%s_comp', compfilename);

outputfile = sprintf('%s_costrap_comp', compfilename); % COSTRAPPED!

fprintf('Saving %s to...\n %s\n', outputfile, PREOUT)
save(fullfile(PREOUT, outputfile), 'comp', 'data');


