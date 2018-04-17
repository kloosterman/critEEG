function rejectEOGandCSD(cfg)
% Do EOG rejection around stim onset
% do threshold artifact rej
% compute csd and save

PREIN = cfg.PREIN;

cd(PREIN)
% inputfile = dir('*data_clean.mat'); % only EOG ic removed
inputfile = dir('*data_costrap_clean.mat'); % only EOG ic removed

fprintf('Loading %s from...\n %s\n', inputfile.name, PREIN)
load(inputfile.name);

%EOG artifact rejection
fprintf('\n\nLooking for EOG artifacts . . .\n')
cfg     = [];
% cfg.padding = 0;
% cfg.continuous = 'no';
% cutoff and padding
% select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
cfg.artfctdef.zvalue.channel = {'EOGV'};
cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.trlpadding  = 0;
%     cfg.artfctdef.zvalue.artpadding  = 0.1;
%     cfg.artfctdef.zvalue.fltpadding  = 0.2;
cfg.artfctdef.zvalue.artpadding  = 0;
cfg.artfctdef.zvalue.fltpadding  = 0;
% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';
% feedback
cfg.artfctdef.zvalue.interactive = 'no';
if isfield(data, 'sampleinfo')
    data = rmfield(data, 'sampleinfo'); % sampleinfo is bogus if data was appended
end
[cfgout, artfctdef.eog_ver.artifact] = ft_artifact_zvalue(cfg, data);
fprintf('%d vertical EOG artifacts found\n', length(artfctdef.eog_ver.artifact))

% horizontal movements (saccades)
cfg.artfctdef.zvalue.channel = {'EOGH'};
cfg.artfctdef.zvalue.cutoff      = 6;

cfg.artfctdef.zvalue.interactive = 'no';

[cfg, artfctdef.eog_hor.artifact] = ft_artifact_zvalue(cfg, data);
fprintf('%d horizontal EOG artifacts found\n', length(artfctdef.eog_hor.artifact))

cfg_rej.artfctdef.eog_ver.artifact = artfctdef.eog_ver.artifact;
cfg_rej.artfctdef.eog_hor.artifact = artfctdef.eog_hor.artifact;

% only reject if around (0.1 s) target onset,
cfg_rej.artfctdef.crittoilim = [zeros(size( data.trialinfo,1),1)-0.1  zeros(size( data.trialinfo,1),1)+0.1];

data  = ft_rejectartifact(cfg_rej,data);

%         cfg=[];
%         cfg.layout = 'elec1010.lay';
%         cfg.channel = 1:48;
%         ft_rejectvisual(cfg,data)
%         data

% do threshold art rejection
range=200; 
% MIBexp_artefact_threshold
fprintf('\n\nLooking for transients artifacts . . .\n')
cfg     = [];
cfg.trl = ft_findcfg(data.cfg, 'trl');
cfg.artfctdef.threshold.channel = {'EEG'};
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.range = range; % cf JJ PNAS

% cfg.showcallinfo = 'no';

artfctdef = [];
[~, artfctdef.threshold.artifact ] = ft_artifact_threshold(cfg, data);

fprintf('%d threshold artifacts found\n', size(artfctdef.threshold.artifact,1))

cfg_rej = [];
cfg_rej.artfctdef.threshold.artifact = artfctdef.threshold.artifact;
data  = ft_rejectartifact(cfg_rej,data);

% compute scd maps AS LAST STEP
cfg=[];
if ismac
    cfg.elecfile = 'standard_1005.elc'; % for topo take elec1010.lay
else
    cfg.elecfile = '/path/fieldtrip-20161220/template/electrode/standard_1005.elc'; % for topo take elec1010.lay
end
cfg.method = 'spline';
data = ft_scalpcurrentdensity(cfg, data); %DEBUGGED: return data instead of scd

[~, datafilename]= fileparts(inputfile.name);
outputfile = sprintf('%s_costrap_CSD', datafilename);

fprintf('Saving %s to...\n %s\n', outputfile, PREIN)
save(outputfile, 'data');


