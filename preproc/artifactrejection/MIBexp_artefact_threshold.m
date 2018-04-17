% MIBexp_artefact_threshold
fprintf('\n\nLooking for transients artifacts . . .\n')
cfg     = [];
% cfg.datafile = cfg1.datafile;
% cfg.headerfile = cfg1.datafile;
% cfg.headerformat = cfg1.headerformat;
cfg.trl = ft_findcfg(data.cfg, 'trl');
% cfg.continuous = 'yes';

cfg.artfctdef.threshold.channel = {'EEG'};

cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.range = range; % cf JJ PNAS

[cfg, artfctdef.threshold.artifact ] = ft_artifact_threshold(cfg, data);

fprintf('%d threshold artifacts found\n', length(artfctdef.threshold.artifact))
