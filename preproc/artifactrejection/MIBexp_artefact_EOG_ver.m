% MIBexp_artefact_EOG
fprintf('\n\nLooking for EOG artifacts . . .\n')
cfg     = [];
cfg.datafile = cfg1.datafile;
cfg.headerfile = cfg1.datafile;
cfg.headerformat = cfg1.headerformat;
cfg.trl = ft_findcfg(data.cfg, 'trl');
% cfg.padding = 1;
cfg.padding = 0;
cfg.continuous = 'yes';
% cutoff and padding
% select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
cfg.artfctdef.zvalue.channel = {'EXG3', 'EXG4'};
cfg.artfctdef.zvalue.cutoff      = cfg1.eogverthr;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.fltpadding  = 0.2;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';
% feedback
cfg.artfctdef.zvalue.interactive = cfg1.artf_feedback;
[cfg, artfctdef.eog_ver.artifact] = ft_artifact_zvalue(cfg);

fprintf('%d EOG artifacts found\n', length(artfctdef.eog_ver.artifact))
