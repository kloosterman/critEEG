% MIBexp_artefact_jump
fprintf('\n\nLooking for JUMP artifacts . . .\n')
cfg     = [];
cfg.datafile = cfg1.datafile;
cfg.headerfile = cfg1.datafile;
cfg.headerformat = cfg1.headerformat;
cfg.trl = data.cfg.trl;

cfg.padding = 0;
cfg.continuous = 'yes';
% cutoff and padding
% select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
% cfg.artfctdef.zvalue.channel= [cfg1.channel {'-EOG'}];
cfg.artfctdef.zvalue.channel = {'MEG'}; %, '-M*1'};  %grads only
cfg.artfctdef.zvalue.cutoff  = cfg1.jumpthr;
cfg.artfctdef.zvalue.trlpadding= 0.5*cfg.padding;
cfg.artfctdef.zvalue.artpadding= 0.5*cfg.padding;
cfg.artfctdef.zvalue.fltpadding= 0.1;
% algorithmic parameters
cfg.artfctdef.zvalue.cumulative= 'yes';
cfg.artfctdef.zvalue.medianfilter= 'yes';
cfg.artfctdef.zvalue.medianfiltord= 9;
cfg.artfctdef.zvalue.absdiff= 'yes';
% feedback (artifact viewer)
cfg.artfctdef.zvalue.interactive= cfg1.artf_feedback;
[cfg, artfctdef.jump.artifact] = ft_artifact_zvalue(cfg);
% filesave = sprintf('%s_artf_jump', outputfile);
% fprintf('Saving %s . . .\n',filesave);
% save(filesave,'artf_jump');

fprintf('%d jump artifacts found\n', length(artfctdef.jump.artifact))