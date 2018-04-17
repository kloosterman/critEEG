% MIBexp_artefact_muscle
fprintf('\n\nLooking for MUSCLE artifacts . . .\n')
cfg     = [];
% cfg.datafile = cfg1.datafile;
% cfg.headerfile = cfg1.datafile;
% cfg.headerformat = cfg1.headerformat;
cfg.datafile = ft_findcfg(data.cfg, 'datafile');
cfg.headerfile = cfg.datafile;
cfg.headerformat = ft_findcfg(data.cfg, 'headerformat');
cfg.trl = ft_findcfg(data.cfg, 'trl');
% cfg.trl = cfg.trl(1:30,:);
cfg.continuous = 'yes';

% cfg.artfctdef.zvalue.channel = {'EEG'};
% % excl_chans = [];
% % for i = 1:length(badchannels)
% %     excl_chans{i} = ['-' badchannels{i}];
% % end
% % cfg.artfctdef.zvalue.channel = [{'EEG'} excl_chans];
%
% % cfg.artfctdef.zvalue.baselinewindow = [-0.45 0]
%
% %            cfg.artfctdef.zvalue.detrend = 'yes'
% %            cfg.artfctdef.zvalue.demean = 'yes'
% % % % %            cfg.artfctdef.zvalue.hpfilter      = 'yes'
% % % % %            cfg.artfctdef.zvalue.hpfreq = 5
%
% % cutoff and padding
% cfg.artfctdef.zvalue.cutoff      = cfg1.musclethr;
% cfg.artfctdef.zvalue.trlpadding  = 0;
% cfg.artfctdef.zvalue.fltpadding  = 0.2;
% cfg.artfctdef.zvalue.artpadding  = 0.1;
% % algorithmic parameters
% cfg.artfctdef.zvalue.bpfilter    = 'yes';
% % cfg.artfctdef.zvalue.bpfreq      = [110 140];
% cfg.artfctdef.zvalue.bpfreq      = [110 127]; % Nyquist: 128 Hz
% cfg.artfctdef.zvalue.bpfiltord   = 3; %9
% cfg.artfctdef.zvalue.bpfilttype  = 'but';
% cfg.artfctdef.zvalue.hilbert     = 'yes';
% cfg.artfctdef.zvalue.boxcar      = 0.2;
% % feedback
% cfg.artfctdef.zvalue.interactive = cfg1.artf_feedback;
% % cfg.artfctdef.zvalue.interactive = 'yes';
%
%%%%
% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel = 'EEG';
%   cfg.artfctdef.zvalue.channel = 'I*';
cfg.artfctdef.zvalue.cutoff      = -1; %% JJ: threshold automatically
%   cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.fltpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter    = 'yes';
cfg.artfctdef.zvalue.bpfreq      = [110 127]; % 110 - 140
cfg.artfctdef.zvalue.bpfiltord   = 6;
cfg.artfctdef.zvalue.bpfilttype  = 'but';
cfg.artfctdef.zvalue.hilbert     = 'yes';
cfg.artfctdef.zvalue.boxcar      = 0.2;

% make the process interactive
% cfg.artfctdef.zvalue.interactive = cfg1.artf_feedback;
cfg.artfctdef.zvalue.interactive = 'no';


[cfg, artfctdef.muscle.artifact ] = ft_artifact_zvalue(cfg);

fprintf('%d muscle artifacts found\n', length(artfctdef.muscle.artifact))

