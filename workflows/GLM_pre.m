cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/R117/R117-2007-06-01');
cfg.output_dir = '/Users/manishm/Work/vanDerMeerLab/GLM/temp';
please.rats = {'R117'};
[cfg.fd,cfg.fd_extra] = getDataPath(please);
cfg.iS = 1;
cfg.nMinSpikes = 100;
% cfg.writeFullError = 1;
cfg.writeOutput = 1;
cfg.plotOutput = 1;
cfg.output_dir = '/Users/manishm/Work/vanDerMeerLab/GLM/temp'; % store files here
cfg.output_prefix = 'S2Pre_'; % prefix filenames with this (identify runs)
cfg.nShuf = 10;
sd = rhythmGLMfit_shufPRE(cfg);