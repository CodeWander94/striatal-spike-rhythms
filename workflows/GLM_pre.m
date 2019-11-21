cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/R117/R117-2007-06-02');
cfg.output_dir = '/Users/manishm/Work/vanDerMeerLab/GLM/temp';
please.rats = {'R117'};
[cfg.fd,cfg.fd_extra] = getDataPath(please);
cfg.iS = 1;
cfg.nMinSpikes = 10;
% cfg.writeFullError = 1;
% cfg.writeOutput = 1;
% cfg.plotOutput = 0;
cfg.output_dir = '/Users/manishm/Work/vanDerMeerLab/GLM/temp'; % store files here
cfg.output_prefix = 'S2_'; % prefix filenames with this (identify runs)
cfg.nShuf = 10;
sd = rhythmGLMfit_shufPRE(cfg);