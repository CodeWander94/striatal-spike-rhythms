cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/R117/R117-2007-06-06');
cfg.output_dir = '/Users/manishm/Work/vanDerMeerLab/GLM/temp';
please.rats = {'R117'};
[cfg.fd,cfg.fd_extra] = getDataPath(please);
cfg.iS = 1;
cfg.nMinSpikes = 100;
%cfg.nShuf = 10;
sd = rhythmGLMfit_shufTask(cfg); 