cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/R119-2007-07-14');
cfg.output_dir = '/Users/manishm/Work/vanDerMeerLab/GLM/temp';
please.rats = {'R119'};
[cfg.fd,cfg.fd_extra] = getDataPath(please);
cfg.iS = 1;
sd = rhythmGLMfit_shuf(cfg);