% cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/R119/R119-2007-07-14');
cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/R117/R117-2007-06-08');
cfg.output_dir = '/Users/manishm/Work/vanDerMeerLab/GLM/temp';
please.rats = {'R117'};
[cfg.fd,cfg.fd_extra] = getDataPath(please);
cfg.iS = 1;
sd = rhythmGLMfit_v1(cfg);