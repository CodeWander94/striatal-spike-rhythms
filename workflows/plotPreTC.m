%%Plot the prediction improvement of models across days%%
cd('/Users/manishm/Work/vanDerMeerLab/GLM/temp/')
rats = {'R117'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat('S2_',curRat,'*sd.mat');
    ofiles = dir(searchString);
    paxis = linspace(-3.14,3.14,36);
    for jdx = 1:length(ofiles)
        figure;
        graphTitle = strcat(curRat,'-Day-',num2str(jdx),'(Before Task)');
        load(ofiles(jdx).name);
        iN = find(sd.cellType==2);
        tp = 0;
        for iM = 1:length(iN)
            if (tp == 0)
                g = subplot(length(iN)+1,5,1:5);
                t = title(graphTitle,'FontSize',40);
                g.Visible = 'off';
                t.Visible = 'on';
                tp = 6;
            else
                tp = tp+5;
            end
            subplot(length(iN)+1,5,tp);
            plot(paxis,sd.beta_phase{iM});
            set(gca,'LineWidth',1,'FontSize',18); box off;
            subplot(length(iN)+1,5,tp+1);
            plot(paxis,sd.theta_phase{iM});
            set(gca,'LineWidth',1,'FontSize',18); box off;
            subplot(length(iN)+1,5,tp+2);
            plot(paxis,sd.delta_phase{iM});
            set(gca,'LineWidth',1,'FontSize',18); box off;
            subplot(length(iN)+1,5,tp+3);
            plot(paxis,sd.lowGamma_phase{iM});
            set(gca,'LineWidth',1,'FontSize',18); box off;
            subplot(length(iN)+1,5,tp+4);
            plot(paxis,sd.highGamma_phase{iM});
            set(gca,'LineWidth',1,'FontSize',18); box off;
        end
    end
end