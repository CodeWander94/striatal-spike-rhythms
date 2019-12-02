%%Plot the prediction improvement of models across days%%
cd('/Users/manishm/Work/vanDerMeerLab/GLM/temp/PreTrial')
rats = {'R119'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat('*',curRat,'*sd.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        figure;
        graphTitle = strcat(curRat,'-Day-',num2str(jdx));
        load(ofiles(jdx).name);
        mn = fieldnames(sd.m);
        tp = 0;
        for iM = 1:length(mn)
            if (tp == 0)
                g = subplot(length(mn),1,1);
                t = title(graphTitle,'FontSize',40);
                g.Visible = 'off';
                t.Visible = 'on';
                tp = iM+1;
            else
               tp = tp+1;
            end
            if strcmp(mn{iM},'baseline')
                tp = tp-1;
                continue;
            end
            this_err = sd.m.(mn{iM}).err;
   
            % mean error over cells
            celldiffmean = nanmean(this_err,2);

            subplot(length(mn),1,tp);
            plot(celldiffmean,'LineWidth',1); hold on;
            keep = find(celldiffmean > 0); plot(keep,celldiffmean(keep),'.g','MarkerSize',20);
            keep = find(celldiffmean < 0); plot(keep,celldiffmean(keep),'.r','MarkerSize',20);
            set(gca,'LineWidth',1,'FontSize',18); box off;
%             ylabel('Prediction improvement'); xlabel('cell #');
            title(mn{iM});
        end
    end
end