for kIdx=1:length(kGroups)
    k=kGroups(kIdx,1);
    for p=1:8
        idx= vertcat(qIdxHold{k,p,:});
        try
            R=corr(allD.data(idx,:)',(mean(allD.data(vertcat(qIdxHold{k,:,:}),:))'))
            [tmp1,tmp]=sort(R,'descend');
            idx=idx(tmp);
        end
        display(p).ch{kIdx}=allD.ch(idx);
    end
end
 display(8).ch{7}=[140 124 123 62 109]
%%
% PLOT ALIGNED BY PERC AND PROD
load([mainPath filesep 'redblackwhite'])
load([mainPath filesep 'blueblackwhite'])
set(gcf,'renderer','painters')
clf
set(gcf,'Color','w')
samps=1:900;
eventSamp={210,300};
colorcell={'blue','red','green'};
set(gcf,'Color','w')
col{2}=[7 9]
col{4}=[11 13]
col{5}=[13 15]
sigPos=[-.6 -.9]
fontsize=8;
rows=length(kGroups);
for p=1:8
    %%
    clf
    for i=1:rows
        for eidx=2
            e=find([2 5]==eidx);
            AllP{p}.Tall{eidx}.Data.segmentedEcog.event(find(cellfun(@isempty,AllP{p}.Tall{eidx}.Data.segmentedEcog.event)))={NaN}
            %%
            for j=length(display(p).ch{i}):-1:1
                if j>length(display(p).ch{i})
                    continue
                end
                ch=display(p).ch{i}(j);
                
                %SHOW BRAIN
                pos=subplotMinGray(rows,9,i,0);
                pos(1)=pos(1)-.02
                pos(3)=pos(3)*1;
                pos(4)=pos(4)*1;
                subplot('Position',pos);
                if j==1%length(display(p).ch{i})
                    flag=1;
                else
                    flag=0;
                end
                try
                    if j==1
                        a=imread(['E:\DelayWord\allBrainPics\' 'EC16' 'scaled.jpg']);
                        imshow(a)
                        hold on
                    end
                    scatter(BrainCoord(p).newXY(1,ch),BrainCoord(p).newXY(2,ch),30,rgb(colorcell{j}),'fill')
                    
                catch
                    colorcell{j}='black';
                end
                axis([0 1000 10 620])
                set(gca,'Visible','off','XTickLabel',[],'YTickLabel',[])
                
                
                %GET DATA
                badtr=unique([find(cellfun(@isnan,AllP{p}.Tall{2}.Data.segmentedEcog(1).event(:,15)))']);
                usetrials1=setdiff(1:size(AllP{p}.Tall{2}.Data.segmentedEcog.event,1),[AllP{p}.Tall{2}.Data.Artifacts.badTrials badtr]);%AllP{p}.Tall{eidx}.Data.Params.usetr;
                [s,sidx]=sort(cell2mat(AllP{p}.Tall{2}.Data.segmentedEcog(1).event(usetrials1,9))- cell2mat(AllP{p}.Tall{eidx}.Data.segmentedEcog(1).event(usetrials1,7)));
                usetrials1=usetrials1(sidx);
                data=squeeze(AllP{p}.Tall{eidx}.Data.segmentedEcog.smoothed100(ch,:,:,usetrials1));
                data1=data(90:590,:);
                
                badtr=unique([find(cellfun(@isnan,AllP{p}.Tall{5}.Data.segmentedEcog(1).event(:,15)))'])
                usetrials2=setdiff(1:size(AllP{p}.Tall{5}.Data.segmentedEcog.event,1),[AllP{p}.Tall{5}.Data.Artifacts.badTrials badtr]);%AllP{p}.Tall{eidx}.Data.Params.usetr;                [s,sidx]=sort(cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog(1).event(usetrials,9))- cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog(1).event(usetrials,7)));
                [s,sidx]=sort(cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog(1).event(usetrials2,9))...
                    - cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog(1).event(usetrials2,7)));
                usetrials2=usetrials2(sidx);
                data2=squeeze(AllP{p}.Tall{5}.Data.segmentedEcog(1).smoothed100(ch,1:500,:,usetrials2));
                
                
                %PLOT AVE PRES
                pos=subplotMinGray(rows,9,i,e);
                pos(3)=pos(3)*2;
                pos(4)=pos(4)*.9;
                pos(1)=pos(1)-.025;
                subplot('position',pos);
                try
                    plot(mean(data1,2),'Color',rgb(colorcell{j}),'LineWidth',2)
                catch
                    plot(mean(data1,2),'Color','k','LineWidth',2)
                end

                hold on
                [~,pval]=ttest2(data1',repmat(reshape(data1([70:150],:),[],1),[1,size(data1,1)]'),0.01,'Right');
                [pval,sig]=MT_FDR_PRDS(pval,0.01);
                try
                    plot(find(sig),6+j*.5,'.','Color',rgb(colorcell{j}),'LineWidth',3,'MarkerSize',3)
                end
                
                res=cell2mat(AllP{p}.Tall{eidx}.Data.segmentedEcog.event(usetrials1,col{2}(2)))-cell2mat(AllP{p}.Tall{eidx}.Data.segmentedEcog.event(usetrials1,col{2}(1)));
                res(find(res<=0))=NaN;
                hold on
                plot(repmat(eventSamp{1}+nanmean(res)*100,[1 10]) ,[-1:8],'Color',rgb('darkblue'),'LineWidth',2)
                plot(repmat(eventSamp{1}-200,[1 10]) ,[-1:8],'g','LineWidth',2)


                axis([0 500 -1 8])
                hl=line([eventSamp{1} eventSamp{1}],[-1 8]);
                set(hl,'Color','k','LineWidth',2,'LineStyle','--')
                set(gca,'XTick',10:100:500,'XTickLabel',-2:1:3,'FontSize',fontsize)
                ylabel('Zscore')
                if i==9
                    xlabel('Time (s) around Presentation Onset')
                end
                set(gca,'XGrid','on','YGrid','on')
                
                if j==1
                    %PLOT STACK PRES
                    pos=subplotMinGray(rows,18,i,1*2+6)
                    pos(4)=pos(4)*.9;
                    pos(3)=pos(3)*1.3;
                    pos(1)=pos(1)-.04;
                    subplot('position',pos);

                    imagesc(data1')
                    set(gca,'clim',[0 3])
                    hl=line([eventSamp{1} eventSamp{1}],[1 length(usetrials1)]);
                    set(hl,'Color','k','LineWidth',2,'LineStyle','--')
                    res=cell2mat(AllP{p}.Tall{eidx}.Data.segmentedEcog.event(usetrials1,col{2}(2)))...
                        -cell2mat(AllP{p}.Tall{eidx}.Data.segmentedEcog.event(usetrials1,col{2}(1)));
                    res(find(res<=0))=NaN;
                    hold on
                    plot(eventSamp{1}+(res)*100,1:length(usetrials1),'k','LineWidth',2)
                    plot(eventSamp{1}-200,1:length(usetrials1),'g','LineWidth',2)

                    set(gca,'XTick',10:100:500,'XTickLabel',-2:1:3,'FontSize',fontsize)
                    if i~=rows
                        set(gca,'XTickLabel',[]);
                    end
                    set(gca,'YTickLabel',[])
                    if j==1
                        ylabel('Sorted Trials')
                    end
                end
                
                %PLOT AVE PROD
                pos=subplotMinGray(rows,9,i,3);
                pos(1)=pos(1)-.02;
                pos(3)=pos(3)*2;
                pos(4)=pos(4)*.9;
                pos(1)=pos(1)-.040;

                subplot('position',pos);
                try
                    plot(mean(data2,2),'Color',rgb(colorcell{j}),'LineWidth',2,'LineWidth',2)
                catch
                    plot(mean(data2,2),'Color','k','LineWidth',2)
                end

                hold on
                
                [~,pval]=ttest2(data2',repmat(reshape(data1(70:150,:),[],1),[1,size(data2,1)]'),0.01,'Right');
                [pval,sig]=MT_FDR_PRDS(pval,0.01);
                try
                    plot(find(sig),6+j*.5,'.','Color',rgb(colorcell{j}),'LineWidth',3,'MarkerSize',3)
                    hold on
                    hl=line([eventSamp{2} eventSamp{2}],[-1 8]);
                end
                res=cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog.event(usetrials2,col{5}(2)-2))...
                    -cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog.event(usetrials2,col{5}(1)-2));
                res(find(res<=0))=NaN;
                hold on
                plot(repmat(eventSamp{2}-nanmean(res)*100,[1 10]) ,[-1:8],'g','LineWidth',2)
                
                
                res=cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog.event(usetrials2,col{5}(2)))...
                    -cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog.event(usetrials2,col{5}(1)));
                res(find(res<=0))=NaN;
                hold on
                plot(repmat(eventSamp{2}+nanmean(res)*100,[1 10]) ,[-1:8],'Color',rgb('darkblue'),'LineWidth',2)
                
                set(hl,'Color','k','LineWidth',2,'LineWidth',2,'LineStyle','--')
                axis([0 500 -1 8])
                set(gca,'XTick',0:100:500,'XTickLabel',-3:1:2,'FontSize',fontsize)
                set(gca,'YTickLabel',[])
                if i==9
                    xlabel('Time (s) around Production Onset')
                end
                set(gca,'XGrid','on','YGrid','on')
                
                if j==1
                    %PLOT STACK PROD
                    pos=subplotMinGray(rows,18,i,1*2+7)
                    pos(1)=pos(1)-.017;
                    pos(4)=pos(4)*.9;
                    pos(1)=pos(1)-.025;
                    pos(3)=pos(3)*1.25;

                    subplot('position',pos);
                    %pcolor((data2'))
                    %shading interp
                    imagesc(data2')

                    set(gca,'clim',[0 3])
                    hl=line([eventSamp{2} eventSamp{2}],[1 length(usetrials2)]);
                    set(hl,'Color','k','LineWidth',2,'LineStyle','--')
                    res=cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog.event(usetrials2,col{5}(2)))...
                        -cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog.event(usetrials2,col{5}(1)));
                    res(find(res<=0))=NaN;
                    hold on
                    plot(eventSamp{2}+(res)*100,1:length(usetrials2),'k','LineWidth',2)

                    res=cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog.event(usetrials2,col{5}(2)-2))...
                        -cell2mat(AllP{p}.Tall{5}.Data.segmentedEcog.event(usetrials2,col{5}(1)-2));
                    res(find(res<=0))=NaN;
                    hold on
                    plot(eventSamp{2}-(res)*100,1:length(usetrials2),'LineWidth',2,'Color','g','LineStyle','-')

                    set(gca,'XTick',0:100:500,'XTickLabel',-3:1:2,'FontSize',fontsize)
                    set(gca,'YTickLabel',[])
                    if i~=rows
                        set(gca,'XTickLabel',[]);
                    end

                    if eidx==5
                        colormap(redblackwhite)
                    else
                        colormap(blueblackwhite)
                    end
                    set(gca,'XLim',[0 500])
                    freezeColors
                end
            end            
        end
    end
    exportfig(gcf,['C:\Users\Angela_2\Dropbox\AngelaUCSFFiles\DelayWordBackup\ERPs\p' int2str(p)],         'Color','rgb','FontMode','fixed','FontSize',12,'preview','tiff','SeparateText',0);
    %input('b')
end