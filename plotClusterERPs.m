col={[7 9],[11 13],[13 15]}
for p=1:8
    for eidx=[2 4 5]
        e=find([2 4 5]==eidx);
        ind=AllP{p}.Tall{eidx}.Data.findTrials('1','n','n');
        ind2=AllP{p}.Tall{eidx}.Data.findTrials('2','y','n');
        ind.cond2=setdiff(ind.cond2,ind2.cond2);
        for i=1:2
            rt=cell2mat(AllP{p}.Tall{eidx}.Data.segmentedEcog.event(ind.(['cond' int2str(i)]),col{e}(2)))-...
                cell2mat(AllP{p}.Tall{eidx}.Data.segmentedEcog.event(ind.(['cond' int2str(i)]),col{e}(1)));
            pr=prctile(rt,99);
            idx=find(rt<p & rt>-pr);
            responseTime{p,e,i}=mean(rt(idx));
            
            hist(rt(idx));
            %keyboard
        end
    end
end
clf
%%
%useQ=[2 3]
useQ=[2 3 4 5]
clf
eventSamps={210,300};
colorjet=jet;
opengl software
set(gcf,'Color','w')
load('C:\Users\Angela_2\Dropbox\AngelaUCSFFiles\AngelaSVN\Basic\redblackwhite')
%set(gcf,'renderer','painters')
set(gcf,'renderer','painters')

%samps={[1:400],[1601:2000]}
%samps2={[1:400],[801:1200]}

samps={[90:600],[2400:2900]}
colorcell={'darkred','darkblue'}
useeidx=[1 3]
for kIdx=1:size(kGroups,1)    
    k=kGroups(kIdx,1);
    pos=subplotMinGray(size(kGroups,1) ,2,kIdx,1);
    pos(1)=pos(1)+.25;
    pos(2)=pos(2)+.01
    pos(3)=pos(3)*.7
    ha=axes('Position',pos)
    a=imread(['E:\General Patient Info\EC16\brain5.jpg']);
    imshow(repmat(a,[1 1 3]))
    hold on
    for p=[1:8]
        idx=vertcat(qIdxHold{k,p,kGroups(kIdx,useQ)});
        if isempty(idx)
            continue
        end
        ch=allD.ch(idx);
        tmp=find(~ismember(ch,[6 9 70 73 134 137 192 198]));
        ch=ch(tmp);
        patientIdx=intersect(idx,find(allD.p==p));
        
        R=corr(allD.data(patientIdx,:)',nanmean(allD.data(idx(tmp),:),1)')
        %R=S(patientIdx);
        
        %color=colorjet(round(length(colorjet)*R),:);
        color=rgb('red')
        color=colorjet(findNearest(floor((p-1)*length(colorjet)/(8-1)),1:length(colorjet)),:);
        plotManyPolygons(BrainCoord(p).newXY(:,ch)',100,color,R.^2,12)
    end
    
    
    
    %clf
    idx=vertcat(qIdxHold{k,:,kGroups(kIdx,[2:end])});    
    for e=1:2
        pos=subplotMinGray(size(kGroups,1) ,4,kIdx,e-1);
        pos(4)=pos(4)*.2;
        pos(2)=pos(2)+.02;
        pos(1)=pos(1)*1.55+.03
        pos(3)=pos(3)*1.5;
        axes('Position',pos)
        idx=vertcat(qIdxHold{k,:,kGroups(kIdx,[4:5])});
        idx=sort(idx);
        tmp=find(~ismember(allD.ch(idx),[6 9 70 73 134 137 192 198]));        
        data1=allD.dataAll(idx,samps{e});
        
        
        idx=vertcat(qIdxHold{k,:,kGroups(kIdx,[2:3])});
        idx=sort(idx);
        data2=allD.dataAll(idx,samps{e}+600);
        color=colorjet(findNearest(floor((kIdx-1)*length(colorjet)/(size(kGroups,1)-1)),1:64),:);
        %pcolor(flipud(vertcat(zscore(allD.dataAll(idx,samps{e}),[],2), zscore(allD.dataAll(idx,samps{e}+399),[],2))));
        %shading interp
        %set(gca,'CLim',[0 3])
        %imagesc((vertcat((allD.dataAll(idx,samps{e}),[],2), zscore(allD.dataAll(idx,samps{e}+399),[],2))),[0 3])
        
            imagesc(vertcat(data1,data2),[0 2])
        
%         if abs(size(data1,1)-size(data2,1))<20
%             imagesc(vertcat(data1,data2),[0 2])
%         else
%             idx=vertcat(qIdxHold{k,:,kGroups(kIdx,[2:5])});
%             idx=sort(idx);
%             data2=allD.dataAll(idx,samps{e}+399);
%             imagesc(data2,[0 2])
%         end
        %set(gca,'HandleVisibility','off')
        
        %axes('Position',pos)
        %im=pcolor(zscore(allD.data(idx,samps{e}),[],2))
        hl=line([eventSamps{e} eventSamps{e}],[-1 600]);
        set(hl,'Color','k','LineStyle','--','LineWidth',2)%         set(hl,'Color','k')
%         axis([1 401 -1 3])
%         set(gca,'XTick',1:100:401,'XTickLabel',-2:2)
         set(gca,'YTick',[])
%         set(gca,'XGrid','on','YGrid','on')
        %xlabel('Time (s)')
        %ylabel(['Stacked Ave ERPs'])
                set(gca,'XLim',[0 500])
        if e==1
            set(gca,'XTick',11:100:600,'XTickLabel',-2:23)
        else
            set(gca,'XTick',1:100:600,'XTickLabel',-3:23)
        end
        set(gca,'YTick',20:20:300)
        set(gca,'YTickLabel',['20';repmat('  ',20,1)])
        %set(gca,'XTickLabel',[])
        %xlabel('Time (s)')
        colormap(redblackwhite)

        hold on
%         rt1=[];
%         for i=1:length(idx)            
%             rt1(i)=responseTime{allD.p(idx(i)),useeidx(e),1}
%         end
% 
%         
%         for i=1:length(idx)            
%             rt1(i+length(idx))=responseTime{allD.p(idx(i)),useeidx(e),2}
%         end
%         
%         %plot([nanmean(rt1(1:length(idx))*100+200) mean(rt1(1:length(idx))*100+200)],[-1 500],'-','LineWidth',2,'Color',rgb(colorcell{1}))
%         plot([nanmean(rt1(length(idx)+1:end)*100+200) mean(rt1(length(idx)+1:end)*100+200)]+1,[-1 500],'-','LineWidth',2,'Color',rgb(colorcell{2}))
%                
        if e==2
             rt2=[];
%             for i=1:length(idx)            
%                 rt2(i)=responseTime{allD.p(idx(i)),2,1}
%             end


%             for i=1:length(idx)            
%                 rt2(i+length(idx))=responseTime{allD.p(idx(i)),2,2}
%             end
            rt2=[responseTime{:,2,:}];
            plot([eventSamps{e}-median(rt2)*100 eventSamps{e}-median(rt2)*100],[1 500],'g','LineWidth',2','LineStyle','-')
        else
            plot([eventSamps{e}-200 eventSamps{e}-200],[-1 500],'g-','LineWidth',1)
        end
        
%          rt1=[];
%         for i=1:length(idx)            
%             rt1(i)=responseTime{allD.p(idx(i)),useeidx(e),1}
%         end
%         for i=1:length(idx)            
%             rt1(i+length(idx))=responseTime{allD.p(idx(i)),useeidx(e),2}
%         end
        
        %plot([mean(rt1(1:length(idx))*100+200) mean(rt1(1:length(idx))*100+200)],[-1 500],'-','LineWidth',2,'Color',rgb(colorcell{1}))
        
        rt1=[responseTime{:,useeidx(e),2}];
        plot([median(rt1*100+eventSamps{e}) median(rt1*100+eventSamps{e})],[-1 500],'-','LineWidth',2,'Color',rgb(colorcell{2}))
      
         
        
        plot([0 600],[length(idx)+1 length(idx)+1],'k','LineWidth',1)
        set(gca,'FontSize',10)
        set(gca,'XLim',[0 510])
    end    
    
    for eidx=[2 5]
        e=find([2 5]==eidx);
        idx=vertcat(qIdxHold{k,:,kGroups(kIdx,useQ)});
        idx=sort(idx); 
        data1=allD.dataAll(idx,samps{e});
        data2=allD.dataAll(idx,samps{e}+600);
        clear data
        baselinedata=vertcat(allD.dataAll(idx,samps{1}), allD.dataAll(idx,samps{1}+600));
        data{1}=data1;
        data{2}=data2;
        pos=subplotMinGray(size(kGroups,1) ,4,kIdx,e-1);
        pos(1)=pos(1)*1.55+.03
        pos(2)=pos(2)+pos(4)*.2+.02
        %pos(2)=pos(2)*1.3
        pos(3)=pos(3)*1.5;
        pos(4)=pos(4)*.7;
        axes('Position',pos)
        idx=vertcat(qIdxHold{k,:,kGroups(kIdx,useQ)});
        idx=sort(idx);
        idx=idx(find(~ismember(allD.ch(idx),[6 9 70 73 134 137 192 198])));        
        %idx=idx(find(~ismember(allD.p(idx),[7])));        
        %color=colorjet(findNearest(floor((kIdx-1)*length(colorjet)/(size(kGroups,1)-1)),1:64),:);
        color=rgb(colorcell{1});
        
         %PLOT SIG
        for sl=1:2
            if eidx==2
                shift=0;
            else
                shift=0;
            end
            
            try
                alpha=.01;
                %tmp=ps_raw{k,sl}(samps2{e});
                %psidx=find(~isnan(tmp));
                %[ps_FRD,sig]=MT_FDR_PRDS(tmp(psidx),alpha);
                %[~,pval]=ttest2(data{sl},repmat(reshape(data{1}(:,[70:150]+10)',[],1),1,size(samps{e},2)),alpha,'right');

                [~,pval]=ttest2(data{sl},repmat(reshape(baselinedata(:,[100:150]+10)',[],1),1,size(samps{e},2)),alpha,'right');
                [ps_FDR,sigHold{sl}]=MT_FDR_PRDS(pval,alpha);
                psidx=1:length(ps_FDR)
                plot(psidx(find(sigHold{sl}))+shift,sl*.4+2,'.','Color',rgb(colorcell{sl}))
                hold on
                sigHold2{kIdx,sl,e}=sigHold{sl};
            end
            try
                %tmp=ps_raw_2cond{k}(samps2{e});
                %psidx=find(~isnan(tmp));
                %[ps_FRD,sig]=MT_FDR_PRDS(tmp(psidx),alpha);
                %plot(psidx(find(sig))+shift,3,'r.')
                useidx=find(sigHold{1} | sigHold{2});
                 [~,pval]=ttest2(data1(:,useidx),data2(:,useidx),alpha);
                [ps_FDR,sig]=MT_FDR_PRDS(pval,alpha);
                psidx=1:length(ps_FDR);
                idx=vertcat(qIdxHold{k,:,kGroups(kIdx,[2:end])});
                hl=line(vertcat(useidx(find(sig))+shift, useidx(find(sig))+shift),...
                    vertcat(nanmean(data1(:,useidx(find(sig)))),nanmean(data2(:,useidx(find(sig))))));
                set(hl,'Color','m')
                hold on
            end
            
        end
        
        idx=vertcat(qIdxHold{k,:,kGroups(kIdx,[2:end])});
        idx=sort(idx);
        color=rgb(colorcell{1});
        %plot(nanmean(allD.dataAll(idx,samps{e})),'Color',color)
        hold on
        %[hl,hp]=errorarea(nanmean(allD.dataAll(idx,samps{e})),2*nansem(allD.dataAll(idx,samps{e})))
        %set(hl,'Color','magenta');
        %set(hp,'FaceColor',rgb(['light',colorcell{1}]),'FaceAlpha',1);
        plot(nanmean(data1),'Color',color,'LineWidth',3)

        color=rgb(colorcell{2});
        %plot(nanmean(allD.dataAll(idx,samps{e}+399)),'Color',color)
        %[hl,hp]=errorarea(nanmean(allD.dataAll(idx,samps{e}+399)),2*nansem(allD.dataAll(idx,samps{e})))
        %set(hl,'Color',rgb(['',colorcell{2}]));
        %set(hp,'FaceColor',rgb(['light',colorcell{2}]),'FaceAlpha',1);
        plot(nanmean(data2),'Color',color,'LineWidth',3)

        
        
        hold on
        hl=line([eventSamps{e} eventSamps{e}],[-1 100]);
        set(hl,'Color','k','LineStyle','--','LineWidth',2)
        axis([0 510 -.2 3])
        set(gca,'YTick',-0:3)
        if e==1
            set(gca,'XTick',11:100:500,'XTickLabel',[])
        else
            set(gca,'XTick',1:100:500,'XTickLabel',[])
        end

        %set(gca,'XTick',0:100:400,'XTickLabel',[])
        set(gca,'XGrid','on','YGrid','on')
        %ylabel('Zscore')
        set(gca,'FontSize',10)

        %set(gca,'YTickLabel',[],'XTickLabel',[])       
        
%          rt1=[];
%         for i=1:length(idx)            
%             rt1(i)=responseTime{allD.p(idx(i)),useeidx(e),1}
%         end
%         for i=1:length(idx)            
%             rt1(i+length(idx))=responseTime{allD.p(idx(i)),useeidx(e),2}
%         end
%         
%         plot([median(rt1(1:length(idx))*100+eventSamps) median(rt1(1:length(idx))*100+eventSamps)],[-1 500],'-','LineWidth',2,'Color',rgb(colorcell{1}))
%         plot([median(rt1(length(idx)+1:end)*100+eventSamps) median(rt1(length(idx)+1:end)*100+eventSamps)],[-1 500],'-','LineWidth',2,'Color',rgb(colorcell{2}))
%       
        rt1=[responseTime{:,useeidx(e),1}];
        plot([median(rt1*100+eventSamps{e}) median(rt1*100+eventSamps{e})],[-1 500],'-','LineWidth',2,'Color',rgb(colorcell{1}))
      
        
        rt1=[responseTime{:,useeidx(e),2}];
        plot([median(rt1*100+eventSamps{e}) median(rt1*100+eventSamps{e})],[-1 500],'-','LineWidth',2,'Color',rgb(colorcell{2}))
      
       
        
        if eidx==5
            plot([eventSamps{e}-median(rt2(1:end)*100) eventSamps{e}-median(rt2(1:end)*100)],[-1 5],'g-','LineWidth',1)
        else
            plot([eventSamps{e}-200 eventSamps{e}-200],[-1 5],'g-','LineWidth',1)
        end

    end        
    
    
    
  

%     idx=vertcat(qIdxHold{k,:,kGroups(kIdx,2:end)});    
% 
%     pos=subplotMinGray(1,2,1,1);
% %     pos(3)=pos(3)*1.2;
% %     pos(4)=pos(4)*1.2;
%     
%     ha=subplot('Position',pos);
%     try
%         imshow(repmat(a,[1 1 3]))
%     catch
%         imshow(a)
%     end
%     hold on
%     set(ha,'handlevisibility','off', ...
%         'visible','off')
%     ha=axes('position',get(ha,'Position'));
%     freezeColors
%     VqAve=mean(cat(3,Vq{kIdx,unique(allD.p(idx))}),3);
%     imshow(VqAve)
%     alpha(VqAve./max(max(VqAve)))
%     colormap(flipud(autumn))
    %title(num2str([k mean(R)]))
    set(gcf,'PaperPositionMode','auto')
    %print(gcf,'-depsc',['C:\Users\Angela_2\Dropbox\AngelaUCSFFiles\DelayWordBackup\clusters\painter_c',int2str(kIdx)])
    %%set(gcf,'WindowStyle','modal')
    %{
     exportfig(gcf,['C:\Users\Angela_2\Dropbox\AngelaUCSFFiles\DelayWordBackup\clusters2\painter_c',int2str(kIdx)],...
          'Color','rgb','FontMode','fixed','FontSize',12,'preview','tiff','SeparateText',0);
    %}
     %input('n')
end