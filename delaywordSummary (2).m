%{ Documentation of Delay Word Repetition Project }
%% PREPROCESSED FILES
mainPath='E:\DelayWord\Summary\Files'
%load([mainPath filesep 'AllP']);%Segmented data
load([mainPath filesep 'BrainCoord2']);%Sulci, grid, normalized grid coordinates
load([mainPath filesep 'allConditions']);%Sulci, grid, normalized grid coordinates




%%
% STEP 1: GET EVENT TIMES
% ALREADY CREATED EVENT FILES AVAILABLE AT:
% TO CREATE EVENT FILES FROM SCRATCH:

%Paths to blocks of raw data
mainDataPath=['E:\DelayWord'];
blocks = {
    'EC18\EC18_B1'
    'EC16\EC16_B1'
    'EC22\EC22_B1'
    'EC23\EC23_B1'
    'EC24\EC24_B1'
    'EC29\EC29_B2'
    'EC30\EC30_B1'
    'EC31\EC31_B1'
    'EC18\EC18_B2'
    'EC23\EC23_B3'
    'EC24\EC24_B3'}

filepath=cellfun(@(x) [mainDataPath filesep x],blocks,'UniformOutput',0)
%%
wpath=[mainPath filesep 'brocaSoundFiles\'];
load([wpath filesep 'BrocaWords.mat'])

for p=1:length(filepath)
    %%
    [dtpath,expt]=fileparts(filepath{p});
    getAnalogEvents(dtpath, wpath, expt, names) 
    save([dtpath filesep expt filesep 'Analog/allConditions'],'allConditions')
end

%%
% STEP 2: LOAD AND SEGMENT DATA
% ALREADY SEGMENTED DATA FILES AVAILABLE IN OBJECTS AT:
% TO SEGMENT AND LOAD DATA FROM SCRATCH:



%Segment around events 1-5 for patients p1-8
for p=1:8
    %%
    for event=[2 4 ]
        
        cd([filepath{p} filesep 'RawHTK'])
        numCh=length(dir)-2;
        
        %Initialize object
        data=SegmentedData([filepath{p} filesep 'HilbAA_70to150_8band'],[],1:numCh);

        %conditions types (all words)
        lh=1:40;
        %conditions to segment by
        seg={[repmat([41],[1 length(lh)]);lh],[lh;repmat([42],[1 length(lh)])],...
            [repmat([42],[1 length(lh)]);lh],[repmat([43],[1 length(lh)]);lh],...
            [repmat([44],[1 length(lh)]);lh],[45;41],[1:50;1:50]};

        %load data segments (HG)
        data.segmentedDataEvents40band(seg(:,event),{[3000 3000]},[],'aa',31:38,0,'DelayRep')
        
        %store in structure
        Tall{event}.Data=data;      
    end
    
    %Calculate baseline
    Tall{2}.Data.ecogBaseline.data=reshape(Tall{2}.Data.segmentedEcog.data(:,[500:800]+400,:,:),Tall{2}.Data.channelsTot,[]);
    Tall{2}.Data.ecogBaseline.std(:,1,1)=std(Tall{2}.Data.ecogBaseline.data,[],2);
    Tall{2}.Data.ecogBaseline.mean(:,1,1)=mean(Tall{2}.Data.ecogBaseline.data,2);
    
    Tall{2}.Data.ecogBaseline.zscore_separate=zscore(Tall{2}.Data.ecogBaseline.data,[],2)
    [a,b]=find(Tall{2}.Data.ecogBaseline.zscore_separate(setdiff(1:Tall{2}.Data.channelsTot,Tall{2}.Data.Artifacts.badChannels),:)>6);
    Tall{2}.Data.ecogBaseline.zscore_separate(:,unique(b))=[];
    try
        Tall{5}.Data.BaselineChoice='rest';
        Tall{5}.Data.ecogBaseline=Tall{2}.Data.ecogBaseline;
        Tall{5}.Data.calcZscore(1,1);
    catch
    end
    Tall{4}.Data.BaselineChoice='rest';
    Tall{2}.Data.BaselineChoice='rest';
    
    Tall{4}.Data.ecogBaseline=Tall{2}.Data.ecogBaseline;
    Tall{4}.Data.calcZscore(1,1);
    Tall{2}.Data.calcZscore(1,1);
    
    for event=[2 4 ]
         %smooth data
        Tall{event}.Data.ecogBaseline.smoothed100=[];
        for c=1:Tall{event}.Data.channelsTot
            for tr=1:size(Tall{event}.Data.segmentedEcog.zscore_separate,4)
                Tall{event}.Data.segmentedEcog.smoothed100(c,:,:,tr)=(resample(squeeze(Tall{event}.Data.segmentedEcog.zscore_separate(c,:,:,tr)),1,4));
            end
            Tall{event}.Data.ecogBaseline.smoothed100(c,:)=(resample(Tall{event}.Data.ecogBaseline.zscore_separate(c,:),1,4));
        end
        Tall{event}.Data.segmentedEcog.data=[];
        Tall{event}.Data.segmentedEcog.zscore_separate=[];
        
    end
    AllP{p}.Tall=Tall;
end
save('AllP_nosmooth','AllP','-v7.3')
%% FIND BAD TRIALS
clear d
for p=1:8
    for eidx=[2 4 5]
        AllP{p}.Tall{eidx}.Data.Artifacts.badTrials=[];
        ch=setdiff(1:AllP{p}.Tall{eidx}.Data.channelsTot,AllP{p}.Tall{eidx}.Data.Artifacts.badChannels);
        d=zscore(squeeze(mean(squeeze(AllP{p}.Tall{eidx}.Data.segmentedEcog.smoothed100(ch,:,:,:)),1)),[],1);
        [a,b]=find(zscore(mean(d,1))>5);
        AllP{p}.Tall{eidx}.Data.Artifacts.badTrials=unique([AllP{p}.Tall{eidx}.Data.Artifacts.badTrials b']);
        imagesc(d')
    end
    %input('b')
end
%% SET EMPTY EVENTS TO NAN
for p=1:8
    for eidx=[2 4 5]
        AllP{p}.Tall{eidx}.Data.segmentedEcog.event(find(cellfun(@isempty,AllP{p}.Tall{eidx}.Data.segmentedEcog.event)))={NaN}
    end
end
%% FIGURE 1: CLUSTER ERPS

%% MAKE DATA STRUCTURE OF AVE WAVEFORM FOR SHORT AND LONG WORDS: ONLY SAMPLES USED FOR CLUSTERING
allD.data=[];
allD.p=[];
allD.ch=[];
samps={[1:400]+100, [1:250]+100, [101:400]+100}

i=1;
for p=1:8
    ch=1:AllP{p}.Tall{2}.Data.channelsTot;
    ch=setdiff(ch,AllP{p}.Tall{2}.Data.Artifacts.badChannels);
    d=[];
    ehold=[];
    for eidx=[2 4 5]
        d{eidx}=[];
        e=find([2 4 5]==eidx);
        ind=AllP{p}.Tall{eidx}.Data.findTrials('1','n','n');
        for condIdx=1:2
            useind=setdiff(ind.(['cond' int2str(condIdx)]),AllP{p}.Tall{eidx}.Data.Artifacts.badTrials);            
            d{eidx}=[d{eidx} mean(AllP{p}.Tall{eidx}.Data.segmentedEcog.smoothed100(ch,samps{e},:,useind),4)];
            ehold=[ehold repmat(eidx,1,length(samps{find([2 4 5]==eidx)}))];
        end
    end
    dcat=horzcat(d{2},d{4},d{5});
    allD.data=vertcat(allD.data,dcat);
    allD.p=vertcat(allD.p,repmat(p,size(dcat,1),1));
    allD.ch=vertcat(allD.ch,ch');
    i=i+1;
end

%% MAKE DATA STRUCTURE OF AVE WAVEFORM FOR SHORT AND LONG WORDS:ALL DATA
allD.dataAll=[];
i=1;
for p=1:8
    ch=1:AllP{p}.Tall{2}.Data.channelsTot;
    ch=setdiff(ch,AllP{p}.Tall{2}.Data.Artifacts.badChannels);
    d=[];
    for eidx=[2 4 5]
        d{eidx}=[];
        e=find([2 4 5]==eidx);
        ind=AllP{p}.Tall{eidx}.Data.findTrials('1','','n');
        ind2=AllP{p}.Tall{eidx}.Data.findTrials('2','y','n');
        ind.cond2=setdiff(ind.cond2,ind2.cond2)
        for condIdx=1:2
            useind=setdiff(ind.(['cond' int2str(condIdx)]),AllP{p}.Tall{eidx}.Data.Artifacts.badTrials);
            d{eidx}=[d{eidx} mean(AllP{p}.Tall{eidx}.Data.segmentedEcog.smoothed100(ch,1:600,:,useind),4)];         
        end
    end
    dcat=horzcat(d{2},d{4},d{5});
    allD.dataAll=vertcat(allD.dataAll,dcat);
end
%% CLUSTER
for knum=3:30
    %%
    [kgroup,centroid]=kmeans(allD.data(:,:),knum,'Distance','correlation','EmptyAction','drop','Replicates',10);
    allD.kgroup=kgroup;
    allD.centroid=centroid;
    S=[];
    for k=1:knum
        idx=find(allD.kgroup==k);
        R=corr(allD.data(idx,:)',centroid');
        S=[S (R(:,k)./max(R(:,setdiff(1:knum,k)),[],2))'];
    end
    hist(S)
    Smean(knum)=mean(S(:));
end
    title(mean(S))


%% KEEP ONLY THOSE WITH HIGH CORRELATION
patients={'EC18','EC16','EC22','EC23','EC24','EC29','EC30','EC31','EC28'}

brainPath='E:\General Patient Info\'
figure(1)
clear qIdxHold
kGroups=[];
for k=1:knum%[3 5:10]
    figure(1)
    idx=find(allD.kgroup==k);   
    if isempty(idx)
        %keyboar
        continue
    end
    allD.centroid(k,:)=mean(allD.data(idx,:),1);
    idx=idx(find(corr(allD.data(idx,:)',(mean(allD.data(idx,:))'))>.7));
    allD.centroid(k,:)=mean(allD.data(idx,:),1);
    if isempty(idx)     
        continue
    end
    if length(unique(allD.p(idx)))>4
        kGroups=[kGroups k];
    end   
    
    pCur=unique(allD.p(idx,:))';
    for p=pCur
        chCur=allD.ch(idx);
        ch=chCur(find(allD.p(idx)==p));
        tmp=find(~ismember(ch,[6 9 70 73 134 137 192 198]));
        ch=ch(tmp);
        idx2=idx(find(allD.p(idx)==p));
        idx2=idx2(tmp);
        q=getBrainQuadrant(BrainCoord(p).xySF,BrainCoord(p).xyCS,BrainCoord(p).xy(:,ch))
        subplot(3,4,p)        
        a=imread([mainPath filesep 'allBrainPics' filesep patients{p} 'scaled.jpg']);
        imshow(a)
        hold on
      
        for qIdx=1:4
            kQ{k}=1:4;
            scatter(BrainCoord(p).xy(1,ch(find(q==qIdx))),BrainCoord(p).xy(2,ch(find(q==qIdx))),50,'b','fill')
            figure(1)
            qIdxHold{k,p,qIdx}=idx2(find(q==qIdx));
        end
    end
    tmp=[];
    l=[];
    for p=1:8
        try
            tmp=vertcat(tmp,allD.dataAll(vertcat(qIdxHold{k,p,:}),:));
            l=[l;size(tmp,1)];
        end
    end          
    
    for e=1:3
        subplot(3,3,e+6)
        imagesc(tmp(:,(e-1)*400+1:(e-1)*400+399),[-1 5])
        hl=line([0 1200],[l l]);        
        set(hl,'Color','k')
    end
    title(k)
    %input('n')
    clf
end


[a,b]=max(allD.centroid(kGroups,210:300),[],2)
[~,idx]=sort(b);
kGroups=kGroups(idx);

kGroups=horzcat(kGroups',repmat([1:4],length(kGroups),1));


%%
idxHold=[]
for k=kGroups(:,1)'
    idxHold=vertcat(idxHold,vertcat(qIdxHold{k,:,:}));
end
silhouette(allD.data(idxHold,:),allD.kgroup(idxHold))

%%
for kIdx=1:length(kGroups)
    tmp=mean([sigHold2{kIdx,1,1};sigHold2{kIdx,2,1}])
    [~,idx]=find(tmp(210:end),1,'first');
    firstSig(kIdx)=idx;
end
[~,sortIdx]=sort(firstSig);
kGroups=kGroups(sortIdx,:)
%% GET GRID POINTS ON COMMON BRAIN

% GET SULCI COORDINATES FROM IMAGE
cd([mailPath filesep 'allBrainPics'])
getSulciFromPic

% GET NORMALIZED COORDINATES
[~,BrainCoord]=gridToNormalizedBrain(BrainCoord,patients)

%% PLOT DATA
plotClusterERPs

%% FIGURE 2: INDIVIDUAL ERPS
load([mainPath filesep 'redblackwhite'])
load([mainPath filesep 'blueblackwhite'])
plotIndividualERPs
%% FIGURE 3: STIMULATION ERROR SITES
load([mainPath filesep 'stimsites'])
load([mainPath filesep 'siteCoord'])
load([mainPath filesep 'stimTaskCoord.mat'])

% READ FROM CELL ARRAY INTO STRUCTURE
subj=unique(cellstr(char(stimsites{:,1})));
sites=unique(cellstr(char(stimsites{:,3})));
sites=sites([3 4 5 2]);
errors={'per','ph','ms','off','nr'}

for i=1:length(subj)
    for j=1:length(sites)
        stimResults(i).subj=subj{i}
        stimIdx=intersect(find(strcmp(stimsites(:,1),subj{i})),find(strcmp(stimsites(:,3),sites{j})));
        if ~isempty(stimIdx)        
            stimResults(i).(sites{j}).total=length(stimIdx);
            stimResults(i).(sites{j}).allError=length(find(~cell2mat(stimsites(stimIdx,7))))
            for k=1:length(errors)
                if ~isempty(errors{k})
                    stimResults(i).(sites{j}).(errors{k})=length(find(strcmp(stimsites(stimIdx,9),errors{k})));
                end
            end
        end        
    end
end

for i=1:length(stimTaskCoord)
    if isnumeric(stimTaskCoord{i,4})
        idx=intersect(find(strcmp(stimsites(:,1),stimTaskCoord{i,1})),...
            find(strcmp(stimsites(:,2),int2str(stimTaskCoord{i,4}))))
    else
        idx=intersect(find(strcmp(stimsites(:,1),stimTaskCoord{i,1})),...
            find(strcmp(stimsites(:,2),stimTaskCoord{i,4})))'
    end
    stimTaskCoord{i,18}=length(idx);
    stimNum=stimsites(idx,9);
    for k=1:length(errors)
        stimTaskCoord{i,18+k}=length(find(strcmp(stimsites(idx,9),errors{k})))/length(stimNum);
    end
end

% CONVERT NUMBERS TO STRINGS
idx=find(cellfun(@isnumeric,siteCoord(:,5)))
siteCoord(idx,5)=cellfun(@int2str,siteCoord(idx,5),'UniformOutput',0)

idx=find(cellfun(@isnumeric,stimsites(:,2)))
stimsites(idx,2)=cellfun(@int2str,stimsites(idx,2),'UniformOutput',0)


%% GET TASK POINTS ON NORMALIZED BRAIN
brainfile=[mainPath filesep 'brain_3Drecon_coordinates2.jpg'];
a=imread(brainfile);
imshow(a)
hold on
xysites=plotPointsBrainNormalized3([stimTaskCoord{:,14}],[stimTaskCoord{:,15}],xySF,xyCS,a);
scatter(xysites(:,1),xysites(:,2))
stimTaskCoord(:,16)=num2cell(double(xysites(:,1)))
stimTaskCoord(:,17)=num2cell(double(xysites(:,2)))

%% PLOT TASK POINTS ON BRAIN
load([mainPath filesep 'hot_adjusted'])
a=imread([mainPath filesep 'brain2.jpg']);
colorjet=flipud(hot);
taskType=unique(stimTaskCoord(:,3))
for k=1:length(taskType)
    idx=find(strcmp(stimTaskCoord(:,3),taskType{k}))
    G(k+length(errors)).points2=[cell2mat(stimTaskCoord(idx,16)) cell2mat(stimTaskCoord(idx,17))];
    G(k+length(errors)).val=ones(1,length(idx))*k;
end
clf
opengl software
ha = axes('units','normalized', ...
            'position',[0 0 1 1]);
imshow(a)
hold on
strings={'S','N','R','C'}
for k=(1:length(taskType))+length(errors)
    hold on    
    text(G(k).points2(:,1), G(k).points2(:,2),strings{k-5},'Color',colorjet(ceil((k-5)*(length(colorjet)/4)),:,:),'FontWeight','bold','FontSize',12)
end

%% FIGURE 4: STIMULATION PROBABILITY MAPS
% GET ERROR POINTS ON NORMALIZED BRAIN
load([mainPath filesep 'sulciCoord'])
clf
brainfile=[mainPath filesep 'brain_3Drecon_coordinates2.jpg'];
a=imread(brainfile);
imshow(a)
hold on
xySF=sulciCoord.xySF;
xyCS=sulciCoord.xyCS;
clear xysites
xysites=plotPointsBrainNormalized3([stimTaskCoord{:,14}],[stimTaskCoord{:,15}],xySF,xyCS,a);

scatter(xysites(:,1),xysites(:,2))
stimTaskCoord(:,16)=num2cell(double(xysites(:,1)))
stimTaskCoord(:,17)=num2cell(double(xysites(:,2)))

%% PLOT GRIDDATA POINTS ON BRAIN: ERROR TYPE
im=imread([mainPath filesep 'brain2.jpg']);
for k=1:length(errors) 
    figure
    interpDataBrain(im,cell2mat(stimTaskCoord(useidx,16:17))',cell2mat(stimTaskCoord(useidx,k+18))',hot_adjusted)
    input('n')
    close
end

%% SUPPLEMENTAL: ACTIVATION MOVIES

%% MAKE MOVIE PERC AND PROD ON ONE PLOT

opengl software
%figure
colorjet=flipud(hot)
stepSize=4
maxZ=3;
for p=8
        close all
        figure('position',[ 100 31 800 800],'units','Normalize','Color','w')
        for s=1:400/stepSize
            clf
            for eidx=[2 5]
                %%
                r=rem(s,10);
                if r==0
                    r=10;
                end
                pos=subplotMinGray(1,1,1,0);
                pos(4)=pos(4)/2;
                if eidx==2
                    pos(2)=pos(2)+pos(4);
                end
                samps=[s*stepSize:(s+1)*stepSize];
                usetr=1:size(AllP{p}.Tall{eidx}.Data.segmentedEcog.event,1);
                usech=setdiff(1:AllP{p}.Tall{2}.Data.channelsTot,AllP{p}.Tall{2}.Data.Artifacts.badChannels);
                d=mean(mean(AllP{p}.Tall{eidx}.Data.segmentedEcog.smoothed100(usech,samps,:,usetr),4),2);
                d(find(d<0))=0;
                im=imread([mainPath filesep 'allBrainPics' filesep AllP{p}.Tall{5}.Data.patientID 'scaled.jpg']);
                interpDataBrain(im,BrainCoord(p).xy(:,usech),d',hot_adjusted,20,30,20,pos,0,maxZ)
                hold on                
                if s>200
                    stext=-(samps(1)-100)*10;
                else
                    stext=(samps(1)-200)*10;
                end
                
                if eidx==2
                    h=text(pos(1),pos(2)+50,[num2str(stext) ' ms'],'Fontsize',16);    
                    if stext>=0
                        set(h,'BackgroundColor','y')
                    end
                end
                if eidx==2
                    text(-20,pos(2)+200,['Perception'],'Fontsize',12);
                else
                    text(-20,pos(2)+200,['Production'],'Fontsize',12);
                end
            end
            M(s)=getframe(gcf);
        end
        input('n')
        MovieHold{p,flag,eidx}.M=M;    
end

