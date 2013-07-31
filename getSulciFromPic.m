
for p=1:8
    clf
    b=imread([patients{p} 'CS.jpg']);
    imagesc(b)
    hold on
    bcopy=b(:,:,1);
    bcopy(find(ceil(bcopy/10)))=255;
    %bcopy(find(floor(bcopy/10)))=0;
    imagesc(bcopy)
    [r,c]=find(diff(bcopy')'>.001)
    plot(c,r,'.')
    %%
    imagesc(b)
    hold on
    edgeC=[]
    uniqueR=unique(r)
    for idx=1:length(uniqueR)
        u=uniqueR(idx)
        i=find(r==u)
        edgeC(idx)=max(c(i));
    end
    
    plot(smooth(edgeC,15),uniqueR,'r')
    BrainCoord(p).xyCS=vertcat(smooth(edgeC,10)',uniqueR');
    input('n')
    %%
    
      b=imread([patients{p} 'SF.jpg']);      
      clf
    bcopy=b(:,:,1);
    bcopy(find(ceil(bcopy/10)))=255;
    bcopy(find(floor(bcopy/10)))=0;
    [r,c]=find(diff(bcopy)>.1)
    plot(c,r,'.')
    %%
    imagesc(b)
    hold on
    edgeC=[]
    uniqueR=unique(c)
    for idx=1:length(uniqueR)
        u=uniqueR(idx)
        i=find(c==u)
        edgeC(idx)=max(r(i));
    end
    plot(uniqueR,smooth(edgeC,20)','r')
    BrainCoord(p).xySF=vertcat(uniqueR',smooth(edgeC,10)');
    input('n')

    
end
