function interpDataBrain(im, coord, data,colormat,grain,grain2,rad,pos,flag,maxVq)

limx=size(im,2);
limy=size(im,1);
%grain=1;
gx=1:grain:limx;
gy=1:grain:limy;
[X,Y] = MESHGRID(gx,gy);
Z=zeros(size(X));

%grain2=10;
gx2=1:grain2:limx;
gy2=1:grain2:limy;
[X2,Y2] = MESHGRID(gx2,gy2);
Z2=zeros(size(X2));

clf
%figure
colorjet=flipud(hot);
%idx=ceil(rand(1,500)*size(X2,1)*size(X2,2));
idx=1:size(X2,1)*size(X2,2);
X3=reshape(X2,1,[]);
Y3=reshape(Y2,1,[]);

XY2=vertcat(X3,Y3);
set(gcf,'Color','w')
nidx=[];

for i=1:size(coord,2)
    d=squareform(pdist([coord(:,i) XY2]'));
    d2=d(1,2:end);
    
    [aidx,bidx]=find(d2<rad);
    nidx=[nidx bidx];
end
idx=setdiff(1:length(X3),unique(nidx));

[Xq,Yq,Vq]=griddata([X3(idx) coord(1,:)],[ Y3(idx) coord(2,:)],[reshape(Z2(idx),1,[]) data],X,Y);


grain=1
Vq(find(isnan(Vq)))=0;
tmp2=imresize(Vq,[size(im,1)*10 size(im,2)*10]);
tmp2(isnan(tmp2))=0;
tmp3=imresize(tmp2,[size(im,1) size(im,2)]);
tmp1=tmp3;
tmp3(find(isnan(tmp3)))=0;
tmp3(find(tmp3<0))=0;
%tmp4=tmp3./maxVq;
tmp4=smooth2(tmp3,10);
tmp4=tmp4./maxVq;
ha = axes('units','normalized', ...
    'position',pos);
try
    imshow(repmat(im,[1 1 3]))
catch
    imshow(im)
end
%imshow(a)
hold on
set(ha,'handlevisibility','off', ...
    'visible','off')
ha = axes('units','normalized', ...
    'position',pos);
imshow(tmp4)

if flag==1
    [Xq,Yq,Vq]=griddata([X3(idx) coord(1,:)],[ Y3(idx) coord(2,:)],[reshape(Z2(idx),1,[]) data+1],X,Y);
    Vq(find(Vq>.4))=1;
    tmp6=smooth2(Vq,5)./(max(max(smooth2(Vq,5))));
    tmp6(find(tmp6>.8))=.8;
    alpha(tmp6)
else
    tmp4(find(tmp4>1))=1;
    alpha(tmp4)
end

if nargin>3
    colormap(flipud(colormat))
else
    colormap(flipud(hot))
end


