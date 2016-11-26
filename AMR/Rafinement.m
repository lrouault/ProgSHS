function [ Tnext ] = Rafinement(T,Num,nivRaf,dt,dx)

% TR : zoom sur la maille a raffiner

%% Calcul de la taille de TR
nbpoint = 2;
for i = 1:nivRaf
    nbpoint = nbpoint + 2^(i-1);
end

TR = zeros (nbpoint,1);

%% Remplissage de TR
TR(1)=(T(Num-1)+T(Num))/2;
TR(nbpoint)=(T(Num)+T(Num+1))/2;

% a=1;
% b=nbpoint;
% d=b;
% n=ceil(b/2);
% c=n;
% for i=1:nivRaf
%     for j=1:2^(i-1)
%         TR(c)=(TR(a)+TR(b))/2;
%         a=b;
%         c=c+d-1;
%         b=b+d-1;
%     end
%     a=1;
%     b=n;
%     d=n;
%     n=ceil(b/2);
%     c=n;
% end

a=1;
b=ceil(nbpoint/2);
TR(b)=T(Num);
d=b;
n=ceil(b/2);
c=n;
for i=2:nivRaf
    for j=1:2^(i-1)
        TR(c)=(TR(a)+TR(b))/2;
        a=b;
        c=c+d-1;
        b=b+d-1;
    end
    a=1;
    b=n;
    d=n;
    n=ceil(b/2);
    c=n;
end


%% DF
% TR=[T(Num-1);TR;T(Num+2)];
% vectRaf=ones(size(TR,1),1)*nivRaf;
% vectRaf(1)=0;
% vectRaf(size(TR,1),1)=0;
% TRnext=TR;
% 
% for i=2:nbpoint+1
%     dxm=dx/(2^vectRaf(i-1));
%     dxp=dx/(2^vectRaf(i+1));
%     TRnext(i)=DF(TR(i-1),TR(i),TR(i+1),dt,dxp,dxm);  %pb  dt/(2^vectRaf(i))
% end

TRnext=TR;
dx=dx/(2^nivRaf);
for i=2:nbpoint-1
    TRnext(i)=DF(TR(i-1),TR(i),TR(i+1),dt,dx,dx);
end

%% Choix de ce que l'on renvoie
%cas=1;
% if cas==1
%     Tnextk1=TRnext(2);
%     Tnextk2=TRnext(size(TRnext,1)-1);
% elseif cas==2
%     poids=2^nivRaf;
%     Tnextk1=((TRnext(1)+poids*TRnext(3))/(poids+1)+TRnext(2))/2;
%     Tnextk2=((TRnext(size(TRnext,1)-2)*poids+TRnext(size(TRnext,1)))/(poids+1)+TRnext(size(TRnext,1)-1))/2;
% end




Tnext=TRnext(ceil(nbpoint/2));

end

