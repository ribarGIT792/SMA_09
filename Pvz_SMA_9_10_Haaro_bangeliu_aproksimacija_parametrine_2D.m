%
%   Parametrine Haro bangeliu transformacija
%

function main
clc;close all;clear all;
spalvos={'r-','g-','m-','c-','k-','y-','r.','g.','m.','c.','k.','y.'};
% Is failu ivedami duomenys:6
n=6   
nnn=2^n;
fclose all; fhx=fopen('chairx.txt','r'); fhy=fopen('chairy.txt','r');
figure(1); axis equal,hold on,grid on
SX=fscanf(fhx,'%g '); SY=fscanf(fhy,'%g ');  
% pervedimas i paramatrini pavidala:
nP=length(SX);
t(1)=0; for i=2:nP, t(i)=t(i-1)+norm([SX(i) SY(i)]-[SX(i-1) SY(i-1)]);  end
fclose all; plot(SX,SY,'c');

a=min(t),b=max(t),t1=[a:(b-a)/(nnn-1):b];
ts=interp1(t,SX,t1); clear SX,  SX=ts;
ts=interp1(t,SY,t1); clear SY,  SY=ts; t=t1;
plot(SX,SY,'k.'); plot(t,SX,'b.');plot(t,SY,'g.');
title(sprintf('duota funkcija, tasku skaicius 2^%d',n));
xmin=min(SX);xmax=max(SX);
ymin=min(SY);ymax=max(SY);


m=6  % detalumo lygiu skaicius
[smoothx,detailsx]=Haar_wavelet_approximation(t,SX,n,m);
[smoothy,detailsy]=Haar_wavelet_approximation(t,SY,n,m);

smoothx
smoothy

% Funkcijos rekonstrukcija:

hx=zeros(1,nnn); hy=zeros(1,nnn);
for k=0:2^(n-m)-1   % suglodinta funkcija
    hx=hx+smoothx(k+1)*Haar_scaling(t,n-m,k,a,b);
    hy=hy+smoothy(k+1)*Haar_scaling(t,n-m,k,a,b); 
end  
leg={sprintf('suglodinta funkcija, detalumo lygmuo %d',n-m)};
figure(2);subplot(m+1,1,1), axis equal,axis([xmin xmax ymin ymax]),hold on,grid on, plot(hx,hy,'.','Linewidth',2);title(sprintf('lygyje %d suglodinta funkcija',0));

for i=0:m-1 %detalumo didinimo ciklas
    % apskaiciuojamos funkcijos detales:  
    h1x=zeros(1,nnn); h1y=zeros(1,nnn); 
    for k=0:2^(n-m+i)-1, 
        h1x=h1x+detailsx{m-i}(k+1)*Haar_wavelet(t,n-m+i,k,a,b); 
        h1y=h1y+detailsy{m-i}(k+1)*Haar_wavelet(t,n-m+i,k,a,b);
    end
    figure(3),subplot(m,1,i+1), axis equal,hold on,grid on 
    xshift=(xmin+xmax)/2; yshift=(ymin+ymax)/2; axis([xmin-xshift xmax-xshift ymin-yshift ymax-yshift])
    plot(h1x,h1y,'b-','Linewidth',2);title(sprintf('%d lygio detales',i));
    leg={leg{1:end},sprintf('lygmens %d detales',n-m+i)};
    hx=hx+h1x; hy=hy+h1y; % detales pridedamos prie ankstesnio suglodinto vaizdo
    figure(2);subplot(m+1,1,i+2), axis equal,axis([xmin xmax ymin ymax]),hold on,grid on, plot(hx,hy,'Linewidth',2);title(sprintf('lygyje %d suglodinta funkcija' ,i+1));
end

return
end


function h=Haar_scaling(x,j,k,a,b)   % ***********************************************************
eps=1e-9; % Parametrines aproksimacijos atveju tektu grizti Haro mastelio funkcijos priekinis ir 
          % galinis frontai nevaizduojami
xtld=(x-a)/(b-a); % (a,b) intervale duota kintamojo reiksme perskaiciuojama i "standartini" 
                        % intervala (0,1), kuriame uzrasyta bangeles formule  
xx=2^j*xtld-k; h=2^(j/2)*(sign(xx+eps)-sign(xx-1-eps))/(2*(b-a));
return
end

function h=Haar_wavelet(x,j,k,a,b)   % ************************************************************
eps=1e-9;
xtld=(x-a)/(b-a); % (a,b) intervale duota kintamojo reiksme perskaiciuojama i "standartini" 
                        % intervala (0,1), kuriame uzrasyta bangeles formule  
xx=2^j*xtld-k; h=2^(j/2)*(sign(xx+eps)-2*sign(xx-0.5)+sign(xx-1-eps))/(2*(b-a));
return
end

function [smooth,details]=Haar_wavelet_approximation(SX,SY,n,m)
% Aproksimavimas Haro bangelemis:
a=min(SX);b=max(SX);
nnn=2^n;
smooth=(b-a)*SY*2^(-n/2); % auksciausio detalumo suglodinimas (pagal duota funkcija)

for i=1:m
    smooth1=(smooth(1:2:end)+smooth(2:2:end))/sqrt(2);
    details{i}=(smooth(1:2:end)-smooth(2:2:end))/sqrt(2);
    fprintf(1,'\n details %d :  ',i);fprintf('%g ', details{i});
    smooth=smooth1;
end
fprintf(1,'\n smooth  %d :  ',i);fprintf('%g ', smooth);fprintf('\n');

return,end
