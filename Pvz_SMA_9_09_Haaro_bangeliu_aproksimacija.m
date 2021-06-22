%
%   Haro bangeliu aproksimacija
%

function main
clc;close all;clear all;
spalvos={'r-','g-','m-','c-','k-','y-','r.','g.','m.','c.','k.','y.'};

% Is failu ivedami duomenys: 
n=8   
nnn=2^n;
fclose all; fhx=fopen('carx.txt','r'); fhy=fopen('cary.txt','r');
figure(1); axis equal,hold on,grid on
SX=fscanf(fhx,'%g '); SY=fscanf(fhy,'%g '); fclose all; plot(SX,SY); 
a=min(SX),b=max(SX),t=[a:(b-a)/(nnn-1):b];
ts=interp1(SX,SY,t);
clear SX SY, SX=t;SY=ts;plot(SX,SY,'r.');   
title(sprintf('duota funkcija, tasku skaicius 2^%d',n));
xmin=min(SX);xmax=max(SX);
ymin=min(SY);ymax=max(SY);


% Aproksimavimas Haro bangelemis:
m=n  % detalumo lygiu skaicius
smooth=(b-a)*SY*2^(-n/2); % auksciausio detalumo suglodinimas (pagal duota funkcija)

for i=1:m
    smooth1=(smooth(1:2:end)+smooth(2:2:end))/sqrt(2);
    details{i}=(smooth(1:2:end)-smooth(2:2:end))/sqrt(2);
    fprintf(1,'\n details %d :  ',i);fprintf('%g ', details{i});
    smooth=smooth1;
end
fprintf(1,'\n smooth  %d :  ',i);fprintf('%g ', smooth);fprintf('\n');
% Funkcijos rekonstrukcija:

h=zeros(1,nnn); for k=0:2^(n-m)-1, h=h+smooth(k+1)*Haar_scaling(SX,n-m,k,a,b); end  % suglodinta funkcija
leg={sprintf('aproksimuotas signalas, detalumo lygmuo %d',n-m)};
figure(2);subplot(m+1,1,1),axis equal,axis([xmin xmax ymin ymax]); hold on,grid on, plot(SX,h,'Linewidth',2);title(sprintf('lygyje %d aproksimuotas signalas',0));

for i=0:m-1 %detalumo didinimo ciklas
    % apskaiciuojamos funkcijos detales:  
    h1=zeros(1,nnn); for k=0:2^(n-m+i)-1, h1=h1+details{m-i}(k+1)*Haar_wavelet(SX,n-m+i,k,a,b);  end
    figure(3),subplot(m,1,i+1), axis equal,hold on,grid on
    yshift=(ymin+ymax)/2;axis([xmin xmax ymin-yshift ymax-yshift]), plot(SX,h1,'b-','Linewidth',2);title(sprintf('%d lygio detales',i));
    leg={leg{1:end},sprintf('lygmens %d detales',n-m+i)};
    h=h+h1; % detales pridedamos prie ankstesnio suglodinto vaizdo
    figure(2);subplot(m+1,1,i+2),axis equal,axis([xmin xmax ymin ymax]), hold on,grid on, plot(SX,h,'Linewidth',2);title(sprintf('lygyje %d aproksimuotas signalas' ,i+1));
end

return
end


function h=Haar_scaling(x,j,k,a,b)   % ***********************************************************
eps=1e-9;
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

