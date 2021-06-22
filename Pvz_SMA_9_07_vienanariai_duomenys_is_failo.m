%
%   Aproksimuojama vienanariu bazeje. Duomenys imami is failo
%   Parodoma, kokia yra aproksimuojancios funkcijos elgsena 
%   uz duomenu intervalo ribu   

function main
clc,close all,clear all

% Is failu ivedami duomenys:
npower=9   
n=2^npower-1;
fclose all; fhx=fopen('carx.txt','r'); fhy=fopen('cary.txt','r');
figure(100); axis equal,hold on,grid on
SX=fscanf(fhx,'%g '); SY=fscanf(fhy,'%g '); fclose all; 
ppp=0*60,   % sukuriamas "atstumas" tarp signalu  %********************************** 
if ppp~= 0   % sukuriamas vaizdas "ilgame periode"
    SX(2:end+1)=SX(1:end);SX(1)=SX(2)-ppp;SX(end+1)=SX(end)+ppp;
    SY(2:end+1)=SY(1:end);SY(1)=0;SY(end+1)=0;
end
plot(SX,SY); 
a=min(SX);b=max(SX);t=[a:(b-a)/n:b];
fff=interp1(SX,SY,t);   % perskaiciuojama i naujas abscises

plot(t,fff,'r.');  
title(sprintf('duota funkcija, tasku skaicius 2^%d',npower));

pause
% Maziausiu kvadratu metodo lygciu sistema:
m= 5;
G=base(m,t)
c=(G'*G)\(G'*fff')
sss=sprintf('%5.2g',c(1));
for i=1:m-1
    sss=[sss,sprintf('+%5.2gx^%1d',c(i+1),i)];
end
sss=strrep(sss,'+-','-');
% Aproksimuojanti funkcija:
nnn=200; %vaizdavimo tasku skaicius
tmin=min(t);tmax=max(t);
% tmin=-14;tmax=28;
ttt=[tmin:(tmax-tmin)/(nnn-1):tmax]; %vaizdavimo taskai 
Gv=base(m,ttt);
fff1=Gv*c;
plot(ttt,fff1,'b.');
legend({'duoti taskai',sprintf('f(x)=%s',sss)})
title(sprintf('aproksimavimas maziausiu kvadratu metodu \n  tasku skaicius %d,  funkciju skaicius  %d',n,m));

end

function G=base(m,x)
     for i=1:m,  G(:,i)=x.^(i-1); end
return
end