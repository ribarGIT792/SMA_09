%
%   Viename periode aprasytai funkcijai atliekama diskrecioji Furje
%   transformacija. Duomenys imami is failo
%

function main
clc,close all,clear all

% Is failu ivedami duomenys:
npower=9   
n=2^npower-1;
fclose all; fhx=fopen('carx.txt','r'); fhy=fopen('cary.txt','r');
figure(100); axis equal,hold on,grid on
SX=fscanf(fhx,'%g '); SY=fscanf(fhy,'%g '); fclose all; 
ppp=1*60,   % sukuriamas "atstumas" tarp signalu  %********************************** 
if ppp~= 0   % sukuriamas vaizdas "ilgame periode"
    SX(2:end+1)=SX(1:end);SX(1)=SX(2)-ppp;SX(end+1)=SX(end)+ppp;
    SY(2:end+1)=SY(1:end);SY(1)=0;SY(end+1)=0;
end
plot(SX,SY); 
a=min(SX);b=max(SX);t=[a:(b-a)/n:b];
fff=interp1(SX,SY,t);   % perskaiciuojama i naujas abscises

plot(t,fff,'r.');  
title(sprintf('duota funkcija, tasku skaicius 2^%d',npower));

m=floor((n+1)/2)  % m - harmoniku skaicius
% m=32   % apribotas harmoniku skaicius   %**********************************
 
T=b-a;
slenkstis=0 ; % harmoniku amplitudziu slenkstis triuksmu filtravimui
dt=T/n
N=1000 % vaizdavimo tasku skaicius
dttt=T/N


ttt=[a-T:dttt:b+T];

% disp('kontrole:'),disp(sum(fC(3,T,t).*fC(0,T,t)))

figure(1),hold on,grid on,axis equal, plot(t,fff,'b.-','MarkerSize',8);
legend(sprintf('n=%d tasku, m=%d harmoniku',n,m))

ac0=dot(fff,fC(0,T,t))/n;
for i=1:m-1
    ac(i)=dot(fff,fC(i,T,t))*2/n;
    as(i)=dot(fff,fS(i,T,t))*2/n;
end
ac,as

figure(2),hold on
bar(0:m-1,[ac0,sqrt(ac.^2+as.^2)],0.01)
xx=axis; plot([xx(1),xx(2)],slenkstis*[1 1],'m--','LineWidth',3); % braizo slenkscio linija
legend(sprintf('n=%d tasku, m=%d harmoniku, slenkstis=%g ',n,m,slenkstis))


fffz=ac0*fC(0,T,ttt)
frequencies=[1:m-1];
% frequencies=[1:3,5,6,10,30,50];
for i=frequencies
    if sqrt(ac(i)^2+as(i)^2) > slenkstis
        fffz=fffz+ac(i)*fC(i,T,ttt)+as(i)*fS(i,T,ttt);    
    end
end

figure(3),hold on,grid on,axis equal, plot(ttt,fffz,'r','LineWidth',2);plot(t,fff,'b-','LineWidth',1);
legend(sprintf('n=%d tasku, m=%d harmoniku, slenkstis=%g ',n,m,slenkstis))

return
end

function c=fC(i,T,t), if i==0,c=1*cos(0*t); else, c=cos(2*pi*i/T*t); end, return, end
function s=fS(i,T,t), s=sin(2*pi*i/T*t); return, end

