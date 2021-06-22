%
%   Viename periode aprasytai funkcijai atliekama diskrecioji Furje
%   transformacija
%

function main
clc,close all,clear all
n=1000 ; n=round(n/2)*2+1 % tasku skaicius, nelyginis
m=(n+1)/2  % m - harmoniku skaicius
 m=50
T=1;
dazniu_slenkstis= m  %50 %m/2    % harmoniku numeriu slenkstis triuksmu filtravimui 
ampl_slenkstis=0 %0.12   %0.6  % harmoniku amplitudziu slenkstis triuksmu filtravimui
dt=T/n
N=10000  % vaizdavimo tasku skaicius
dttt=T/N

t=[0:dt:T-dt];
ttt=[-T:dttt:2*T];

% disp('kontrole:'),disp(sum(fC(3,T,t).*fC(0,T,t)))

fff=fnk(T,t); % apskaiciuojame ir pavaizduojame duota tasku seka 
figure(1),hold on,grid on,plot(t,fff,'b.-','MarkerSize',8);
legend(sprintf('Vienas funkcijos periodas n=%d tasku',n))
title('Duotoji funkcija')

ac0=dot(fff,fC(0,T,t))/n;
for i=1:m-1
    ac(i)=dot(fff,fC(i,T,t))*2/n;
    as(i)=dot(fff,fS(i,T,t))*2/n;
end
ac,as

figure(2),hold on
bar(0:m-1,[ac0,sqrt(ac.^2+as.^2)],0.8)
xx=axis; 
plot(xx(1:2),ampl_slenkstis*[1 1],'m--','LineWidth',3); % braizo ampl slenkscio linija
plot(dazniu_slenkstis*[1 1],xx(3:4),'g--','LineWidth',3); % braizo dazniu slenkscio linija
title('Furje amplitudziu spektras pagal kompleksinio skaiciaus moduli ')
legend({'amplitudziu spektras';'amplitudes slenkstis';'dazniu slenkstis'})


fffz=ac0*fC(0,T,ttt)
frequencies=[1:m-1];
frequencies=frequencies(find(frequencies < dazniu_slenkstis))
for i=frequencies
    if sqrt(ac(i)^2+as(i)^2) > ampl_slenkstis
        fffz=fffz+ac(i)*fC(i,T,ttt)+as(i)*fS(i,T,ttt);    
    end
end

figure(3),hold on,grid on, plot(ttt,fffz,'r'); 
% plot(t,fff,'b-','LineWidth',1);
plot(t,fff,'b.-','MarkerSize',8);
legend({sprintf('n=%d tasku, m=%d harmoniku, a-slenkstis=%g  d-slenkstis=%g',n,m,ampl_slenkstis,dazniu_slenkstis),'iseities signalas'})
title(sprintf('Pagal Furje spektra atkurta funkcija'))

return
end

function c=fC(i,T,t), if i==0,c=1*cos(0*t); else, c=cos(2*pi*i/T*t); end, return, end
function s=fS(i,T,t), s=sin(2*pi*i/T*t); return, end

%------------------------------------------------------------------------------------------------------
% function rez=fnk(T,t), rez=sin(2*pi*t/T); return, end
% function rez=fnk(T,t), rez=sin(2*pi*t/T)+cos(2*pi*2*t/T)+sin(2*pi*10*t/T); return, end
% function rez=fnk(T,t), rez=sign(sin(2*pi*t/T))+sign(cos(2*pi*2*t/T))+sign(cos(2*pi*10*t/T)); return, end
function rez=fnk(T,t), rez=sin(2*pi*t/T)+cos(2*pi*2*t/T)+0.05*cos(2*pi*100*t/T)+0.1*cos(2*pi*120*t/T)+0.8; return, end
