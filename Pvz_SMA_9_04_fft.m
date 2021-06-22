%
%   Viename periode aprasytai funkcijai atliekama diskrecioji Furje
%   transformacija
%
%   Sulyginimas su fft rezultatais 

function main
clc,close all,clear all
n=1000; n=round(n/2)*2+1 % tasku skaicius, visuomet nelyginis
m=(n+1)/2  % m - harmoniku skaicius interpoliuojanciai Furje aproksimacijai

T=10;   % periodas
dazniu_slenkstis=m  %50;    % harmoniku numeriu slenkstis triuksmu filtravimui 
ampl_slenkstis=0.2*0  % harmoniku amplitudziu slenkstis triuksmu filtravimui
dt=T/n  % zingsnis 
N=1000  % vaizdavimo tasku skaicius
dttt=T/N  % vaizdavimo zingsnis

t=[0:dt:T-dt];
ttt=[-T:dttt:2*T];

% disp('kontrole:'),disp(sum(fC(3,T,t).*fC(0,T,t)))

fff=fnk(T,t); % apskaiciuojame ir pavaizduojame duota tasku seka 
figure(1),hold on,grid on,plot(t,fff,'b.-','MarkerSize',8);
legend(sprintf('n=%d tasku, m=%d harmoniku',n,m))
title('vienas duotosios funkcijos periodas')

% Spektras apskaiciuojamas naudojant fft:  ************************
yyy=fft(fff)/n;
spektras=abs(2*yyy(1:m));    % harmoniku amplitudes 
spektras(1)=spektras(1)/2;   % pastovi dedamoji 

spektras_c0=real(yyy(1));   % pastovi dedamoji
spektras_c=real(2*yyy(2:m)); % cos amplitudes 
spektras_s=-imag(2*yyy(2:m)); % sin amplitudes

% ******************************************************************

figure(2),hold on,grid on,plot([0:m-1],spektras)
xx=axis; 
plot(xx(1:2),ampl_slenkstis*[1 1],'m--','LineWidth',3); % braizo ampl slenkscio linija
plot(dazniu_slenkstis*[1 1],xx(3:4),'g--','LineWidth',3); % braizo dazniu slenkscio linija
title('Furje amplitudziu spektras pagal kompleksinio skaiciaus moduli ')
legend({'amplitudziu spektras';'amplitudes slenkstis';'dazniu slenkstis'})

% apskaiciuojame cos ir sin koeficientus:
ac0=dot(fff,fC(0,T,t))/n;
for i=1:m-1
    ac(i)=dot(fff,fC(i,T,t))*2/n;
    as(i)=dot(fff,fS(i,T,t))*2/n;
end
ac,as

figure(3),hold on
bar(0:m-1,[ac0,sqrt(ac.^2+as.^2)],0.01)
legend(sprintf('n=%d tasku, m=%d harmoniku, a-slenkstis=%g  d-slenkstis=%g',n,m,ampl_slenkstis,dazniu_slenkstis))
title('Furje amplitudziu spektras pagal DFT sin ir cos amplitudes')

fffz=ac0*fC(0,T,ttt);
frequencies=[1:m-1];
% frequencies=[1:3,5,6,10,30,50];
frequencies=frequencies(find(frequencies < dazniu_slenkstis)) 

for i=frequencies
    if sqrt(ac(i)^2+as(i)^2) > ampl_slenkstis  %  filtravimas pagal amplitude
        fffz=fffz+ac(i)*fC(i,T,ttt)+as(i)*fS(i,T,ttt);    
    end
end

figure(4),hold on,grid on, plot(ttt,fffz,'r');plot(t,fff,'b-','LineWidth',2);
legend(sprintf('n=%d tasku, m=%d harmoniku, a-slenkstis=%g  d-slenkstis=%g',n,m,ampl_slenkstis,dazniu_slenkstis))
title('Pagal harmonikas atstatytas funkcijos grafikas')


% Sulyginimui pavaizduojame cos ir sin amplitudes gautas fft ir DFT:
figure(5),hold on,grid on
bar(0:m-1,[ac0,ac],0.01);plot(0:m-1,[spektras_c0,spektras_c])
title('FFT ir DFT gautu cos amplitudziu sulyginimas');
legend({'FFT','DFT'})
figure(6),hold on,grid 
figure(6),hold on,grid on
bar(1:m-1,as,0.01);plot(1:m-1,spektras_s)
title('FFT ir DFT gautu sin amplitudziu sulyginimas');
legend({'FFT','DFT'})
figure(6),hold on,grid 

return
end

function c=fC(i,T,t), if i==0,c=1*cos(0*t); else, c=cos(2*pi*i/T*t); end, return, end
function s=fS(i,T,t), s=sin(2*pi*i/T*t); return, end
% function rez=fnk(T,t), rez=sign(sin(2*pi*t/T))+sign(cos(2*pi*2*t/T))+sign(cos(2*pi*10*t/T)); return, end
% function rez=fnk(T,t), rez=sin(2*pi*t/T)+cos(2*pi*2*t/T)+cos(2*pi*10*t/T); return, end
function rez=fnk(T,t), rez=sin(2*pi*t/T)+cos(2*pi*2*t/T)+0.1*cos(2*pi*100*t/T)+0.05*cos(2*pi*120*t/T)+0.8; return, end
