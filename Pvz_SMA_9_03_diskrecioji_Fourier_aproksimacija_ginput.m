%
%   Viename periode aprasytai funkcijai atliekama diskrecioji Furje
%   transformacija
%   Funkcija ivedama taskais i brezini, ginput
%

function main
clc,close all,clear all


T=10
n=1000;n=round(n/2)*2+1; % n visuomet nelyginis
dt=T/n

xmin=0;ymin=-10;xmax=T;ymax=10; % koordinaciu sistemos ribos
figure(1),hold on, axis([xmin,xmax,ymin,ymax]);grid on
title('vienas duotosios funkcijos periodas')

% Pele ivedami taskai. Baigiama, kai taskas parenkamas uz koord. sistemos ribu
ibreak=0;icount=1;
x1=0;y1=0;plot(x1,y1,'*');tinp(1)=x1;fffinp(1)=y1;
while 1     % pele ivedami taskai
    [x1,y1]=ginput(1); icount=icount+1;  
    if x1 < xmin || x1 > xmax || y1 < ymin || y1 > ymax,  x1=T; y1=0; ibreak=1; end
    tinp(icount)=x1; fffinp(icount)=y1;
    plot(x1,y1,'*');
    if ibreak, break; end
end
t=[0:dt:T-dt]; fff=interp1(tinp,fffinp,t); cla, plot(t,fff,'k.');

m=(n+1)/2  % m - harmoniku skaicius interpoliuojanciai Furje aproksimacijai, 


N=10000 % vaizdavimo tasku skaicius
dttt=T/N
ttt=[-T:dttt:2*T];

% disp('kontrole:'),disp(sum(fC(3,T,t).*fC(0,T,t)))

figure(2),hold on,grid on,plot(t,fff,'b.-','MarkerSize',8);
legend(sprintf('n=%d tasku, m=%d harmoniku',n,m))

ac0=dot(fff,fC(0,T,t))/n;
for i=1:m-1
    ac(i)=dot(fff,fC(i,T,t))*2/n;
    as(i)=dot(fff,fS(i,T,t))*2/n;
end
ac,as

figure(3),hold on
bar(0:m-1,[ac0,sqrt(ac.^2+as.^2)],0.7);
legend(sprintf('n=%d tasku, m=%d harmoniku ',n,m))


fffz=ac0*fC(0,T,ttt);
frequencies=[1:m-1];
% frequencies=[1:3,5,6,10,30,50];
for i=frequencies
        fffz=fffz+ac(i)*fC(i,T,ttt)+as(i)*fS(i,T,ttt);    
end

figure(4),hold on,grid on, plot(ttt,fffz,'r');plot(t,fff,'b.','LineWidth',2);
legend(sprintf('n=%d tasku, m=%d harmoniku',n,m))

return
end

function c=fC(i,T,t), if i==0,c=1*cos(0*t); else, c=cos(2*pi*i/T*t); end, return, end
function s=fS(i,T,t), s=sin(2*pi*i/T*t); return, end

