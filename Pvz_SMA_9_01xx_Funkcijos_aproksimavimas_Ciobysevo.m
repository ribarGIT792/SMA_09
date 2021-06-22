
% Lagranzo_interpoliavimas_1D_pagal_duota_funkcija
% Programa demonstruoja, kaip skiriasi interpoliavimo kokybe, 
% kai interpoliavimui parenkame N tolygiai paskirstytu tasku, 
% ir kai interpoliavimo taskais parenkame "Ciobysevo abscises"

function pagrindine
clc,close all

xmin=-1.2;xmax=1.2;  % duotas funkcijos apibrezimo intervalas 
N=5;   % interpoliavimo tasku skaicius
X=[xmin:(xmax-xmin)/(N-1):xmax];  % tolygiai paskirstytu interpoliavimo tasku abscises
k=[0:N-1];
XC=(xmax+xmin)/2+(xmax-xmin)/2*cos((2*k+1)*pi/(2*N)); % "Ciobysevo abscises"

Y=funkcija(X);     %  tolygiai paskirstytu interpoliavimo tasku ordinates
YC=funkcija(XC);   %  ordinates "Ciobysevo abscisiu" taskuose
x=min(X):(max(X)-min(X))/1000:max(X);   %x reiksmes vaizdavimui

leg={'duota funkcija',... 
    'tolygiai isdestyti mazgai',...
    'interpoliavimas per tolygiai isdestytus mazgus',...
    'netiktis interpoliuojant per tolygiai isdestytus mazgus',...
    'Ciobysevo abscises',...
    'interpoliavimas per Ciobysevo mazgus',...
    'netiktis interpoliuojant per tolygiai isdestytus mazgus'};

figure(1), hold on, grid on,box on,axis equal, set(gcf,'Color','w'); 
plot(x,funkcija(x),'b-')   % vaizduojama duotoji funkcija (t.y. pagal kuria interpoliuojama) 
hg=legend(leg{1});pause
plot(X,Y,'ro','MarkerFaceColor','r','MarkerSize',8) % vaizduojami tolygiai isdestyti interpoliavimo taskai 
delete(hg);hg=legend(leg{1:2});pause
F=0;
FC=0;
for j=1:N
    L=Lagranzo_daugianaris(X,j,x);  % Lagranzo daugianariaipagal tolygiai paskirstytas abscises
    LC=Lagranzo_daugianaris(XC,j,x);% Lagranzo daugianariai pagal "Ciobysevo abscises"
    F=F+L*Y(j);                     % kaupiamos sumos interpoliuojanciu funkciju vaizdavimui
    FC=FC+LC*YC(j);
end

plot(x,F,'r-')      % vaizduojama funkcija, interpoliuojanti tolygiai paskirstytuose mazguose 
plot(x,funkcija(x)-F,'r--'),      % vaizduojama netiktis duotos funkcijos atzvilgiu 
delete(hg);hg=legend(leg{1:4});pause

plot(XC,YC,'go')   % vaizduojami interpoliavimo mazgai ties Ciobysevo abscisemis
delete(hg);hg=legend(leg{1:5});pause
plot(x,FC,'g-')   % vaizduojama funkcija, interpoliuojanti Ciobysevo mazguose 
plot(x,funkcija(x)-FC,'g--'),  % vaizduojama netiktis duotos funkcijos atzvilgiu
delete(hg);hg=legend(leg);



return
end

function L=Lagranzo_daugianaris(X,j,x)
% X - interpoliavimo tasku abscises
% j - Lagranzo daugianario numeris (atitinka interp.tasko numeriui)
% x - abscises, kuriose apskaiciuojama daugianario reiksme
    n=length(X);
    L=1;
    for k=1:n, if k ~= j, L=L.*(x-X(k))/(X(j)-X(k)); end, end
    % daugianario reiksmes apskaiciuojamos visuose vaizdavimo taskuose,
    % kuriu abscises yra masyve x
return
end



function fnk=funkcija(x)
% apskaiciuoja interpoliuojamos funkcijos reiksmes taskuose x

% fnk=sin(5*x)+x.^2/10;
% fnk=exp(-10*x.^2);
fnk=1./(1+3*x.^2);
return
end




function G=TschebysevBase(m,x)
     G(:,1)=ones(size(x));
     G(:,2)=x;
     for i=3:m,  G(:,i)=2*x.*G(:,i-1)-G(:,i-2); end
return
end