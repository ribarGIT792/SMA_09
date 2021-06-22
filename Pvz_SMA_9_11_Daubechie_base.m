%
%   Debauchy bangeliu baze
%

function main
clc,clear all,close all

N=5     % Daubechie bangeles eile

if 1,  % sprendziama koeficientu lygciu sistema
    
    c0=ones(1,N+1);
    options=optimoptions('fsolve','TolFun',1e-16);
    c=fsolve(@fun,c0,options)
    
    else    % koeficientai duoti
    
    switch N
        case 1,  c=[1,1]   
        case 3,  c=[1+sqrt(3),3+sqrt(3),3-sqrt(3),1-sqrt(3)]/4; 
        case 5,  c=[  0.332670552950,   0.806891509311,  0.459877502118,...
            -0.135011020010,  -0.085441273882,  0.035226291885]*sqrt(2)
        case 7 , c= [0.230377813309,  0.714846570553,  0.630880767929, -0.027983769417, ...
        -0.187034811719,  0.030841381835,  0.032883011667, -0.010597401785]*sqrt(2)
    end
end

if abs(c(end)) > abs(c(1)), c=c(end:-1:1),end  % parenkamas mf plëtinio sprendinys ið dviejø galimø
g=c(end:-1:1); % bangeliu koeficientai 
g(2:2:end)=-g(2:2:end); 

N=numel(c)-1;  % pagrindo plotis
xmin=-0.5;xmax=N+0.5;n=10000;
dx=(xmax-xmin)/(n-1); xxx=[xmin:dx:xmax];

titc=sprintf('N=%d,  c=%s',N,sprintf(' %g  ',c));
titg=sprintf('N=%d,  c=%s',N,sprintf(' %g  ',g));
figure(1),hold on,grid on, set(gcf,'Color','w'), axis([-0.5 N+0.5 -2  2]), title(titc); 
figure(2);set(gcf,'Position',[200 200 560 420] ),hold on,grid on, set(gcf,'Color','w'),axis([-0.5 N+0.5 -2  2]),
figure(3);set(gcf,'Position',[200 10 560 420] ),hold on,grid on, set(gcf,'Color','w'),axis([-0.5 N+0.5 -2  2]) 

xs=[0:0.5:N], fs=zeros(size(xs(1:end-1)));fs(1:2)=[1,1]; % pradinis artinys, Haro mastelio funkcija
figure(1), for i=1:length(xxx),fff(i)=step_function(xxx(i),xs,fs);end, plot(xxx,fff,'b-','Linewidth',2);

% pradinio artinio vaizdavimas
    npoints=length(fs); x=linspace(0,N,npoints);
    gs=fs;
    figure(2), cla,for i=1:length(xxx),fff(i)=step_function(xxx(i),x,fs);end, plot(xxx,fff,'r-','Linewidth',2); title([titc,sprintf('   iteracija  0')]);
    figure(3), cla,for i=1:length(xxx),fff(i)=step_function(xxx(i),x,gs);end, plot(xxx,fff,'b-','Linewidth',2); title([titg,sprintf('   iteracija  0')]);
    pause

for iit=1:10 % ------------- iteracijos mastelio funkcijos ir bangeles apskaiciavimui
    iit
    fscomp=fs;
    fs=zeros(1,2*length(fscomp)); nf=length(fscomp); % dvigubinamas artinio reiksmiu skaicius
    gs=zeros(1,2*length(fscomp)); 
    fs(1:nf)=fs(1:nf)+fscomp*c(1);  % sekantis artinys gaunamas, sumuojant 2k suspaustus ir perstumtus buvusius artinius 
    gs(1:nf)=gs(1:nf)+fscomp*g(1);
    nshift=length(fs)/(2*N); % per kiek poziciju perstumiama 
    for iii=1:N  % pletinio suma
        ishift=iii*nshift;
        fs(ishift+1:ishift+nf)=fs(ishift+1:ishift+nf)+fscomp*c(iii+1);
        gs(ishift+1:ishift+nf)=gs(ishift+1:ishift+nf)+fscomp*g(iii+1);
    end
    
% artinio vaizdavimas
    npoints=length(fs);
    x=linspace(0,N,npoints);
    figure(2), cla,for i=1:length(xxx),fff(i)=step_function(xxx(i),x,fs);end, plot(xxx,fff,'r-','Linewidth',2); title([titc,sprintf('   iteracija  %d',iit)]);
    figure(3), cla,for i=1:length(xxx),fff(i)=step_function(xxx(i),x,gs);end, plot(xxx,fff,'b-','Linewidth',2); title([titg,sprintf('   iteracija  %d',iit)]);
%     figure(2), cla, plot(x,fs,'r-');
%     figure(3), cla, plot(x,gs,'g-');
%     pause  
deubechie_int=interp1(x,fs,[0:N])
pause
end         % -----------------
figure(1); x=linspace(0,N,length(fs)); plot(x,fs,'k-','Linewidth',2);plot(x,gs,'m-','Linewidth',2);
legend({'pradinis artinys (Haro mf)'  'Daubechie mf'  'Daubechie bangele'})
deubechie_int=interp1(x,fs,[0:N])

% scaling functions values at integer points:
%N=3  deubechie_int =    0.0005    1.3654   -0.3659         0
%N=5  deubechie_int =    0.0005    1.2861   -0.3859    0.0952    0.0042         0
%N=7  deubechie_int =    0.0000    1.0086   -0.0357    0.0403   -0.0118   -0.0012    0.0000         0
%N=9  deubechie_int =    0.0000    0.6978    0.4467   -0.1812    0.0370    0.0015   -0.0017    0.0000    0.0000         0

return
end

function rez=add_step_functions(x1,f1,x2,f2)
% apskaiciuoja laiptuotos funkcijos reiksme taske x
% x - funkcijos suoliu abscises 
% f - funkcijos reiksmes intervaluose, length(fs) =  length(xs-1)

if x<xs(1) | x>xs(end), f=0; return,end
v=find((xs-x) > 0);f=fs(v(1)-1);
return
end

function f=step_function(x,xs,fs)
% apskaiciuoja laiptuotos funkcijos reiksme taske x
% xs - funkcijos suoliu abscises 
% fs - funkcijos reiksmes intervaluose, length(fs) =  length(xs-1)

if x<xs(1) | x>xs(end), f=0; return,end
v=find((xs-x) > 0);f=fs(v(1)-1);
return
end

function f=fun(c)
        N=length(c)-1;
        f(1)=sum(c.^2)-2;
        for L=1:(N+1)/2-1
            cc=[c,zeros(1,2*L)];
            f(L+1)=dot(c(1:N+1),cc(2*L+1:2*L+N+1));
        end      
        for i=0:(N+1)/2-1
            coef=[0:N].^i; cc=c(end:-1:1).*coef;
            f((N+1)/2+i+1)=sum(cc(1:2:end))-sum(cc(2:2:end));    
        end
return
end
