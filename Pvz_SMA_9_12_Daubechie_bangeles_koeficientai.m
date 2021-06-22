%
%   Bendroji Daubechie bangeliu aproksimacija
%
function main
clc;close all;clear all;

n=5             % aproksimavimo lygiu skaicius
nnn=2^n        % tasku skaicius smulkiausiame mastelyje
N=3             % Daubechie bangeles pagrindas
m=n-1             % atkuriamu lygiu skaicius
     
% mastelio funkcijos reiksmes sveikaskaitiniuose taskuose:
switch N
    case 1,  deubechie_int = [   1         0]
    case 3,  deubechie_int = [   0.0005    1.3654   -0.3659         0]
    case 5,  deubechie_int = [   0.0005    1.2861   -0.3859    0.0952    0.0042         0]
    case 7,  deubechie_int = [   0.0000    1.0086   -0.0357    0.0403   -0.0118   -0.0012    0.0000         0]
    case 9,  deubechie_int = [   0.0000    0.6978    0.4467   -0.1812    0.0370    0.0015   -0.0017    0.0000    0.0000         0]
    otherwise, '********* nera mastelio funkciju reiksmiu signalo rekonstravimui',return
end

% Is failu ivedami funkcijos taskai: 
fclose all; fhx=fopen('carx.txt','r'); fhy=fopen('cary.txt','r'); % signalo reiksmiu failai
figure(1); axis equal,hold on,grid on
SX=fscanf(fhx,'%g '); SY=fscanf(fhy,'%g '); fclose all; % signalas skaitomas is failo
% SY=SY(end:-1:1);SX=SX(end)-SX(end:-1:1);  % paveikslas "apsukamas"
plot(SX,SY); 

% interpoliuojama i smulkiausia glodinimo tinkleli:
a=min(SX),b=max(SX),t=[a:(b-a)/(nnn-1):b]; ts=interp1(SX,SY,t);clear SX SY, SX=t;SY=ts;
%--------------- sukuriamas "triuksmas" ------------------
% ind=[floor(length(SY)*0.4)];SY(ind:ind+5)=SY(ind:ind+5)*1.5;   
% ind=[floor(length(SY)*0.3)];SY(ind:ind+5)=SY(ind:ind+5)*1.3;
%---------------------------------------------------------
plot(SX,SY,'r.'); title(sprintf('duotas signalas, tasku skaicius 2^%d',n));
xmin=min(SX);xmax=max(SX); ymin=min(SY);ymax=max(SY); 



if 1,  % sprendziama mf pletinio koeficientu lygciu sistema
    
    options=optimoptions('fsolve','TolFun',1e-20);
    c0=ones(1,N+1);
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
g=c(end:-1:1);  g(2:2:end)=-g(2:2:end);        % bangeliu koeficientai 
g

smooth=[SY*2^(-n/2),zeros(1,N-1)]; % auksciausio detalumo glodinimas (pagal duota funkcija)ir nulinis pratesimas
fprintf(1,'\n smooth  ');fprintf('%g ', smooth);fprintf('\n');fprintf('\n');

for i=1:m
        smooth(end+1:end+N-1)=0;  % nulinis pratesimas del nevienetinio bangeles pagrindo
        smooth1=[];
    % pagalbiniai koeficientai srities pradzioje:    
    nstpnt=(N+1)/2-1; smooth_start{i}(1:nstpnt)=0;details_start{i}(1:nstpnt)=0;
    for j=1:nstpnt
        smooth_start{i}(j)=dot(smooth(1:N+1-2*j),c(2*j+1:end))/sqrt(2); details_start{i}(j)=dot(smooth(1:N+1-2*j),g(2*j+1:end))/sqrt(2); 
    end
    % bangeliu ir detaliu koeficientai:   
        for j=1:2^(n-i)
            smooth1(j)=dot(smooth(2*(j-1)+1:2*(j-1)+N+1),c)/sqrt(2);
            details{i}(j)=dot(smooth(2*(j-1)+1:2*(j-1)+N+1),g)/sqrt(2);
        end
        fprintf(1,'\n details %d :  ',i);fprintf('%12.5g ', details{i});
        smooth=smooth1;
        fprintf(1,'\n smooth  %d :  ',i);fprintf('%12.5g ', smooth);fprintf('\n');
end

% Funkcijos rekonstrukcija:
fprintf(1,'\n rekonstrukcija:\n ');

for i=m+1:-1:1 %detalumo didinimo ciklas
    if i < m+1 %---------------------                
        smooth1=[smooth,zeros(1,N-1)]; % nulinis pratesimas i desine del nevienetinio bangeles pagrindo
        smooth=zeros(1,2^(n-i+1)+N);   % kaupiama mf koeficientu suma smulkesnio detalumo lygyje 
    % pagalbiniu koeficientu srities pradzioje panaudojimas:    
    for j=1:nstpnt
        smooth(1:N+1-2*j)=smooth(1:N+1-2*j)+(smooth_start{i}(j)*c(2*j+1:end)+details_start{i}(j)*g(2*j+1:end))/sqrt(2);
    end 
    % bangeliu ir detaliu koeficientu panaudojimas:  
        for j=1:2^(n-i)
            smooth((j-1)*2+1:(j-1)*2+N+1)=smooth((j-1)*2+1:(j-1)*2+N+1)+(smooth1(j)*c+details{i}(j)*g)/sqrt(2);
        end
        smooth=smooth(1:2^(n-i+1)); 
        fprintf(1,'\n smooth  %d :  ',i);fprintf('%12.5g ', smooth);fprintf('\n');
    end     %------------------------

    % kreives braizymas pagal mf koeficientus  
    tstep=(b-a)/(length(smooth));
    t=[a+tstep/2:tstep:b-tstep/2];
    curve=smooth*2^((n-i+1 )/2); 
    figure(2);subplot(m+1,1,m-i+2),axis equal,axis([xmin xmax ymin-0.2 ymax]); hold on,grid on, plot(SX,SY,'Linewidth',1);
    title(sprintf('lygyje %d aproksimuotas signalas   n=%d      N=%d',n-i+1,n,N));
    if i~=1 ,plot(t,curve,'r-*','Linewidth',1.5), else, plot(t,curve,'k-','Linewidth',1.5);end
    
        % pagal MF reiksmes atkurtos kreives braizymas 
    lsmooth=length(smooth);    
    tstep=(b-a)/(lsmooth);
    t=[a:tstep:b];
    curve=zeros(1,lsmooth+1); 
    for i1=1:lsmooth
        N1=N; if i1+N > lsmooth, N1=lsmooth-i1;end
        curve(i1:i1+N1)= curve(i1:i1+N1)+smooth(i1)*deubechie_int(1:N1+1)*2^((n-i+1 )/2);
    end
    % korekcija pagal pagalbinius aproksimavimo koeficientus srities pradzioje
    for i1=1:nstpnt
        start=smooth_start{i1}(1:nstpnt)
        cr=start(i1)*deubechie_int*2^((n-i+1 )/2);
        curve(1:N-nstpnt-1)= curve(1:N-nstpnt-1)+cr(end-(N-nstpnt)+2:end);      
    end
    
    figure(2);subplot(m+1,1,m-i+2),axis equal,axis([xmin xmax ymin-0.2 ymax]); hold on,grid on, plot(SX,SY,'Linewidth',1);
    title(sprintf('lygyje %d aproksimuotas signalas   n=%d      N=%d',n-i+1,n,N));
    if i~=1 ,plot(t,curve,'g-o','Linewidth',1.5), else, plot(t,curve,'m-','Linewidth',1.5);end
end
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
            f((N+1)/2+i+1)=sum(cc (1:2:end))-sum(cc(2:2:end));    
        end

return
end