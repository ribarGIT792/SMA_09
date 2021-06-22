%
%   Bendroji Daubechie bangeliu aproksimacija
%
function main
clc;close all;clear all;

n=6        % glodinimo lygiu skaicius
nnn=2^n;    % tasku skaicius smulkiausiame mastelyje
N=3    % Daubechie bangeles pagrindas
wavelet=sprintf('db%1d',floor(N/2+1))

% Is failu ivedami funkcijos taskai: 
fclose all; fhx=fopen('carx.txt','r'); fhy=fopen('cary.txt','r'); % signalo reiksmiu failai
figure(1); axis equal,hold on,grid on
SX=fscanf(fhx,'%g '); SY=fscanf(fhy,'%g '); fclose all; % signalas skaitomas is failo
% SY=SY(end:-1:1);SX=SX(end)-SX(end:-1:1);  % paveikslas "apsukamas"
plot(SX,SY);

% signalas interpoliuojamas i smulkiausia tinkleli:
a=min(SX),b=max(SX),dt=(b-a)/(nnn-1);t=[a:dt:b]; ts=interp1(SX,SY,t);clear SX SY, SX=t;SY=ts;

%--------------- sukuriamas "triuksmas" ------------------
% ind=[floor(length(SY)*0.4)];SY(ind:ind+5)=SY(ind:ind+5)*1.5;   
% ind=[floor(length(SY)*0.3)];SY(ind:ind+5)=SY(ind:ind+5)*1.3;
%---------------------------------------------------------
plot(SX,SY,'r.'); title(sprintf('duota funkcija, tasku skaicius 2^%d',n));
xmin=min(SX);xmax=max(SX); ymin=min(SY);ymax=max(SY); 

% sprendziama mf pletinio koeficientu lygciu sistema --   
    options=optimoptions('fsolve','TolFun',1e-20); c0=ones(1,N+1);
    c=fsolve(@fun,c0,options)
if abs(c(end)) > abs(c(1)), c=c(end:-1:1),end  % parenkamas mf plëtinio sprendinys iğ dviejø galimø
g=c(end:-1:1);  g(2:2:end)=-g(2:2:end);        % bangeliu koeficientai 
% -----------------------------------------------------

m=n   % nagrinejamu detalumo lygiu skaicius   (<=n)

fprintf(1,'\n original:\n');fprintf('%12.5g ', SY*2^(-n/2));fprintf('\n');fprintf('\n');

%************* wavelet decomposition in MATLAB:
fprintf(1,'*****************************************************************************\n\n');
[C,L]=wavedec(SY*2^(-n/2),n,wavelet);
for i=n:-1:1
    fprintf(1,'\nMATLAB(wavedec)    details at level %d,  %d coefficients:\n',i-1,L(i+1)); 
    fprintf(1,'%g12.5   ',C(sum(L(1:i))+1:sum(L(1:i))+L(i+1)));fprintf(1,'\n\n'); 
end
fprintf(1,'\nMATLAB(wavedec)    approx at level %d:',n-n);fprintf(1,'%g12.5   ',C(1:L(1)));fprintf(1,'\n\n');
fprintf(1,'*****************************************************************************\n\n');
%************

%************* wavelet decomposition script:
fprintf(1,'-----------------------------------------------------------------------------\n\n');
smooth=[SY*2^(-n/2),zeros(1,N-1)]; % auksciausio detalumo aproksimacija (pagal duota funkcija)ir nulinis pratesimas
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
        fprintf(1,'\n details %d (level %d),%d+%d coefficients:\n',i,n-i,length(details_start{i}),length(details{i}));
        fprintf('%12.5g ', details_start{i},details{i});
        smooth=smooth1;
        fprintf(1,'\n smooth  %d (level %d):\n',i,n-i);fprintf('%12.5g ',smooth_start{i}, smooth);fprintf('\n');
        fprintf(1,'\n\n');
end 
fprintf(1,'-----------------------------------------------------------------------------\n\n');
fprintf(1,'-----------------------------------------------------------------------------\n\n');

%************* wavelet reconstruction in MATLAB:
fprintf(1,'*****************************************************************************\n\n');
X = waverec(C,L,wavelet);
fprintf(1,'\n reconstructed in MATLAB:\n');fprintf('%12.5g ',X);fprintf('\n');fprintf('\n');
fprintf(1,'*****************************************************************************\n\n');
%*************

% Funkcijos rekonstrukcija:
fprintf(1,'-----------------------------------------------------------------------------\n\n');
fprintf(1,'\n------------------------- rekonstrukcija:\n ');
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
        fprintf(1,'\n smooth  %d (level %d):  ',i,n-i+1);fprintf('%12.5g ', smooth);fprintf('\n');
    end     %------------------------

    % musu pavyzdyje apskaiciuotu aproksimacijos koeficientu vaizdavimas  
    tstep=(b-a)/(length(smooth)); t=[a+tstep/2:tstep:b-tstep/2]
    curve=smooth*2^((n-i+1 )/2) 
    figure(2);subplot(m+1,1,m-i+2),axis equal,axis([xmin xmax ymin ymax]); hold on,grid on, plot(SX,SY,'b.-','Linewidth',1);
    title(sprintf('lygyje %d suglodinta funkcija   n=%d      N=%d',n-i+1,n,N));
    if i~=1 ,plot(t,curve,'r-*','Linewidth',1.5), else, plot(t,curve,'k-','Linewidth',1.5);end
    %legend(leg);   leg={'suglodinta funkcija' 'duota funkcja'};
    
    % *** MATLAB gautu aproksimacijos koeficientu apskaiciavimas ir vaizdavimas
    A = appcoef(C,L,wavelet,i-1) 
    plot(t,A(end-length(smooth)+1:end)*2^((n-i+1 )/2),'b-*','Linewidth',1.5);
 
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