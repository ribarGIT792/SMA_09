%
%   Bendroji Daubechie bangeliu aproksimacija MATLAB Wavelet TB
%
function main
clc;close all;clear all;

n=5        % aproksimavimo lygiu skaicius
nnn=2^n;    % tasku skaicius smulkiausiame mastelyje
N=3    % Daubechie bangeles pagrindas
wavelet=sprintf('db%1d',floor(N/2+1))
m=n   % nagrinejamu detalumo lygiu skaicius   (<=n)

% Is failu ivedami funkcijos taskai: 
fclose all; fhx=fopen('carx.txt','r'); fhy=fopen('cary.txt','r'); % signalo reiksmiu failai
figure(1); axis equal,hold on,grid on
SX=fscanf(fhx,'%g '); SY=fscanf(fhy,'%g '); fclose all; % signalas skaitomas is failo
plot(SX,SY); 
% interpoliuojama i smulkiausia tinkleli:
a=min(SX),b=max(SX),dt=(b-a)/(nnn-1);t=[a:dt:b]; ts=interp1(SX,SY,t);clear SX SY, SX=t;SY=ts;
% %--------------- sukuriamas "triuksmas" ------------------
% ind=[floor(length(SY)*0.4)];SY(ind:ind+5)=SY(ind:ind+5)*1.5;   
% ind=[floor(length(SY)*0.3)];SY(ind:ind+5)=SY(ind:ind+5)*1.3;
%---------------------------------------------------------
plot(SX,SY,'r.'); title(sprintf('duota funkcija, tasku skaicius 2^%d',n));
xmin=min(SX);xmax=max(SX); ymin=min(SY);ymax=max(SY); 
fprintf(1,'\n original:\n');fprintf('%12.5g ', SY*2^(-n/2));fprintf('\n');fprintf('\n');

%************* wavelet decomposition in MATLAB:
fprintf(1,'*****************************************************************************\n\n');
[C,L]=wavedec(SY*2^(-n/2),n,wavelet);
L
for i=n:-1:1
    fprintf(1,'\nMATLAB(wavedec)    details at level %d,  %d coefficients:\n',i-1,L(i+1)); 
    fprintf(1,'%g12.5   ',C(sum(L(1:i))+1:sum(L(1:i))+L(i+1)));fprintf(1,'\n\n'); 
end
fprintf(1,'\nMATLAB(wavedec)    approximation coefficients at level %d:',n-n);fprintf(1,'%g12.5   ',C(1:L(1)));fprintf(1,'\n\n');
fprintf(1,'*****************************************************************************\n\n');
%************

%************* wavelet reconstruction in MATLAB:
fprintf(1,'*****************************************************************************\n\n');
X = waverec(C,L,wavelet);
fprintf(1,'\n reconstructed in MATLAB:\n');fprintf('%12.5g ',X);fprintf('\n');fprintf('\n');
fprintf(1,'*****************************************************************************\n\n');
%*************

% Funkcijos rekonstrukcija kiekviename lygyje:
fprintf(1,'-----------------------------------------------------------------------------\n\n');
fprintf(1,'\n------------------------- rekonstrukcija:\n ');
for i=m+1:-1:1 %detalumo didinimo ciklas
    figure(2);subplot(m+1,1,m-i+2),axis equal,axis([xmin xmax ymin ymax]); hold on,grid on, plot(SX,SY,'Linewidth',1);
    title(sprintf('lygyje %d  aproksimuota funkcija   n=%d      N=%d',n-i+1,n,N));

    % ****daline rekonstrukcija:
C1=C; C1([sum(L(1:n-i+2))+1:end])=0; % atsisakoma smulkesniu detaliu
tstep=(b-a)/(2^n-1);t1=[a:tstep:b];
X1 = waverec(C1,L,wavelet);
plot(t1,X1*2^(n/2),'go','Linewidth',1.5);    
   
end

end
