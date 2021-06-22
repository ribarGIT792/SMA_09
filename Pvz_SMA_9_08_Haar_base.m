%
%   Haro bangeliu baze
%

function main
clc,close all,clear all

a=0;b=1.;
xmin=a-0.5;xmax=b+0.5;n=1000;
dx=(xmax-xmin)/(n-1);
xxx=[xmin:dx:xmax];
figure(1),hold on,grid on
subplot(4,1,1);hold on,grid on;plot(xxx,Haar_wavelet(xxx,0,0,a,b),'b-','LineWidth',4);title('Haro bangeles');
subplot(4,1,2); hold on,grid on;plot(xxx,Haar_wavelet(xxx,1,0,a,b),'b-','LineWidth',4);plot(xxx,Haar_wavelet(xxx,1,1,a,b),'r-','LineWidth',3);
subplot(4,1,3); hold on,grid on;plot(xxx,Haar_wavelet(xxx,2,0,a,b),'b-','LineWidth',4);plot(xxx,Haar_wavelet(xxx,2,1,a,b),'r-','LineWidth',3);
                        plot(xxx,Haar_wavelet(xxx,2,2,a,b),'g-','LineWidth',2);plot(xxx,Haar_wavelet(xxx,2,3,a,b),'m-','LineWidth',1);  
subplot(4,1,4); hold on,grid on;plot(xxx,Haar_wavelet(xxx,3,0,a,b),'b-','LineWidth',4);plot(xxx,Haar_wavelet(xxx,3,1,a,b),'r-','LineWidth',3);
                        plot(xxx,Haar_wavelet(xxx,3,2,a,b),'g-','LineWidth',2);plot(xxx,Haar_wavelet(xxx,3,3,a,b),'m-','LineWidth',1);
                        plot(xxx,Haar_wavelet(xxx,3,4,a,b),'g-','LineWidth',2);plot(xxx,Haar_wavelet(xxx,3,5,a,b),'m-','LineWidth',1);
                        plot(xxx,Haar_wavelet(xxx,3,6,a,b),'g-','LineWidth',2);plot(xxx,Haar_wavelet(xxx,3,7,a,b),'m-','LineWidth',1);
                        

figure(2),hold on,grid on,
Haar_scaling(xxx,0,0,a,b)
subplot(4,1,1);hold on,grid on;plot(xxx,Haar_scaling(xxx,0,0,a,b),'b-','LineWidth',4);


title('Haro mastelio funkcijos')
subplot(4,1,2); hold on,grid on;plot(xxx,Haar_scaling(xxx,1,0,a,b),'b-','LineWidth',4);plot(xxx,Haar_scaling(xxx,1,1,a,b),'r-','LineWidth',3);
subplot(4,1,3); hold on,grid on;plot(xxx,Haar_scaling(xxx,2,0,a,b),'b-','LineWidth',4);plot(xxx,Haar_scaling(xxx,2,1,a,b),'r-','LineWidth',3);
                        plot(xxx,Haar_scaling(xxx,2,2,a,b),'g-','LineWidth',2);plot(xxx,Haar_scaling(xxx,2,3,a,b),'m-','LineWidth',1);
subplot(4,1,4); hold on,grid on;plot(xxx,Haar_scaling(xxx,3,0,a,b),'b-','LineWidth',4);plot(xxx,Haar_scaling(xxx,3,1,a,b),'r-','LineWidth',3);
                        plot(xxx,Haar_scaling(xxx,3,2,a,b),'g-','LineWidth',2);plot(xxx,Haar_scaling(xxx,3,3,a,b),'m-','LineWidth',1);   
                        plot(xxx,Haar_scaling(xxx,3,4,a,b),'g-','LineWidth',2);plot(xxx,Haar_scaling(xxx,3,5,a,b),'m-','LineWidth',1);
                        plot(xxx,Haar_scaling(xxx,3,6,a,b),'g-','LineWidth',2);plot(xxx,Haar_scaling(xxx,3,7,a,b),'m-','LineWidth',1);                        
                        
    title('Haro mastelio funkcijos')
return





return
end

function h=Haar_scaling(x,j,k,a,b)
xtld=1/(b-a)*x+a/(a-b); % (a,b) intervale duota kintamojo reiksme perskaiciuojama i "standartini" 
                        % intervala (0,1), kuriame uzrasyta bangeles formule  
xx=2^j*xtld-k;
h=2^(j/2)*(sign(xx)-sign(xx-1))/(2*(b-a));
end
function h=Haar_wavelet(x,j,k,a,b)
xtld=1/(b-a)*x+a/(a-b); % (a,b) intervale duota kintamojo reiksme perskaiciuojama i "standartini" 
                        % intervala (0,1), kuriame uzrasyta bangeles formule  
xx=2^j*xtld-k;
h=2^(j/2)*(sign(xx)-2*sign(xx-0.5)+sign(xx-1))/(2*(b-a));
return
end


