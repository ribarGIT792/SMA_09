function main
clc,close all,clear all

xmin=-10;ymin=-10;xmax=10;ymax=10; % koordinaciu sistemos ribos
figure(1),hold on, axis([xmin,xmax,ymin,ymax]);grid on
m=6; %aproksimuojanciu funkciju skaicius bazeje;

% Pele ivedami taskai. Baigiama, kai taskas parenkamas uz koord. sistemos ribu
X=[];Y=[];
while 1
    [X(end+1,1),Y(end+1,1)]=ginput(1); % ,1 rasome, kad gautume vektorius-stulpelius
    if X(end) < xmin || X(end) > xmax || Y(end) < ymin || Y(end) > ymax, 
        X(end)=[];Y(end)=[]; break;
    end
     plot(X(end),Y(end),'ko');
end
cla, plot(X,Y,'ko');
X,Y
number=floor((ymax-ymin)/2);
for i=1:number 
    h(i)=text(xmax+1,ymin+2*i-1,sprintf('m=%d',i),'BackgroundColor','c');
end % sukuriami aproksimavimo eiles valdikliai 
i=number+1;
h(i)=text(xmax,ymin+2*i-1,'Pabaiga','BackgroundColor',[0.5 1 0.5]);


% Pele parenkama aproksimavimo eile

hcurve=[];
iprint=1;m1=0;
while 1 %--------------- aproksimavimo eiliu ciklas -------------------
    pause(0.1) 
    m=[];pt=[];
    while 1 % pele nustatomas valdiklio(aproksimavimo eiles) numeris    
%         pt=get(gca,'Currentpoint');  % perskaitoma peles padetis
        pause(0.01)
        m=find(gco == h);
        if m~=m1, m1=m;iprint=1;end
        if ~isempty(m), set(h(m),'BackgroundColor','r'); pause(0.5), set(h(m),'BackgroundColor','c'); break,end
    end
    if m == number+1, break, end
    n=length(X);   % tasku skaicius

    % Maziausiu kvadratu metodo lygciu sistema:
    G=MonomialBase(m,X);
%         G=TschebysevBase(m,X);
    if iprint, fprintf(1,"m=%d, \nG'*G=",m);G'*G, iprint=0;end
    c=(G'*G)\(G'*Y);
    sss=sprintf('%5.2g',c(1));
    for i=1:m-1,  sss=[sss,sprintf('+%5.2gx^%1d',c(i+1),i)]; end
    sss=strrep(sss,'+-','-');
    
    % Aproksimuojanti funkcija:
    nnn=200; %vaizdavimo tasku skaicius
    xxx=[xmin:(xmax-xmin)/(nnn-1):xmax]'; %vaizdavimo taskai 
    Gv=MonomialBase(m,xxx);
%     Gv=TschebysevBase(m,xxx);
    fff=Gv*c;
    if ~isempty(hcurve), delete(hcurve);end
    hcurve=plot(xxx,fff,'r-');
    legend({'duoti taskai',sprintf('f(x)=%s',sss)})
    title(sprintf('aproksimavimas maziausiu kvadratu metodu \n  tasku skaicius %d,  funkciju skaicius  %d',n,m));

end %------------------------------------------------------------
close all
return
end


function G=MonomialBase(m,x)
     for i=1:m,  G(:,i)=x.^(i-1); end
return
end


function G=TschebysevBase(m,x)
     G(:,1)=ones(size(x));
     G(:,2)=x;
     for i=3:m,  G(:,i)=2*x.*G(:,i-1)-G(:,i-2); end
return
end