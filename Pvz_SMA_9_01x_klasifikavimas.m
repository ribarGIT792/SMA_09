function main
clc,close all,clear all

xmin=-10;ymin=-10;xmax=10;ymax=10; % koordinaciu sistemos ribos
figure(1),hold on, axis([xmin,xmax,ymin,ymax]);grid on

% Pele ivedamos dvi tasku grupes. Baigiama, kai taskas parenkamas uz koord. sistemos ribu
hhh(1)=text(xmin-3,ymin+1,'Vesti A','BackgroundColor',[0.5 1 0.5]);
hhh(2)=text(xmin-3,ymin+3,'Vesti o','BackgroundColor',[0.5 1 0.5]);
hhh(3)=text(xmin-3,ymin+5,'Skaiciuoti','BackgroundColor',[0.5 1 0.5]);


fprintf(1,'Pazymekite ivedamu tasku grupe: \n');
mmm=[]; % vedamu tasku grupes numeris 
while 1 % pele nustatoma vedamu tasku grupe    
      pause(0.01), if ~isempty(gco),mmm=find(gco == hhh);end
      if ~isempty(mmm), if mmm>=3, break,end, set(hhh(mmm),'BackgroundColor','r'); pause(0.5); break, end
end
fprintf(1,'Ivedamu tasku grupes nr  %d: \n',mmm);
if mmm < 3 %----------------------------------------------------------
    X={};Y={};m1=mmm;
    count(1)=0;count(2)=0;
    while 1 
        count(mmm)=count(mmm)+1;[X{mmm}(count(mmm)),Y{mmm}(count(mmm))]=ginput(1); 
        if X{mmm}(end) < xmin || X{mmm}(end) > xmax || Y{mmm}(end) < ymin || Y{mmm}(end) > ymax, 
            X{mmm}(end)=[];Y{mmm}(end)=[];count(mmm)=count(mmm)-1; set(hhh(mmm),'BackgroundColor',[0.5 1 0.5]);
            fprintf(1,'Pazymekite ivedamu tasku grupe: \n');
            while 1 % pele nustatoma vedamu tasku grupe    
                pause(0.01), if ~isempty(gco),mmm=find(gco == hhh);end
                if ~isempty(mmm), if mmm>=3,fprintf(1,'Skaiciuojama \n',mmm); break,end, 
                    set(hhh(mmm),'BackgroundColor','r');
                    fprintf(1,'Ivedamu tasku grupes nr  %d: \n',mmm);
                    break, 
                end
            end
            
        end
        if mmm ~= m1 & mmm ~= 3, m1=mmm; continue,end
        switch mmm
            case 1, plot(X{mmm}(end),Y{mmm}(end),'k^');
            case 2, plot(X{mmm}(end),Y{mmm}(end),'mo');  
            otherwise, break
        end
    end
end   % ---------------------------------------------------------------

X{1},Y{1}
X{2},Y{2}
delete(hhh)

% pirmosios grupes tasku svorio centras turi buti "auksciau", nei antrosios 
c1=sum(Y{1});c2=sum(Y{2}); if c1 <  c2,'keitimas', yyy=Y{2};Y{2}=Y{1};Y{2}=yyy; xxx=X{2};X{2}=X{1};X{2}=xxx;end

% Sudalijama i intervalus pagal Ox asi ir randami skiriamieji taskai
ninterv= 10;
minX=xmin;maxX=xmax;dx=(maxX-minX)/ninterv;
for i=1:ninterv
     x1=minX+dx*(i-1);x2=minX+i*dx;
     plot([x2,x2],[ymin,ymax],'--c');
     ind1=find((X{1}(:)>=x1 & X{1}(:)< x2));
     ind2=find((X{2}(:)>=x1 & X{2}(:)< x2));
     ycc=takoskyra(Y{1}(ind1),Y{2}(ind2),ymin,ymax);
     if ~isempty(ycc), yc(i)=ycc; else, continue,end
     xc(i)=x1+dx/2;
     plot(xc(i),yc(i),'r*');
%      pause
end

    n=length(xc);   % tasku skaicius
    
hcurve=[];first=1;iend=0;
disp('Press mouse button');

while 1  %--skirtingu aproksimavimo eiliu ciklas---------------------------  
waitforbuttonpress;
m=[];
if first
    number=floor((ymax-ymin)/2);
    for i=1:number, h(i)=text(xmax+1,ymin+2*i-1,sprintf('m=%d',i-1),'BackgroundColor','c'); 
        i=number+1;h(i)=text(xmax,ymin+2*i-1,'Pabaiga','BackgroundColor',[0.5 1 0.5]);
    end % sukuriami aproksimavimo eiles valdikliai 
end
first=0;
fprintf(1,'Pazymekite aproksimavmo eile:\n');

while 1 % pele nustatomas valdiklio(aproksimavimo eiles) numeris    
%         pt=get(gca,'Currentpoint');  % perskaitoma peles padetis
    pause(0.01)
    if ~isempty(gco),m=find(gco == h);end
    if m == number+1; iend=1;break,end
    if ~isempty(m), set(h(m),'BackgroundColor','r');
        if ~isempty(mmm),set(h(mmm),'BackgroundColor','c');end, mmm=m;break,
    end
end
if iend, break,end

fprintf(1,'aproksimavmo eile:  %d\n',m-1);
      
    % Maziausiu kvadratu metodo lygciu sistema:
    G=base(m,xc);
    c=(G'*G)\(G'*yc');
    sss=sprintf('%5.2g',c(1));
    for i=1:m-1,  sss=[sss,sprintf('+%5.2gx^%1d',c(i+1),i)]; end
    sss=strrep(sss,'+-','-');
    
    % Aproksimuojanti funkcija:
    nnn=200; %vaizdavimo tasku skaicius
    xxx=[xmin:(xmax-xmin)/(nnn-1):xmax]; %vaizdavimo taskai 
    Gv=base(m,xxx);
    fff=Gv*c;
    if ~isempty(hcurve), delete(hcurve);end
    hcurve=plot(xxx,fff,'r-');
    
end  %--skirtingu aproksimavimo eiliu ciklo pabaiga------------------------    
close all
disp('Pabaiga')
    return
end


function yc=takoskyra(y1,y2,ymin,ymax)
%
% funkcija randa takoskyros taska tarp vektoriuose pateiktui skaiciu grupiu
% Taskas yra "svorio centro koordinate, imant taskus persikirtimo zonoje 
%
    if 0,     
         if isempty(y1) | isempty(y2),yc=[]; return, end
    else
         if isempty(y1) & isempty(y2),yc=[]; return, end
        if isempty(y1), y1=max(y2); end
        if isempty(y2), y2=min(y1); end 
    end
    
     Y1=min(y1);Y2=max(y2);

 ind1=find(y1 <= Y2 | y1 == Y1)
 ind2=find(y2 >= Y1 | y2 == Y2)
     
     s1=length(y1(ind1));
     s2=length(y2(ind2));
     
     m1=median(y1(ind1));
     m2=median(y2(ind2));
     yc=(m1*s1+m2*s2)/(s1+s2);

return
end

function G=base(m,x)
     for i=1:m,  G(:,i)=x.^(i-1); end
return
end
