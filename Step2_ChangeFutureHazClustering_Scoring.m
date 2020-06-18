clear all; close all; clc

%% Code developed by Cristina Izaguirre (izaguirrelasac@gmail.com/izaguirrec@unican.es) and Paula Camus (P.Camus-Brana@soton.ac.uk/camusp@unican.es)
%% within the works of the paper "Climate change risk to global port operations". June 2020
%% IHCantabria-Instituto de Hidraulica Ambiental de la Universidad de Cantabria, Santander, Spain

addpath([pwd '\Clustering'])

path_data=[pwd '\Data\'];
path_colors=[pwd '\Colormaps\'];
path_res=[pwd '\Results\'];

% Load data
portsWPI=load([path_data 'SelectionPortsWPI.dat']); 
CB=load([path_data 'CBposition.dat']);

load([path_data 'ClimaticDrivers.mat']); % climatic drivers dataset
CF=load([path_data 'CoastalFlooding.dat']);
OV=load([path_data 'OvertoppingRM.dat']); % overtopping for Rubble Mound breakwater

colores1=load([path_colors 'colormapChangeHazardGroups.dat']);
colores=flipud(colores1);

cmap=load([path_colors 'cmapChangeHazard.dat']);
rankH=load([path_res 'rankH_tam8.dat']); % output of Step1

colores2=load([path_colors 'colormapScore.dat']);
colores22=flipud(colores2);

% Obtaining normalized changes for every climate impact
% extreme wind (daily resolution)
vWPmax=max(VV(:,1,1)); 
vWFmax=max(VV(:,1,2));
vW=((VV(:,1,2)-VV(:,1,1))./vWPmax).*(VV(:,1,2)./vWFmax);

% high temperature (daily resolution)
vTPmax=max(VV(:,2,1));
vTFmax=max(VV(:,2,2));
vT=((VV(:,2,2)-VV(:,2,1))./vTPmax).*(VV(:,2,2)./vTFmax);

% heavy precipitation (daily resolution)
v20Pmax=max(VV(:,3,1));
v20Fmax=max(VV(:,3,2));
v20=((VV(:,3,2)-VV(:,3,1))./v20Pmax).*(VV(:,3,2)./v20Fmax);

% wave agitation in the Navigation Zone (hourly resolution)
vNZPmax=max(VV(:,5,1))./24; % 
vNZFmax=max(VV(:,5,2))./24;
vNZ=((VV(:,5,2)-VV(:,5,1))./(vNZPmax.*24)).*(VV(:,5,2)./(24*vNZFmax));

% overtopping (hourly resolution)
vOvPmax=max(OV(:,1))./24;
vOvFmax=max(OV(:,2))./24;
vOv=((OV(:,2)./24-OV(:,1)./24)./vOvPmax).*((OV(:,2)./24)./vOvFmax);

% coastal flooding (hourly resolution)
vCFPmax=max(CF(:,1))./24;
vCFFmax=max(CF(:,2))./24;
vCF=((CF(:,2)./24-CF(:,1)./24)./vCFPmax).*((CF(:,2)./24)./vCFFmax);

% TCs (exceedance probability)
vTCPmax=max(VV(:,9,1));
vTCFmax=max(VV(:,9,2));
vTC=((VV(:,9,2)-VV(:,9,1))./vTCPmax).*(VV(:,9,2)./vTCFmax);

% TCs >cat3 (exceedance probability)
vMTCPmax=max(VV(:,10,1));
vMTCFmax=max(VV(:,10,2));
vMTC=((VV(:,10,2)-VV(:,10,1))./vMTCPmax).*(VV(:,10,2)./vMTCFmax);


vP=[vWPmax vTPmax v20Pmax vNZPmax vOvPmax vCFPmax vTCPmax vMTCPmax]; % maximum present values
vF=[vWFmax vTFmax v20Fmax vNZFmax vOvFmax vCFFmax vTCFmax vMTCFmax]; % maximum future values

data=[vW vT v20 vNZ vOv vCF vTC vMTC]; % matrix with the 8 normalized changes for each climate impact

% Data to classify
dumNaN=find(isnan(data(:,1)));
data(dumNaN,:)=[];
[N,dim]=size(data);
novalidos=dumNaN';

names(1)={'DaysW15mph'};
names(2)={'DaysT40Deg'};
names(3)={'DaysP20mm'};
names(4)={'DaysHs2.5m'};
names(5)={'DaysOv0.1l'};
names(6)={'DaysCF'};
names(7)={'ProbTC'};
names(8)={'ProbMajorTC'};

% Scalar and directional parameters
escalar=1:dim;
direccional=[];

% Data normalization
minimos = zeros(dim,1);
maximos = zeros(dim,1);
maxiabs = zeros(dim,1);
for i=1:dim
    minimos(i)=min(data(:,i));
    maximos(i)=max(data(:,i));
    maxiabs(i)=max(abs(minimos(i)),abs(maximos(i)));
end
datos_n = zeros(N,dim);
for i=1:dim
    datos_n(:,i)=data(:,i)./(2*maxiabs(i)); % values ranging from -0.5 to 0.5, taking into account the positive or negative relative change
end

% Clustering using K-means
mascara=ones(dim,1);
% Inicialization using MDA, Camus et al. (2011)
tam=9;
[bmus,centers_n]=kmeans_modificado_initMDA(datos_n,tam,escalar,direccional,mascara);
centers=zeros(tam,3);
for i=1:dim
    centers(:,i)=centers_n(:,i)*2*maxiabs(i);
end

poOv=find(bmus==6 & CB==1); % We consider only changes in overtopping to those ports featured as coastal breakwater
ponoOv=find(bmus==6 & CB==0);

bmus_CB=bmus; % Ports without coastal breakwater are assigned to the next most similar cluster
for i=1:length(ponoOv)
    dis=sqrt((centers_n(:,1)-datos_n(ponoOv(i),1)).^2+(centers_n(:,2)-datos_n(ponoOv(i),2)).^2+(centers_n(:,3)-datos_n(ponoOv(i),3)).^2+(centers_n(:,4)-datos_n(ponoOv(i),4)).^2+(centers_n(:,5)-datos_n(ponoOv(i),5)).^2+(centers_n(:,6)-datos_n(ponoOv(i),6)).^2+(centers_n(:,7)-datos_n(ponoOv(i),7)).^2+(centers_n(:,8)-datos_n(ponoOv(i),8)).^2);
    [mini pomini]=sort(dis);
    bmus_CB(ponoOv(i))=pomini(2);
    
    clear dis mini pomini
end

save([path_res 'bmus_tam' num2str(tam) '_CB.dat'],'bmus_CB','-ascii')

% We consider the maximum value of the centroid but in the case of cluster
% 6 and non breakwater existence, where we consider de following maximum
% value of the centroid
maxicentroid=ones(tam,2).*NaN;
poma=ones(tam,2).*NaN;
for i=1:tam
    [maxicentroid_n(i,1), poma(i,1)]=max(abs(centers_n(i,:)));
    if i==9 && poma(i)==5
        [sorti, posorti]=sort(abs(centers_n(i,:)),'descend');
        maxicentroid_n(i,1:2)=sorti(1:2);
        poma(i,1:2)=posorti(1:2);
    end
    
    maxicentroid(i,1)=centers_n(i,poma(i,1))*2*maxiabs(poma(i,1));
    if isnan(poma(i,2))==0
        maxicentroid(i,2)=centers_n(i,poma(i,2))*2*maxiabs(poma(i,2));
    end
end

% Estimation of changes in the time of stoppage
timeStop=ones(tam,2).*NaN;
for i=1:tam
    
    if i==9
        
        for j=1:2
            
            pp=poma(i,j);
            
            if pp<=3
                dum=find(bmus==i & CB==0);
                VVF=mean(VV(dum,pp,2));
                timeStop(i,j)=maxicentroid(i,j)*vP(pp)*vF(pp)./VVF;
            elseif pp==5
                dum=find(bmus==i & CB==1);
                VVF=mean(OV(dum,2)/24);
                timeStop(i,j)=maxicentroid(i,j)*vP(pp)*vF(pp)./VVF;
            end
            clear dum VVF pp
        end
        
    else
        pp=poma(i,1);
        dum=find(bmus==i);
        if pp<=3
            VVF=mean(VV(dum,pp,2));
            timeStop(i,1)=maxicentroid(i)*vP(pp)*vF(pp)./VVF;
        elseif pp==4
            VVF=mean(VV(dum,5,2)/24);
            timeStop(i,1)=maxicentroid(i)*vP(pp)*vF(pp)./VVF;
        elseif pp==5
            VVF=mean(OV(dum,2)/24);
            timeStop(i,1)=maxicentroid(i)*vP(pp)*vF(pp)./VVF;
        elseif pp==6
            VVF=mean(CF(dum,2)./24);
            timeStop(i,1)=maxicentroid(i)*vP(pp)*vF(pp)./VVF;
        elseif pp==7
            VVF=mean(VV(dum,9,2));
            timeStop(i,1)=maxicentroid(i)*vP(pp)*vF(pp)./VVF;
        elseif pp==8
            VVF=mean(VV(dum,10,2));
            timeStop(i,1)=maxicentroid(i)*vP(pp)*vF(pp)./VVF;
        end
        clear dum pp VVF
    end
end

TS=timeStop;
TS=double(TS);

diasStop=20;
TS(3)=(48+diasStop).*timeStop(3); % Change in the days of stoppage due to TC assuming every port closes 1 day before the hit of the TC and remains closed during its passage (considering 1 day) and 20 days afterwards 
TS(5)=timeStop(5)./2.5; % Representative days/year of stoppage of cluster 5 have been corrected by dividing them by 2.5

% Mean values, 5 and 95 percentiles of the data of each cluster
medios = zeros(tam,dim);
P75 = zeros(tam,dim);
P25 = zeros(tam,dim);
P95 = zeros(tam,dim);
P5 = zeros(tam,dim);
for i=1:tam
        dum=find(bmus_CB==i);
        for d=1:dim
            medios(i,d)=mean(datos_n(dum,d));
            P95(i,d)=prctile(datos_n(dum,d),95);
            P5(i,d)=prctile(datos_n(dum,d),5);
            P75(i,d)=prctile(datos_n(dum,d),75);
            P25(i,d)=prctile(datos_n(dum,d),25);
        end
end

% Number of data in each cluster
acierto=zeros(tam,1);
for t=1:tam    
    dum=find(bmus_CB==t);
    acierto(t)=length(dum);
end
acierto_ = acierto/length(bmus_CB)*100;
acierto_ = round(acierto_*100)/100;
[racierto, pacierto]=sort(acierto_);
posis=flipud(pacierto);


%% Spatial distribution of the clusters. Figure SM3a

figure(1);
ax = worldmap('world');
setm(ax,'FontSize',16)
load coastlines
geoshow(ax, coastlat, coastlon,...
    'DisplayType', 'polygon', 'FaceColor', [.8 .8 .8])
scatterm(portsWPI(:,2),portsWPI(:,1),30,bmus_CB,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001) 
colormap(colores)
colorbar
caxis([1 9])
colorbar('Ticks',[1.5:0.889:9],'TickLabels',[1:9],'Fontsize',16)

file1='FigureSM3a.eps';
print(file1,'-depsc2')

%% Cluster characteristics. Figure SM3b

figure(2)  
for t=1:tam
   subplot(3,3,t)

    for d=1:dim
        
        plot([P25(posis(t),d) P75(posis(t),d)],[d d],'color',colores(posis(t),:),'LineWidth',5)
        hold on
        plot([centers_n(posis(t),d) centers_n(posis(t),d)],[d d],'+k')
        hold on
        plot([P5(posis(t),d) P95(posis(t),d)],[d d],'--o','MarkerSize',2,'color','k')
    end
    axis([-0.5 0.5 0 dim+1])
    if t==1 || t==4 || t==7
        set(gca,'Ytick',1:1:dim,'YtickLabel',names,'FontSize',7.5)
    else
        set(gca,'Ytick',1:1:dim,'YtickLabel',[],'FontSize',7.5)
    end

    grid on
    title([num2str(acierto_(posis(t))) '% ports. Cluster # ' num2str(posis(t))],'FontSize',8.5)
end

file2='FigureSM3b.eps';
print(file2,'-depsc2')

%% Scoring changes

for i=1:length(TS)
    
    score(i,:)=round(TS(i,:)./5.25*10)/10;
    
end

score(score>5)=5; %We established an upper boundary of change in the days of stoppage, considering that more than 23 days/year results in a maximum score

rankchangeH=ones(length(bmus),1).*NaN;

for i=1:tam
    
    if i==9
        po1=find(bmus_CB==i & CB==1);
        rankchangeH(po1)=score(i,1);
        po2=find(bmus_CB==i & CB==0);
        rankchangeH(po2)=score(i,2);
        clear po1 po2
    else
        pob=find(bmus_CB==i);
        rankchangeH(pob)=score(i,1);
        clear pob
    end
    
end


save([path_res 'rankChangeH_tam' num2str(tam) '.dat'],'rankchangeH','-ascii')

%% Multihazard change severity for each port in the year 2100 under RCP8.5 in terms of the improvement and worsening of present multihazard conditions

figure(3);
ax = worldmap('world');
setm(ax,'FontSize',14)
load coastlines
geoshow(ax, coastlat, coastlon,...
    'DisplayType', 'polygon', 'FaceColor', [.8 .8 .8])
scatterm(portsWPI(:,2),portsWPI(:,1),30,rankchangeH,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001) % 
colorbar
colormap(cmap(1:end-1,:))
caxis([-1.5 5.5])
c=colorbar('Ticks',[-1:1:5],'TickLabels',{'Very Low Improvement','No Change','Very Low Worsening','Low Worsening','Medium Worsening','High Worsening','Very High Worsening'});
c.FontSize=11;

file3='FigureSM4.eps';
print(file3,'-depsc2')


%% Estimation of future multi-hazard impact level 

rankFuture=rankH+rankchangeH;

% Scores to plot

poH1=find(rankFuture<=1.5);
poH2=find(rankFuture>1.5 & rankFuture<=2.5);
poH3=find(rankFuture>2.5 & rankFuture<=3.5);
poH4=find(rankFuture>3.5 & rankFuture<=4.5);
poH5=find(rankFuture>4.5 & rankFuture<=5.5);
poH6=find(rankFuture>5.5);

save([path_res 'rankFutureH.dat'],'rankFuture','-ascii')

rankFH=ones(length(rankFuture),1).*NaN;
rankFH(poH1)=1;
rankFH(poH2)=2;
rankFH(poH3)=3;
rankFH(poH4)=4;
rankFH(poH5)=5;
rankFH(poH6)=6;


figure(4);
ax = worldmap('world');
setm(ax,'FontSize',16)
load coastlines
geoshow(ax, coastlat, coastlon,...
    'DisplayType', 'polygon', 'FaceColor', [.8 .8 .8])
scatterm(portsWPI(:,2),portsWPI(:,1),30,rankFH,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001) % 
colormap(colores22)
caxis([1 6])
colorbar('Ticks',[1.4165 2.25 3.083 3.916 4.75 5.583],'TickLabels',{'Very Low','Low','Medium','High','Very High','Extremely High'})


file4='Figure2b.eps';
print(file4,'-depsc2')

