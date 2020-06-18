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
load([path_data 'ClimaticDrivers.mat']); % climatic drivers dataset
CF=load([path_data 'CoastalFlooding.dat']);
OV=load([path_data 'OvertoppingRM.dat']); % overtopping for Rubble Mound breakwater

Scores=load([path_data 'ScoringHistImpacts.dat']);

color=load([path_colors 'colormapHazardGroups.dat']);
colores=flipud(color);

colorScore=load([path_colors 'colormapScore.dat']);
colorScore1=flipud(colorScore);
colorScore1=colorScore1(1:end-1,:);


vW=VV(:,1,1); 
vT=VV(:,2,1); 
v20=VV(:,3,1);

vNZ=VV(:,5,1)./24;
vOv=OV(:,1)./24;

vCF=CF(:,1)./24;

vTC=VV(:,9,1);
vMTC=VV(:,10,1);

datos=[vW vT v20 vNZ vOv vCF vTC vMTC]; % Seven impacts due to: extreme wind, high temperature, heavy precipitation, wave agitation, overtopping, coastal flooding, tropical cyclones (cat 1 to 5) and major TC (>cat5).

% Data to classify
dumNaN=find(isnan(datos(:,1)));
datos(dumNaN,:)=[];
[N,dim]=size(datos);
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
for i=1:dim
    minimos(i)=min(datos(:,i));
    maximos(i)=max(datos(:,i));
    if i==5 % normalization of overtopping with the maximum value of the future
       minimos(i)=min(OV(:,1));
       maximos(i)=max(OV(:,2)./24); % maximun value in the future
    end
end

datos_n = zeros(N,dim);
for i=1:dim
    datos_n(:,i)=(datos(:,i)-minimos(i))./(maximos(i)-minimos(i));
end

% Clustering using K-means
mascara=ones(dim,1);
% Inicialization using MDA, Camus et al. (2011)
tam=8;
[bmus,centers_n]=kmeans_modificado_initMDA(datos_n,tam,escalar,direccional,mascara);

centers=zeros(tam,dim);
for i=1:dim
    centers(:,i)=centers_n(:,i)*(maximos(i)-minimos(i))+minimos(i);
end

% Mean values, 5 and 95 percentiles of the data of each cluster
medios = zeros(tam,dim);
P75 = zeros(tam,dim);
P25 = zeros(tam,dim);
P95 = zeros(tam,dim);
P5 = zeros(tam,dim);
for i=1:tam
        dum=find(bmus==i);
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
    dum=find(bmus==t);
    acierto(t)=length(dum);
end
acierto_ = acierto/length(bmus)*100;
acierto_ = round(acierto_*100)/100;
[racierto, pacierto]=sort(acierto_);
posis=flipud(pacierto);

for i=1:dim
        [maxicentroid_n(i,1), poma(i,1)]=max(centers_n(i,:));
        maxicentroid(i,1)=maxicentroid_n(i,1)*(maximos(poma(i))-minimos(poma(i)))+minimos(poma(i));
end

% Calculating days of stoppage
timeStop(1)=maxicentroid(1);
timeStop(3)=maxicentroid(3);
timeStop(6)=maxicentroid(6);
timeStop(7)=maxicentroid(7);


% Days of stoppage due to TC assuming every port closes 1 day before the hit of the TC and remains closed during its passage (considering 1 day) and 20 days afterwards 
diasStop=20;
timeStop(2)=(48+diasStop).*maxicentroid(2);
timeStop(8)=(48+diasStop).*maxicentroid(8);

timeStop(4)=maxicentroid(4)/2.5; % We divided the centroid value by 2.5, considering that the maximum daily temperature occurs during less than half of the day
timeStop(5)=20; % We adopt 20 days/year, which is the maximum closure period accepted by the Spanish Recommendations for Maritime Works (ROM3.1-99) in port design 

save([path_res 'TS.dat'],'timeStop','-ascii')

%% Spatial distribution of the clusters. Figure 1a

figure(1)
ax = worldmap('world');
setm(ax,'FontSize',12)
load coastlines
geoshow(ax, coastlat, coastlon,...
    'DisplayType', 'polygon', 'FaceColor', [.8 .8 .8])
scatterm(portsWPI(:,2),portsWPI(:,1),30,bmus,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001) % 
colormap(colores)
colorbar
caxis([0.5 8.5])
colorbar('Ticks',[1 2 3 4 5 6 7 8],'TickLabels',[1 2 3 4 5 6 7 8],'Fontsize',12)

file1='Figure1a.eps';
print(file1,'-depsc2')

%% Cluster characteristics. Figure 1b

figure(2)  
for t=1:tam
    subplot(2,4,t)
    for d=1:dim
        
        plot([P25(posis(t),d) P75(posis(t),d)],[d d],'color',colores(posis(t),:),'LineWidth',5)
        hold on
        plot([centers_n(posis(t),d) centers_n(posis(t),d)],[d d],'+k')
        hold on
        plot([P5(posis(t),d) P95(posis(t),d)],[d d],'--o','MarkerSize',2,'color','k')
    end
    axis([0 1 0 dim+1])
    if t==1 || t==5
        set(gca,'Ytick',1:1:dim,'YtickLabel',names,'FontSize',8)
    else
        set(gca,'Ytick',1:1:dim,'YtickLabel',[],'FontSize',8)
    end
    grid on
    title([num2str(acierto_(posis(t))) '% ports. Cluster # ' num2str(posis(t))],'FontSize',7.5)
end

file2='Figure1b.eps';
print(file2,'-depsc2')

%% Pie chart

tikilabel={'Wind','Major TC','Coastal Flooding','Temperature','Waves','Precipitation','Multi-hazard','TC'};

explode=[1 1 1 1 1 1 1 1];

figure(3)
p=pie(double(acierto_),explode,tikilabel);
colormap(colores)
set(findobj(p,'type','text'),'fontsize',10)

file3='Figure1c.eps';
print(file3,'-depsc2')

%% Clustering rank

rankH=ones(length(bmus),1).*NaN;

for i=1:8
    
    posi=bmus==i;
    posiScore=Scores(:,1)==i;
    rankH(posi)=Scores(posiScore,2);
    
    clear posi posiScore
end

save([path_res 'rankH_tam' num2str(tam) '.dat'],'rankH','-ascii')

%% Spatial disribution of present multihazard impact level
figure(4);

ax = worldmap('world');
setm(ax,'FontSize',16)
load coastlines
geoshow(ax, coastlat, coastlon,...
    'DisplayType', 'polygon', 'FaceColor', [.8 .8 .8])
scatterm(portsWPI(:,2),portsWPI(:,1),30,rankH,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001) % 
colormap(colorScore1)
colorbar('Ticks',[1.4 2.2 3 3.8 4.6],'TickLabels',{'Very Low','Low','Medium','High','Very High'})

file4='Figure2a.eps';
print(file4,'-depsc2')

