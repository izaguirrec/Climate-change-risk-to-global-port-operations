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
load([path_data 'Exposure_Vulnerability_WPI.mat'])% port characteristics from the WPI database

rankH=load([path_res 'rankH_tam8.dat']);
V=load([path_res 'rankV_tam6.dat']);
E=load([path_res 'ExposureIndicator.dat']);

color=load([path_colors 'colormapRiskGroups.dat']);
colores=flipud(color);

colores1=load([path_colors 'colormapScore.dat']);
colores11=flipud(colores1);
colores11=colores11(1:end-1,:);

% Clustering of risk using the three components: hazard, exposure and
% vulnerability
datos=[rankH E V]; 
dumNaN=find(isnan(datos(:,1)));
datos(dumNaN,:)=[];
[N,dim]=size(datos);
novalidos=dumNaN';

nombres(1)={'Hazard'};
nombres(2)={'Exposure'};
nombres(3)={'Vulnerability'};

escalar=1:dim;
direccional=[];

% Data normalization
minimos = zeros(dim,1);
maximos = zeros(dim,1);
for i=1:dim
    minimos(i)=min(datos(:,i));
    maximos(i)=max(datos(:,i));
end
datos_n = zeros(N,dim);
for i=1:dim
    datos_n(:,i)=(datos(:,i)-minimos(i))./(maximos(i)-minimos(i));
end

% Clustering using K-means
mascara=ones(dim,1);
% Inicialization with MDA, Camus et al. (2011)
tam=10;
[bmus,centers_n]=kmeans_modificado_initMDA(datos_n,tam,escalar,direccional,mascara);
centers=zeros(tam,dim);
for i=1:dim
    centers(:,i)=centers_n(:,i)*(maximos(i)-minimos(i))+minimos(i);
end

for i=1:tam
    podum=find(bmus==i);
    dist=sqrt((centers_n(i,1)-datos_n(podum,1)).^2+(centers_n(i,2)-datos_n(podum,2)).^2+(centers_n(i,3)-datos_n(podum,3)).^2);
    meandist(i,1)=mean(abs(dist));
    meandist(i,2)=std(abs(dist));
    meandist(i,3)=max(abs(dist));
    meandist(i,4)=min(abs(dist));
    
    clear podum dist   
end

save([path_res 'bmusR_tam' num2str(tam) '.dat'],'bmus','-ascii')
save([path_res 'centersR_n_tam' num2str(tam) '.dat'],'centers_n','-ascii')
save([path_res 'minimosR.dat'],'minimos','-ascii')
save([path_res 'maximosR.dat'],'maximos','-ascii')
save([path_res 'esta_distancias_tam' num2str(tam) '.dat'],'meandist','-ascii')

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

%% Spatial distribution of vulnerability clusters. Figure SM8a

figure(1);
ax = worldmap('world');
setm(ax,'FontSize',16)
load coastlines
geoshow(ax, coastlat, coastlon,...
    'DisplayType', 'polygon', 'FaceColor', [.8 .8 .8])
scatterm(portsWPI(:,2),portsWPI(:,1),30,bmus,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001) % 
colormap(colores)
colorbar
caxis([1 10])
colorbar('Ticks',[1.45:0.9:10],'TickLabels',[1:10],'Fontsize',16)

file1='FigureSM8a.eps';
print(file1,'-depsc2')

%% Cluster characteristics. Figure SM8b

figure(2)
for t=1:tam
    subplot(3,4,t)
    
    for d=1:dim
        
        plot([P25(posis(t),d) P75(posis(t),d)],[d d],'color',colores(posis(t),:),'LineWidth',5)
        hold on
        plot([centers_n(posis(t),d) centers_n(posis(t),d)],[d d],'+k')
        hold on
        plot([P5(posis(t),d) P95(posis(t),d)],[d d],'--o','MarkerSize',2,'color','k')
    end
    axis([0 1 0 dim+1])
    if t==1 || t==5 || t==9
        set(gca,'Ytick',1:1:dim,'YtickLabel',nombres,'FontSize',7.5)
    else
        set(gca,'Ytick',1:1:dim,'YtickLabel',[],'FontSize',7.5)
    end
    
    grid on
    title([num2str(acierto_(posis(t))) '% ports. Cluster # ' num2str(posis(t))],'FontSize',8)
end

file2='FigureSM8b.eps';
print(file2,'-depsc2')


%% Scoring risk

sum_cen=sum(centers_n,2);

[SC poSC]=sort(sum_cen);

rangos=linspace(0,3,6);

rankR=ones(length(bmus),1).*NaN;
grupran=[];
for i=1:length(rangos)-1
    poran=find(sum_cen>rangos(i) & sum_cen<=rangos(i+1));
    for j=1:length(poran)
        po=find(bmus==poran(j));
        rankR(po)=i;
        clear po
    end
    pupi=[poran ones(length(poran),1)*i];
    grupran=[grupran;pupi];
    clear poran pupi
end

save([path_res 'rankRisk' num2str(tam) '.dat'],'rankR','-ascii')
save([path_res 'rankgruposRiskP' num2str(tam) '.dat'],'grupran','-ascii')

figure(3);
ax = worldmap('world');
setm(ax,'FontSize',14)
load coastlines
geoshow(ax, coastlat, coastlon,...
    'DisplayType', 'polygon', 'FaceColor', [.8 .8 .8])
scatterm(portsWPI(:,2),portsWPI(:,1),30,rankR,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.1) % 
colormap(colores11)
colorbar
c=colorbar('Ticks',[1.4 2.2 3 3.8 4.6],'TickLabels',{'Very Low','Low','Medium','High','Very High'});
c.FontSize=12;

file3='FigureSM9.eps';
print(file3,'-depsc2')




