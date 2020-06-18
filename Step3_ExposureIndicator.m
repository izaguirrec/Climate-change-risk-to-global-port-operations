clear all; close all; clc

%% Code developed by Cristina Izaguirre (izaguirrelasac@gmail.com/izaguirrec@unican.es) and Paula Camus (P.Camus-Brana@soton.ac.uk/camusp@unican.es)
%% within the works of the paper "Climate change risk to global port operations". June 2020
%% IHCantabria-Instituto de Hidraulica Ambiental de la Universidad de Cantabria, Santander, Spain

path_data=[pwd '\Data\'];
path_colors=[pwd '\Colormaps\'];
path_res=[pwd '\Results\'];

% Load data
portsWPI=load([path_data 'SelectionPortsWPI.dat']); 
load([path_data 'Exposure_Vulnerability_WPI.mat'])% port characteristics from the WPI database

colores=load([path_colors 'colormapScore.dat']);
colores1=flipud(colores);
colores1=colores1(2:end-1,:);


%% Harbour Type
cod_HT={'CN','CB','CT','TH','OR'};
punt_HT=[2 3 1 3 4]; % Score from table SM6

HT=PORTS_WPI.HarborType;
HT_punt=ones(length(HT),1).*NaN;
for i=1:5
    
    poT=find(strcmp(cod_HT{i},HT));
    HT_punt(poT)=punt_HT(i);
    
    clear poT
end

save([path_res 'ExposureIndicator.dat'],'HT_punt','-ascii')

figure(1);
ax = worldmap('world');
setm(ax,'FontSize',14)
load coastlines
geoshow(ax, coastlat, coastlon,...
    'DisplayType', 'polygon', 'FaceColor', [.8 .8 .8])
scatterm(portsWPI(:,2),portsWPI(:,1),30,HT_punt,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.1) % 
colormap(colores1)
colorbar
c=colorbar('Ticks',[1.375 2.125 2.875 3.625],'TickLabels',{'Low','Medium','High','Very High'});
c.FontSize=12;

file='FigureSM5.png';
print(file,'-depsc2')




