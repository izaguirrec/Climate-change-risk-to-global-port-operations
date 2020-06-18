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

colores=load([path_colors 'colormapScore.dat']);
colores1=flipud(colores);

rankFH=load([path_res 'rankFutureH.dat']);
V=load([path_res 'rankV_tam6.dat']);
E=load([path_res 'ExposureIndicator.dat']);

numF=[rankFH E V]; 

tamr=10;
rankRisk=load([path_res 'rankRisk' num2str(tamr) '.dat']);
centers_n=load([path_res 'centersR_n_tam' num2str(tamr) '.dat']);
bmus=load([path_res 'bmusR_tam' num2str(tamr) '.dat']);
distp=load([path_res 'esta_distancias_tam' num2str(tamr) '.dat']);
grupran=load([path_res  'rankgruposRiskP' num2str(tamr) '.dat']);

% Data normalization
minimosR = load([path_res 'minimosR.dat']);
maximosR = load([path_res 'maximosR.dat']);
[N,dim]=size(numF);
datos_n = zeros(N,dim);
for i=1:dim
    datos_n(:,i)=(numF(:,i)-minimosR(i))./(maximosR(i)-minimosR(i));
end

% Check how future triplets of risk change present clusters
rankFRisk=ones(length(E),1).*NaN;
grupo=ones(length(E),1).*NaN;
j=1;
for i=1:length(E)
    dis=sqrt((centers_n(:,1)-datos_n(i,1)).^2+(centers_n(:,2)-datos_n(i,2)).^2+(centers_n(:,3)-datos_n(i,3)).^2);
    
    if bmus(i)==2 || bmus(i)==3 || bmus(i)==4 || bmus(i)==8 % will be able to shift into groups 5, 10, 7 and 1, respectively

        [mini, pomini]=min(dis);
        grupo(i,1)=pomini;
        
    else
        if dis(bmus(i))<=distp(bmus(i),3) % Euclidean distance to the centroid of the original cluster is less than the historical maximum distance so no need to change cluster
            grupo(i,1)=bmus(i);
            
        else % Euclidean distance to the centroid of the original cluster is higher than the historical maximum distance
            sumcen=sum(centers_n(bmus(i),:));
            sumport=sum(datos_n(i,:));
            cambio=sumport-sumcen;
            pasoScore=round(cambio/0.6);% increase in the score  considering 0.6 for the score step in the risk scoring range
            rankFRisk(i)=rankRisk(i)+pasoScore;

            if rankFRisk(i)==6
                
                if datos_n(i,2)==1 && datos_n(i,3)==1
                    rankFRisk(i)=6;
                    
                elseif datos_n(i,1)>=2 && datos_n(i,2)==1
                    rankFRisk(i)=6;
                
                elseif datos_n(i,1)>=2 && datos_n(i,3)==1
                    rankFRisk(i)=6;
                    
                else
                    rankFRisk(i)=5; % only ports with very high exposure and vulnerability or those with very high future hazard and either the exposure or the vulnerability very high will reach this risk category

                end
                
            end
                
            
            clear cambio pasoScore
        end
        
        clear mini pomini
    end
    clear dis
end

for i=1:length(centers_n(:,1))
    porgrupo(i)=length(find(grupo==i))./length(grupo)*100;  
end


%% Future risk scoring

for i=1:tamr
    pogrupo=find(grupo==i);
    pp=find(grupran(:,1)==i);
    rankFRisk(pogrupo)=grupran(pp,2);
    clear pogrupo pp
end

difRisk=rankFRisk-rankRisk;

posm1=find(difRisk==-1);
pos1=find(difRisk==1);
pos2=find(difRisk==2);
pos3=find(difRisk==3);
pos4=find(difRisk==4);

pos0=find(difRisk==0);


figure(1);
ax = worldmap('world');
setm(ax,'FontSize',16)
load coastlines
geoshow(ax, coastlat, coastlon,...
    'DisplayType', 'polygon', 'FaceColor', [.8 .8 .8])
scatterm(portsWPI(pos4,2),portsWPI(pos4,1),250,rankFRisk(pos4),'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001);%,'MarkerEdgeColor',[1 1 1]) %
hold on
scatterm(portsWPI(pos3,2),portsWPI(pos3,1),175,rankFRisk(pos3),'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001);%,'MarkerEdgeColor',[1 1 1]) % 
scatterm(portsWPI(pos2,2),portsWPI(pos2,1),100,rankFRisk(pos2),'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001);%,'MarkerEdgeColor',[1 1 1]) % 
scatterm(portsWPI(pos1,2),portsWPI(pos1,1),50,rankFRisk(pos1),'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001);%,'MarkerEdgeColor',[1 1 1]) % 
scatterm(portsWPI(pos0,2),portsWPI(pos0,1),25,rankFRisk(pos0),'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001);%,'MarkerEdgeColor',[1 1 1]) % 
scatterm(portsWPI(posm1,2),portsWPI(posm1,1),15,rankFRisk(posm1),'filled','MarkerEdgeColor',[0 0 0],'LineWidth',0.001);%,'MarkerEdgeColor',[1 1 1]) % 
colormap(colores1)
caxis([1 6])
colorbar('Ticks',[1.4165 2.25 3.083 3.916 4.75 5.583],'TickLabels',{'Very Low','Low','Medium','High','Very High','Extremely High'})
file1='Figure3.eps';
print(file1,'-depsc2')



