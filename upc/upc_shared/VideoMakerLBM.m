clc
clear all
close all

%% TELL ME FROM WHERE TO READ!
cd ..
location = 'upc_shared/'

%% Options file , BCells file , InitialData file
options       = importdata([location,'SetUpData.ini']);
% BoundaryCells = importdata([location,'/Results/MyBoundaryCells.csv']);
fname         = [location,'Results/InitialData.csv'];

%% Read options from SetUpData.ini file
iterstart = options(11)+options(10);
dt        = options(10);
iterend   = options(9);

%% Resolution
nx = 200*1;
ny = 200*1;

%% Choose what to plot ::: XX
% 1  2  3  4     5     6    7     8    9
% x  y  u  v  vel_mag rho press fluid  TH#
XX = 5;

%% Read data
adat = importdata([fname]);
init = adat.data;
vars = adat.textdata;

%% Find Solid and Fluid cells, separate them
j = 1;
k = 1;
for i = 1:length(init)
    if init(i,8) == 0
        solid(j,:) = init(i,[1 2 8]);
        j = j+1;
    elseif init(i,8) == 1
        fluid(k,:) = init(i,[1 2 8]);
        k = k+1;
    end
end

%% Plot solid and fluid cells
figure(1)
plot3(solid(:,1),solid(:,2),solid(:,3),'k.'), hold on
plot3(fluid(:,1),fluid(:,2),fluid(:,3),'b.'), grid on
view(0,90)
axis equal tight
pause(0.5)

%% Interpolate, create grid
xlin  =  linspace(min(init(:,1)), max(init(:,1)), nx) ;
ylin  =  linspace(min(init(:,2)), max(init(:,2)), ny) ;

[X,Y] =  meshgrid(xlin, ylin);

%% Plot chosen variable
fig100 = figure(11);
TriU = scatteredInterpolant(init(:,1), init(:,2), init(:,XX), 'linear');
Lin = TriU(X,Y) ;
surf(Y, X, Lin, 'EdgeColor', 'none','FaceColor','interp'), grid on
view(90,-90)
    axis equal tight

xlabel('y')
ylabel('x')
title(char(vars(XX)))

pause(0.1)

%% Plot boundary cells and Load distribution
% figure(2)
% for i = 1:length(init)
%  plot3(init(:,1),init(:,2),init(:,9),'bo'), grid on, hold on
 %drawnow
 
% end
% plot3(BoundaryCells.data(:,1),BoundaryCells.data(:,2),BoundaryCells.data(:,9),'r*')


%% Start iteration ::: Create video
XX = 5;
k = 1;
for i = iterstart:dt:iterend
    
    fname = sprintf('%sResults/autosave_iter%05d.csv',location,i);
    %fname = sprintf('%sResults/FinalData.csv',location);
    adat = importdata([fname]);
    init = adat.data;
    fig100 = figure(100);
    hold off
    TriU = scatteredInterpolant(init(:,1), init(:,2),  init(:,XX), 'linear');
    Lin = TriU(X,Y) ;
    surf(Y, X, Lin, 'EdgeColor', 'none')%,'FaceColor','interp'), grid on
    hold on, plot3(solid(:,2),solid(:,1),solid(:,3),'ko', 'MarkerFaceColor',[0 0 0])
    grid on
    view(90,-90)
    axis equal tight
    %pbaspect([2 2 1])
    xlabel('y')
    ylabel('x')
    title(['Velocity magnitude, step: ', num2str(i)])
    drawnow%, pause
%     saveas(fig100,sprintf('pics/pic%05d',k),'png'), k = k+1;
    %pause(0.01)
% plot3(BoundaryCells.data(:,2),BoundaryCells.data(:,1),BoundaryCells.data(:,5),'w*'), pause
end
% figure(100)
% hold on, plot3(BoundaryCells.data(:,2),BoundaryCells.data(:,1),BoundaryCells.data(:,5),'w*')