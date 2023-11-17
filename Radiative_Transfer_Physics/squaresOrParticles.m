%% ---- Random Squares in a Larger Square -----

% How much light gets through if I scatter random squares within a square?
% What numbers cause the light to change significantly, allowing a random
% position? What numbers cause the light to change negligably despite
% random position?

% By Andrew J. Buggee

function [pgon] =  squaresOrParticles(edgeLength,numParticles,plotFlag)



% ---- Lets make a single square with no overlapping small squares ---

% for there to be overlap between squares, the resolution of our grid must
% be different from the resolution of our squares. Lets set the resolution
% of our grid to be an order of magnitude less the the length of our
% squares

% edgeLength = 2;
dx = 0.0001;
x = -(edgeLength/2):dx:(edgeLength/2);
y = x;
[XX,YY] = meshgrid(x,y);


particleSize = dx/10;


% numParticles = 100;

% for all particles to be within the boundary of the box, the CG must be
% atleast 1/2 * particleSize away from the edges. To do this we can shrink
% to boundaries by particleSize/2.

indexCG_X = x>(min(x)+particleSize/2) & x<(max(x)-particleSize/2);
indexCG_Y = y>(min(y)+particleSize/2) & y<(max(y)-particleSize/2);

CGx = randsample(x(indexCG_X),numParticles,true);
CGy = randsample(y(indexCG_Y),numParticles,true);

%% --- Can they be at the same location? Can they Overlap? ---

%% ---- Plot our Boxed Particles ----

% polyshape() creates a 2D polygon using [x,y] vertices defined by the
% user. polyshap objects are powerful, because we can preform many
% opertaions on them like union(), subtract(), overlap() etc.

newx = CGx-particleSize/2;
newy = CGy - particleSize/2;
cornersX = [newx',newx'+particleSize,newx'+particleSize,newx'];
cornersY = [newy',newy',newy'+particleSize,newy'+particleSize];

pgon = polyshape();
for ii = 1:numParticles
    
    pgon(ii) = polyshape(cornersX(ii,:),cornersY(ii,:));
    
end

pgonUnion = union(pgon);

%% ---- PLOTS ----

if plotFlag == true
    
    % plots the outer square
    w = max(x)-min(x);
    h = max(y)-min(y);
    figure; rectangle('position',[min(x),min(y),w,h],'EdgeColor','w','LineWidth',2)
    axis([min(x)-(0.1*w),max(x)+0.1*w,min(y)-(0.1*h),max(y)+0.1*h])
    hold on
    
    % plot the polygons
    plot(pgon,'EdgeColor','w','FaceColor','w','FaceAlpha',0.6)
    title('All Squares, including Overlap','color','white')
    
    % now plot only the union of each - meaning we dont double count
    % overlapping squares
    
    % plots the outer square
    figure; rectangle('position',[min(x),min(y),w,h],'EdgeColor','w','LineWidth',2)
    axis([min(x)-(0.1*w),max(x)+0.1*w,min(y)-(0.1*h),max(y)+0.1*h])
    hold on
    
    % plot the polygons
    plot(pgonUnion,'EdgeColor','w','FaceColor','w','FaceAlpha',0.6)
    title('Total Area - Overlaping Area','color','white')
    
    
else
    
end


end


