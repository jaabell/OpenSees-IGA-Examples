function [] = refinePatch(controlPts, noPtsX, noPtsY, uKnot, vKnot, refineLevel, refU, refV)

controlPts = controlPts;

% controlPts = reshape(controlPts,4,2,2);
% controlPts(1)=transpose(controlPts(1));
% controlPts(2)=transpose(controlPts(2));
uKnot = uKnot;
vKnot = vKnot;


controlPts_def  = zeros(4,noPtsX,noPtsY);
for i=1:noPtsX
    controlPts_def(:,:,i) = transpose(controlPts(:,:,i));
end
% controlPts_def(:,:,1) = transpose(controlPts(:,:,1));
% controlPts_def(:,:,2) = transpose(controlPts(:,:,2));


controlPts = controlPts_def;




%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

%% p-refinment

solid = nrbdegelev(solid,[refU refV]); % to cubic-linear NURBS

%% h-refinement

for i=1:refineLevel
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX newKnotsY};
    solid     = nrbkntins(solid,newKnots);
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end


% figure
% hold on
% nrbkntplot(solid)
% nrbctrlplot(solid)
% %% 


p          = solid.order(1)-1;
q          = solid.order(2)-1;
uKnot      = cell2mat(solid.knots(1));
vKnot      = cell2mat(solid.knots(2));
noPtsX     = length(uKnot)-p-1;
noPtsY     = length(vKnot)-q-1;
weights    = reshape(solid.coefs(4,:,:),noPtsX*noPtsY,1);

controlPts = [];

for iy=1:noPtsY
    controlPts = [controlPts; solid.coefs(1:3,:,iy)'];
end

% our controlPts only stores (x,y,z) not (w*x,w*y,w*z)

controlPts(:,1) = controlPts(:,1)./weights;
controlPts(:,2) = controlPts(:,2)./weights;
controlPts(:,3) = controlPts(:,3)./weights; % for shell problems, z-coord exists
controlPts(:,4) = weights; % for shell problems, z-coord exists


save("shell_refined.mat", "p", "q", "uKnot", "vKnot", "noPtsX", "noPtsY", "weights", "controlPts",'-v7')