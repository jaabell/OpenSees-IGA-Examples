%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for Kirchoff-Love shell problems.
%
% Rotation-free thin shells. Fully clamped or simply supported 
% rectangular plates.
%
% Vinh Phu Nguyen,
% Cardiff University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [] = makeCantileverShell(a,b, refineLevel, refU, refV)

addpath ../migfem/fem_util/
addpath ../migfem/C_files/
addpath ../migfem/data/
addpath ../migfem/meshing/
addpath ../migfem/post-processing/
addpath ../migfem/fem-functions/
addpath ../migfem/meshing/
addpath ../migfem/nurbs-util/
addpath ../migfem/nurbs-geopdes/inst/

%% Geometry data


% plate dimensions
% a = 54.0;
% b = 30.0;


% knots
uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

% control points
controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [a;0;];

controlPts(1:2,1,2) = [0;b];
controlPts(1:2,2,2) = [a;b];

% weights
controlPts(4,:,:)   = 1;


% controlPts


%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});


%% p-refinement

solid = nrbdegelev(solid,[refU refV]); % to quadratic-quadratic NURBS


%% h-refinement

% refineLevel = 4;
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



save("shell.mat", "p", "q", "uKnot", "vKnot", "noPtsX", "noPtsY", "weights", "controlPts",'-v7')