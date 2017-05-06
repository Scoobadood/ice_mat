% ICP Testing

% Start tidy
clc;
clear;

% Make some points as a model
Nx = 100;
X = abs( rand( 3, Nx ) ) * 100;

% Make a Rotation and translation
t = rand( 3, 1 ) * 25;

% Select random rotation angles
theta = rand() * 20;
phi   = rand() * 20;
psi   = rand() * 20;

R = eul2rotm( deg2rad( [ theta, phi, psi] ), 'ZYX' );
  
% Select a subset of the points
Np = 40;
P = X( :, uint32( ceil( rand( 1, Np ) * Nx ) ) );

% Transform points
P = R' *  (P - repmat( t, 1, Np ));

%% 
% Recover transformation
[Rr,tr, k] = icp( P, X, 0.5 );

fprintf( 'Iterations : %d\n', k );

% Compare to the originals
eul = rotm2eul( Rr, 'ZYX' );
eul = rad2deg( eul );

fprintf( 'Theta : %3.f, estimate : %3.f\n', theta, eul(1) );
fprintf( '  Phi : %3.f, estimate : %3.f\n', phi,   eul(2) );
fprintf( '  Psi : %3.f, estimate : %3.f\n', psi,   eul(3) );

fprintf( '   tx : %3.f, estimate : %3.f\n', t(1),   tr(1) );
fprintf( '   ty : %3.f, estimate : %3.f\n', t(2),   tr(2) );
fprintf( '   tz : %3.f, estimate : %3.f\n', t(3),   tr(3) );
    
%% Plot the originals and recovered points
% Model
scatter3( X(1,:), X(2,:) , X(3,:) , 'go' );
axis vis3d
hold

% Original P
scatter3( P(1,:), P(2,:) , P(3,:) , 'b+');

% P fitted to X
np = Rr * P + repmat( tr, 1, Np );    
scatter3( np(1,:), np(2,:) , np(3,:) , 'rx');