function [R, t, k] = icp( srcPoints, targetPoints, tau )
% ICP  Compute rigid transformation between two point sets.
%   ICP is an implementation of the Iterative Closest Point algorithm
%   described in Besl, P.J. & McKay, N.D. 1992, 'A Method for Registration 
%   of 3-D Shapes', IEEE Transactions on Pattern Analysis and Machine 
%   Intelligence, vol. 14, no. 2,  IEEE Computer Society.
%
%   [R, t] = ICP( srcPoints, targetPoints, tau ) iteratively seeks the
%   rotation R and translation t which align srcPoints most closely with
%   targetPoints. srcPoints and targetPoints are matrices of 3D points with
%   a column for each point. Note that there is no requirement for 
%   srcPoints to have the same number of points as targetPoints.
%
%   tau is a threshold such that when the RMSE between transformed 
%   srcPoints does not change by more than this threshold between
%   iterations the algorithm is deemed to have converged
%
%   [R, t, k] = ICP( srcPoints, targetPoints, tau ) returns k, the number 
%   of iterations required to complete the match.
%
%
%   For more information see <a href="matlab:
%   web('http://graphics.stanford.edu/courses/cs164-09-spring/Handouts/paper_icp.pdf')">The ICP paper at ResearchGate</a>.


    % Number of args
    if nargin ~= 3
        error( 'icp expects 3 arguments' );
    end
    
    % Get number of points in each and align
    if ~ismatrix( srcPoints) || ~ismatrix( targetPoints)
        error( 'icp only handles 2D matrices' );
    end
    
    % Check they're both 3x something matrices
    if size( srcPoints, 1 ) ~= 3 || size( targetPoints, 1 ) ~= 3
        error( 'icp expects both source and target points to be 3x something matrices' );
    end
    
    % Iterate
    k = 0;
    P_k = srcPoints;
    last_rmse = 0;
    while true
        closestPoints = icp_find_closest_points( P_k, targetPoints );
        [R, t]        = icp_compute_registration ( srcPoints, closestPoints );
        P_k           = icp_apply_registration( srcPoints, R, t );
        rmse          = icp_compute_rmse( closestPoints, P_k );
        
        % Check variation in RMSE
        if abs(rmse - last_rmse) < tau
            break
        end
        last_rmse = rmse;
        k = k + 1;
    end
end    
%%
function closestPoints = icp_find_closest_points( srcPoints, targetPoints )
    % Find the closest point in targetPoints to each point in srcPoints
    % and return them in closestPoints
    %
    % closestPoints = icp_find_closest_points( srcPoints, targetPoints )
    % srcPoints, targetPoints must both be 3x something matrices
    
    % Check for both args present
    if nargin ~= 2 
        error( 'icp_find_closest_point expects two arguments' );
    end
    
    % Assure they're both 2-dimensional
    if ~ismatrix( srcPoints )  || ~ismatrix( targetPoints )
        error( 'icp_find_closest_point expects 2D matrices as input' );
    end
    
    % Assure they're both 3x N
    [srcPointRows,~] = size( srcPoints );
    [targetPointRows, ~] = size( targetPoints);
    
    if srcPointRows ~= 3 || targetPointRows ~= 3
        error( 'icp_find_closest_point expects points to be in columns' );
    end

    % We flip these to use pdist2
    srcPoints = srcPoints';
    targetPoints = targetPoints';
    
    % Compute nearest points
    pointDistances = pdist2( srcPoints, targetPoints, 'euclidean' );
    [~,i] = min( pointDistances, [], 2 );
    % i is now the index of the nearest point in target to each point in
    % src
    closestPoints = targetPoints( i, : );
    closestPoints = closestPoints';
end
%%
function [R, t] = icp_compute_registration ( srcPoints, targetPoints )
% Compute the rigid registration which best maps points P to X
% [R, t] = compute_registration( srcPoints, targetPoints )
%
% sourcePoints is the set of source points (3xN)
% targetPoints is the set of target points (3xN)
% source and target points must be the same size
%
% R is the 3x3 rotation matrix
% t is the 3x1 translation


    % Check for both args present
    if nargin ~= 2 
        error( 'icp_compute_registration expects two arguments' );
    end
    
    % Assure they're both 2-dimensional
    if ~ismatrix( srcPoints )  || ~ismatrix( targetPoints )
        error( 'icp_compute_registration expects 2D matrices as input' );
    end
    
    [srcPointRows,srcPointCols] = size( srcPoints );
    [targetPointRows, targetPointCols] = size( targetPoints);
    
    if srcPointRows ~= 3 || targetPointRows ~=3
        error( 'icp_compute_registration expects points to be in rows' )
    end

    if srcPointCols == targetPointCols
        nPoints = srcPointCols;
    else
        error( 'icp_compute_registration Expects the same number fo source and target points' );
    end
    
    %  Compute centre of mass for both sets of points
    sourceCentreOfMass = sum( srcPoints, 2 ) / nPoints;
    targetCentreOfMass = sum( targetPoints, 2 ) / nPoints;

    % Compute Cross Co-variance matrix
    Sigma_px = zeros(3);
    for i=1:nPoints
        Sigma_px = Sigma_px + srcPoints(:,i) * targetPoints(:,i)';
    end
    Sigma_px = Sigma_px / nPoints;
    Sigma_px = Sigma_px - sourceCentreOfMass * targetCentreOfMass';

    A = Sigma_px - Sigma_px';

    delta = [A(2,3); A(3,1); A(1,2)];

    % Build Q
    Q = [ trace( Sigma_px )  delta';
           delta, Sigma_px + Sigma_px' - (trace( Sigma_px ) * eye(3) )];

    [vec, d] = eigs( Q );
    [~,ci] = max( sum( d ) );

    % Extract unit quaternion from evec corresponding to max eval
    quat = vec(:,ci);

    R = quat2rotm(quat');

    % Compute translation
    t = targetCentreOfMass - R * sourceCentreOfMass;
end
%%
function transformedPoints = icp_apply_registration( srcPoints, rotation, translation ) 
    % Apply the registration in R and t to the points in P
    % returning the transformed points
    
    % Check for all args present
    if nargin ~= 3 
        error( 'icp_apply_registration expects three arguments' );
    end
    
    % dimensionally correct?
    if size( srcPoints, 1 ) ~= 3 
       error( 'icp_apply_registration expects srcPoints is 3xN' );
    else
        nPoints = size( srcPoints, 2 );
    end

    transformedPoints = rotation * srcPoints + repmat( translation, 1, nPoints );
end
%%
function rmse = icp_compute_rmse( points1, points2 )
% Computes the root mean squared error between the sets of 3D points.
%
% rmse = compute_rmse( Pk, Yk )
% Pk is the first set of points and should be 3xN
% Yk is the second set of points and should be 3xN

    % Assert both args are present
    if nargin ~= 2 
        error( 'icp_compute_rmse expects two arguments. only got %d', nargin );
    end

    d = points1 - points2;
    d = d .* d;
    d = sum( d, 1 );
    d = sqrt( d );
    rmse = sum( d, 2 );
end
