% Lab 6, Q3
% Author: Chng Eng Siong
% 13 Oct 2020

clear all; close all;
A = [1 1 1; 1 1 1; 1 1 -1; 1 -1 -1];
[Q,R] = qr(A,0)   % economy size Q

%sanityCheck
err = A-Q*R


% Q3b:  Q'*Q = = I is orthogonal matrix
% BUT Q*Q' IS NOT I, -> Q is not orthogonal matrix
% bcos it is economy size decomposition
% Q*Q' is the projection matrix with the column space of A
P1 = Q'*Q  % is orthogonal
P2 = Q*Q'  % is the projection matrix


%Q3c) 
% manually find means COUNT the number of pivots
% after rref (reduced row echelon form) of A
r = rank(A)


%Q3d)
r1 = rank([A Q])
% shows that Q and A has the same column space

% Q3e)


%Q3f
y = [1 2 3 4]'
[ProjY_ontoA_Method1, PyW] = my_orthogonalProjection_Method1(Q,y)
ProjY_ontoA_Method2  = (Q*Q')*y    % Method 2 looks trivial 
err1 = y-ProjY_ontoA_Method1 
err2 = y-ProjY_ontoA_Method2
disp('end');



% Q3fe
function [retPyW, PyW] = my_orthogonalProjection_Method1(Q,y)
    [nr,nc] = size(Q);
    %assert(Q'*Q = I) else die
    
    PyW = zeros(nr,nc);
    retPyW = zeros(nr,1);
    for i=1:nc
        PyW(:,i)  = y'*Q(:,i)*(Q(:,i));
        retPyW = retPyW +PyW(:,i);
    end
    
end

