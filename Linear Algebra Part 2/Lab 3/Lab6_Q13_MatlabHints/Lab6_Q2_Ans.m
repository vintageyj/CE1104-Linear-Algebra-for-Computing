% Lab Chapter 6:
% Name: Chng Eng Siong
% Date: 12 Oct 2020

% Q2 
A = [1 1 1; 1 1 1; 1 1 -1; 1 -1 -1];
B = [ 3 -1 -1/2; 1 2 -2; 1 1 7/2; 0 1 0];

%Q2a
disp('\n\n\n=============Q2a\n')
A'*A   % since it is not diagonal => columns are not orthgogonal
B'*B   % shows that B is orthogonal matrix
       % and the diahonal of B is (norm of column vector) squared

%Q2b)
disp('\n\n\n=============Q2b\n')
y1 = [2.5 11 -7.5 2]'
% Q2c - brute force, we can use Theorem 5, bcos B HAS orthogonal column
% the weights are dot product(col and y)/(u_i'*u_i)
x = zeros(3,1);
for (i=1:3)
    u_i = B(:,i);
    x(i) = (y1'*u_i)/(u_i'*u_i);
end
y1_tst = B*x
disp('x estimated by Theorem 5:'); 
x
check_approx = y1-(B*x)  % showing that x is found correctly

% actually there is AN easier way to find x_hat
x_hat=pinv(B)*y1  
% the above is my way of getting y1
% This is how I will do it, BUT
% they are not taught this yet, so they have to do it by
% orthgonal decomposition


%Q2c
disp('\n\n\n=============Q2c\n')
a1 = A(:,1); a2 = A(:,2); a3=A(:,3);
y2 = 1*a1+2*a2-3*a3
xA_hat=pinv(A)*y2   % they are not taught this yet,
                    % BUT if they know pinv is a great way to solve this!

% O2c - using theorem 5 WRONGLY to get xA_hat
% BCOS A IS NOT orthogonal
[nr,nc] = size(A);
xA = zeros(nc,1);
for i=1:nc
    a_i = A(:,i);
    xA(i) = y2'*a_i/(a_i'*a_i);
end
xA
check_approx = y2-(A*xA)  % showing that it is NOT equals to y2
% Ans: bcos A is NOT an orthogonal matrix, we cannot use theorem 5

%Q2d
disp('\n\n\n=============Q2d\n')
Bnorm = zeros(nr,nc);
length_u_i = zeros(3,1);
for i=1:nc
    u_i = B(:,i);
    length_u_i(i) = sqrt(u_i'*u_i);
    Bnorm(:,i) = u_i/length_u_i(i);
end
Bnorm'*Bnorm  % showing that it is an orthonormal set
% Q2d - brute force, we can use Theorem 5, bcos Bnorm HAS orthoNORMALcolumn
% the weights are dot product(col and y) and denominator == 1
% Q2d: bcos columns of B are now normalized,
% using theorem 5, the denominator can be ignored.
xN = zeros(3,1);
for (i=1:3)
    u_i = Bnorm(:,i);
    xN(i) = (y1'*u_i);
end
y1_tst = Bnorm*xN
disp('x estimated by Theorem 5:'); 
xN
check_approx = y1-(Bnorm*xN)  % showing that x is found correctly
% sanity Check, using pinv to calculate
xN_hat=pinv(Bnorm)*y1

% Lets find out HOW is xN related to x
% x is calculated by   y1'*(col i of B)/(length(col i of B)^2)
% xN is calculated by  y1'*(col i of Bnormalized)  
%                          // each col of  Bnormalized == 1
%
length_u_i = zeros(3,1);
for i=1:nc
    u_i = B(:,i);
    length_u_i(i) = sqrt(u_i'*u_i);
end
check_xN = xN - (x.*length_u_i)
% xN weight  =  x*length(B col) //  (x = unNormalized B weight) 
% the ABOVE is the relatiobship between x and xN
% the weight found by unNormalized B vs normalizedB 

%Q2e
disp('\n\n\n=============Q2e\n')
% We are projecting ONTO Bn (the normalized vesion of B)
Py_Bn = zeros(nr,nc);
Residual_PyBn = zeros(nr,nc);
err_Residual = zeros(nc,1);
err_Residual2 = zeros(nc,1);
xN = zeros(3,1);
for (i=1:3)
    u_i = Bnorm(:,i);
    xN(i) = (y1'*u_i);
    Py_Bn(:,i) = xN(i)*Bnorm(:,i); 
    Residual_PyBn(:,i) = y1-Py_Bn(:,i);
    err_Residual(i)  = norm(Residual_PyBn(:,i))
    err_Residual2(i) = sqrt(y1'*y1 - Py_Bn(:,i)'*Py_Bn(:,i));
        % err_Residual2 calculated AS scalar (length) of each side ^2
        % Pythagoras theorem!
end

Py_Bn
Residual_PyBn
err_Residual   % we see that there are 2 ways to calculate residual
err_Residual2  % the 2 vectors err_Residual and err_Residual2 are SAME
% the order of importance is the smaller residual error first
% In this example col 3 is the most important col to use!!!


%Q2f
disp('\n\n\n=============Q2f\n')
% exhaustively form 2 columns for tmpB from columns of B
% then find the error! from using tmmB to approx y2
listCombi = [1 2; 1 3; 2 3];  % LIST the combinations of columns
nCombi = length(listCombi);
err_Q2f=zeros(nCombi,1);
tmpB = zeros(nr,2);           % tmpB is holding 2 columns of B
for i=1:nCombi
    idx1 = listCombi(i,1);
    idx2 = listCombi(i,2);
    
    tmpB(:,1) = B(:,idx1);
    tmpB(:,2) = B(:,idx2);
    % Copying the respective 2 columns from B to tmpB

    xhat(1) = y1'*tmpB(:,1)/(tmpB(:,1)'*tmpB(:,1));
    xhat(2) = y1'*tmpB(:,2)/(tmpB(:,2)'*tmpB(:,2));
    xhat = xhat(:);
    residualErr = norm(y1- (tmpB*xhat));
    err_Q2f(i) = residualErr;
    tmpStr = sprintf('i=%d, j=%d, err= %.3e', idx1,idx2, residualErr);
    disp(tmpStr);
end

% We found that using column 2, and 3 PRODUCED the smallest error!!!


