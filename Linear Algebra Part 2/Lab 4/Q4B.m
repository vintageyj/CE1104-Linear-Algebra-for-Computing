% Q4B Matrix Algebra
%
A = [1 -2 3 ; 5 -6 7 ; 9 10 -11]
B = [-1 -2 -3 ; 5 -6 7 ; 9 -10 11]
C = [-1 2 -3 ; -5 -6 7 ; -9 10 -11]
% given above A,B,C find X

disp('Q4B ==============');
m = 3;
X = inv(inv(A)*C*inv(B)-eye(m))*(-eye(m))
tst_Left  = inv(A*X-A)
tst_Right = inv(X)*B*inv(C)






