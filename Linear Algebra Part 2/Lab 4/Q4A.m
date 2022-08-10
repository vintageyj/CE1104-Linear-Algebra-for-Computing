% Lab Chapter 7

% Q4A)
% 
% this question gives a A matrix4x2, b = 4x1, and solve for x
% and hence estimate hat{b} = A*x_est
% solve for x_est using two methods:
%    a) pseudo inverse
%    b) normal equation
%    c) Comment if x their results are the same?

% d) find \hat{b} using pseudo inverse vs normal_equation solution, 
%    find hat_b1a = A*x1a (pseudo inverse),  hat_b1b = A*x1b (normal equation)

% e) what is the norm of error of approximating b? which solution is better
%    if x1a and x1b are different?


% f) question: is b in the column space of A? % How many ways to show this ?
% g) Is the error vector  (produced by normal eq) is orthogonal to column space of A
% h) Can we use the normal equation to solve for x?

% i) Construct Q with orthogonal basis spanning the column space of A
%    i1) how can you show that columns of Q is in the column space of A?
%      - answer use rank([A q_j]) == rank([A]) for j=1..n (number of col in Q)
%    i2) solve Qx = b, what is the error of approximating b? is it the same
%    as Q1e? Why?
%
%  j) Evaluate N point (e.g, 10) around the least squares solution (normal equation)
%     when it exist
%     to show that it gives the least error squares in approximating b.
%     Make delta_x has norm(1),


% lets create different A,b to test OUT
%[A,b] = makeAb_Ex1([1 2 3 4])

% This is full rank example
[A,b] = makeAb_Ex2([1 2 3 4])
check_Ab(A,b)
noisy_b  = b+rand(length(b),1);
check_Ab(A,noisy_b)


% This is a rank 1 matrix A, ONLY 1 col is independent
% and b is in the space of A
function [A,b] = makeAb_Ex1(ipCol)
ipCol=ipCol(:);
A = [ipCol 2*ipCol 3*ipCol];
b = A(:,1)-2*A(:,2)+3*A(:,3);
end

% This is a full col rank  matrix A, 
% and b is in the space of A
function [A,b] = makeAb_Ex2(ipCol)
ipCol=ipCol(:);
A = [ipCol ipCol ipCol];
A(2,2) = 0; A(3,3) = 0;
b = A(:,1)-2*A(:,2)+3*A(:,3);
end



function Res_Str = check_Ab(A,b)
   
% solve for x_est using two methods:
%    a) pseudo inverse
%    b) normal equation
%    find hat_b1a = A*x1a (pseudo inverse),  hat_b1b = A*x1b (normal equation)
%    c) Comment if their results are the same?

opStr=sprintf('\n\n\n\n==========================\n\n');
disp(opStr);

[nr nc] = size(A);
disp('1a,b)Ques: Solve for x_est using two methods:');
disp('1a)Ans pseudo inverse:');  x1a = pinv(A)*b
disp('1b)Ans normal equation:'); x1b = inv(A'*A)*A'*b
flag_x1a_x1b_SAME = 0;
disp('1c)Ques: are the solution of pseudo inverse and normal equation same?');
if (norm(x1a-x1b)<1e-6)
    disp('1c)Ans Results of a,b are the same');
    flag_x1a_x1b_SAME = 1;
else
    disp('1c)Ans: Results of a,b are the DIFFERENT');
end

disp('d)Ques: find \hat{b} using pseudo inverse vs normal_equation solution');
disp('d)Ans \hat{b} using pseudo inverse:');
hat_b_x1a = A*x1a
disp('d)Ans \hat{b} using normal eqn:');
hat_b_x1b = A*x1b

% e) what is the norm of error of approximating b? which solution is better
%    if x1a and x1b are different?
disp('1e)Ques: what is the norm of error of approximating b? which solution is better if x1a and x1b are different');
disp('1e)Ans (pseudoInverse error) = ');
err_b_x1a = b-hat_b_x1a;  norm_err_b_x1a = norm(err_b_x1a)
disp('1e)Ans (normal equation) = ');
err_b_x1b = b-hat_b_x1b;  norm_err_b_x1b = norm(err_b_x1b)
if (abs(norm_err_b_x1a- norm_err_b_x1b)<1e-6)
    disp('1e)Ans both solutions produced the same error to approximate b');
else
    opStr = sprintf('Error using pseudoInv = %.3e  vs normalEq =%.3e',norm_err_b_x1a, norm_err_b_x1b);
    disp(opStr)
    if (norm_err_b_x1a < norm_err_b_x1b)
        disp('e) USING pseudo inverse is better');
    else
        disp('e) USING normal equation is better');
    end
end



% f) question: is b in the column space of A? How many ways to show this ?
disp('1f)Quest, is b in the column space of A? How many ways to show this ?');

disp('1f)Ans method 1: check rank([A b]), if rank is unchanged than col b to col space A');
flag_b_inColSpaceA_method1 = 0;
if (rank([A b],1e-6) == rank(A,1e-6))
    % this shows that number of pivot == number of col of A
    % this implies that the col b is in the col space of A
    flag_b_inColSpaceA_method1 = 1;
end

disp('1f)Ans method 2: check if norm(Ax-b) is close to zero');
flag_b_inColSpaceA_method2 = 0;
if (norm_err_b_x1a < 1e-6)  % some very small value
    flag_b_inColSpaceA_method2 = 1;
end

opStr = sprintf('1f)Ans 2 methods to show that col b is IN COL SPACE of A : %d,%d,%d',  flag_b_inColSpaceA_method1,flag_b_inColSpaceA_method2);
disp(opStr)
    


% g) Is the error vector (produced by normal eq) is orthogonal to column space of A
disp('1g)Quest Is the error vector  (produced by normal eq)  orthogonal to column space of A:');
if (norm_err_b_x1b < 1e-6)
    opStr = sprintf('1g)Ans: error vector is orthogonal to col space of A BCOS erro vector is ZERO');  
    disp(opStr)
else
   flag_dotProductErrVecToColA = 0;
   for i=1:nc
      col_Ai = A(:,i);
      if norm(col_Ai'*err_b_x1b) > 1e-6
           flag_dotProductErrVecToColA = flag_dotProductErrVecToColA +1
      end % of if
   end  %of for i

   if (flag_dotProductErrVecToColA == 0)    
      opStr = sprintf('1g)Ans: error vector is orthogonal to col space of A : CHECK dotProduct colA to errVector');  
   else
      opStr = sprintf('1g)Ans: error vector is NOT-ORTHOGONAL to col space of A: CHECK dotProduct colA to errVector'); 
   end %of if
   disp(opStr)
end % of  if norm_err_b_x1b


disp('1h)Quest Can we use the normal equation to solve for x?');
if (nc == rank(A'*A,1e-6))
    opStr = sprintf('1h)Ans: Yes because rank AtA is rank %d for %d variables - inverse AtA exits', rank(A'*A,1e-6), nc);
else
    opStr = sprintf('1h)Ans: NO because AtA is rank %d for %d variables - AtA SINGULAR ', rank(A'*A, 1e-6), nc);
end
    disp(opStr);

% i) Construct Q with orthogonal basis spanning the column space of A
%    i1) Show that the columns of Q is in the column space of A?
%      - answer use rank([A q_j]) == rank([A]) for j=1..n (number of col in Q)
%    i2) solve Qx = b, what is the error of approximating b? is it the same as Q2.1c?

disp('i1:Quest Show that the columns of Q is in the column space of A?')
[Q,R] = qr(A,0); % using economy QR
tmpCount1 = 0;
for i=1:nc
    qi = Q(:,i);
    if (rank([A qi],1e-6)~=rank(A,1e-6))
        tmpCount1=tmpCount1+1;
    end
end
if (tmpCount1  == 0)
    opStr = sprintf('i1)Ans: TRUE Q col are in col space of A - using Rank([A qi])');
else 
    opStr = sprintf('i1)Ans: FALSE Q col are NOT in col space of A - using Rank([A qi])');
end
disp(opStr);


disp('i2) solve Qx = b, what is the error of approximating b? is it the same as Q1c?')
est_xQ = pinv(Q)*b;
err_b_xQ = b- (Q*est_xQ);
if (norm(err_b_xQ-err_b_x1b) < 1e-6)
    opStr = sprintf('i1)Ans: TRUE , error approximating b is the same  (%.3e  (QR) vs (%.3e (normalEqn) vs (%.3e (pinv)', norm(err_b_xQ), norm(err_b_x1b), norm(err_b_x1a));
else
    opStr = sprintf('i1)Ans: FALSE, error approximating b is the DIFFERENT  (%.3e  (QR) vs (%.3e (normalEqn)  vs (%.3e (pinv)', norm(err_b_xQ), norm(err_b_x1b), norm(err_b_x1a));
end
disp(opStr);


disp('j) Evaluate N point (e.g, 10) around the least squares solution (normal equation)when it exist to show that it gives the least error squares in approximating b. Make delta_x has norm(1)')
N=10;
errExample = zeros(N,1);
for i=1:N
    tmp_delta_x = rand(nc,1);
    tmp_delta_x = tmp_delta_x/norm(tmp_delta_x);
    tmp_x = x1b + tmp_delta_x;
    tmp_hat_b = A*tmp_x;
    errExample(i) = norm(b-tmp_hat_b);
end
mean_errExample = mean(abs(errExample));
opStr = sprintf('norm approximation by LS = %.3e, by %d examples (minErr=%.3e, meanErr=%.3e)',norm(err_b_x1b), N, min(errExample), mean_errExample);
opStr


end  % of function
