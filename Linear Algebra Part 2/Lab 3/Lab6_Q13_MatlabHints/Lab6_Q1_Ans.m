% Lab 6, Q1
% Author: Chng Eng Siong

% Q1a
% Each row represents 1 person, to measure cosine similarity, 
% the one with the largest cosine similarity is the most similar.
% Q1a)
%Answer: the temperature should be ignored as it flutuates daily.

% Q1b)
A1 =   [62.2, 171.1, 17, 36.8; 
        52.2, 162.6, 19, 36.5;
        72.3, 178.2, 22, 36.7;
        80.8, 185.2, 24, 37.9;
        73.5, 178.3, 22, 37.8];

A1(:,4)=0.*A1(:,4);
% we will ignore temperature

x = A1(1,:)';   % 1st row
y = A1(2,:)';   % second row

% Q1d)
sprintf('comparing my_norm (%f) vs norm (%f)  on pax1 of A1: ',   my_norm(x), norm(x))
sprintf('comparing pax1-pax2 : my_dot (%f) vs dot(%f): ',   my_dot(x,y), x'*y)
sprintf('comparing pax1-pax2 : my_cosSimilarity(%f)',   my_cosSimilarity(x,y))
    

% Q1e)
A1_norm =  my_normalizeMatrix(A1)
x = A1_norm(1,:)';   % 1st row
y = A1_norm(2,:)';   % second row
% Q1f)
sprintf('comparing my_norm (%f) vs norm (%f)  on pax1 of A1: ',   my_norm(x), norm(x))
sprintf('comparing pax1-pax2 : my_dot (%f) vs dot(%f): ',   my_dot(x,y), x'*y)
sprintf('comparing pax1-pax2 : my_cosSimilarity(%f)',   my_cosSimilarity(x,y))


% Q1g)
cosineSimilarityVec = zeros(nr,1);
x = A1_norm(5,1);
for (i=1



% Q1c    
function ret_val = my_norm(x)
    [nr,nc] = size(x);
    assert(nc==1)
    ret_val = 0;
    for i=1:nr
        ret_val = ret_val+x(i)*x(i);
    end
    ret_val = sqrt(ret_val);
end
   
% Q1c
function ret_val = my_dot(x,y)
    [nr,nc] = size(x);
    assert(nc==1)
    ret_val = 0;
    for i=1:nr
        ret_val = ret_val+x(i)*y(i);
    end
end
   
% Q1c
function ret_val = my_cosSimilarity(x,y)
    ret_val = x'*y/(my_norm(x)*my_norm(y));
end
   

%Q1d    
function retMat = my_normalizeMatrix(X)
[nr,nc] = size(X)
retMat = zeros(nr,nc);
    for i=1:nr
       retMat(i,:) =  X(i,:)./norm(X(i,:));
    end
end    
        