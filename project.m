%Project by Students:
%
%	Παναγιώτης Πράττης / Panagiotis Prattis
%


%Code for:
%	Isolate columns, join tables and create a new matrix 
%	containing the common elements of the u, u1 matrixes
%
u = u(:, [1 2]);
u1 = u1(:, [1 6:24]);
[~, length_i] = ismember(u(:, 2), u1(:, 1));
Matrix = [u, u1(length_i, 2:20)];
Matrix = Matrix(:, [1 3:21]);

%Code for:
%	Create a new matrix whose each line is a 
%	vector of attributes for each user
%
[C, ~, Matrix_id] = unique(Matrix(:, 1));
N = numel(C);
M = zeros(N, size(Matrix, 2) -1);
for length_i = 1:N
	M(length_i, :) = sum(Matrix(Matrix_id == length_i, [2:20]));
end
Matrix = M;

%Code for:
%	Create a new matrix whose each line of vectors is 
%	calculated according to the set of vectors
%
S = size(Matrix);
for length_i = 1:S(1)
    Summary = sum(Matrix(length_i, :));
    for length_j = 1:S(2)
        Matrix_Normalized(length_i,length_j) = Matrix(length_i,length_j)/Summary;
    end
end
Matrix = Matrix_Normalized;

%Code for:
%	Clustering vectors of features, according to a representative vector
%	where the distance between each vector and the representative group vector  
%	is compared with to the threshold of dissimilarity and accordingly 
%	the vector joins the group or a new group is created
%
Matrix = Matrix';

[l,N] = size(Matrix);
D = zeros(N, N);
for i=1:N
    for j = i+1:N
        D(i, j) = sqrt(sum((Matrix(:, i)-Matrix(:, j)).^2));
        D(j, i) = D(i, j);
    end
end
max_i = max(max(D));
min_i = min(D(~logical(eye(N))));
mean_i = (min_i + max_i)/2;
min_theta = 0.25 * mean_i;
max_theta = 1.75 * mean_i;
theta = 50;
s = (max_theta - min_theta)/(theta - 1);
q = N + 1;
nr = 10;
total = [];
for theta = min_theta:s:max_theta
    list = zeros(1,q);
        order = randperm(N);
		[l, N1] = size(Matrix);
		if(length(order) == 0)
			order = 1:1:N1;
		end
		clust = 1;
		[l, N1]=size(Matrix);
		cl = zeros(1, N1);
		cl(order(1)) = clust;
		m = Matrix(:, order(1));
		for i=2:N1
		   [m1, m2] = size(m);
		   [s1, s2] = min(sqrt(sum((m-Matrix(:, order(i))*ones(1, m2)).^2)));
		   if(s1 > theta) && (clust < q)
			   clust = clust + 1;
			   cl(order(i)) = clust;
			   m = [m Matrix(:, order(i))];
		   else
			   cl(order(i)) = s2;
			   m(:, s2)=((sum(cl == s2) -1)*m(:, s2) + Matrix(:, order(i)))/sum(cl == s2);
		   end
		end
        list(size(m, 2)) = list(size(m, 2))+1;
    end
    [q1,si] = max(list);
    total = [total si];
