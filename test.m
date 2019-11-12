[idx, centers] = kmeans (u, 4);
figure;
plot (u (id==1, 1), u (id==1, 4), 'ro');
hold on;
plot (u (id==2, 1), u (id==2, 4), 'bs');
hold on;
plot (u (id==3, 1), u (id==3, 4), 'y*');
hold on;
plot (u (id==4, 1), u (id==4, 4), 'gx');
plot (centers (:, 1), centers (:, 4), 'kv', 'markersize', 10);
hold off;
%--------------------------------------------------%
u = u(:, [1 2]);
u1 = u1(:, [1 6:24]);
[~, length_i] = ismember(u(:, 2), u1(:, 1));
Matrix = [u, u1(length_i, 2:20)];
Matrix = Matrix(:, [1 3:21]);

%--------------------------------------------------%
%subsum function. A function to calculate sum of submatrix defined by an id
%   A matrix is the input. Id must be on the first column
%   The function calculates the unique occurences of id in the given matrix. 
%   An empty matrix is also created, which size is number of unique ids x range of
%   data(vector).
%   The sum for every data of the vector is then calculated seperately and
%   then added to the new matrix, when their id are the same
%   the new matrix is the output
[C, ~, Matrix_id] = unique(Matrix(:, 1));
N = numel(C);
M = zeros(N, size(Matrix, 2) -1);
for length_i = 1:N
	M(length_i, :) = sum(Matrix(Matrix_id == length_i, [2:20]));
end
Matrix = M;
%--------------------------------------------------%
%UNTITLED5 Given a matrix(Set of Vectors), normalize given data
%For the given matrix:
% For each vector (row)
%Divide each value with the sum value of vector
S = size(Matrix);
for length_i = 1:S(1)
    Summary = sum(Matrix(length_i, :));
    for length_j = 1:S(2)
        Matrix_Normalized(length_i,length_j) = Matrix(length_i,length_j)/Summary;
    end
end
Matrix = Matrix_Normalized;
%--------------------------------------------------%
[coeff,score,latent,tsquared,explained] = pca(Matrix)
%--------------------------------------------------%
[~, score] = pca(Matrix, 'NumComponents', 5);
%--------------------------------------------------%
[coeff,score,latent,tsquared,explained] = pca(Matrix)
%--------------------------------------------------%
[~, score] = pca(Matrix, 'NumComponents', 3);
Matrix = score;
%--------------------------------------------------%
[idx, centers] = kmeans (u, 2);
figure;
plot (u (idx==1, 1), u (idx==1, 2), 'ro');
hold on;
plot (u (idx==2, 1), u (idx==2, 2), 'bs');
plot (centers (:, 1), centers (:, 2), 'kv', 'markersize', 10);
hold off;
%--------------------------------------------------%
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

		% Cluster determination phase
		clust = 1;  % no. of clusters
		[l, N1]=size(Matrix);
		cl = zeros(1, N1);
		cl(order(1)) = clust;
		m = Matrix(:, order(1));
		for i=2:N1
		   [m1, m2] = size(m);
		% Determining the closest cluster msentative}}
		   [s1, s2] = min(sqrt(sum((m-Matrix(:, order(i))*ones(1, m2)).^2)));
		   if(s1 > theta) && (clust < q)
			   clust = clust + 1;
			   cl(order(i)) = clust;
			   m = [m Matrix(:, order(i))];
		   else

		% Pattern classification phase(*4)}}
			   cl(order(i)) = s2;
			   m(:, s2)=((sum(cl == s2) -1)*m(:, s2) + Matrix(:, order(i)))/sum(cl == s2);
		   end
		end
        list(size(m, 2)) = list(size(m, 2))+1;
    end
    [q1,si] = max(list);
    total = [total si];
%--------------------------------------------------%
[theta, cl] = kmeans(Matrix,3);
count = histcounts(theta);
scatter3(Matrix(:,1),Matrix(:,2),Matrix(:,3),10,theta);
%--------------------------------------------------%
Z = linkage(Matrix,'complete','squaredeuclidean');
%--------------------------------------------------%
c = cluster(Z,'maxclust',4);
count2 = histcounts(c);
%--------------------------------------------------%
dendrogram(Z,30);
%--------------------------------------------------%
scatter3(Matrix(:,1),Matrix(:,2),Matrix(:,3),10,c);
