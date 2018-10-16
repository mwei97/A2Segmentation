function [norm_pdm,mean_pdm,V, moment_mat] = normalize(pdm)

% Param
% pdm: 128*1 column vector of correct pdm

pdm = reshape(pdm, 64, 2);


% Mean of x column and mean of y column
mean_pdm = mean(pdm, 1);
% Subtract center of mass
pdm = pdm - mean_pdm;

% 2x2 2nd-moment matrix
moment_mat = zeros(2);
% Upperleft
moment_mat(1,1) = dot(pdm(:,1), pdm(:,1));
% Lowerleft and Upperright
moment_mat(1,2) = dot(pdm(:,1), pdm(:,2));
moment_mat(2,1) = dot(pdm(:,1), pdm(:,2));
% Lowerright
moment_mat(2,2) = dot(pdm(:,2), pdm(:,2));

% Eigenvector and eigenvalue
[V,D] = eig(moment_mat);
if (D(1,1) < D(2,2))
    tmp = V(:,1);
    V(:,1) = V(:,2);
    V(:,2) = tmp;
end

% V is now the rotation matrix
% Rotate with its inverse (transpose)
rotation = transpose(V);

pdm = pdm*rotation;

% Change dim back to 128*1
norm_pdm = zeros(128,1);
i = 1;
while i<128
    norm_pdm(i) = pdm((i+1)/2,1);
    norm_pdm(i+1) = pdm((i+1)/2,2);
    i = i+2;
end


end