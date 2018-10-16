% Load the cpdm mat
cpdms = load('correctpdms.mat');
cpdms = cpdms.correctpdms;

% Load the greyimages
greyimages = load('greyimages.mat');
greyimages = greyimages.greyimages;


% Normalization and get factor for initialization
center_mass = 0;
rotation = zeros(2);
norm_cpdms = cpdms;
for i=1:32
    cpdm = cpdms(:,i);
    [norm_pdm, mean_pdm, V] = normalize(cpdm);
    norm_cpdms(:,i) = norm_pdm;
    center_mass = center_mass + mean_pdm;
    rotation = rotation+V;
end

center_mass = center_mass/i;
rotation = rotation/i;


% Randomly divide pdm into training and testing sets
idx = randperm(32);
train_idx = (idx<=24);
test_idx = (idx>24);

% train and test set
train = norm_cpdms(:,train_idx);
test = norm_cpdms(:,test_idx);
save('test.mat', 'test');
test_ori = cpdms(:,test_idx);
save('test_ori.mat', 'test_ori');
%train_img = greyimages(:,train_idx);
test_img = greyimages(:,test_idx);


% Calculate the mean object PDM
mean_obj_pdm = mean(train, 2);
% Subtract mean_obj_pdm
train_norm = train - mean_obj_pdm;


% Do PCA, get eigenvector & eigenvalue
% coeff_k is the chosen eigenvectors
[coeff, scr, var] = pca(transpose(train_norm));
k = 6;
[~,id_k]=sort(var, 'descend');
coeff_k = coeff(:,id_k(1:k));
% Get the standard deviation rather than variance
sd = sqrt(var(1:k));


% Apply to test image
% Get the image correct pdm
for im=1:8
    pdm = test(:,im);
    greyimage = reshape(test_img(:,im),256,256);
    [derx, dery] = Derivative(double(greyimage), 1);

    % Loop to get new_pdm and new coef factor
    step = 0;
    while step<12
        % Fix up points (the ?normal? and directional derivative part)
        new_pdm = normal(pdm, derx, dery);

        % Subtract mean object pdm
        new_pdm = new_pdm - mean_obj_pdm;

        % Compute ai coefficient factor
        a = zeros(1,k);
        for i=1:k
            a(i) = dot(new_pdm, coeff_k(:,i)) / sd(i);

            if a(i) > 2.5
                a(i) = 2.5;
            elseif a(i) < -2.5
                a(i) = -2.5;
            end
        end

        % Projected pdm
        sum = zeros(128,1);
        for i = 1:k
            sum = sum + a(i)*var(i)*coeff_k(:,i);
        end
        proj_pdm = mean_obj_pdm + sum;
        pdm = proj_pdm;

        step = step+1;
    end

    img = test_img(:,im);
    name = strcat('img ',string(im),'.mat');
    save(name, 'img');
    name = strcat('pdm ',string(im),'.mat');
    save(name, 'pdm');

end
