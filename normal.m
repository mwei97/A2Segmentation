function new_pdm = normal(pdm, derx, dery)

% pdm: pdm of one image
%      a column vector of size 128*1
% greyimage: 256*256 greyimage

pdm_rs = reshape(pdm, 64,2);

new_pdm = zeros(64,2);

%[derx, dery] = Derivative(double(greyimage), 1);
[X,Y] = meshgrid(1:256,1:256);

for i=1:64
    % Get the normal
    p = pdm_rs(i,:);
    
    if i==1
        prev = pdm_rs(64,:);
    else
        prev = pdm_rs(i-1,:);
    end
    
    if i==64
        next = pdm_rs(1,:);
    else
        next = pdm_rs(i+1,:);
    end
    
    diff1 = p-prev;
    diff2 = p-next;
    
    a = diff1(1);
    b = diff1(2);
    div = sqrt(a^2+b^2);
    n1 = [-b/div a/div];
    
    a = diff2(1);
    b = diff2(2);
    div = sqrt(a^2+b^2);
    n2 = [-b/div a/div];
    
    n = (n1+n2)/2;
    
    % Bilinear interpolation for p
    
    grad = [interp2(X,Y,derx,p(1),p(2)) interp2(X,Y,dery,p(1),p(2))];
    
    max = dot(-n, grad);
    max_p = p;
    
    % Traverse through 11 points
    for j = -5:5
        if j==0
            continue
        end
        
        new_p = p + j*n;
        % Bilinear interpolation of new_p gradient
        new_grad = [interp2(X,Y,derx,new_p(1),new_p(2)) interp2(X,Y,dery,new_p(1),new_p(2))];
        if dot(-n, new_grad) > max
            max = dot(-n, new_grad);
            max_p = new_p;
        end
    end
    
    new_pdm(i,:) = max_p; 
end


tmp = zeros(128,1);
i = 1;
while i<128
    tmp(i) = new_pdm((i+1)/2,1);
    tmp(i+1) = new_pdm((i+1)/2,2);
    i = i+2;
end
new_pdm = tmp;


end