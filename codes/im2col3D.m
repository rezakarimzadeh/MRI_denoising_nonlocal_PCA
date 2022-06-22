function Y = im2col3D(im,r,step)
% extract patches with desired size and step
[m,n,w] = size(im);
counter = 1;
for i = 1:step:m-r+1
    for j = 1:step:n-r+1
        for q = 1:step:w-r+1
            temp = im(i:i+r-1,j:j+r-1,q:q+r-1);
            Y(:,counter) = temp(:);
            counter = counter +1;
        end
    end
end
end