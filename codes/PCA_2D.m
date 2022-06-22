clc;clear;close all
 s = 0.1             % noise standard deviation
 u = im2double(rgb2gray(imread('peppers.png'))); 
 v = u+s*randn(size(u)); 
 
 f=5;
 Y = im2col(v,[2*f+1 2*f+1],'sliding');
%  Y = Y';
 %%
 [nr, nc] = size(u);
 [mY,C,Yc]=moyCov(Y);
 [X,D] = eig(C);
 %%
 % plot eigenvalues in decreasing order
[D,I] = sort(diag(D), 'descend'); plot(D);                       
X = X(:,I); 
% plot first eigenvectors
figure;colormap(gray);
for i=1:25
    subplot(5,5,i);imagesc(reshape(X(:,i),2*f+1,2*f+1));axis image;
end
%%
 Yproj=X'*Yc;                           % computation of all scalar products X_i^t(Y_k-m_Y)
 eta = 2.1*s;            
 Yproj(abs(Yproj)<eta)=0;               % hard thresholding  
 Z = X*Yproj + mY*ones(1,size(Y,2));    % reconstruction of the patches after hard thresholding
%%
tmp = zeros(nr,nc,(2*f+1)^2);
nb = zeros(nr,nc);
for x = 1:2*f+1
  for y = 1:2*f+1
    i = (2*f+1) *(y-1)+x;
    w = reshape(Z(i,:),nr-2*f,nc-2*f);           % use matrix instead of reshape in scilab
    tmp(x:nr-2*f+x-1,y:nc-2*f+y-1,i) = w; 
    nb(x:nr-2*f+x-1,y:nc-2*f+y-1) = nb(x:nr-2*f+x-1,y:nc-2*f+y-1) +1;    
  end
end
vdenoised = sum(tmp,3)./nb;
%%
figure;imshow([u, v, vdenoised],[])
 %%
function [mY,C,Yc]=moyCov(Y)
 mn = size(Y,2)
 mY  = (sum(Y,2)/mn);
 C = cov(Y');  %sum((Y-mY)*(Y-mY)')/mn
 Yc = minus(Y , mY);
end
