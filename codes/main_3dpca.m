clc; clear; close all

% read volume
name ='t1_icbm_normal_1mm_pn0_rf0.rawb';
fid = fopen(name,'r');    
s=[181,217,181];
ima=zeros(s(1:3));
for z=1:s(3),    
  ima(:,:,z) = fread(fid,s(1:2),'uchar');
end;
fclose(fid);
ima=double(ima);



% sigma = 20;
Rician = 1 %rician distribution 1 & gaussian 0
variable = 1 %variable=0(Homogeneus noise)  variable=1(spatially variable noise)

% subvolume (do a test with a smaller volume)
ima=ima(50:80,50:80,50:60);
s=size(ima);
for i=1:2:9
  i  
sigma=i*max(ima(:))/100;
randn('seed',0)
if(variable)
  map = ones(3,3,3);
  map(2,2,2)=3;
  [x1,y1,z1] = meshgrid(1:3,1:3,1:3);
  [x2,y2,z2] = meshgrid(1:2/(s(2)-1):3,1:2/(s(1)-1):3,1:2/(s(3)-1):3);
  map = sigma*interp3(x1,y1,z1,map,x2,y2,z2,'cubic'); 
  if(Rician) rima=sqrt((ima+randn(size(ima)).*map).^2+(randn(size(ima)).*map).^2);
  else       rima=ima+randn(size(ima)).*map;
  end
else
  if(Rician) rima=sqrt((ima+randn(size(ima))*sigma).^2+(randn(size(ima))*sigma).^2);
  else       rima=ima+randn(size(ima))*sigma;
  end 
  map=ones(s)*sigma;
end


%% Median Filter
med_filt = medfilt3(rima,[3 3 3]);
%% Non Local PCA

t = 3; %block half size (2*t+1)
r = 4;  % patch size
step = 1;
[denoised_step1, map_est] = NL_PCA(med_filt,r,t,step,Rician);
%% Non Local mean
v = 5 %half size of search window
t = 1  %half size of surrounding window in pixel (i,j,k)
denoised_step2 = PRI_NL_PCA(denoised_step1,v,t,map_est,Rician);
%% show result
figure;
clf
colormap(gray);
n=round(size(denoised_step1,3)/2);
suptitle(['Results for ', num2str(i), '% Noise Level'])
subplot(2,2,1),imagesc(imrotate(ima(:,:,n),90));title('Original Noise free Image')
subplot(2,2,2),imagesc(imrotate(rima(:,:,n),90));;title('Noisey Image')
subplot(2,2,3),imagesc(imrotate(denoised_step1(:,:,n),90));;title('NL-PCA denoised')
subplot(2,2,4),imagesc(imrotate(denoised_step2(:,:,n),90));;title('PRI-NL-PCA denoised')
%% fidelity
% find max intensity of image for psnr
R=max(ima(:));
indi=find(ima>10);
sw = [1 1 1]; 
% for noisy image
oerror0(i)=sqrt(mean((ima(indi)-rima(indi)).^2));
opsnr0(i)=20*log10(R/oerror0(i));
ossim0(i)= ssim_index3d(rima,ima,sw,(indi));
% for NL_PCA
oerror1(i)=sqrt(mean((ima(indi)-denoised_step1(indi)).^2));
opsnr1(i)=20*log10(R/oerror1(i));
ossim1(i)= ssim_index3d(denoised_step1,ima,sw,indi);
% for PRI-NL-PCA
oerror2(i)=sqrt(mean((ima(indi)-denoised_step2(indi)).^2));
opsnr2(i)=20*log10(R/oerror2(i));
ossim2(i)= ssim_index3d(denoised_step2,ima,sw,indi);

ER = abs(1-map_est./map);
MER(i) = mean(ER(:));
end
%% plots
op0=mean(opsnr0(1:2:9))
op1=mean(opsnr1(1:2:9))
op2=mean(opsnr2(1:2:9))
 
figure
title('PSNR for noisy, NL-PCA and PRI-NL-PCA output')
clf
plot(1:2:9,opsnr0(1:2:9),'g')
hold on
plot(1:2:9,opsnr1(1:2:9),'b')
plot(1:2:9,opsnr2(1:2:9),'r')
xlabel('Noise level(%)')
ylabel('PSNR')
legend('noisy', 'NL-PCA', 'PRI-NL-PCA')
os0=mean(ossim0(1:2:9))
os1=mean(ossim1(1:2:9))
os2=mean(ossim2(1:2:9))
 
figure
title('SSIM for noisy, NL-PCA and PRI-NL-PCA output')
clf
plot(1:2:9,ossim0(1:2:9),'g')
hold on
plot(1:2:9,ossim1(1:2:9),'b')
plot(1:2:9,ossim2(1:2:9),'r')
xlabel('Noise level(%)')
ylabel('SSIM')
legend('noisy', 'NL-PCA', 'PRI-NL-PCA')
%% ER
figure
title('MER for noise levels')
clf
plot(1:2:9,MER(1:2:9)/4,'g')
xlabel('Noise level(%)')
ylabel('MER')
 %% functions
function out = etta(phi)
phi24 = (phi^2)/4;
out =( sqrt(pi/2)*exp(-phi24)*((1+2*phi24)*besseli(0,phi24) + (2*phi24)*besseli(1,phi24))^2);
end

function [final_img, map] = NL_PCA(noisy_im,r,t,step,Rician)
s = size(noisy_im);
denoised_img = zeros(s);
normalization = zeros(s);
map = zeros(s);
for ii=1:s(1)-2*t
    for jj=1:s(2)-2*t
        for ll = 1:s(3)-2*t
            noisy_block = noisy_im(ii:ii+2*t,jj:jj+2*t,ll:ll+2*t);
            [denoised_block, block_map] = Nl_PCA_blockwise(noisy_block,r,step,Rician);
            denoised_img(ii:ii+2*t,jj:jj+2*t,ll:ll+2*t) =denoised_block +  denoised_img(ii:ii+2*t,jj:jj+2*t,ll:ll+2*t);
            map(ii:ii+2*t,jj:jj+2*t,ll:ll+2*t) =block_map +  map(ii:ii+2*t,jj:jj+2*t,ll:ll+2*t);
            normalization(ii:ii+2*t,jj:jj+2*t,ll:ll+2*t) =1 +  normalization(ii:ii+2*t,jj:jj+2*t,ll:ll+2*t);
        end
    end
end
final_img = denoised_img./normalization;
map = map./normalization;
if Rician
bias_corrected_img = zeros(s);
mdi = max(final_img(:));
for ii=1:s(1)
    for jj=1:s(2)
        for ll = 1:s(3)
           bias_corrected_img(ii,jj,ll)= map(ii,jj,ll) * etta( (final_img(ii,jj,ll)/(map(ii,jj,ll)*mdi))); 
        end
    end
end
% final_img = final_img - bias_corrected_img;
end
end


function  [denoised_block, block_map] = Nl_PCA_blockwise(noisy_block,r,step,Rician)
%% pca
Y = im2col3D(noisy_block,r,step);
[mY,C,Yc]=moyCov(Y);
 [X,D] = eig(C);
 %% eigenvalues in decreasing order
[D,I] = sort(diag(D), 'descend'); %plot(D);                       
X = X(:,I); 
%% Noise estimation
est_sigma = Noise_std_est(D,Rician, Y);

%% hard thresholding and vector patch reconstruction
Z = hrdThresh_vPatchRec(X,Yc,mY,Y,est_sigma);

 %% reconstructing original block
block_shape = size(noisy_block);
denoised_block = block_reconstruction(Z,block_shape,r,step);
block_map = est_sigma*ones(block_shape);
end

 
 
 function Z = hrdThresh_vPatchRec(X,Yc,mY,Y,est_sigma)
 Yproj=X'*Yc;                           % computation of all scalar products X_i^t(Y_k-m_Y)
 eta = 2.1*est_sigma;            
 Yproj(abs(Yproj)<eta)=0;               % hard thresholding  
 Z = X*Yproj + mY*ones(1,size(Y,2));    % reconstruction of the patches after hard thresholding
 end

function est_sigma = Noise_std_est(D,Rician,Y)
med = median(sqrt(D));
lambda_t = 0;
for i=1:length(D)
    if sqrt(D(i)) < 2*sqrt(med)
        lambda_t = [lambda_t D(i)];
    end
end
lambda_t = lambda_t(2:end);
beta = 12.9;
est_sigma = beta*sqrt(abs(median(lambda_t)));

if Rician
    l_m  = mean(Y(:));
    l_std = std(Y(:));
    gamma = l_m/l_std; 
    if gamma >1.86
        phi_gamma  = (0.9846*(gamma-1.86)+0.1983)/((gamma-1.86)+0.1175);
    else
        phi_gamma  = 0;
    end
    est_sigma = est_sigma*phi_gamma;
end
end

 function vdenoised = block_reconstruction(Z,block_shape,r,step)
 nr = block_shape(1);nc = block_shape(2);nw = block_shape(3);
 tmp = zeros(nr,nc,nw,(nr*nc*nw));
nb = zeros(nr,nc,nw);
i=0;
for x = 1:step:nr-r+1
  for y = 1:step:nc-r+1
      for l =1:step:nw-r+1
    i = i+1;
    w = reshape(Z(:,i),r,r,r);           % use matrix instead of reshape in scilab
    tmp(x:r+x-1,y:r+y-1,l:r+l-1,i) = w; 
    nb(x:r+x-1,y:r+y-1,l:r+l-1) = nb(x:r+x-1,y:r+y-1,l:r+l-1) +1;    
  end
  end
end
vdenoised = (sum(tmp,4))./nb;
 end

