%%                PHOTOMETRIC STEREO
% |===================================================|
% | Nguyen Tuan Nghia                                 |
% | nghianguyentuan@snu.ac.kr                         |
% | Seoul National University                         |
% | Department of Electrical and Computer Engineering |
% |===================================================|
%
%% Preprocessing

close all;
clear;
clc;

size_after_crop_y = 176; % CHANGE THIS: Size y of cropped image
size_after_crop_x = 176; % CHANGE THIS: Size x of cropped image

crop_y = 117; % CHANGE THIS: Start pixel y
crop_x = 257; % CHANGE THIS: Start pixel x

crop_range_y = crop_y:crop_y + size_after_crop_y - 1;
crop_range_x = crop_x:crop_x + size_after_crop_x - 1;

[ambimage, imarray, lightdirs] = LoadFaceImages();

imarray = permute(imarray,[2 3 1]);

imarray = imarray(crop_range_y,crop_range_x,:);
ambimage = ambimage(crop_range_y,crop_range_x);

[width, height, nImages] = size(imarray);

% Subtract the ambient image from other images
% Threshold so that no pixel value is negative
imarray = max(imarray - repmat(ambimage,1,1,nImages),0);

imarray = permute(imarray,[3 1 2]);

%% PART 1: Calibrated Photometric Stereo

%% 3 Images:

% Choose random 3 images to calculate Albedo and Normal
% 3 first indexes will be used
ridx = randperm(nImages);
train = ridx(1:3);
test = ridx(4:end);

im_compute = imarray(train,:,:);
li_compute = lightdirs(train,:);

im_compute = reshape(im_compute,[3, width * height]);

kdn = li_compute \ im_compute;
albedo = zeros(size(kdn,2),1);
normal = zeros(size(kdn));
for i = 1:size(kdn,2)
    albedo(i) = norm(kdn(:,i));
    if (albedo(i) == 0)
        normal(:,i) = kdn(:,i);
    else
        normal(:,i) = kdn(:,i) / albedo(i);
    end
end

% Recover Albedo and Normal maps
alb_map = reshape(albedo,[width, height]);
nrm_map = reshape(normal',[width, height, 3]);
nrm_map = 0.5 * nrm_map + 0.5;

% Show results
figure(1);
subplot(1,2,1);
imshow(alb_map);
title('3 Images: Albedo');
subplot(1,2,2);
imshow(nrm_map);
title('3 Images: Normal');

% Reconstruction Test 1;
li_rcs = lightdirs(test,:);
im_rcs = li_rcs * kdn;
im_rcs = reshape(im_rcs',[width,height,8]);

figure(2);
set(figure(2),'units','normalized','outerposition',[0 0 1 1]);
for i = 1:8
    subplot(2,8,i);
    imshow(reshape(imarray(i,:,:),[width,height]));
    title(sprintf('Ground-truth %d',i))
    subplot(2,8,i+8);
    imshow(im_rcs(:,:,i));
    title(sprintf('Reconstructed %d',i))
end

%% Least Square Method

ridx = randperm(nImages);
n_train = 7; % CHANGE THIS: Number of images will be used to find Normal and Albedo
train = ridx(1:n_train);
test = ridx(n_train+1:end);

im_compute = imarray(train,:,:);
li_compute = lightdirs(train,:);

im_compute = reshape(im_compute,[n_train, width * height]);

kdn = (li_compute' * li_compute) \ (li_compute' * im_compute);
albedo = zeros(size(kdn,2),1);
normal = zeros(size(kdn));
for i = 1:size(kdn,2)
    albedo(i) = norm(kdn(:,i));
    if (albedo(i) == 0)
        normal(:,i) = kdn(:,i);
    else
        normal(:,i) = kdn(:,i) / albedo(i);
    end
end

% Recover Albedo and Normal maps
alb_map = reshape(albedo,[width, height]);
nrm_map = reshape(normal',[width, height, 3]);
nrm_map = 0.5 * nrm_map + 0.5;

% Show results
figure(3);
subplot(1,2,1);
imshow(alb_map);
title(sprintf('%d Images: Albedo',n_train));
subplot(1,2,2);
imshow(nrm_map);
title(sprintf('%d Images: Normal',n_train));

% Reconstruction Test 2;
li_rcs = lightdirs(test,:);
im_rcs = li_rcs * kdn;
im_rcs = reshape(im_rcs',[width,height,nImages - n_train]);

figure(4);
set(figure(4),'units','normalized','outerposition',[0 0 1 1]);
for i = 1:nImages - n_train
    subplot(2,nImages - n_train,i);
    imshow(reshape(imarray(i,:,:),[width,height]));
    title(sprintf('Ground-truth %d',i))
    subplot(2,nImages - n_train,i+nImages - n_train);
    imshow(im_rcs(:,:,i));
    title(sprintf('Reconstructed %d',i))
end

%% PART 2: Uncalibrated Photometric Stereo

% Preprocessing
resize_factor = 0.8; % CHANGE THIS: Resize image to work with SVD
imarray = permute(imarray,[2 3 1]);
imarray = imresize(imarray,resize_factor);
[width, height, nImages] = size(imarray);
imarray = permute(imarray,[3 1 2]);

% We will use all 11 images for svd
im_compute = imarray;
im_compute = reshape(im_compute,[nImages, width * height]);
im_compute = im_compute';

% SVD
[U,S,V] = svd(im_compute);

% Since S is sorted, three best eigenvalues are at index 1, 2, 3
u = U(:,1:3);
s = S(1:3,1:3);
v = V(:,1:3);

lightdirs = v;
kdn = u*s;

kdn = permute(kdn,[2,1]);
albedo = zeros(size(kdn,2),1);
normal = zeros(size(kdn));
for i = 1:size(kdn,2)
    albedo(i) = norm(kdn(:,i));
    if (albedo(i) == 0)
        normal(:,i) = kdn(:,i);
    else
        normal(:,i) = kdn(:,i) / albedo(i);
    end
end

% Recover Albedo and Normal maps
alb_map = reshape(albedo,[width, height]);
nrm_map = reshape(normal',[width, height, 3]);
nrm_map = 0.5 * nrm_map + 0.5;

figure(5);
% subplot(1,2,1);
% imshow(alb_map);
% subplot(1,2,2);
imshow(nrm_map);
title('Estimated Normal by SVD');

% Reconstruction from s
% im_rcs = lightdirs * kdn;
% im_rcs = reshape(im_rcs',[width,height,nImages]);

% Reconstruction from s*
rand_s_star = randperm(nImages-3) + 3;
lightdirs_star = V(:,rand_s_star(1:3));
im_rcs_star = lightdirs_star * kdn;
im_rcs_star = reshape(im_rcs_star',[width,height,nImages]);

% Show results
figure(6);
set(figure(6),'units','normalized','outerposition',[0 0 1 1]);
for i = 1:nImages
    subplot(2,nImages,i);
    imshow(reshape(imarray(i,:,:),[width,height]));
    title(sprintf('GT %d',i))
%     subplot(3,nImages,i+nImages);
%     imshow(im_rcs(:,:,i));
    subplot(2,nImages,i+nImages);
    imshow(im_rcs_star(:,:,i));
    title(sprintf('s* %d',i))
end
