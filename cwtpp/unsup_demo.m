%% A demo script to do unsupervised analysis/visualisation 
clc; close all; clear all
%% load in preprocessed data
[data,dims,mz] = h5toMat(['datacube.h5']);

%% generate a mask using KCA
[kcadata,Imagekca,C] = cluster_analysis((data),dims,3);
mask = (kcadata==2|kcadata==3);% here cluster 2 & 3, check cluster images accordingly
mask = reshape(mask,dims);
% mask = imfill(mask,'holes');
figure,imagesc(mask);axis image

%% filter data matrix with mask & do component analysis
data_r = data;
data_r(~mask,:) = 0;
output = component_analysis((data_r),3,dims);% here set to show first 3 components, can change accordingly
