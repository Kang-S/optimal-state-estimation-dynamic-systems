function v=correlated_noise(r,m)
%function v=correlated_noise(r,m)
%
% This m-file produces an m x n matrix of zero-mean 
% Gaussian noise with covariance r, which includes
% correlated terms.
%
%  The inputs are:
%       r = covariance matrix (nxn)
%       m = number of points
%
%  The output is:
%       v = noise matrix (mxn)

% Written by John L. Crassidis 9/9/04

% Decompose the Covariance Matrix
n=length(r);
[u,r_diag]=eig(r);

% Get Uncorrelated Noise
v_uncorr=randn(m,n)*r_diag.^(0.5);

% Get Correlated Noise
v=(u*v_uncorr')';