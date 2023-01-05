clc;clear;close all;

%% Define input file
inputFile = 'problem4.inp';

%% Solve Model
[xc,yc,xf,yf,A,b,U] = ss2d_conv(inputFile);

%% Write out model
writeModel(A,b);