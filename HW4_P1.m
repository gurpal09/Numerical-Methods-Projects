% Gurpal Singh
% October 25, 2017
% Math 567 Homework 4

% Problem 1 Parts a & b

% Part a
c = [-1/12 4/3 -5/2 4/3 -1/12]; % Stencil Coefficients
x = [-2 -1 0 1 2]; % x values

k = 2; % Order of derivative
p = 4; % Order of Accuracy

fprintf("-------------- Solution --------------\n\n");
fprintf("\n Part a:\n\n");
LeadingOrderTerm(c,x,k,p) % Function to compute and print


% Part b
clear all;

c = [-1/6, 2, -13/2, 28/3, -13,2, 2, -1/6];
x = [-3, -2, -1, 0, 1, 2, 3];

k = 4;
p = 6;

fprintf("\n Part b:\n\n");
LeadingOrderTerm(c,x,k,p)








