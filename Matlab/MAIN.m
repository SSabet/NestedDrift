%==========================================================================
% Supporting code for 'The Nested Drift Algorithm'
% Written by: Soroush Sabet & Patrick Schneider
% Date:       20 November 2024
%
% This code solves the household block of a continuous time model with
% idiosyncratic risk and two-assets with a convex adjustment cost on 
% transfers between the two as in e.g. Kaplan, Moll & Violante (2018)
%
% Key functions:
% - The household HJB updates are computed in the 'updateHousehold' function 
%   with the nested-drift algorithm. 
% - The updating process is guided in the 'solveHousehold' function by the
%   tentative guess algorithm to help when the starting point is too far away 
%   from the solution.
%
%==========================================================================

clear all; close all; clc;

% -------------------------------------------------------------------------
% Read parameters & grids into memory

par   = parameters;
grids = makegrids(par);

% Check the parameter values are valid
assert(par.rho - par.ra > 0, 'ra too large, more than compensates discount rate')
assert(par.rho - par.rb_pos > 0, 'rb_pos too large, more than compensates discount rate')
assert(par.chi0 < 1,'chi0 prohibitvely large, not interesting!');
assert(par.ra*par.chi1 < 1-par.chi0, "ra too large, illiquid wealth accumulates faster than the fastest withdrawal rate that's optimal under adjustment costs -> desired infinite accumulation!");
assert(par.rb_neg * par.bmin + (1-par.xi)*par.w*par.z(1) > 0, "Liquid flow must be positive at all points in the state space! Consider raising w, z(low), or moving the borrowing constraint closer to zero");

% -------------------------------------------------------------------------
% Solve stationary policies, value & measure

sol = solveHousehold(nan,grids,par);

% -------------------------------------------------------------------------
% Make plots

stationaryFigures;
