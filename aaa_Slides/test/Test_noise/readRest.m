function [Rest1, Rest2,Rest11, Rest22] = readRest(filename)
disp(filename)
load(filename,'Rest1','Rest2','Rest11','Rest22');