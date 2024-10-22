function [MVC1, MVC2,MVC11, MVC22] = readMVC(filename)
disp(filename)
load(filename,'MVC1','MVC2','MVC11','MVC22');