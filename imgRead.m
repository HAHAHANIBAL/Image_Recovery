function imgIn = imgRead(fileName)
% load the input image into a matrix;
% INPUTS:
%   fileName: name of the input file
%
% OUTPUTS:
%   imgIn: a matrix of the input image
%
% Examples:  imgIn=imgRead('lena.bmp')
%
% @ 2011 Huapeng Zhou -- huapengz@andrew.cmu.edu

imgIn = imread(fileName);
% imgIn = rgb2gray(imgIn);
imgIn = double(imgIn);

end