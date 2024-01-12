function [ GT ] = read_correspondence( filename )
%READ_CORRESPONDENCE Reads a two column file into 
% a n x 2 vector.

file = fopen(filename,'r');
myformat1 = '%6d %6d\n';
GT = fscanf(file, myformat1);
GT = reshape(GT, 2, size(GT,1)/2)';
fclose(file);

end

