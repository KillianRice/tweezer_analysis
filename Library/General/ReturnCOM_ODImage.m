function [xcom, ycom] = ReturnCOM_ODImage(ImageOD)
    [rows,cols] = size(ImageOD);
    [X,Y] = meshgrid(1:cols,1:rows);
    total_mass = sum(ImageOD(:));
    xcom = sum(X(:).*ImageOD(:)) / total_mass;
    ycom = sum(Y(:).*ImageOD(:)) / total_mass;
end