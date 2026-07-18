function rawImage = read_raw_image(analyVar, baseFilename, suffix)
% Read one full raw camera image.

filePath = [analyVar.dataDir baseFilename suffix];

fid = fopen(filePath, 'rb', 'ieee-be');

if fid < 0
    error('Could not open raw image file:\n%s', filePath);
end

rawImage = fread(fid, analyVar.matrixSize, '*int16');
fclose(fid);

rawImage = double(rawImage);

end