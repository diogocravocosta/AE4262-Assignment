function sizes = showB16(filenames)
% showB16  Display one or more PCO *.b16 images and return their sizes
%
% sizes = showB16(filename)
% sizes = showB16([filename1, filename2, ...])
%
% Output: sizes - Nx2 array where each row is [height, width]
%
% Example:
%   sizes = showB16("image.b16")
%   sizes = showB16(["image1.b16", "image2.b16", "image3.b16"])

sizes = zeros(length(filenames), 2);

for i = 1:length(filenames)
    image = readB16(filenames{i});
    sizes(i,:) = size(image);
    figure;
    imagesc(image);
    colormap(gray);
    colorbar;
    axis image;
    title(filenames{i}, 'Interpreter', 'none');
end
end