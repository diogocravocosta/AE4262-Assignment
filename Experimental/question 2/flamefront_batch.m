function results = flamefront_batch(filenames,baseline,centers,tops,thersholds,verbose)



n = length(filenames);
results = struct('filename', {}, 'flame_left', {},'flame_center', {}, 'flame_right', {});

for i = 1:n

    [fl,fc, fr] = flameFront(filenames{i}, centers(i), thersholds(i), baseline(i), tops(i), verbose);
    results(i).filename = filenames{i};
    results(i).flame_left = fl;
    results(i).flame_center = fc;
    results(i).flame_right = fr;
end


end