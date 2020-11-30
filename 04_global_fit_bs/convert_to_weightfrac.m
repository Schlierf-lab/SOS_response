load('S119A_res.mat');

D = [0.033, 0.205, 0.913, 3.584];
z_all = (allP1/D(1) + allP2/D(2) + allP3/D(3) + (1-allP1-allP2-allP3)/D(4));

allP1_log = (allP1 / D(1)) ./ z_all;
allP2_log = (allP2 / D(2)) ./ z_all;
allP3_log = (allP3 / D(3)) ./ z_all;

z_arr = (arrP1/D(1) + arrP2/D(2) + arrP3/D(3) + (1-arrP1-arrP2-arrP3)/D(4));
z_arr(z_arr < 0) = 1;

arrP1_log = (arrP1 / D(1)) ./ z_arr;
arrP1_log(arrP1_log < 0) = -1;
arrP2_log = (arrP2 / D(2)) ./ z_arr;
arrP2_log(arrP2_log < 0) = -1;
arrP3_log = (arrP3 / D(3)) ./ z_arr;
arrP3_log(arrP3_log < 0) = -1;

save('S119A_log_res.mat', 'allP1_log', 'allP2_log', 'allP3_log', 'arrP1_log', 'arrP2_log', 'arrP3_log');