% with photophysical parameters determined from fixed cell measurements

swift_location = pwd;
dir = uigetdir;

exp_displacement = '164';
max_displacement = '1000';
max_blinking_duration = '10';
tau = '20';
precision = '55';
exp_noise_rate = '10';
p_bleach = '0.2739';
p_blink = '0.2741';
p_reappear = '0.1017';
p_switch = '0.001';



swift_command = ['.\swift_cli.exe ' dir '\*filt.csv -f --exp_displacement ' exp_displacement ' --max_displacement ' max_displacement ' --exp_noise_rate ' exp_noise_rate ' --p_reappear ' p_reappear ' --p_bleach ' p_bleach ' --p_blink ' p_blink ' --p_switch ' p_switch ' --max_blinking_duration ' max_blinking_duration ' --tau ' tau ' --precision ' precision ' --out_values "mjd mjd_std_err mjd_n track_pos D_msd dynamics" --separate_files --splitby "cell_id"'];

system(swift_location(1:2));
system(['cd ' swift_location]);
[status, cmdout] = system(swift_command);
%disp(cmdout);
disp('Done.');