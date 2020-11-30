% running swift with minimal parameter assumptions

swift_location = pwd;
dir = uigetdir;

exp_displacement = '55'; % (nm) expected displacement of a single molecules ~30nm for fixed cells | 164nm old value
max_displacement = '125'; % (nm) maximum displacement
max_blinking_duration = '0'; % (frames) maximum duration of the blinking state
tau = '20'; % (ms) time between consecuitive frames
precision = '55'; % (nm) localization precision
exp_noise_rate = '10'; % (%) false positive localizations

p_bleach = '0.0001'; % probability per frame to bleach
p_blink = '0'; % probability per frame to blink
p_switch = '0'; % probability per frame to switch between diffusion states

swift_command = ['.\swift_cli.exe ' dir '\*filt.csv -f --exp_displacement ' exp_displacement ' --max_displacement ' max_displacement ' --exp_noise_rate ' exp_noise_rate ' --p_bleach ' p_bleach ' --p_blink ' p_blink ' --p_switch ' p_switch ' --max_blinking_duration ' max_blinking_duration ' --tau ' tau ' --precision ' precision ' --out_values "mjd mjd_std_err mjd_n track_pos D_msd dynamics" --separate_files --splitby "cell_id"'];

system(['cd ' swift_location]);
[status, cmdout] = system(swift_command);

disp('Done.');