%% swift to SMAUG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minlength = 4; % localizations
maxTraceLen = 50;
SampleRand = 0; % draw set of random tracks yes (1) or no (0) ?
NumSample = 17000; % number of random tracks to draw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reset(RandStream.getGlobalStream,sum(100*clock));
% -column 1: track number
% -column 2: step number within that track
% -column 3: a placeholder I used for some simulations and testing
% -column 4: x position (in pixels)
% -column 5: y position (in pixels)

[files, path] = uigetfile('*filt.csv','MultiSelect','on');
% if only one file is selected
if ~iscell(files)
    file_count = 1;
    files = {files};
else
    file_count = length(files);
end

datall = struct('tracks', []);
len_init = 0;
for movie = 1:file_count
    in_filename = [path files{movie}];    
    tracks = importdata([in_filename(1:end-4) '.tracked.loc.txt'], ',', 1);    
    datall(movie).tracks = tracks.data;
    len_init = len_init + length(tracks.data);
end

trfile = zeros(len_init, 5);
trackcount = 0;
currpos = 1;
for i = 1:length(datall)
   numtracks = max(datall(i).tracks(:,18));
   for ii = 1:numtracks
       tracklength = sum(datall(i).tracks(:,18)==ii);
       if tracklength >= minlength && tracklength <= maxTraceLen
           trackcount = trackcount + 1;
           trackarray = datall(i).tracks(datall(i).tracks(:,18)==ii,:);
           trfile(currpos:currpos+tracklength-1,:) = [trackcount.*ones(tracklength, 1) cumsum([1; diff(trackarray(:,1))]) zeros(tracklength, 1) trackarray(:,2:3)];
           currpos = currpos+tracklength;
       end
   end
end

trfile(trfile(:,1)==0,:) = [];

if SampleRand
    if NumSample >= trackcount
        error(['Not enough tracks (' num2str(trackcount) ')! Reduce NumSample or switch off random sampling.']);
    end
    trfile_rand = [];
    iter = 0;
    rand_int = datasample(1:trackcount, NumSample, 'replace', false);
    for sel = rand_int
        iter = iter + 1;
        rand_track = trfile(trfile(:,1) == sel, :);
        rand_track(:,1) = iter;
        trfile_rand = [trfile_rand; rand_track];
    end
    
    trfile = trfile_rand;
    disp([num2str(NumSample) ' tracks randomly drawn from total of ' num2str(trackcount) '.']);
else
    disp([num2str(trackcount) ' tracks extracted to "trfile"']);
end