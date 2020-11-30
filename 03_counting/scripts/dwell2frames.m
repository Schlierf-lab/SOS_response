function framesOUT=dwell2frames(actT,Ton,Toff)

    startF=floor(actT);
    endF=ceil(actT+sum(Ton)+sum(Toff));

    frames=(startF:1:endF);
    exposure=zeros(1,length(frames)-1);

    t1=actT;

    for iterT=1:length(Ton)

        t2=t1+Ton(iterT);

        t=[startF t1 t2 endF];
        amplitude=[0 0 (t2-t1) (t2-t1)];

        intAMP=interp1(t,amplitude,frames);

        exposure=exposure+diff(intAMP);

        if length(Ton)>iterT

            t1=t2+Toff(iterT);
        end
    end

    ff=frames([exposure>=0.5 false]);    
    framesOUT=ff(:);
end