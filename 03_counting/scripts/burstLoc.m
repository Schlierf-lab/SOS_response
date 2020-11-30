function [bStartAcc bLengthAcc]=burstLoc(Arr,minDistance)

    Arr=[Arr;Arr(end)];
    bIndex=find(Arr(2:end)-Arr(1:end-1)>1)+1;
    bIndex=[1;bIndex];
    bLength=zeros(size(bIndex,1),1);
    bStart=Arr(bIndex);

    for iter=1:size(bIndex,1)

        length=1;
        index=bIndex(iter);
        while Arr(index+1)==Arr(index)+1 
            length=length+1;
            index=index+1;
        end

        bLength(iter)=length;
    end
    
    bStartAcc=bStart(1);
    bLengthAcc=bLength(1);
    for iter=2:size(bStart)
        
        if (bStart(iter)-(bStart(iter-1)+bLength(iter-1)-1))>minDistance
            
            bStartAcc=[bStartAcc;bStart(iter)];
            bLengthAcc=[bLengthAcc;bLength(iter)];
        end
    end
end