function lmrArray=lrSD(funArray)

    % funArray = [fun1(:) fun2(:) ... funN(:)]

    numPoints=size(funArray,1);
    
    meanValues=mean(funArray,2);
    
    leftSD=zeros(numPoints,1);
    rightSD=zeros(numPoints,1);

    for iter=1:numPoints
        
        indLeft=funArray(iter,:)<meanValues(iter);
        indRight=funArray(iter,:)>meanValues(iter);
        
        if ~isempty(indLeft)
            
            leftSD(iter)=sqrt(sum((meanValues(iter)-funArray(iter,indLeft)).^2)/(sum(indLeft)));
        else
            leftSD(iter)=0;
        end
            
        if ~isempty(indRight)
        
            rightSD(iter)=sqrt(sum((meanValues(iter)-funArray(iter,indRight)).^2)/(sum(indRight)));
        else
            rightSD(iter)=0;
        end
    end
        
%     figure();
%     hold on;
%     
%     for iter=1:size(funArray,2)
%         
%        plot(funArray(:,iter),'-k'); 
%     end
%     
%     plot(meanValues,'-r');
%     hold on;
%     plot(meanValues-leftSD,'--r');
%     plot(meanValues+rightSD,'--r');
    
    lmrArray=[leftSD(:) meanValues(:) rightSD(:)];
end