% Draw Nsamp random numbers from a arbitrary monotonical increasing
% probability density function
%
% pdfX - x array of the probability density function
% pdfY - y array of the probability density function
%
% Andreas Hartmann 2017-04-17

function valueRND=randDISTR(pdfX,pdfY,Nsamp)

    cumY=cumsum(pdfY);
    cumY=cumY./max(cumY);
    
    valueRND=zeros(Nsamp,1);
    randNUM=rand(Nsamp,1);
    
    for iter=1:Nsamp
    
        [~,minIDX]=min(abs(cumY-randNUM(iter)));
        valueRND(iter)=pdfX(minIDX);
    end
end