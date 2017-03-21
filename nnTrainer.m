%   Training Artificial Neural Network.
%   Input Layer: 5, Hidden Layer: 10, Output Layer: 3 neurons. 
%	Returns mean square error between desired and actual outputs.

% function [output] = name(input)
function ObjVal = nnTrainer(Chrom)
[Nind Nvar]=size(Chrom);

%read input values file
dataIn=load('LR1.txt');
%read output values file
dataOu=load('LR22.txt');


for i=1:Nind

	x=Chrom(i,:);
	ObjVal(i)=0;
	mseSum=0;
	
	for n=1:4:60
	
		hidden=10;
		trin=dataIn([n,n+1,n+2,n+3],:);
		trout=dataOu([n,n+1,n+2,n+3],:);
		
		inpRow=size(trin,1);%input rows number, rows of data(4)
		inpCol=size(trin,2);%input colums number, number of inputs(5)
		
		outRow=size(trout,1);%output rows number, rows of data(4)
		outCol=size(trout,2);%output cols number, number of outputs(3)
    
        w1 = reshape(x(1:hidden*inpCol),hidden,inpCol);
        b1 = reshape(x(hidden*inpCol+1:hidden*inpCol+hidden),hidden,1);
        w2 = reshape(x(hidden*inpCol+hidden+1:hidden*inpCol+hidden+hidden*outCol),outCol,hidden);
        b2 = reshape(x(hidden*inpCol+hidden+hidden*outCol+1:hidden*inpCol+hidden+hidden*outCol+outCol),outCol,1);
    
    
        %logsig: activaton function
        y = logsig(logsig(trin*w1'+repmat(b1',size(trin,1),1))*w2'+repmat(b2',size(trin,1),1));
    
    
        output=reshape(y,12,1);
        desire=reshape(trout,12,1);
    
		mseSum=mseSum+mse(output-desire);
        
    end;%for n=1:4:60
	
	ObjVal(i)=mseSum/15;%objval(i) is the mse of 60 rows of data, x(i) as the weights and bias.
   
end;%for i=1:Nind



