
clear all
close all
clc

%/* Control Parameters of ABC algorithm*/
NP=50; %/* The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2; %/*The number of food sources equals the half of the colony size*/
maxCycle=100; %/*The number of cycles for foraging {a stopping criteria}*/


%/* Problem specific variables*/
objfun='nnTrainer'; %cost function to be optimized
D=93; %/*The number of parameters of the problem to be optimized*/
ub=ones(1,D)*0.5; %/*lower bounds of the parameters.ones: generate an array full of 1*/
lb=ones(1,D)*(-0.5);%/*upper bound of the parameters.*/
%/*A food source which could not be improved through "limit" trials is abandoned*/
%limit=FoodNumber*D; 
limit=10; 
runtime=1;%/*Algorithm can be run many times in order to see its robustness*/

%%%%%%%%%%%%%%%%%%%%%% PARAMETERS USED %%%%%%%%%%%%%%%%%%%%
%ObjVal[FoodNumber];  /*f is a vector holding objective function values associated with food sources */
%Fitness[FoodNumber]; /*fitness is a vector holding fitness (quality) values associated with food sources*/
%trial[FoodNumber]; /*trial is a vector holding trial numbers through which solutions can not be improved*/
%prob[FoodNumber]; /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
%solution [D]; /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomlu chosen solution different from i*/
%ObjValSol; /*Objective function value of new solution*/
%FitnessSol; /*Fitness value of new solution*/
%neighbour, param2change; /*param2change corrresponds to j, neighbour corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
%GlobalMin; /*Optimum solution obtained by ABC algorithm*/
%GlobalParams[D]; /*Parameters of the optimum solution*/
%GlobalMins[runtime]; /*GlobalMins holds the GlobalMin of each run in multiple runs*/

GlobalMins=zeros(1,runtime);

for r=1:runtime
  
% /*All food sources are initialized */
%/*Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub */
    
%MAX and MIN

Range = repmat((ub-lb),[FoodNumber 1]);
Lower = repmat(lb, [FoodNumber 1]);
%Foods is the population of food sources. 
%Each row of Foods matrix is a vector holding D parameters to be optimized. 
%The number of rows of Foods matrix equals to the FoodNumber
Foods = rand(FoodNumber,D) .* Range + Lower;
%fprintf('foods(1,1)=%g\n',Foods(1,1));
fprintf('foods size=%d\n',size(Foods));
%fprintf('length(size(Foods))=%d\n',length(size(Foods)));


ObjVal=feval(objfun,Foods);
%fprintf('ObjVal=%d\n',ObjVal);

fprintf('FOODS OVER\n');

Fitness=calculateFitness(ObjVal);

%reset trial counters
trial=zeros(1,FoodNumber);

%/*The best food source is memorized*/
BestInd=find(ObjVal==min(ObjVal));
BestInd=BestInd(end);
GlobalMin=ObjVal(BestInd);
GlobalParams=Foods(BestInd,:);


iter=1;
%write optimized weights and bias to file
fid=fopen('WeightsAndBias.txt','a+');
fprintf(fid,'*********************WeightsAndBias*****************************\r\n');
while ((iter <= maxCycle)),

%%%%%%%%% EMPLOYED BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:(FoodNumber)
        
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
      
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        ind=find(sol<lb);
        sol(ind)=lb(ind);
        ind=find(sol>ub);
        sol(ind)=ub(ind);
        
        %evaluate new solution
        ObjValSol=feval(objfun,sol);
        %fprintf('FOODS2 OVER\n');
        FitnessSol=calculateFitness(ObjValSol);
        
       % /*a greedy selection is applied between the current solution i and its mutant*/
       %/*If the mutant solution is better than the current solution i, 
       %replace the solution with the mutant and reset the trial counter of solution i*/
       if (FitnessSol>Fitness(i)) 
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
         
         
    end;

%%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/* A food source is chosen with the probability which is proportioal to its quality*/
%/*Different schemes can be used to calculate the probability values*/
%/*For example prob(i)=fitness(i)/sum(fitness)*/ 
%/*or in a way used in the below prob(i)=a*fitness(i)/max(fitness)+b*/

prob=(0.9.*Fitness./max(Fitness))+0.1;
  
%%%%%%%%%%%%%%%%%%%%%%%% ONLOOKER BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))
        t=t+1;
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*D)+1;
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
            while(neighbour==i)
                neighbour=fix(rand*(FoodNumber))+1;
            end;
        
       sol=Foods(i,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       %mutant operation
       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2;
        
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        ind=find(sol<lb);
        sol(ind)=lb(ind);
        ind=find(sol>ub);
        sol(ind)=ub(ind);
        
        %evaluate new solution
        ObjValSol=feval(objfun,sol);
       
        FitnessSol=calculateFitness(ObjValSol);
        
       % /*a greedy selection is applied between the current solution i and its mutant*/
       %/*If the mutant solution is better than the current solution i, 
       %replace the solution with the mutant and reset the trial counter of solution i*/
       if (FitnessSol>Fitness(i)) 
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
        else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end;
    end;
    
    i=i+1;
    if (i==(FoodNumber)+1) 
        i=1;
    end;   
end; 


%/*The best food source is memorized*/
         ind=find(ObjVal==min(ObjVal));
         ind=ind(end);
         if (ObjVal(ind)<GlobalMin)
         GlobalMin=ObjVal(ind);
         GlobalParams=Foods(ind,:);
         end;
         
         
%%%%%%%%%%%% SCOUT BEE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%/*determine the food sources whose trial counter exceeds the "limit" value. 
%In Basic ABC, only one scout is allowed to occur in each cycle*/

ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    Bas(ind)=0;
    sol=(ub-lb).*rand(1,D)+lb;
    ObjValSol=feval(objfun,sol);
    FitnessSol=calculateFitness(ObjValSol);
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
end;



fprintf('iter=%d ObjVal=%g\n',iter,GlobalMin);

%fprintf('%g\n',GlobalParams);
  

iter=iter+1;



end % End of ABC

GlobalMins(r)=GlobalMin;
end; %end of runs

fprintf('%g\n',GlobalParams);
fprintf(fid,'%g\r\n',GlobalParams);
fclose(fid);

save all

