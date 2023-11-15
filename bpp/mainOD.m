%Optimization and Decision 2019/20
%Project 8: Bin Packing and Cuting Stock Problem
%Group 10:
%Nº80998  Name: Pedro Miguel Menezes Ramalho
%Nº85183  Name: Ricardo Miguel Diogo de Oliveira Chin

clc
close all
clear all %#ok


fprintf('\nOptimization and Decision 2019/20\n');
fprintf('This program solves the BPP with Genetic Algorithm or Particle Swarm Optimization\n');

%% Model choice

selec = 0;

while(selec==0) 
    fprintf('\nPick an option for the model:\n');
    fprintf('1: Small model (15 bins optimal)      2: Big model (100 bins optimal)\n');
    fprintf('0: Exit program\n');
    fprintf('\n'); 
    option = input('');
    valid=0;
  
    switch(option)
        case 0
            valid=1;
        case 1
            valid=1;
            model = CreateModel(1);
            
        case 2
            valid=1;
            model = CreateModel(2);
    end
    
        if(valid==0)
            fprintf('\nError! Invalid choice!\n');
        end
       
    if valid == 1 && option ~= 0

        
%% Algorithm choice

selec2=0;
while(selec2==0)
    fprintf('\nPick an option for the solution:\n');
    fprintf('1: Perform GA now          2: Check saved GA solution\n');   
    fprintf('3: Perform PSO now         4: Check saved PSO solution\n');   
    fprintf('0: Press 0 if you want to go back\n'); 
    fprintf('\n'); 
    option2 = input('');
    valid2=0;
    close all
    switch(option2)
        case 0
            valid2=1;
            selec2=1;
        
        case 1
            valid2=1;
            tic
            [AvgCost, BestCost, BestSol, pop, it, GAdata] = GA(model); 
            time = toc;
            [Solution] = Plotter(BestCost, BestSol);
            
            disp(['Iterations performed: ' num2str(it) ' , Best Cost = ' ...
                num2str(BestCost(it)) ' , Simulation time: ' num2str(time) ' seconds']); 
        
        case 3
            valid2=1;
            tic
            [AvgCost, BestCost, BestSol, particle, it, PSOdata] = PSO(model);
            time = toc;
            [Solution] = Plotter(BestCost, BestSol);
            
            disp(['Iterations performed: ' num2str(it) ' , Best Cost = ' ...
                num2str(BestCost(it)) ' , Simulation time: ' num2str(time) ' seconds']); 
        
        case 2
            valid2=1;
            if option == 1
                load('GA15.mat')
            else
                load('GA100.mat')
            end
            [Solution] = Plotter(BestCost, BestSol);
            
            disp(['Iterations performed: ' num2str(it) ' , Best Cost = ' ...
                num2str(BestCost(it)) ' , Simulation time: ' num2str(time) ' seconds']); 
        
        case 4
            valid2=1;
            if option == 1
                load('PSO15.mat')
            else
                load('PSO100.mat')
            end
            [Solution] = Plotter(BestCost, BestSol);
            
            disp(['Iterations performed: ' num2str(it) ' , Best Cost = ' ...
                num2str(BestCost(it)) ' , Simulation time: ' num2str(time) ' seconds']); 
    end 
        if(valid2==0)
            fprintf('\nError! Invalid choice!\n');
        end
end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%First Option             
    elseif valid == 1 && option == 0
        selec=1; %Closes the loop
        fprintf('\nProgram ended.\n\n');          
    end  
end

%Delete non important variables
clear option option2 selec selec2 valid valid2