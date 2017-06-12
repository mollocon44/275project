% Main Test Script
%% Clean Stuff
clear
clc
close all


c = 40; %number of cities
salesmen = 3;
min_tour = 0;
max_tour = 2000;
tw = 0; %time window
pop_size = 80; %size of population
num_iter = 50; %number of iterations within genetic algorithm
use_complex = 0;
show_prog = 1;
show_res = 0;
test_iter = 2; %number of times each GA is called
map_iter = 3; %number of maps tested
num_algos = 5; % number of algos
clr = [1 0 0; 0 0 1; 0.67 0 1; 0 1 0; 1 0.5 0];


%% Parameters
for j = 1:map_iter
    
    c = 40; %number of cities
    xy = c*rand(c,2);
    N = size(xy,1); %map of locations
    a = meshgrid(1:N); %make a grid
    dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),N,N);  %creates symetric cost matrix, the diagonol is zeroes

    


    %% Run throught the code a ton of freaking time
    for i = 1:test_iter
        %% Calling the different algorithms
        [opt_rte_td, smd_td, dist_history_td] =  tsp_ga(xy,dmat,pop_size, num_iter);
        [opt_rte_pd, smd_pd, dist_history_pd] =  mtsp_percent_diff(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);
        [opt_rte_std, smd_std, dist_history_std] =  mtsp_std_dev(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);
        [opt_rte_t, smd_t, dist_history_t] = mtsp_tour_2(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);
        [opt_rte_tnt, smd_tnt, dist_history_tnt] = mtsp_tnt(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);

        %[opt_rte_tnt, min_dist_tnt, smd_tnt, dist_history_tnt] =  mtsp_tnt(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);  
        %% Calculating Time/Max Distance traveled by salesman
 
        time(j,i,1) = max(smd_td);
        time(j,i,2) = max(smd_pd);
        time(j,i,3) = max(smd_std);
        time(j,i,4) = max(smd_t);
        time(j,i,5) = max(smd_tnt);


        

        %% MEGA PLOT
     % for the civilized man

        for k = 1:num_algos
            plot(j*ones(1,test_iter),time(j,:,k),'.','Color',clr(k,:));
            hold on;
        end

    end
end    