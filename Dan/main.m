% Main Test Script
%% Clean Stuff
clear
clc
close all

%% Parameters
c = 40; %number of cities
xy = c*rand(c,2);
N = size(xy,1); %map of locations
a = meshgrid(1:N); %make a grid
dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),N,N);  %creates symetric cost matrix, the diagonol is zeroes
salesmen = 3;
min_tour = 5;
max_tour = 100;
tw = 0; %time window
pop_size = 80; %size of population
num_iter = 1000; %number of iterations within genetic algorithm
use_complex = 0;
show_prog = 1;
show_res = 1;
test_iter = 4; %number of times each GA is called


clr = [1 0 0; 0 0 1; 0.67 0 1; 0 1 0; 1 0.5 0];


%% Run throught the code a ton of freaking time
for i = 1:test_iter
    %% Calling the different algorithms
    [opt_rte_td, min_dist_td, dist_history_td] =  tsp_ga(xy,dmat,pop_size, num_iter);
    [opt_rte_pd, min_dist_pd, smd_pd, dist_history_pd, total_dist_pd, cost_pd] =  mtsp_percent_diff(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);
    [opt_rte_std, min_dist_std, smd_std, dist_history_std, total_dist_std, cost_std] =  mtsp_std_dev(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);

    %[opt_rte_tnt, min_dist_tnt, smd_tnt, dist_history_tnt] =  mtsp_tnt(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);  
    %% Calculating Time/Max Distance traveled by salesman
    time(i,1) = min_dist_td;
    time(i,2) = max(smd_pd);
    time(i,3) = max(smd_std);
    
    if i == test_iter
        time
        fprintf("Average time of single salesman = %1.4f \n", mean(time(:,1)))
        fprintf("std deviation of time of single salesman = %1.4f \n", std(time(:,1)))
        fprintf("Average time of multiple salesman with cost function of percent difference = %1.4f \n", mean(time(:,2)))
        fprintf("std deviation of time of cost function of percent difference = %1.4f \n", std(time(:,2)))
        fprintf("Average time of cost function of standard deviation = %1.4f \n", mean(time(:,3)))
        fprintf("std deviation of time of standard deviation = %1.4f \n", std(time(:,3)))
    end

    %% Plots

    %Plot of City Locations
    if i == 5
        figure('Name','City Map','NumberTitle','off','Color','white')
        title('City Locations');
        plot(xy(:,1),xy(:,2),'k.')
    end

    %Plot comparing Max travel distance/time
    figure('Name','Time of Travel','NumberTitle','off','Color','white')
    bar(time);

    %Plot of Salesmen for cost(total distance)
    figure('Name','Salesman Travel','NumberTitle','off','Color','white')
    subplot(3,1,1);
            rte = opt_rte_td([1:c 1]);
            plot(xy(rte,1),xy(rte,2) );
            title(sprintf('Total Distance for cost(total distance) = %1.4f',min_dist_td));
            hold on;       

    %Plot of Salesmen for cost(percent difference)
    subplot(3,1,2)
    for s = 1:salesmen
            rte = [1 opt_rte_pd.ch{s} 1];
            plot(xy(rte,1),xy(rte,2),'.-','Color',clr(s,:));
            title(sprintf('Total Distance for cost(percent difference) = %1.4f',min_dist_pd));
            hold on;
    end

    %Plot of Salesmen for cost(std_dev)
    subplot(3,1,3)
    for s = 1:salesmen
            rte = [1 opt_rte_pd.ch{s} 1];
            plot(xy(rte,1),xy(rte,2),'.-','Color',clr(s,:));
            title(sprintf('Total Distance for cost(standard deviation) = %1.4f',min_dist_std));
            hold on;
    end


    % Best Solution History Travel Cost Total distance
    figure('Name','Best Solution History','NumberTitle','off','Color','white')
        subplot(3,1,1);
        title('Best Solution History 1');
        plot(dist_history_td)
    % Best Solution History Travel cost Percent Difference
        subplot(3,1,2);
        title('Best Solution History 2');
        plot(dist_history_pd)
    % Best Solution History Travel cost standard deviation
        subplot(3,1,3);
        title('Best Solution History 2');
        plot(dist_history_std)
end