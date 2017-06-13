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
num_iter = 80; %number of iterations within genetic algorithm
use_complex = 0;
show_prog = 1;
show_res = 0;
test_iter = 3; %number of times each GA is called
map_iter = 3; %number of maps tested
num_algos = 5; % number of algos
clr = [1 0 0; 0 0 1; 0.67 0 1; 0 1 0; 1 0.5 0];
algo_type = char('Single','Percent Difference', 'Std. Dev.', 'Longest Salesman Route', 'Longest Salesman and Total Route', ...
                 'Single Average', 'Percent Difference Ave.', 'Std. Dev. Ave.', 'Longest Salesman Route Ave.', 'Longest Salesman and Total Route Ave.');


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
        
        std_dev(j,i,1) = std(smd_td);
        std_dev(j,i,2) = std(smd_pd);
        std_dev(j,i,3) = std(smd_std);
        std_dev(j,i,4) = std(smd_t);
        std_dev(j,i,5) = std(smd_tnt);
        
        %molly's old code, its hard to create a legend for, also weird
%         %failure for  increasing test_iter
%         for k = 1:num_algos
%             plot(j*ones(1,test_iter),time(j,:,k),'.','Color',clr(k,:));
%             hold on;
%         end
        
    end
    
    %% MEGA PLOT Time
%             t = j*ones(1,test_iter); %city visited
            t_a1 = time(j,:, 1);
            t_a2 = time(j,:, 2);
            t_a3 = time(j,:, 3);
            t_a4 = time(j,:, 4);
            t_a5 = time(j,:, 5);
            
            %at = linspace(1,map_iter
            ave_time(j,1) = mean(t_a1);
            ave_time(j,2) = mean(t_a2);
            ave_time(j,3) = mean(t_a3);
            ave_time(j,4) = mean(t_a4);
            ave_time(j,5) = mean(t_a5);
            
            %effin colors won't work, no idea why
%             plot (t,t_a1,'mo',...
%                  'MarkerSize', 3,...
%                  'MarkerFaceColor', 'm');
%             hold on;
%             plot (t,t_a2,'ro', ...
%                  'MarkerSize', 3,...
%                  'MarkerFaceColor', 'r');
%             hold on;
%             plot (t,t_a3,'go',...'MarkerEdgeColor',[0 0 0], ...
%                  'MarkerSize', 3,...
%                  'MarkerFaceColor', 'g');
%             hold on;
%             plot (t,t_a4,'bo', ...'Color',clr(4,:), ...
%                  'MarkerSize', 3,...
%                  'MarkerFaceColor', 'b');
%             hold on;
%             plot (t,t_a5,'ko',... % 'Color',clr(5,:));
%                  'MarkerSize', 3,...
%                  'MarkerFaceColor', 'k');
%             set(gca, 'XLim', [0 map_iter+1], 'XTick' , 1:map_iter);
%             xlabel('Map Number');
%             ylabel('Time');
%             
%             title('Time For Each Algorithms solution')
%             hold on;
            %plot the averages as asterisks
            
%             plot (j,ave_time(j,1), '-ms','MarkerSize', aMS,'LineWidth',aLW);
%             plot (j+0.08*sp,ave_time(j,2), '-rs','MarkerSize', aMS,'LineWidth',aLW);
%             plot (j-0.08*sp,ave_time(j,3), '-gs','MarkerSize', aMS,'LineWidth',aLW);
%             plot (j+0.18*sp,ave_time(j,4), '-bs','MarkerSize', aMS,'LineWidth',aLW);
%             plot (j-0.18*sp,ave_time(j,5), '-ks','MarkerSize', aMS,'LineWidth',aLW);
%            legend(algo_type);
%             hold on;
            
            %% Mega Plot Std Dev
            std_a1 = std_dev(j,:, 1);
            std_a2 = std_dev(j,:, 2);
            std_a3 = std_dev(j,:, 3);
            std_a4 = std_dev(j,:, 4);
            std_a5 = std_dev(j,:, 5);
            
            ave_std(j,1) = mean(std_a1);
            ave_std(j,2) = mean(std_a2);
            ave_std(j,3) = mean(std_a3);
            ave_std(j,4) = mean(std_a4);
            ave_std(j,5) = mean(std_a5);
%             
%             %effin colors won't work, no idea why
%             plot (t,std_a1,'mo',...
%                  t,std_a2,'r.', ...
%                  t,std_a3,'g.',...'MarkerEdgeColor',[0 0 0], ...
%                  t,std_a4,'.b', ...'Color',clr(4,:), ...
%                  t,std_a5,'k.'); % 'Color',clr(5,:));
%             set(gca, 'XLim', [0 map_iter+1], 'XTick' , 1:map_iter);
%             xlabel('Map Number');
%             ylabel('Standar Deviation of Time');
%             legend(algo_type);
%             title('Standard Deviaton of Time for Each Algorithm')
%             hold on;
end
% maps = linspace(1,map_iter,map_iter);
% plot (maps,ave_time(:,1), '-ms','MarkerSize', aMS,'LineWidth',aLW);
% plot (maps+0.08*sp,ave_time(:,2), '-rs','MarkerSize', aMS,'LineWidth',aLW);
% plot (maps-0.08*sp,ave_time(:,3), '-gs','MarkerSize', aMS,'LineWidth',aLW);
% plot (maps+0.18*sp,ave_time(:,4), '-bs','MarkerSize', aMS,'LineWidth',aLW);
% plot (maps-0.18*sp,ave_time(:,5), '-ks','MarkerSize', aMS,'LineWidth',aLW);
% legend(algo_type);
% hold on;


%get lowest average time and std deviation overall
for b =1:num_algos
    total_ave_time(b) = mean(ave_time(:,b));
    total_ave_std(b) = mean(ave_std(:,b));
end