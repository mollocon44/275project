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
min_tour = 1;
max_tour = 100;
tw = 0; %time window
pop_size = 80; %size of population
num_iter = 1000; %number of itterations
use_complex = 0;
show_prog = 1;
show_prog_t = 1;
show_prog_t = logical(show_prog_t(1));
change_clrs = 1;
change_clrs = logical(change_clrs(1));
show_res = 1;

clr = [1 0 0; 0 0 1; 0.67 0 1; 0 1 0; 1 0.5 0];

%% Calling the different algorithms
%[opt_rte_1, min_dist_1, smd_1, dist_history_1] =  mtsp_1(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);
%[opt_rte_2, min_dist_2, smd_2, dist_history_2,total_dist, cost] =  mtsp_3(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);
%tour == longest tour length optimization -> still needs parameter correcting
[opt_rte_t, opt_out_t, soln_history_t, history_t, smd_t] = mtsp_tour(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);
%tnt == tour and total optimization - this works now!!!
%[opt_rte_tnt, opt_out_tnt, soln_history_tnt, history_tnt, smd_tnt] = mtsp_tnt(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);
%% Calculating Time/Max Distance traveled by salesman
% time(1) = max(smd_1);
% time(2) = max(smd_2);
% 
% time
%% Unpack Tour, Tour+Total
%this just saves me super disgusting inputs/outputs on functions that would
%take forever to debug

%Longest tour optimization

% opt_out is the final sol'n info
opt_ltour_t = opt_out_t(1);
opt_dist_t = opt_out_t(2);
opt_lt_i_t = opt_out_t(3);
opt_iter_t = opt_out_t(4);
%history is for animations if we want them
longest_tour_history_t = history_t{1};
dist_history_t = history_t{2};
lt_i_history_t = history_t{3};

%Tour+Total

% %final sol'n
% opt_ltour_tnt = opt_out_tnt(1);
% opt_dist_tnt = opt_out_tnt(2);
% opt_lt_i_tnt = opt_out_tnt(3);
% opt_iter_tnt = opt_out_tnt(4);
% %history for animations
% longest_tour_history_tnt = history_tnt{1};
% dist_history_tnt = history_tnt{2};
% lt_i_history_tnt = history_tnt{3};

%% Animation Plots to Show Class
if show_prog_t
    % Plot the Best Route
    for h = 1:num_iter
        opt_rte_h = soln_history_t{h};
        lt_i_t = lt_i_history_t{h};
        longest_tour_t = longest_tour_history_t{h};
        total_dist_t = dist_history_t{h};
        
        figure(pfig);
        for s = 1:salesmen
            rte = [1 opt_rte_h.ch{s} 1];
            %change_clrs determines whether to change colors as we're
            %going, which makes the display a little dizzying
            clrs = [clr_3;clr_2;clr_1];
            if change_clrs
                if lt_i_t(h) == 1
                    clrs = [clr_1;clr_2;clr_3];
                elseif lt_i_t(h) == 2
                    clrs = [clr_2;clr_1;clr_3];
                end
            end
            if dims == 3, 
                plot3(xy(rte,1),xy(rte,2),xy(rte,3),'.-','Color',clrs(s,:));
            else
                plot(xy(rte,1),xy(rte,2),'.-','Color',clrs(s,:));
            end
            title(sprintf('Longest Tour = %1.4f, Distance = %1.4f, Iteration = %d',longest_tour_t(h),total_dist_t(h),h));
            hold on
        end
        if dims == 3,
            plot3(xy(1,1),xy(1,2),xy(1,3),'ko');
        else
            plot(xy(1,1),xy(1,2),'ko'); 
        end
        if mod(num_iter,50) == 0 
            if num_iter > 50
                pause(0.1)
            end
        end
        hold off
    end
end

%% Plots

% %Plot of City Locations
% figure('Name','City Map','NumberTitle','off','Color','white')
% title('City Locations');
% plot(xy(:,1),xy(:,2),'k.')
% 
% 
% %Plot comparing Max travel distance/time
% figure('Name','Time of Travel','NumberTitle','off','Color','white')
% bar(time);
% 
% %Plot of Salesmen Travel 1
% figure('Name','Salesman Travel','NumberTitle','off','Color','white')
% subplot(2,1,1);
% for s = 1:salesmen
%         rte = [1 opt_rte_1.ch{s} 1];
%         plot(xy(rte,1),xy(rte,2),'.-','Color',clr(s,:));
%         title(sprintf('Total Distance 1 = %1.4f',min_dist_1));
%         hold on;
% end
% 
%     
% %Plot of Salesmen Travel 2
% subplot(2,1,2)
% for s = 1:salesmen
%         rte = [1 opt_rte_2.ch{s} 1];
%         plot(xy(rte,1),xy(rte,2),'.-','Color',clr(s,:));
%         title(sprintf('Total Distance 2 = %1.4f',min_dist_2));
%         hold on;
% end
% 
% % Best Solution History Travel 1
% figure('Name','Best Solution History','NumberTitle','off','Color','white')
%     subplot(2,1,1);
%     title('Best Solution History 1');
%     plot(dist_history_1)
% % Best Solution History Travel 2
%     subplot(2,1,2);
%     title('Best Solution History 2');
%     plot(dist_history_2)