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
show_res = 1;

clr = [1 0 0; 0 0 1; 0.67 0 1; 0 1 0; 1 0.5 0];

%% Calling the different algorithms
[opt_rte_1, min_dist_1, smd_1, dist_history_1] =  mtsp_1(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);
[opt_rte_2, min_dist_2, smd_2, dist_history_2,total_dist, cost] =  mtsp_3(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);
[opt_rte_tnt, min_dist_tnt, smd_tnt, dist_history_tnt] =  mtsp_tnt(xy,dmat,salesmen,min_tour,max_tour,tw,pop_size,num_iter,use_complex,show_prog,show_res);

%% Calculating Time/Max Distance traveled by salesman
time(1) = max(smd_1);
time(2) = max(smd_2);

time

%% Plots

%Plot of City Locations
figure('Name','City Map','NumberTitle','off','Color','white')
title('City Locations');
plot(xy(:,1),xy(:,2),'k.')


%Plot comparing Max travel distance/time
figure('Name','Time of Travel','NumberTitle','off','Color','white')
bar(time);

%Plot of Salesmen Travel 1
figure('Name','Salesman Travel','NumberTitle','off','Color','white')
subplot(2,1,1);
for s = 1:salesmen
        rte = [1 opt_rte_1.ch{s} 1];
        plot(xy(rte,1),xy(rte,2),'.-','Color',clr(s,:));
        title(sprintf('Total Distance 1 = %1.4f',min_dist_1));
        hold on;
end

    
%Plot of Salesmen Travel 2
subplot(2,1,2)
for s = 1:salesmen
        rte = [1 opt_rte_2.ch{s} 1];
        plot(xy(rte,1),xy(rte,2),'.-','Color',clr(s,:));
        title(sprintf('Total Distance 2 = %1.4f',min_dist_2));
        hold on;
end

% Best Solution History Travel 1
figure('Name','Best Solution History','NumberTitle','off','Color','white')
    subplot(2,1,1);
    title('Best Solution History 1');
    plot(dist_history_1)
% Best Solution History Travel 2
    subplot(2,1,2);
    title('Best Solution History 2');
    plot(dist_history_2)