aMS = 3; %average Marker Size
aLW = 1.4; %average line width
sp = 4; %spacing factor for averages

maps = linspace(1,map_iter,map_iter);
plot (maps,ave_time(:,1), '-m','LineWidth',aLW);
hold on
plot (maps,ave_time(:,2), '-r','LineWidth',aLW);
hold on
plot (maps,ave_time(:,3), 'b-','LineWidth',aLW);
hold on
plot (maps,ave_time(:,4), '-g','LineWidth',aLW);
hold on
plot (maps,ave_time(:,5), '-k','LineWidth',aLW);
hold on

legend(algo_type);
hold on;

for j = 1:map_iter
    t = j*ones(1,test_iter); %city visited
    t_a1 = time(j,:, 1);
    t_a2 = time(j,:, 2);
    t_a3 = time(j,:, 3);
    t_a4 = time(j,:, 4);
    t_a5 = time(j,:, 5);
    
    plot (t,t_a1,'mo',...
         'MarkerSize', 3,...
         'MarkerFaceColor', 'm');
    hold on;
    plot (t-0.03*sp,t_a2,'ro', ...
         'MarkerSize', 3,...
         'MarkerFaceColor', 'r');
    hold on;
    plot (t-0.01*sp,t_a3,'bo',...'MarkerEdgeColor',[0 0 0], ...
         'MarkerSize', 3,...
         'MarkerFaceColor', 'b');
    hold on;
    plot (t+0.01*sp,t_a4,'go', ...'Color',clr(4,:), ...
         'MarkerSize', 3,...
         'MarkerFaceColor', 'g');
    hold on;
    plot (t+0.03*sp,t_a5,'ko',... % 'Color',clr(5,:));
         'MarkerSize', 3,...
         'MarkerFaceColor', 'k');
    set(gca, 'XLim', [0 map_iter+1], 'XTick' , 1:map_iter);
    xlabel('Map Number');
    ylabel('Time');

    title('Algorithm Solution Times')
    hold on;
end