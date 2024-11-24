function GROUPED_BOXPLOT(x,group,legend_str)
n_OBJ=group(end)/size(legend_str,2);

positions = []; n_shift=0;
for i=1:max(group)/n_OBJ
    temp = [n_OBJ*(i-1)+1:n_OBJ*i] + n_shift;
    positions = [positions temp];
    n_shift = n_shift + 1;
end
boxplot(x,group, 'positions', positions);

mu_pos=[];  color=[];
for i=1:size(legend_str,2)
    mu_pos(1,i)=mean(positions(n_OBJ*(i-1)+1:n_OBJ*i));
    if n_OBJ == 5
        color = [color 'b', 'r', 'y', 'c', 'm'];
    else
        color = [color 'b', 'r', 'y', 'c', 'm', 'g'];
    end
    
end

set(gca,'xtick',mu_pos)
set(gca,'xticklabel',legend_str)
color=[color color];

h = findobj(gca,'Tag','Box');
for j=1:length(h)
%    patch(get(h(j),'XData'),get(h(j),'YData'),color(j));
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.55);
end

end