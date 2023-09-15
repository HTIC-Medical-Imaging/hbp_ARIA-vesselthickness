function plotcontours(fighandle,boundaryptslist,name)

figure(fighandle)
hold on
for k=1:length(boundaryptslist)
    boundpts = boundaryptslist{k};
    plot(boundpts(:,2),boundpts(:,1),'-','linewidth',2,'DisplayName', name + "-" + k)
end
hold off