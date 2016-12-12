figure(3)
subplot(2,2,1)
for i=[1 2 3 7 8 10]
    semilogy(set{i,7}(:,2),set{i,7}(:,1))
    xlim([0 150]);
    hold on
end
hold off
subplot(2,2,2)
for i=[1 2 3 7 8 10]
    normT = (set{i,1}(:,1)-min(set{i,1}(:,1)))/(max(set{i,1}(:,1))-min(set{i,1}(:,1)));
    plot(set{i,1}(:,2),normT)
    xlim([0 150]);
    hold on
end
for i=[1 2 3 7 8 10]
    normT = (set{i,2}(:,1)-min(set{i,2}(:,1)))/(max(set{i,2}(:,1))-min(set{i,2}(:,1)));
    plot(set{i,2}(:,2),normT,'LineStyle', '--')
    xlim([0 150]);
    hold on
end
hold off
subplot(2,2,3)
for i=[1 2 3 7 8 10]
    plot(set{i,3}(:,2),set{i,3}(:,1))
    xlim([0 150]);
    hold on
end
for i=[1 2 3 7 8 10]
    plot(set{i,4}(:,2),set{i,4}(:,1), 'LineStyle','--')
    xlim([0 150]);
    hold on
end
hold off

subplot(2,2,4)
for i=[1 2 3 7 8 10]
    plot(set{i,5}(:,2),set{i,5}(:,1))
    xlim([0 150]);
    hold on
end
for i=[1 2 3 7 8 10]
    plot(set{i,6}(:,2),set{i,6}(:,1), 'LineStyle','--')
    xlim([0 150]);
    hold on
end
hold off
