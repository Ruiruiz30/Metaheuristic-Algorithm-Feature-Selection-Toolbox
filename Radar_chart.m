
data = randi([2,8],[4,7])+rand([4,7]);

h = radarChart(data);
h = h.draw();

labels = {'A','B','C','D','E'};
legend(labels,'Location','northeast','FontSize',12);
for n=1:h.ClassNum
    h.setPatchN(n,'Color',colorList(n,:),'MarkerFaceColor',colorList(n,:))
end