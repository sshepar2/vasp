figure;
load matbands.mat   % [band][kpoint]  energy eigenvalues
load matkpts.mat
load fermi.mat
load symlines.mat
load ispin.mat
load matpbands.mat  % [ion][orbital][band][kpoints] ; [s,py,pz,px,dxy,dyz,dz2,dxz,dx2-y2]  projected weights

% s(1): how many bands, p(2): how many kpoints
s = size(matbands);
p = size(matkpts);

% for no spin leave as spin = 1. For spin, toggle which spin to plot (1 or 2)
spin = 1;
projected = 'y';



hold;

% plotting band structure for spin 1 bands

if spin == 1
    for b=1:s(1)
        band = matbands(b,:);
        plot(matkpts(1:end-1),band(:)-fermi,'ko','MarkerFaceColor','k','MarkerSize',2);
    end
end


% plotting projected band contributions for spin 1 bands

%projection 1
if spin == 1 & projected == 'y'
    for b=1:s(1)
        band = matbands(b,:);
        for k=1:p(2)-1
            plot(matkpts(k),band(k)-fermi,'go','MarkerSize',20*matpbands(1,3,b,k)+0.00001,'MarkerFaceColor','g')
        end
    end
end
 
% if spin == 1 & projected == 'y'
%     for b=1:s(1)
%         band = matbands(b,:);
%         for k=1:p(2)
%             plot(matkpts(k),band(k)-fermi,'go','MarkerSize',40*((matpbands(19,3,b,k)+matpbands(20,3,b,k)+matpbands(21,3,b,k)+matpbands(22,3,b,k)+matpbands(23,3,b,k)+matpbands(24,3,b,k)+matpbands(25,3,b,k)+matpbands(26,3,b,k))/8)+0.00001,'MarkerFaceColor','g')
%         end
%     end
% end

% projection 2
% if spin == 1 & projected == 'y'
%     for b=1:s(1)
%         band = matbands(b,:);
%         for k=1:p(2)
%             plot(matkpts(k),band(k)-fermi,'go','MarkerSize',40*((matpbands(1,3,b,k)+matpbands(2,3,b,k)+matpbands(3,3,b,k)+matpbands(4,3,b,k)+matpbands(5,3,b,k)+matpbands(6,3,b,k)+matpbands(7,3,b,k)+matpbands(8,3,b,k))/8)+0.00001,'MarkerFaceColor','g')
%         end
%     end
% end
 
%projection 3
if spin == 1 & projected == 'y'
    for b=1:s(1)
        band = matbands(b,:);
        for k=1:p(2)-1
            plot(matkpts(k),band(k)-fermi,'bo','MarkerSize',20*matpbands(1,1,b,k)+0.00001,'MarkerFaceColor','b')   
        end
    end
end

%%projection 4
% if spin == 1 & projected == 'y'
%     for b=1:s(1)
%         band = matbands(b,:);
%         for k=1:p(2)
%             plot(matkpts(k),band(k)-fermi,'ro','MarkerSize',20*((matpbands(7,1,b,k)+matpbands(8,1,b,k))/2)+0.00001,'MarkerFaceColor','r')   
%         end
%     end
% end

% spin bands (spin2)

if ispin == 2 & spin == 2
    load matbands1.mat
    load matpbands1.mat
    s1 = size(matbands1);
    for b=1:s1(1);
        band1 = matbands1(b,:);
        plot(matkpts(:),band1(:)-fermi,'ko','MarkerFaceColor','k','Markersize',2);
    end
end


% spin projected bands (spin2)  See above where matpbands is loaded to
% determine how to specify ion/orbital

if ispin == 2 & spin == 2 & projected == 'y'
    for b=1:s(1)
        band1 = matbands1(b,:);
        for k=1:p(2)
            plot(matkpts(k),band1(k)-fermi,'go','MarkerSize',20*((matpbands1(7,3,b,k)+matpbands1(8,3,b,k))/2)+0.00001,'MarkerFaceColor','g')
        end
    end
end











% finding maximum and minimum eigenvalues for plotting range
ymax = max(max(matbands));
ymin = min(min(matbands));

% plotting vertical lines at high symmetry points
for l=1:length(symlines);
    plot([symlines(l) symlines(l)],[ymin-10 ymax+10],'k','Linewidth',1)
end
hold;

%plotting fermi level
hold;
plot([matkpts(1) matkpts(end)],[0 0],'--k')


% plotting details
axis([matkpts(1) matkpts(end) -10 12])
%xticks([symlines(1) symlines(2) symlines(3) symlines(4)])
%xticklabels({'-B','\Gamma','B','AB','A','A-B','-B','-A-B','\Gamma','AB','A','\Gamma','-A'}) % enter high symmetry points
%xticklabels({'\Gamma','M','K','\Gamma'})

%yticks([-2 -1 0 1 2])
set(gca,'YMinorTick','on','Fontsize',30)
%yticks([-25 -20 -15 -10 -5 0 5 10 15])
ylabel('Energy (eV)')

%set(gca,'YMinorTick','on','Fontsize',30)
pbaspect([2 1 1])
box on

%axis([matkpts(1) matkpts(end) -24 2])
xlabel('wavenumber (k)','Fontsize',30)

