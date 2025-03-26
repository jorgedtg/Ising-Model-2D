    clear;clc;close all;

tic
A = -1;

t = 1.5:0.01:3;

Q = length(t);

runs = 1:100;
nr = length(runs);

MNE = zeros(nr, Q);
MNM = zeros(nr, Q);
VRE = zeros(nr, Q);
VRM = zeros(nr, Q);

% wb = waitbar(0, 'Starting');

for r = runs

c = 1;
for T = t

N = 32;       
M = 10000;

STATES = zeros(N, N, M);
ENERGIES = zeros(1, M);
MAGNETS = ENERGIES;

S = sign(rand(N)-0.5);

for m = 1:M
    for n = 1:N^2

        SN = S;

        ni = randi([1, N]);
        nj = randi([1, N]);

        SN(ni, nj) = -SN(ni, nj);

        dE = denergy(S, SN, ni, nj, A);

        if dE < 0
            S = SN;
        else
            P = exp(-dE/T);
            rn = rand();

            if rn < P
                S = SN;
            else
                continue
            end
        end
    end
    surf(1:N, 1:N, S)
    view(0, 90)
    pause(1/144);

    % STATES(:, :, m) = S;
    ENERGIES(m) = energy(S, A);
    MAGNETS(m) = magnetization(S);

end

MNE(r, c) = mean(ENERGIES(1000:M));
MNM(r, c) = mean(MAGNETS(1000:M));
VRE(r, c) = var(ENERGIES(1000:M));
VRM(r, c) = var(MAGNETS(1000:M));
c = c+1;
% T
% toc
end
% waitbar(r/nr, wb, 'Progress: '+string(floor(r/nr*100)) + ' $\%$');
% r

end
toc
%%
MNE2 = mean(MNE(1:50, :));
VRE2 = mean(VRE(1:50, :));
%%
MNE3 = zeros(1, length(MNE2));

for i = 2:length(MNE3)-1

    MNE3(i) = (1/5)*(2*MNE2(i-1) + MNE2(i) + 2*MNE2(i+1));

end
MNE3(1) = MNE2(1);
MNE3(end) = MNE2(end);
%%
figure;
plot(t, MNE3)
ylabel("$\langle E \rangle$", Rotation=0);
xlabel("T")
axis square

figure;
plot(t, MNE2)
ylabel("$\langle E \rangle$", Rotation=0);
xlabel("T")
axis square

figure;
plot(t, VRE2)
ylabel("$\sigma_E^2$", Rotation=0);
xlabel("T")
axis square

%%
dMNE = derOrd4(MNE3, t);
figure;
plot(t, dMNE)
ylabel("$\frac{d\langle E\rangle}{dT}$", Rotation=0, FontSize=30)
xlabel("T")
text(2.29-0.1, dMNE(t==2.29)+0.05, "T = 2.29", FontSize=18)
ylim([0.05, 2.2])
axis square
%%
VREN = VRE2/sum(VRE2*0.01);
dMNEN = abs(dMNE)/sum(abs(dMNE)*0.01);

figure;
plot(t, VREN, t, abs(dMNEN))
set(gca,'XTick',[],'YTick',[])
legend("$\sigma_E^2$", "$\frac{d\langle E\rangle}{dT}$")
axis square


%%
MNM2 = mean(MNM(1:50, :));
VRM2 = mean(VRM(1:50, :));
%
MNM3 = zeros(1, length(MNM2));

for i = 2:length(MNM2)-1

    MNM3(i) = (1/5)*(MNM2(i-1) + 3*MNM2(i) + MNM2(i+1));

end
MNM3(1) = MNM2(1);
MNM3(end) = MNM2(end);

figure;
plot(t, MNM3)
ylabel("$\langle M \rangle$", Rotation=0);
xlabel("T")
axis square

figure;
plot(t, VRM2)
ylabel("$\sigma_M^2$", Rotation=0);
xlabel("T")
axis square
%%
dMNM = derOrd4(MNM2, t);

figure;
plot(t, abs(dMNM'))
ylabel("$\frac{d\langle M\rangle}{dT}$", Rotation=0, FontSize=30)
xlabel("T")
% text(2.3-0.07, abs(dMNM(t==2.33))+0.1, "T = 2.33", FontSize=18)
% ylim([0, 2])
axis square

%%
VRMN = VRM2/sum(VRM2*0.01);
dMNMN = abs(dMNM)/sum(abs(dMNM)*0.01);

figure;
plot(t, VRMN, t, abs(dMNMN))
axis square
set(gca,'XTick',[],'YTick',[])
legend("$\sigma_M^2$", "$\frac{d\langle M\rangle}{dT}$")


%%
figure;
plot(t, MNE2)
ylabel("$\langle E \rangle$", Rotation=0);
xlabel("T")
axis square

figure;
plot(t, MNM2)
ylabel("$\langle M \rangle$", Rotation=0)
xlabel("T")
axis square

%%
figure;
plot(t, MNE)
ylabel("$\langle E \rangle$", Rotation=0);
xlabel("T")
axis square

figure;
plot(t, MNM)
ylabel("$\langle M \rangle$", Rotation=0)
xlabel("T")
axis square
%%
figure;
plot(t, VRE2)
yl = ylabel("$\sigma_E^2$", Rotation=0);
xlabel("T")
axis square



figure;
plot(t, VRM2)
ylabel("$\sigma_M^2$", Rotation=0)
xlabel("T")
axis square
%%
dMNE = derOrd4(MNE2, t);

figure;
plot(t, dMNE)
ylabel("$\frac{d\langle E\rangle}{dT}$", Rotation=0, FontSize=30)
xlabel("T")
text(2.3-0.1, dMNE(t==2.3)+0.06, "T = 2.3", FontSize=18)
axis square

dMNM = derOrd4(MNM2, t);

figure;
plot(t, dMNM)
ylabel("$\frac{d\langle M\rangle}{dT}$", Rotation=0, FontSize=30)
xlabel("T")
text(2.3-0.09, dMNM(t==2.3)-0.06, "T = 2.3", FontSize=18)
axis square
%%
dMNE(1) = 0.1;
dMNM(1) = -0.01;
%%
VRE2N = VRE2/sum(VRE2*0.01);
VRM2N = VRE2/sum(VRM2*0.01);
dMNEN = dMNE/sum(dMNE*0.01);
dMNMN = dMNM/sum(dMNM*0.01);


figure;
plot(t, VRE2N, t, abs(dMNEN))
xlabel("T")
text(2.3-0.09, dMNM(t==2.3)-0.06, "T = 2.3", FontSize=18)
axis square
%%
figure;
plot(t, 5*VRM2N, t, abs(dMNMN))
xlabel("T")
text(2.3-0.09, dMNM(t==2.3)-0.06, "T = 2.3", FontSize=18)
axis square



%%

figure;
bh = bar(MNE, 'FaceColor','flat', 'EdgeColor','none');  
set(gca,'xtick', 1:3, 'xticklabel', {'T=1.5', 'T=2.3', 'T=3'});
text(1:length(MNE),(MNE-0.11),num2str(MNE'),'vert','bottom','horiz','center', fontsize=14); 
ylim([-2.2, 0])
axis square
bh.CData(1, :) = [0.3010 0.7450 0.9330];   %blue
bh.CData(2, :) = [0.4660 0.6740 0.1880]; %light blue
bh.CData(3, :) = [0.6350 0.0780 0.1840]; %pink
%%
figure;
bh1=bar(MNM, 'FaceColor','flat', 'EdgeColor','none');
set(gca,'xtick', 1:3, 'xticklabel', {'T=1.5', 'T=2.3', 'T=3'});
text(1:length(MNM),MNM,num2str(MNM'),'vert','bottom','horiz','center', fontsize=14); 
ylim([0, 1.1])
axis square
bh1.CData(1, :) = [0.3010 0.7450 0.9330];   %blue
bh1.CData(2, :) = [0.4660 0.6740 0.1880]; %light blue
bh1.CData(3, :) = [0.6350 0.0780 0.1840]; %pink
%%
figure;
bh2=bar(VRE, 'FaceColor','flat', 'EdgeColor','none');
set(gca,'xtick', 1:3, 'xticklabel', {'T=1.5', 'T=2.3', 'T=3'});
text(1:length(VRE),VRE,num2str(VRE'),'vert','bottom','horiz','center', fontsize=14); 
 ylim([0, 0.012])
axis square
bh2.CData(1, :) = [0.3010 0.7450 0.9330];   %blue
bh2.CData(2, :) = [0.4660 0.6740 0.1880]; %light blue
bh2.CData(3, :) = [0.6350 0.0780 0.1840]; %pink
%%
figure; 
bh3 = bar(VRM, 'FaceColor','flat', 'EdgeColor','none');
set(gca,'xtick', 1:3, 'xticklabel', {'T=1.5', 'T=2.3', 'T=3'});
text(1:length(VRM),VRM,num2str(VRM'),'vert','bottom','horiz','center', fontsize=14); 
ylim([0, 0.016])
axis square
bh3.CData(1, :) = [0.3010 0.7450 0.9330];   %blue
bh3.CData(2, :) = [0.4660 0.6740 0.1880]; %light blue
bh3.CData(3, :) = [0.6350 0.0780 0.1840]; %pink



% dE = denergy(S, S1, 1, 1, A);

function E = energy(S, A)

E = 0;
ns = length(S);

for i = 1:ns
    for j = 1:ns
        if i == ns
            im = 1;
        else
            im = i + 1;
        end
        E = E + S(i, j)*S(im, j);
    end
end

for i = 1:ns
    for j = 1:ns
        if j == ns
            jm = 1;
        else
            jm = j + 1;
        end
        E = E + S(i, j)*S(i, jm);
    end
end

E = A*E/ns^2;

end

function M = magnetization(S)

M = abs(sum(S, "all"))/length(S)^2;

end

function dE = denergy(S, S1, i, j, A)

ns = length(S);

if i == ns
    ip = 1;
    im = ns-1;
elseif i == 1
    im = ns;
    ip = 2;
else
    ip = i+1;
    im = i-1;
end

if j == ns
    jp = 1;
    jm = ns-1;
elseif j == 1
    jm = ns;
    jp = 2;
else
    jp = j+1;
    jm = j-1;
end

ES = (S(i, j)*S(ip, j) + S(i, j)*S(im, j) + S(i, j)*S(i, jp) + S(i, j)*S(i, jm));
ES1 = (S1(i, j)*S1(ip, j) + S1(i, j)*S1(im, j) + S1(i, j)*S1(i, jp) + S1(i, j)*S1(i, jm));


dE = A*(ES1-ES);


end

function dV = derOrd4(V, t)

N = length(V);
dx = abs(t(1)-t(2));

D = zeros(N, N);

val = [1/12, -2/3, 0, 2/3, -1/12];
% val = [-1/2, 0, 1/2];

nv = length(val);
nm = floor(nv/2);

for i = -nm:nm

    D = D + diag(ones(N-abs(i), 1)*val(i+nm+1), i);

end

r1 = zeros(1, length(D));
rf = r1;

r1(1) = -25/12;
r1(2) = 4;
r1(3) = -3;
r1(4) = 4/3;
r1(5) = -1/4;

rf(end) = 25/12;
rf(end-1) = -4;
rf(end-2) = 3;
rf(end-3) = -4/3;
rf(end-4) = 1/4;

r2 = [0,    r1];
r2(end) = [];

rf1 = [rf, 0];
rf1(1) = [];

D(1, :) = r1;
D(2, :) = r2;
D(end, :) = rf;
D(end-1, :) = rf1;

dV = (D*V')/dx;

end


