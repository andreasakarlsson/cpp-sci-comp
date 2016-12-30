%% Plot algebraic expressions

% Algebraic expressions of first and second partial derivatives in x

x = linspace(-10,5,50);
x0 = x;

figID = figure(30001);

ux = sin((x./10).^2).*cos(x./10);  % x-component of u(x,y)

dxu = 2/100.*x.*cos((x./10).^2).*cos(x./10) ... % first parital derivative
      - 1/10*sin((x./10).^2).*sin(x./10);

dxxu = 2/100*cos((x./10).^2).*cos(x./10) ... % second partial derivative
       - (2/100.*x).^2.*sin((x./10).^2).*cos(x./10) ...
       - 1/250*cos((x./10).^2).*x.*sin(x./10) ...
       - 1/100*sin((x./10).^2).*cos(x./10);

subplot(1,3,1); plot(x,ux,'linewidth',2); title('(a)');   
subplot(1,3,2); plot(x,dxu,'linewidth',2); title('(b)');
subplot(1,3,3); plot(x,dxxu,'linewidth',2); title('(c)');


% 
% set(figID,'Units','centimeters');
% pos = get(figID,'Position');
% set(figID,'PaperPositionMode','Auto','PaperUnits','centimeters',...
%     'PaperSize',[pos(3), pos(4)])
% 
% print(figID,'algebraicDeriv1','-dpdf')
% close(figID)


%% Load and plot data from bin files

fid = fopen('task1.bin','r');
c = fread(fid,'double');
fclose(fid);

x = c(1:length(c)/2);
y = c(length(c)/2+1:end);

%figure(101)
%plot(x,y,'.')
%axis([-12 7 -1 4])


fid = fopen('task3-4.bin','r');
u = fread(fid,'double');
fclose(fid);
A = vec2mat(u,50);
%figure(102)
%imagesc(A)
figure(103)
X = vec2mat(x,20)';
Y = vec2mat(y,20)';
Y(end:-1:1) = Y;
surf(X,Y,A)
figure(104)
plot(X(end-5,:),A(end-5,:))


%% Compare algebraic expressions and C++ result

figID = figure(105); 

algeb = dxxu;  %  choose dxu for 1st deriv and dxxu for 2nd

hold on
yyaxis left
plot(x0,squeeze(A(end-5,:)),'g','linewidth',2)
plot(x0,algeb,'r','linewidth',2)
ax = gca;
ax.YColor = [1.0 0.0 0.0];
%axis([-10.5 5.5 -0.04 0.023])
%xlim([-10.5 5.5])
  % Choose ylabel depending on 1st or 2nd deriv
%ylabel('$\frac{\partial}{ \partial x} u(x,y)$','Interpreter','latex','FontSize',14)
ylabel('$\Delta u(x,y)$','Interpreter','latex','FontSize',14)
yyaxis right

erd =(algeb-squeeze(A(end-1,:)));
ern = (algeb+squeeze(A(end-1,:)));
plot(x0,abs( erd ),'b','linewidth',2);
ax = gca;
ax.YColor = [0.0 0.0 1.0];
%ylim([-0.1e-4 4.5e-4])
ylabel('error','FontSize',14)



% set(figID,'Units','centimeters');
% pos = get(figID,'Position');
% set(figID,'PaperPositionMode','Auto','PaperUnits','centimeters',...
%     'PaperSize',[pos(3), pos(4)])
% 
% print(figID,'comparison-xx','-dpdf')
% close(figID)
