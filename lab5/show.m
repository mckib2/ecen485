function show(n,N,e,in,out,theta_hat)
    % Plot phase error and signals
    subplot(2,1,1);
    plot(e(2:n),'k-');
    grid on;
    xlim([ 0 N ]);
    title('Phase Error');
    xlabel(sprintf('$ e(n) $ = %g',e(n)),'Interpreter','latex');
    ylabel('$ \theta_e(n) $','Interpreter','latex');
    
    subplot(2,1,2);
    plot(real(in(2:n)),'k-');
    grid on;
    hold on;
    plot(real(out(2:n)),'k--');
    xlabel([ '$ \hat{\theta} $ = ' num2str(theta_hat) ],'Interpreter','latex');
    xlim([ 0 N ]);
    ylim([ -1 1 ]);
    ylabel('$ Re\{e^{j(\cdot)}\} $','Interpreter','latex');
    hold off;
    
    legend('in','out');
    drawnow;
end