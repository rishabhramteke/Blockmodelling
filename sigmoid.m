function sigmoid(gamma, beta_mu, tau)

    vX = 0.0:0.001:1;
    vSigmoid = 1 ./ (1 + gamma * exp(-beta_mu * (vX - tau)));
    
    vGauss = gaussmf(vX, [0.1 0.9]);
%     vGauss = gauss2mf(vX, [0.1 0.5]);
    
    hold on;
    plot(vX, vSigmoid, vX, vX, vX, abs(vSigmoid - vX));
%     plot(vX, vSigmoid);

    
    
%     plot(vX, vGauss);
    
    hold off;
    
end