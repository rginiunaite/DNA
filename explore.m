load u_mutated_polyd_h.mat

load ../optimisation/data/refintervalsMinMax4.mat

% plot(D(1).u(2:3:end))
% plot(D(1).u(3:3:end))
% plot(D(1).u(4:3:end))

figure; 
hold on

for ii = 1:18
    
umax = zeros(147,1);
j = 0;

for k = 2:3:442
    
   %[ui,i] = max(abs(D(ii).u(k:k+2) - D(ii).u(1)));
   j = j+1; 
   %D(ii).S(i,:)
   %umax(j) = D(ii).u(k+i-1);
   umax(j) = max((D(ii).u(k:k+2) - D(ii).u(1)));
   %uu(ii,j) = max(abs(D(ii).u(k:k+2) - D(ii).u(1))); 
end    

plot(umax);
% [i j ] = max(umax);
% j

umaxall(ii,:)=umax;

end

for k = [ref.indc ref.indw]
    xline(k,'--k');
end

ylabel('Nucleosome wrapping energy')
xlabel('Position of a point mutation in the sequence')
