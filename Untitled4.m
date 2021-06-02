% Temperature
for i=2:length(Int)
    T_raw(i) = 1/(1/T_raw(1)- k*a1{1,i}/q);
    T_plus_raw(i)=1/(1/T_raw(1)- k*a1_plus{1,i}/q);
    T_min_raw(i)=1/(1/T_raw(1)- k*a1_min{1,i}/q);
    delta_T_interval_raw(i)=max(abs(T_plus_raw(i)-T_raw(i)),abs(T_raw(i)-T_min_raw(i)));
    delta_T_slope_raw(i)=T_raw(i)^2*sqrt((delta_T_raw(1)/T_raw(1)^2)^2+(k/q)^2*delta_a1{1,i}^2);
    delta_T_raw(i)=delta_T_interval_raw(i)+delta_T_slope_raw(i);
    
    
    mu_raw(i) = T_raw(i)/T_raw(1)*mu_raw(1)+a0{1,i}*T_raw(i)*k;
    mu_plus_raw(i)=T_raw(i)/T_raw(1)*mu_raw(1)+a0_plus{1,i}*T_raw(i)*k;
    mu_min_raw(i)=T_raw(i)/T_raw(1)*mu_raw(1)+a0_min{1,i}*T_raw(i)*k;
    delta_mu_interval_raw(i)=max(abs(mu_plus_raw(i)-mu_raw(i)),abs(mu_raw(i)-mu_min_raw(i)));
    delta_mu_slope_raw(i)=mu_raw(i)^2*sqrt((delta_T_raw(i)/T_raw(i))^2+k^2*delta_a0{1,i}^2/(mu_ref/T_raw(1)+k*a0{1,i})^2);
    delta_mu_raw(i)=delta_mu_interval_raw(i)+delta_mu_slope_raw(i);
end
