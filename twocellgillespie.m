% 
clear all

for i=1:10
    
    
    tot_t = 250;
    N = 4;
    R = 5;
    t = 0;

    beta = 1e-1/2 4;
    k = 1/4.2;
    p = 240;
    delta = 1/2.9;
    c = 1/29;
    n = 2;
    m = 2;
    f = 0.5;
    r_b = 1e-3;
    r_p = 1e-3;  % r_p takes [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1] 

    T = 500;
    f = 0.5;
    T_c = T*f;
    T_n = T*(1-f);   
    E_c = 0;
    E_n = 0;
    I_c = 0;
    I_n = 0;
    V = 10;

    % Initilalizing counters
    
    T_c_count = [T_c t];
    E_c_count = [E_c t];
    I_c_count = [I_c t];
    T_n_count = [T_n t];
    E_n_count = [E_n t];
    I_n_count = [I_n t];
    V_count = [V t];
    a = zeros(1,9);

    %%  propensities

    a(1)    =   beta*T_c*V;
    a(2)    =   r_b*beta*T_n*V;
    a(3)    =   k*E_c;
    a(4)    =   k*E_n;
    a(5)    =   p*I_c;
    a(6)    =   r_p*p*I_n;
    a(7)    =   delta*I_c;
    a(8)    =   delta*I_n;
    a(9)    =   c*V;
    prop    =   a;
    A = sum(a);

   %% SSA

    while (t<=tot_t)

        r = rand(1,2);
        tau = (1/A)*log(1/r(1));
        mu = r(2)*A;
       if (mu>0 && mu<=a(1))
           T_c = T_c-1;
           E_c = E_c+1;
           a(1)    =   beta*T_c*V;
           a(3)    =   k*E_c;
       end
        if (mu>a(1) && mu<=sum(a(1:2)))
           T_n = T_n-1;
           E_n = E_n+1;
           a(2)    =   r_b*beta*T_n*V;
           a(4)    =   k*E_n;
       end
       if (mu>sum(a(1:2)) && mu<=sum(a(1:3)))
           E_c = E_c-1;
           I_c = I_c+1;
           a(3)    =   k*E_c;
           a(5)    =   p*I_c; 
           a(7)    =   delta*I_c;
       end

       if (mu>sum(a(1:3)) && mu<=sum(a(1:4)))
           E_n = E_n-1;
           I_n = I_n+1;
           a(4)    =   k*E_n;
           a(6)    =   r_p*p*I_n; 
           a(8)    =   delta*I_n;
       end
       if (mu>sum(a(1:4)) && mu<=sum(a(1:5)))
           V = V+n;
           a(1)    =   beta*T_c*V;
           a(2)    =   r_b*beta*T_n*V;
           a(5)    =   c*V;

       end

       if (mu>sum(a(1:5)) && mu<=sum(a(1:6)))
           V = V+n;
           a(1)    =   beta*T_c*V;
           a(2)    =   r_b*beta*T_n*V;
           a(5)    =   c*V;
       end

       if (mu>sum(a(1:6)) && mu<=sum(a(1:7)))
           I_c = I_c-1;
           a(5)    =   p*I_c;
           a(7)    =   delta*I_c;
       end
          if (mu>sum(a(1:7)) && mu<=sum(a(1:8)))
           I_n = I_n-1;
           a(6)    =   r_p*p*I_n;
           a(8)    =   delta*I_n;
       end
       if (mu>sum(a(1:8)) && mu<=sum(a(1:9)))
           V = V-m;
           a(1)    =   beta*T_c*V;
           a(2)    =   r_b*beta*T_c*V;
           a(9)    =   c*V;
       end

       A = sum(a);

       t = t+tau;

           T_c_count = [T_c_count; T_c t];
           E_c_count = [E_c_count; E_c t];
           I_c_count = [I_c_count; I_c t];
           T_n_count = [T_n_count; T_n t];
           E_n_count = [E_n_count; E_n t];
           I_n_count = [I_n_count; I_n t];
           V_count = [V_count; V t];

           prop = [prop; a];

    end
    set(i,:) = {T_c_count,T_n_count, E_c_count,E_n_count,I_c_count,...
        I_n_count,V_count, prop};
end



