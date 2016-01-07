within ;
model Woodpecker
  Real z(start=5);
  Real phi_s(start=0);
  Real phi_b(start=0);
  Real lambda1(start=0);
  Real lambda2(start=0);
  Real z_prime(start=0);
  Real phi_s_prime(start=1);
  Real phi_b_prime(start=1.0);
  Real z_bis;
  Real phi_s_bis;
  Real phi_b_bis;
  Integer state(start=1);
  Integer hits(start=0);
  //Real lambda1prev(start=0);
  Real I;

  constant Real m_s=3e-4;
  constant Real J_s=5e-9;
  constant Real m_b=4.5e-3;
  constant Real J_b=7e-7;
  constant Real r_0=2.5e-3;
  constant Real r_s=3.1e-3;
  constant Real h_s=2e-2;
  constant Real l_s=1e-2;
  constant Real l_g=1.5e-2;
  constant Real l_b=2.01e-2;
  constant Real h_b=2e-2;
  constant Real c_p=5.6e-3;
  constant Real g=9.81;
equation
  z_prime=der(z);
  phi_s_prime=der(phi_s);
  phi_b_prime=der(phi_b);
  z_bis=der(z_prime);
  phi_s_bis=der(phi_s_prime);
  phi_b_bis=der(phi_b_prime);

  if state == 1 then
    (m_s+m_b)*z_bis+m_b*l_s*phi_s_bis+m_b*l_g*phi_b_bis = -(m_s+m_b)*g;
    (m_b*l_s)*z_bis+(J_s+m_b*l_s*l_s)*phi_s_bis+(m_b*l_s*l_g)*phi_b_bis = c_p*(phi_b-phi_s)-m_b*l_s*g-lambda1;
    m_b*l_g*z_bis+(m_b*l_s*l_g)*phi_s_bis+(J_b+m_b*l_g*l_g)*phi_b_bis = c_p*(phi_s-phi_b)-m_b*l_g*g-lambda2;
    0 = lambda1;
    0 = lambda2;
  elseif state == 2 then
    (m_s+m_b)*z_bis+m_b*l_s*phi_s_bis+m_b*l_g*phi_b_bis = -(m_s+m_b)*g-lambda2;
    (m_b*l_s)*z_bis+(J_s+m_b*l_s*l_s)*phi_s_bis+(m_b*l_s*l_g)*phi_b_bis = c_p*(phi_b-phi_s)-m_b*l_s*g-h_s*lambda1-r_s*lambda2;
    m_b*l_g*z_bis+(m_b*l_s*l_g)*phi_s_bis+(J_b+m_b*l_g*l_g)*phi_b_bis = c_p*(phi_s-phi_b)-m_b*l_g*g;
    //0 = (r_s-r_0)+h_s*phi_s;
    0 = phi_s_bis;
    //0 = z_prime+r_s*phi_s_prime;
    0 = z_bis+r_s*phi_s_bis;
  elseif state == 3 then
    (m_s+m_b)*z_bis+m_b*l_s*phi_s_bis+m_b*l_g*phi_b_bis = -(m_s+m_b)*g-lambda2;
    (m_b*l_s)*z_bis+(J_s+m_b*l_s*l_s)*phi_s_bis+(m_b*l_s*l_g)*phi_b_bis = c_p*(phi_b-phi_s)-m_b*l_s*g+h_s*lambda1-r_s*lambda2;
    m_b*l_g*z_bis+(m_b*l_s*l_g)*phi_s_bis+(J_b+m_b*l_g*l_g)*phi_b_bis = c_p*(phi_s-phi_b)-m_b*l_g*g;
    //0 = (r_s-r_0)-h_s*phi_s;
    0 = phi_s_bis;
    //0 = z_prime+r_s*phi_s_prime;
    0 = z_bis+r_s*phi_s_bis;
    //when phi_b_prime > 0 and h_b*phi_b >= l_s + l_g - l_b - r_0 then
    //reinit(phi_b_prime,-pre(phi_b_prime));
     //end when;
  else
    (m_s+m_b)*z_bis+m_b*l_s*phi_s_bis+m_b*l_g*phi_b_bis = -(m_s+m_b)*g-lambda2;
    (m_b*l_s)*z_bis+(J_s+m_b*l_s*l_s)*phi_s_bis+(m_b*l_s*l_g)*phi_b_bis = c_p*(phi_b-phi_s)-m_b*l_s*g+h_s*lambda1-r_s*lambda2;
    m_b*l_g*z_bis+(m_b*l_s*l_g)*phi_s_bis+(J_b+m_b*l_g*l_g)*phi_b_bis = c_p*(phi_s-phi_b)-m_b*l_g*g;
    //0 = (r_s-r_0)-h_s*phi_s;
    0 = phi_s_bis;
    //0 = z_prime+r_s*phi_s_prime;
    0 = z_bis+r_s*phi_s_bis;
  end if;

algorithm
  when state == 1 and phi_b_prime < 0 and h_s*phi_s <= -(r_s-r_0) then
    I :=m_b*l_g*z_prime + m_b*l_s*l_g*phi_s_prime + (J_b + m_b*l_g*l_g)*phi_b_prime;
    reinit(z_prime, 0);
    reinit(phi_s_prime,0);
    reinit(phi_b_prime,I/(J_b+m_b*l_g*l_g));
    reinit(lambda1, -1e-12);
    state :=2;
  end when;
  when state == 1 and phi_b_prime > 0 and h_s*phi_s >= (r_s-r_0) then
    I :=m_b*l_g*z_prime + m_b*l_s*l_g*phi_s_prime + (J_b + m_b*l_g*l_g)*phi_b_prime;
    reinit(z_prime, 0);
    reinit(phi_s_prime,0);
    reinit(phi_b_prime,I/(J_b+m_b*l_g*l_g));
    state :=3;
    hits := hits + 1;
    //inte omöjligt att det skulle kunna fungera att räkna hits här
  end when;
  when state == 2 and lambda1 >= 0 then
    state :=1;
  end when;
  when state == 3 and phi_b_prime < 0 and lambda1 >= 0 then
    state :=1;
  end when;
  when state == 3 and phi_b_prime > 0 and h_b*phi_b >= l_s + l_g - l_b - r_0 then
    hits := hits + 1;
    state :=4;
  end when;
  when state == 4 then
    reinit(phi_b_prime,-pre(phi_b_prime));
    state :=3;
  end when;

  annotation (uses(Modelica(version="3.2.1")));
end Woodpecker;
