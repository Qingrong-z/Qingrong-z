function ratio=mu_T_ratio(T, mu, mu_ref,E,T_L,m_e,m_h,Eg,D)

% Code that computes the ratio of PL spectra accounting for carrier heating
% and band filling
%
%@input: T              the carrier temperature at high-intensity PL (in K)
%@input: mu             the chemical potential at high-intensity PL (in eV)
%@input: mu_ref         the chemical potential at reference (low) intensity PL (in eV)
%@input: E              the energy vector (in eV)
%@input: T_L            the lattice temperature (in K)
%@input: m_e            the effective mass of electrons (as a fraction of the mass of free electrons)
%@input: m_h            the effective mass of holes (as a fraction of the mass of free electrons)
%@input: Eg             the band gap of the absorber material
%@input: D              the dimensionality of the absorber (3 for 3D etc.)
%
%@output:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% k = 1.381e-23; %Boltzmann's constant
% q=1.602e-19; % Elementary charge
% 
% BF=sinh(q*(E-mu)/(2*k*T))./(cosh(q*(E-mu)/(2*k*T))+cosh(((m_h-m_e)/(m_h+m_e))*(q*(E-Eg)/(2*k*T))-(D/4)*log(m_h/m_e)));
% BF_ref=sinh(q*(E-mu_ref)/(2*k*T_L))./(cosh(q*(E-mu_ref)/(2*k*T_L))+cosh(((m_h-m_e)/(m_h+m_e))*(q*(E-Eg)/(2*k*T_L))-(D/4)*log(m_h/m_e)));
% 
% ratio=(BF./BF_ref).*(exp(q*(E-mu_ref)/(k*T_L))-1)./(exp(q*(E-mu)/(k*T))-1);
k = 1.381e-23; %Boltzmann's constant
q=1.602e-19; % Elementary charge

BF=sinh(q.*(E-mu)./(2.*k.*T))./(cosh(q.*(E-mu)./(2.*k.*T))+cosh(((m_h-m_e)./(m_h+m_e)).*(q.*(E-Eg)./(2.*k.*T))-(D./4).*log(m_h./m_e)));
BF_ref=sinh(q.*(E-mu_ref)./(2.*k.*T_L))./(cosh(q.*(E-mu_ref)./(2.*k.*T_L))+cosh(((m_h-m_e)./(m_h+m_e)).*(q.*(E-Eg)./(2.*k.*T_L))-(D./4).*log(m_h./m_e)));

ratio=(BF./BF_ref).*(exp(q.*(E-mu_ref)./(k.*T_L))-1)./(exp(q.*(E-mu)./(k.*T))-1);