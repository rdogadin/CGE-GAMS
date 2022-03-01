Sets
i                commodities /c1 agri, c2 mining, c3 mfg, c4 services/
h                households  /h1, h2, h3/
dmx              domestic imports exports /d, m, x/
dm(dmx)          domestic and imports /d, m/
f                factors     /l, k/
r                rates       /td, ts, tm, tl, tk, ty, ttr/
ti(r)            taxes       /td, ts, tm/
tp(r)            taxes transfers /ty, ttr/
tf(r)            factor taxes    /tl, tk/
iuva             intermediate use and value added /ciu, cva/
tx               export taxes /fob, t_exp/
;
alias (i,j);
alias (dm,idm);

parameter
param_IO_Table(*,*,i)
param_exports_Table(*,i)
param_factor_use_Table(*,i)
param_Tax_Rate_Table(*,*,i)
param_pfce_Table(*,*,h)
param_endowment_Table(*,h)
param_hh_tax_rates_Table(*,h)
param_sigma_dm_io_Table(*,i)
param_sigma_dm_fu_Table(*,h)
param_sigma_kl_Table(*,i)
param_sigma_iu_va_Table(*,i)
param_exports_Tax_Table(*,i)
;

*load data from gdx files
$include gdx_load.gms

*** GAMS: Parameters to Link Data
parameters
IO_Use_DM(i,dm,j)        intermediate use of D and M in production of commodity j
exports0(i)              exports of commodity j
exports_eta(i)           elasticity of exports
fob0(i)                  free on board duties costs etc.
t_exp0(i)                export subsidies taxes etc
factor_use(f,i)          factor use in industry i
tax_rates(r,dm,j)        tax rates on indirect use and factors in production of commodity j
pfce0(i,dm,h)            PFCE of commodity i by household h
endowment(f,h)           endowment by households
hh_tax_rates(tp,h)       personal tax rates for households
sigma_dm_io(i,j)         elasticity of substitution between d&m in io intermediate use
sigma_dm_fu(i,h)         elasticity of substitution between d&m in final use by household h
sigma_kl(i)              elasticity of substitution between l&k in value added
sigma_iu(i)              elasticity of substitution between composite IU
sigma_iu_va(i)           elasticity of substitution between IU and VA
;

*** GAMS: Linking Data to Parameters
IO_Use_DM(i,dm,j)        = param_IO_Table(i,dm,j);
exports0(i)              = param_exports_Table('exports',i);
exports_eta(i)           = param_exports_Table('eta_x',i);
fob0(i)                  = param_exports_Tax_Table('fob',i);
t_exp0(i)                = param_exports_Tax_Table('t_exp',i);
factor_use(f,i)          = param_factor_use_Table(f,i) ;
tax_rates(r,dm,j)        = param_Tax_Rate_Table(r,dm,j);
pfce0(i,dm,h)            = param_pfce_Table(i,dm,h);
endowment(f,h)           = param_endowment_Table(f,h);
hh_tax_rates(tp,h)       = param_hh_tax_rates_Table(tp,h);
sigma_dm_io(i,j)         = param_sigma_dm_io_Table(i,j);
sigma_dm_fu(i,h)         = param_sigma_dm_fu_Table(i,h);
sigma_kl(i)              = param_sigma_kl_Table('f',i);
sigma_iu(i)              = param_sigma_iu_va_Table('iu',i);
sigma_iu_va(i)           = param_sigma_iu_va_Table('iuva',i);
*display exports0;
*display param_sigma_kl_Table, sigma_kl, sigma_iu, sigma_iu_va;

***  GAMS Parameters for Calibration
parameters
IU_M(i,j)                intermediate use of imported inputs i in industry j
IU_D(i,j)                intermediate use of domestic inputs i in industry j
FU_M(i,h)                final use of imported goods by household h
FU_D(i,h)                final use of domestic goods by household h
factor_tax(f,i)          factor tax paid on use of factor f in industry i
rti(ti,i)                rates of excise sales and tariffs on use of commodity i
rtp(tp,h)                rates of income taxes and transfers
rtf(f,i)                 rate of factor taxes
tm_io(i,j)               total tax rate on use of imported goods i in sector j
td_io(i,j)               total tax rate on use of domestics goods i in sector j
tm_hh(i,h)               total tax rate on use of imported goods i on household h
td_hh(i,h)               total tax rate on use of domestics goods i on household h
pm0(i)                   domestic price of imported good i
pd0(i)                   domestic price of domestic good i
pf0(f)                   factor price
tiu_imp0(i)              taxes paid by industry i on intermediate imported inputs j
tiu_dom0(i)              taxes paid by industry i on intermediate domestic inputs j
tfu_imp0(i)              taxes paid by all households h on imported good i
tfu_dom0(i)              taxes paid by all households h on domestic good i
nit0(i)                  net indirect tax paid by industry i
factor_tax0(i)           factor tax paid by industry i
income_tax0(h)           income tax paid by household h
va0(i)                   value added in industry i
total_tax_revenue0       total tax revenue
gross_income0(h)         gross income of household h
transfer_hh0(h)          transfer to household h
disposable_income0(h)    disposable income of household h
gross_expenditure0(h)    gross expenditure of household h
;

***
IU_M(i,j)        = IO_Use_DM(i,'m',j);
IU_D(i,j)        = IO_Use_DM(i,'d',j);
FU_M(i,h)        = pfce0(i,'m',h);
FU_D(i,h)        = pfce0(i,'d',h);

rti('td',j)      = tax_rates('td','d',j);
rti('ts',j)      = tax_rates('ts','d',j);
rti('tm',j)      = tax_rates('tm','m',j);
rtp(tp,h)        = hh_tax_rates(tp,h);
rtf('L',i)       = tax_rates('tl','d',i);
rtf('K',i)       = tax_rates('tk','d',i);

tm_io(j,i)        = rti('ts',j) + rti('tm',j);
td_io(j,i)        = rti('ts',j) + rti('td',j);
td_io(i,i)        = 0;
tm_hh(j,h)        = rti('ts',j) + rti('tm',j);
td_hh(j,h)        = rti('ts',j) + rti('td',j);

pm0(i)        = 1;
pd0(i)        = 1;
pf0(f)        = 1;

display rti, rtf;
*rti('ts','c1') = 0;

tiu_imp0(i)        = sum(j,pm0(j)*tm_io(j,i)*IU_M(j,i));
tiu_dom0(i)        = sum(j,pd0(j)*td_io(j,i)*IU_D(j,i));
tfu_imp0(i)        = sum(h,pm0(i)*(rti('ts',i) + rti('tm',i))*pfce0(i,'m',h));
tfu_dom0(i)        = sum(h,pd0(i)*(rti('ts',i) + rti('td',i))*pfce0(i,'d',h));

nit0(i)                 = tiu_imp0(i) + tiu_dom0(i);
factor_tax0(i)          = sum(f, pf0(f)*factor_use(f,i)*rtf(f,i));
va0(i)                  = factor_tax0(i) + sum(f,pf0(f)*factor_use(f,i));
income_tax0(h)          = rtp('ty',h)*sum(f,endowment(f,h)*pf0(f));
total_tax_revenue0      = sum(i,nit0(i) + factor_tax0(i) + tfu_imp0(i) + tfu_dom0(i)) + sum(h,income_tax0(h));
gross_income0(h)        = sum(f,endowment(f,h)*pf0(f));
disposable_income0(h)   =  (1-rtp('ty',h))*gross_income0(h) + rtp('ttr',h)*total_tax_revenue0;
gross_expenditure0(h)   = sum(i,pm0(i)*(1+(rti('ts',i) + rti('tm',i)))*pfce0(i,'m',h))
                        + sum(i,pd0(i)*(1+(rti('ts',i) + rti('td',i)))*pfce0(i,'d',h));
 
*display IU_M, IU_D, FU_M, FU_D, rti, rtp, rtf, tm_io, tm_hh, td_io, tiu_imp0, tiu_dom0, nit0, factor_tax0;
*display income_tax0, tfu_imp0, tfu_dom0, va0, total_tax_revenue0, disposable_income0, gross_expenditure0;

***
loop(h,
if ((abs(disposable_income0(h)-gross_expenditure0(h)) gt 1E-3),
        abort "the SAMS are NOT balanced "
));

** MACROS begins from here
$macro get_sum(i,x) sum(i,x)
$macro get_prod(i,q) prod(i,q)

$macro get_total_cost(i,p,q,sigma) (get_sum(i,p*q))$(sigma eq 1) + (get_sum(i,p*q))$(sigma ne 1)
$macro get_total_share_contribution(i,p,q,sigma) (get_sum(i,p*q))$(sigma eq 1) + (get_sum(i,p*q**(1/sigma)))$(sigma ne 1)

$macro get_factor_share(p,q,sigma,total_factor_expense) ((p*q)/total_factor_expense)$(sigma eq 1) + ((p*q**(1/sigma))/total_factor_expense)$(sigma ne 1)
$macro get_factor_contribution(i,share,q,sigma,rho) (get_prod(i,q**share ))$(sigma eq 1) + ((get_sum(i,share*q**rho))**(1/rho))$(sigma ne 1)

$macro get_scale_factor(output,factor_contribution) output/factor_contribution
$macro compute_production(i,phi,share,q,sigma,rho) (phi*get_prod(i,q**share))$(sigma eq 1) + (phi*get_sum(i,share*q**rho)**(1/rho))$(sigma ne 1)

$macro get_phi(use_cd,use_ces,rho,total,sigma) (total/use_cd)$(sigma eq 1) + (total/use_ces)$(sigma ne 1)
$macro check_calibration(phi,use_cd,use_ces,sigma) (phi*use_cd)$(sigma eq 1) + (phi*use_ces)$(sigma ne 1)

$macro unit_cost(i,p,share,s,phi) ((get_prod(i,(p/share)**share))/phi)$(s eq 1) + (((get_sum(i,(share**s)*(p**(1-s))))**(1/(1-s)))/phi)$(s ne 1)

$macro fac_dd_const(i,p,share,s) (get_prod(i,(share/p)**share))$(s eq 1) + ((get_sum(i,share*(p/share)**(1-s)))**(s/(s-1)))$(s ne 1)
$macro factor_demand(phi,Q,share,p,s,fac_dd_const) (((Q/phi)*(share/p))/fac_dd_const)$(s eq 1) + ((Q/phi)*(share/p)**(s)/fac_dd_const)$(s ne 1)

$macro utility(i,a,s,x) (get_sum(i,a**(1/s)*x**((s-1)/s)))**(s/(s-1))$(s ne 1) + get_prod(i,x**a)$(s eq 1)
*** GAMS Code: Calibration
*** #################  CALIBRATION WILL BEGIN FROM HERE ################## ***
$ontext
We have the following nests
N0P_IU        : CES/CD/Leontief of domestic goods
N0P_VA        : CES/CD between L & K in VA
N1P_IU        : CES/CD/Leontief of IU composite and VA
N1C_FC        : CES/CD in demand function

Each nest will have a price, quantity, scale and share parameter
P_N1P_IU        : is the composite price of the commodity at N1P_IU
Q_N1P_IU        : is the quantity of the composite good at N1P_IU
share_N1P_IU        : share of D & M in the production function at N1P_IU
scale_N1P_IU        : scale factor in the production function at N1P_IU
$offtext

*** GAMS Code: Calibration computing rho=(s-1)/s
parameters
rho_dm_io(i,j)        (sigma_dm_io - 1) : sigma_dm_io
rho_kl(i)             (sigma_kl - 1) : sigma_kl
rho_iu(i)             (sigma_iu - 1) : sigma_iu
rho_iu_va(i)          (sigma_iu_va - 1) : sigma_iu_va
rho_dm_fu(i,h)        (sigma_dm_fu - 1) : sigma_dm_fu
;

rho_dm_io(i,j)        = ((sigma_dm_io(i,j)-1)/sigma_dm_io(i,j))$(sigma_dm_io(i,j) ne 1);
rho_kl(j)             = ((sigma_kl(j)-1)/sigma_kl(j))$(sigma_kl(j) ne 1);
rho_iu(j)             = ((sigma_iu(j)-1)/sigma_iu(j))$(sigma_iu(j) ne 1);
rho_iu_va(j)          = ((sigma_iu_va(j)-1)/sigma_iu_va(j))$(sigma_iu_va(j) ne 1);
rho_dm_fu(i,h)        = ((sigma_dm_fu(i,h)-1)/sigma_dm_fu(i,h))$(sigma_dm_fu(i,h) ne 1);

*** Calibration � Initial parameters
parameters
PMBAR0             benchmark foreign currency value in $
ER0                benchmark exchange rate in Currency Units per $ (ex: � per $)
cif0(i)            cif on import of commodity i
pmd0(i)            domestic price of imports of commodity i
tariff_m0(i,j)     tariff on import of commodity i in sector j
pmd_iu0(i,j)       domestic price of import of commodity i in sector j
pdd0(i)            domestic price of domestic good
pdd_iu0(i,j)       domestic price of domestic good used in sector j
aij_m_iu0(i,j)     use of imported good i in sector j
aij_d_iu0(i,j)     use of domestic good i in sector j
;

*** Calibration - Assigning Data to parameters

PMBAR0         = 1;
ER0            = 1;
cif0(i)        = 0;
pmd0(i)        = PMBAR0*ER0*(1+cif0(i));
tariff_m0(i,j) = tm_io(i,j);
pmd_iu0(i,j)   = pmd0(i)*(1+tariff_m0(i,j))   ;
pdd0(i)        = 1;
pdd_iu0(i,j)   = pdd0(i)*(1+td_io(i,j));

*display tm_io, tm_hh, rti;

*** Calibrating share_N0P_IU0(dm,i,j)  scale_N0P_IU0(i,j)
*** parameters for calibrating the DM composite in A(i,j)
parameters
p_N0P_IU0(dm,i,j)        base price @ nest 0 for D & M of good i in sector j
q_N0P_IU0(dm,i,j)        base qty used @ nest 0 for D & M of good i in sector j
v_N0P_IU0(i,j)           expense on composite DM good i in sector j
sc_N0P_IU0(i,j)          total share contribution
fc_N0P_IU0(i,j)          factor contribution on composite DM good i in sector j
share_N0P_IU0(dm,i,j)    share parameter @ nest 0 for D & M of good i in sector j
scale_N0P_IU0(i,j)       share parameter @ nest 0 for good i in sector j
chk_N0P_IU0(i,j)         check calibration of composite DM good i in sector j
pi_N0P_IU0(i,j)          unit cost of composite DM good i in sector j
qi_N0P_IU0(i,j)          quantity of composite DM good i in sector j
chk_calib_N0P_IU0(i,j)   checking if calibration @ nest 0 is correct
;

*** Calibration � Production Nest 0
aij_m_iu0(i,j)            = IU_M(i,j);
aij_d_iu0(i,j)            = IU_D(i,j);
p_N0P_IU0('d',i,j)        = pdd_iu0(i,j);
p_N0P_IU0('m',i,j)        = pmd_iu0(i,j);
q_N0P_IU0('d',i,j)        = aij_d_iu0(i,j);
q_N0P_IU0('m',i,j)        = aij_m_iu0(i,j);
v_N0P_IU0(i,j)            = get_total_cost(dm,p_N0P_IU0(dm,i,j),q_N0P_IU0(dm,i,j),sigma_dm_io(i,j))  ;
sc_N0P_IU0(i,j)           = get_total_share_contribution(dm,p_N0P_IU0(dm,i,j),q_N0P_IU0(dm,i,j),sigma_dm_io(i,j));
share_N0P_IU0(dm,i,j)     = get_factor_share(p_N0P_IU0(dm,i,j),q_N0P_IU0(dm,i,j),sigma_dm_io(i,j),v_N0P_IU0(i,j));

fc_N0P_IU0(i,j)           = get_factor_contribution(dm,share_N0P_IU0(dm,i,j),q_N0P_IU0(dm,i,j),sigma_dm_io(i,j),rho_dm_io(i,j));
scale_N0P_IU0(i,j)        = get_scale_factor(v_N0P_IU0(i,j),fc_N0P_IU0(i,j));
chk_N0P_IU0(i,j)          = compute_production(dm,scale_N0P_IU0(i,j),share_N0P_IU0(dm,i,j),q_N0P_IU0(dm,i,j),sigma_dm_io(i,j),rho_dm_io(i,j));
pi_N0P_IU0(i,j)           = unit_cost(dm,p_N0P_IU0(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j),scale_N0P_IU0(i,j));
qi_N0P_IU0(i,j)           = compute_production(dm,scale_N0P_IU0(i,j),share_N0P_IU0(dm,i,j),q_N0P_IU0(dm,i,j),sigma_dm_io(i,j),rho_dm_io(i,j)) ;
chk_calib_N0P_IU0(i,j)    = pi_N0P_IU0(i,j)*qi_N0P_IU0(i,j);

*display p_N0P_IU0, rho_dm_io, sigma_dm_io, fc_N0P_IU0, sc_N0P_IU0;
*display chk_N0P_IU0, v_N0P_IU0, chk_calib_N0P_IU0, pi_N0P_IU0, qi_N0P_IU0 ;
*display share_N0P_IU0, scale_N0P_IU0;

*** Calibrating � Production Nest 1 share_N1P_IU0(i,j), scale_N1P_IU0(j)
*** parameters for calibrating the composite intermediate inputs
parameters
p_N1P_IU0(i,j)                base price @ nest 1 for good i in sector j
q_N1P_IU0(i,j)                base quantity @ nest 1 for good i in sector j
v_N1P_IU0(j)                  expense on composite DM good i in sector j
sc_N1P_IU0(j)                 intermediate computation parameter
fc_N1P_IU0(j)                 factor contribution on composite DM good i in sector j
share_N1P_IU0(i,j)            share parameter @ nest 1 for good i in sector j
scale_N1P_IU0(j)              scale parameter @ nest 1 in sector j
chk_N1P_IU0(j)                check calibration of composite DM good i in sector j
pi_N1P_IU0(j)                 unit cost of composite DM good i in sector j
qi_N1P_IU0(j)                 quantity of composite DM good i in sector j
chk_calib_N1P_IU0(j)          check calibration @ nest 1 for IUse
;

*** Calibrating � Production Nest 1 share_N1P_IU0(i,j)
p_N1P_IU0(i,j)       = pi_N0P_IU0(i,j);
q_N1P_IU0(i,j)       = qi_N0P_IU0(i,j);
v_N1P_IU0(j)         = get_total_cost(i,p_N1P_IU0(i,j),q_N1P_IU0(i,j),sigma_iu(j))  ;
sc_N1P_IU0(j)        = get_total_share_contribution(i,p_N1P_IU0(i,j),q_N1P_IU0(i,j),sigma_iu(j));
share_N1P_IU0(i,j)   = get_factor_share(p_N1P_IU0(i,j),q_N1P_IU0(i,j),sigma_iu(j),sc_N1P_IU0(j));
fc_N1P_IU0(j)        = get_factor_contribution(i,share_N1P_IU0(i,j),q_N1P_IU0(i,j),sigma_iu(j),rho_iu(j));
scale_N1P_IU0(j)     = get_scale_factor(v_N1P_IU0(j),fc_N1P_IU0(j));
chk_N1P_IU0(j)       = compute_production(i,scale_N1P_IU0(j),share_N1P_IU0(i,j),q_N1P_IU0(i,j),sigma_iu(j),rho_iu(j));
pi_N1P_IU0(j)        = unit_cost(i,p_N1P_IU0(i,j),share_N1P_IU0(i,j),sigma_iu(j),scale_N1P_IU0(j));
qi_N1P_IU0(j)        = compute_production(i,scale_N1P_IU0(j),share_N1P_IU0(i,j),q_N1P_IU0(i,j),sigma_iu(j),rho_iu(j)) ;
chk_calib_N1P_IU0(j) = pi_N1P_IU0(j)*qi_N1P_IU0(j);

*display chk_N1P_IU0, v_N1P_IU0, chk_calib_N1P_IU0, pi_N1P_IU0, qi_N1P_IU0;
*display share_N1P_IU0, scale_N1P_IU0;

*** Calibrating share_N0P_VA0(f,j), scale_N0P_VA0(j)
parameters
*** parameters for calibrating the LK in Value added
p_N0P_VA0(f,j)              base price @ nest 0 for factor f in VA in sector j
q_N0P_VA0(f,j)              base quantity @ nest 0 for factor f in VA in sector j
v_N0P_VA0(j)                expense on value added in sector j
sc_N0P_VA0(j)               intermediate computation parameter
fc_N0P_VA0(j)               factor contribution of factors L & K  in sector j
share_N0P_VA0(f,j)          share parameter @ nest 0 for factor f in VA in sector j
scale_N0P_VA0(j)            scale parameter @ nest 0 for VA in sector j
chk_N0P_VA0(j)              check calibration of Value Added in sector j
pi_N0P_VA0(j)               unit cost of composite DM good i in sector j
qi_N0P_VA0(j)               quantity of composite DM good i in sector j
chk_calib_N0P_VA0(j)        check if calibration is correct
;

*** GAMS: Calibrating share_N0P_VA0(f,j)
p_N0P_VA0(f,j)        = pf0(f)*(1+rtf(f,j));
q_N0P_VA0(f,j)        = factor_use(f,j);
*display  factor_ue, p_N0P_VA0, q_N0P_VA0;

v_N0P_VA0(j)         = get_total_cost(f,p_N0P_VA0(f,j),q_N0P_VA0(f,j),sigma_kl(j))  ;
sc_N0P_VA0(j)        = get_total_share_contribution(f,p_N0P_VA0(f,j),q_N0P_VA0(f,j),sigma_kl(j));
share_N0P_VA0(f,j)   = get_factor_share(p_N0P_VA0(f,j),q_N0P_VA0(f,j),sigma_kl(j),sc_N0P_VA0(j));

*** GAMS: Calibrating scale_N0P_VA0(j)
fc_N0P_VA0(j)        = get_factor_contribution(f,share_N0P_VA0(f,j),q_N0P_VA0(f,j),sigma_kl(j),rho_kl(j));
scale_N0P_VA0(j)     = get_scale_factor(v_N0P_VA0(j),fc_N0P_VA0(j));
chk_N0P_VA0(j)       = compute_production(f,scale_N0P_VA0(j),share_N0P_VA0(f,j),q_N0P_VA0(f,j),sigma_kl(j),rho_kl(j));
pi_N0P_VA0(j)        = unit_cost(f,p_N0P_VA0(f,j),share_N0P_VA0(f,j),sigma_kl(j),scale_N0P_VA0(j));
qi_N0P_VA0(j)        = compute_production(f,scale_N0P_VA0(j),share_N0P_VA0(f,j),q_N0P_VA0(f,j),sigma_kl(j),rho_kl(j)) ;
chk_calib_N0P_VA0(j) = pi_N0P_VA0(j)*qi_N0P_VA0(j);

*display chk_N0P_VA0, v_N0P_VA0, pi_N0P_VA0, qi_N0P_VA0, chk_calib_N0P_VA0;
*display share_N0P_VA0, scale_N0P_VA0;

*** Calibrating share_N2P_OP0(iuva,j), scale_N2P_OP0(j)
*** parameters for calibrating the Output as a CES/CD/Leontief of IU and VA
parameters
p_N2P_OP0(iuva,j)        base price of IUse & VA in sector j
q_N2P_OP0(iuva,j)        base quantity of IUse & VA in sector j
v_N2P_OP0(j)             expense on composite DM good i in sector j
sc_N2P_OP0(j)            intermediate computation parameter
fc_N2P_OP0(j)            factor contribution on composite DM good i in sector j
share_N2P_OP0(iuva,j)    share parameter @ nest 2 for IUse & VA in sector j
scale_N2P_OP0(j)         scale parameter @ nest 2 for sector j
chk_N2P_OP0(j)           check calibration of composite DM good i in sector j
pi_N2P_OP0(j)            unit cost of composite DM good i in sector j
qi_N2P_OP0(j)            quantity of composite DM good i in sector j
chk_calib_N2P_OP0(j)
;

*** GAMS: Calibrating share_N2P_OP0(iuva,j)
p_N2P_OP0('ciu',j)        = pi_N1P_IU0(j);
p_N2P_OP0('cva',j)        = pi_N0P_VA0(j);
q_N2P_OP0('ciu',j)        = qi_N1P_IU0(j);
q_N2P_OP0('cva',j)        = qi_N0P_VA0(j);

*display p_N2P_OP0, q_N2P_OP0;

v_N2P_OP0(j)                = get_total_cost(iuva,p_N2P_OP0(iuva,j),q_N2P_OP0(iuva,j),sigma_iu_va(j))  ;
sc_N2P_OP0(j)               = get_total_share_contribution(iuva,p_N2P_OP0(iuva,j),q_N2P_OP0(iuva,j),sigma_iu_va(j));
share_N2P_OP0(iuva,j)       = get_factor_share(p_N2P_OP0(iuva,j),q_N2P_OP0(iuva,j),sigma_iu_va(j),sc_N2P_OP0(j));

*** GAMS: Calibrating scale_N2P_OP0(j)
fc_N2P_OP0(j)         = get_factor_contribution(iuva,share_N2P_OP0(iuva,j),q_N2P_OP0(iuva,j),sigma_iu_va(j),rho_iu_va(j));
scale_N2P_OP0(j)      = get_scale_factor(v_N2P_OP0(j),fc_N2P_OP0(j));
chk_N2P_OP0(j)        = compute_production(iuva,scale_N2P_OP0(j),share_N2P_OP0(iuva,j),q_N2P_OP0(iuva,j),sigma_iu_va(j),rho_iu_va(j));
pi_N2P_OP0(j)         = unit_cost(iuva,p_N2P_OP0(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j),scale_N2P_OP0(j));
qi_N2P_OP0(j)         = compute_production(iuva,scale_N2P_OP0(j),share_N2P_OP0(iuva,j),q_N2P_OP0(iuva,j),sigma_iu_va(j),rho_iu_va(j)) ;
chk_calib_N2P_OP0(j)  = pi_N2P_OP0(j)*qi_N2P_OP0(j);

*display chk_N2P_OP0, v_N2P_OP0, pi_N2P_OP0, qi_N2P_OP0,chk_calib_N2P_OP0;
*display scale_N2P_OP0, share_N2P_OP0;

***Calibrating Demand System alpha(i,h)Nest 0 (D,M)
*** parameters for calibrating the DM composite PFCE(i,h)
parameters
pmd_hh0(i,h)                base domestic price of imported good i used by household h
pdd_hh0(i,h)                base domestic price of domestic good i used by household h
p_N0C_HH0(dm,i,h)           base domestic price of domestic & imported good i used by household h
q_N0C_HH0(dm,i,h)           base quantity of domestic & imported good i used by household h
v_N0C_HH0(i,h)              expense on composite DM good i for household h
sc_N0C_HH0(i,h)             intermediate computation parameter
fc_N0C_HH0(i,h)             factor contribution on composite DM good i for household h
share_N0C_HH0(dm,i,h)       share parameter @ nest 0 for D & M good i consumed by household h
scale_N0C_HH0(i,h)          scale parameter @ nest 0 for good i consumed by household h
chk_N0C_HH0(i,h)            check calibration of composite DM good i for household h
pi_N0C_HH0(i,h)             unit cost of composite DM good i for household h
qi_N0C_HH0(i,h)             quantity of composite DM good i for household h
chk_calib_N0C_HH0(i,h)
;

*** Calibrating Demand System alpha(i,h)Nest 0 (D,M)
pmd_hh0(i,h)        = pmd0(i)*(1+tm_hh(i,h));
pdd_hh0(i,h)        = pdd0(i)*(1+td_hh(i,h));

p_N0C_HH0('d',i,h)        = pdd_hh0(i,h);
p_N0C_HH0('m',i,h)        = pmd_hh0(i,h);
q_N0C_HH0('d',i,h)        = FU_D(i,h);
q_N0C_HH0('m',i,h)        = FU_M(i,h);
v_N0C_HH0(i,h)            = get_total_cost(dm,p_N0C_HH0(dm,i,h),q_N0C_HH0(dm,i,h),sigma_dm_fu(i,h))  ;
sc_N0C_HH0(i,h)           = get_total_share_contribution(dm,p_N0C_HH0(dm,i,h),q_N0C_HH0(dm,i,h),sigma_dm_fu(i,h));
share_N0C_HH0(dm,i,h)     = get_factor_share(p_N0C_HH0(dm,i,h),q_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),sc_N0C_HH0(i,h));

fc_N0C_HH0(i,h)        = get_factor_contribution(dm,share_N0C_HH0(dm,i,h),q_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),rho_dm_fu(i,h));
scale_N0C_HH0(i,h)     = get_scale_factor(v_N0C_HH0(i,h),fc_N0C_HH0(i,h));
chk_N0C_HH0(i,h)       = compute_production(dm,scale_N0C_HH0(i,h),share_N0C_HH0(dm,i,h),q_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),rho_dm_fu(i,h));
pi_N0C_HH0(i,h)        = unit_cost(dm,p_N0C_HH0(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),scale_N0C_HH0(i,h));
qi_N0C_HH0(i,h)        = compute_production(dm,scale_N0C_HH0(i,h),share_N0C_HH0(dm,i,h),q_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),rho_dm_fu(i,h)) ;
chk_calib_N0C_HH0(i,h) = pi_N0C_HH0(i,h)*qi_N0C_HH0(i,h);

display chk_N0C_HH0, v_N0C_HH0, chk_calib_N0C_HH0, pi_N0C_HH0, qi_N0C_HH0 ;
display scale_N0C_HH0, share_N0C_HH0;

*** Calibrating Demand System alpha(i,h)Nest 1 (i)
parameters
chk_expenditure0(h)
net_income0(h)
;
chk_expenditure0(h)   = sum(i,chk_calib_N0C_HH0(i,h));
net_income0(h)        = (1-rtp('ty',h))*gross_income0(h);
display gross_income0, net_income0, disposable_income0,chk_expenditure0;

parameter
alpha(i,h)        cobb-douglas demand system parameter
;
alpha(i,h)        =  v_N0C_HH0(i,h)/disposable_income0(h);
display alpha;

*** GAMS Variables 1: Income & Demand Block
tm_io('c1','c1') = 0;
tm_io('c1','c2') = 0;
tm_io('c1','c3') = 0;
tm_io('c1','c4') = 0;

tm_hh('c1','h1') = 0;
tm_hh('c1','h2') = 0;
tm_hh('c1','h3') = 0;

scalar tdbar /0/;
scalar pmbar /1/;

Variables
Z                                dummy optimization variable
TD                               trade deficit
gross_income(h)                  gross income
net_income(h)                    net income
disposable_income(h)             disposable income
dd_N1C_HH(i,h)                   demand for composite commodity i by household h
dd_N0C_HH(dm,i,h)                demand for D and M commodity i by household h
dd_const_N0C_HH(i,h)             computing constant for demand of commodity i by household h

*** GAMS Variables 2: Price Bloc
ER                               exchange rate
pdd(i)                           domestic price of domestic goods
pf(f)                            factor price
pmd(i)                           domestic price of imported good
pmd_iu(i,j)                      domestic price of imported good i in sector j
pdd_iu(i,j)                      domestic price of domestic good i in sector j
pdd_hh(i,h)                      domestic price of domestic good i for household h
pmd_hh(i,h)                      domestic price of imported good i for household h
p_N0P_IU(dm,i,j)                 price at nest 0 of D or M commodity i used in sector j
p_N0C_HH(dm,i,h)              price at nest 0 of D or M commodity i used by household h
pi_N0P_IU(i,j)                price index of composite commodity i used in sector j
p_N1P_IU(i,j)                 price at nest 1 of IU(i j)
pi_N1P_IU(j)                  price index at nest 1 for sector j
p_N0P_VA(f,j)                 price @ nest 0 in value added for factor f in sector j
pi_N0P_VA(j)                  price index of VA in sector j at nest 0
p_N2P_OP(iuva,j)              price of IU and VA @ nest 2 in sector j
pi_N2P_OP(j)                  price index @ nest 2 for sector j
pi_N0C_HH(i,h)                price index at nest 0 for commodity i used by household h
pX_OP(i)                      domestic price of exports

*** GAMS Variables 4: Quantity Block
q_const_N2P_OP(j)        computing constant @ nest 2 for output of sector j
q_N2P_OP(iuva,j)         quantity of IU & VA used in sector j
q_N2P_IU(j)              quantity of IU used in sector j
q_N1P_VA(j)              quantity of VA used in sector j

q_const_N1P_IU(j)        computing constant @ nest 1 for IUse in sector j
q_N1P_IU(i,j)            quantity @ nest 1 of commodity i used by sector j

q_const_N0P_IU(i,j)      computing constant @ nest 0 for use of commodity i in sector j
q_N0P_IU(dm,i,j)         quantity @ nest 0 of D & M commodity i used by sector j

q_const_N0P_VA(j)        computing constant @ nest 0 for VA in sector j
q_N0P_VA(f,j)            quantity of factor f used by setor j @ nest 0

tot_q_OP(dm,i)              total quantity of D & M in sector i
tot_qD_OP(i)                total quantity of domestic good i
tot_qM_OP(i)                total quantity of imported good i
tot_qX_OP(i)                total quantity of exports of good i
tot_qf_SS(f)                total supply of factors f
tot_qf_DD(f)                total demand for fator f
tot_qI_SS(i)                total supply of good i
tot_qI_DD(i)                total demand for good i
qi_N2P_OP(j)                output of sector j

tot_val_M                   total value of imports
tot_val_X                   total value of exports
tax_M_IU                    tax on imported goods in IUse
tax_D_IU                    tax on domestic goods in IUse
tax_M_HH                    tax on imported goods used by households
tax_D_HH                    tax on domestic goods used by households
tax_F_VA                    factor tax in Value added
tax_I_HH                    total income tax
total_tax_Rev               total tax revenue
U0_D(h)                     total utility per household from consuming domestic goods
U0_M(h)                     total utility per household from consuming imported goods
U0_total                    total utility
;

equations
eq_obj                      maximize fictitious obj z=0
eq_pmd(i)                   domestic price of imported good
eq_pmd_iu(i,j)              domestic price of imported good i in sector j
eq_pdd_iu(i,j)              domestic price of domestic good i in sector j
eq_pdd_hh(i,h)              domestic price of domestic good i for household h
eq_pmd_hh(i,h)              domestic price of imported good i for household h
eq_pd_N0P_IU(dm,i,j)        price at nest 0 of D commodity i used in sector j
eq_pm_N0P_IU(dm,i,j)        price at nest 0 of M commodity i used in sector j
eq_pd_N0C_HH(dm,i,h)        price at nest 0 of D commodity i used by household h
eq_pm_N0C_HH(dm,i,h)        price at nest 0 of M commodity i used by household h
eq_pi_N0P_IU(i,j)           price index of composite commodity i used in sector j @ nest 0
eq_p_N1P_IU(i,j)            price at nest 1 of IU(i j)
eq_pi_N1P_IU(j)             price index at nest 1 for sector j
eq_p_N0P_VA(f,j)            price @ nest 0 in VA for factor f in sector j
eq_pi_N0P_VA(j)             price index of VA in sector j at nest 0
eq_p_cva_N2P_OP(iuva,j)     price of VA @ nest 2 in sector j
eq_p_ciu_N2P_OP(iuva,j)     price of IU @ nest 2 in sector j
eq_pi_N2P_OP(j)             price index @ nest 2 for sector j
eq_pi_N0C_HH(i,h)           price index at nest 0 for commodity i used by household h
eq_gross_income(h)          gross income of household h
eq_net_income(h)            net income f household h
eq_disposable_income(h)     disposable income of household h
eq_dd_N1C_HH(i,h)           demand of composite good i by household h

eq_dd_const_N0C_HH(i,h)     computing constant for demand of commodity i by houseold h
eq_dd_N0C_HH(dm,i,h)        demand for D & M good i by household h

eq_q_const_N2P_OP(j)        computing constant @ nest 2 for output of sector j
eq_q_N2P_OP(iuva,j)         quantity of IU & VA used in sector j

eq_q_N2P_IU(j)              quantity of IU used in sector j @ nest 2
eq_q_N1P_VA(j)              quantity of VA used in sector j @ nest 1

eq_q_const_N1P_IU(j)        computing constant @ nest 1 for IUse in sector j
eq_q_N1P_IU(i,j)            quantity @ nest 1 of commodity i used by sector j

eq_q_const_N0P_IU(i,j)      computing constant @ nest 0 for use of commodity i in sector j
eq_q_N0P_IU(dm,i,j)         quantity @ nest 0 of D & M commodity i used by sector j

eq_q_const_N0P_VA(j)        computing constant @ nest 0 for use of VA in sector j
eq_q_N0P_VA(f,j)            quantity of factor f used by sector j @ nest 0

eq_pX_OP(i)                 domestic price of exports
eq_tot_qX_OP(i)             total quantity of exports of good i
eq_tot_q_OP(dm,i)           total quantity of D & M in sector i
eq_tot_qD_OP(i)        total quantity of domestic good i
eq_tot_qM_OP(i)        total quantity of imported good i
eq_tot_qf_SS(f)        total supply of factor f
eq_tot_qf_DD(f)        total demand for factor f
eq_f_eqbm_OP(f)        equilibrium condition for factor f SS >= DD
eq_ZP_N2P_OP(j)        zero profit condition

eq_tot_qI_SS(i)        total supply of good i
eq_tot_qI_DD(i)        total demand for good i
eq_I_eqbm_OP(i)        equilibrium condition for good i SS >= DD

eq_tot_val_M           total value of imports
eq_tot_val_X           total value of exports
eq_TD                  computing trade deficit
eq_zero_td             fixing trade deficit

eq_tax_M_IU                tax on imported goods in IUse
eq_tax_D_IU                tax on domestic goods in IUse
eq_tax_M_HH                tax on imported goods used by households
eq_tax_D_HH                tax on domestic goods used by households
eq_tax_F_VA                factor tax in Value added
eq_tax_I_HH                total income tax
eq_total_tax_Rev           total tax revenue
eq_utility_d(h)            total utility per household from consuming domestic goods
eq_utility_m(h)            total utility per household from consuming imported goods
eq_U0_total                total utility
;

eq_obj..                z =e= 0;
eq_pmd(i)..             pmd(i)             =e= PMBAR*ER*(1+cif0(i));
eq_pmd_iu(i,j)..        pmd_iu(i,j)        =e= pmd(i)*(1+tm_io(i,j));
eq_pdd_iu(i,j)..        pdd_iu(i,j)        =e= pdd(i)*(1+td_io(i,j));
eq_pmd_hh(i,h)..        pmd_hh(i,h)        =e= pmd(i)*(1+tm_hh(i,h));
eq_pdd_hh(i,h)..        pdd_hh(i,h)        =e= pdd(i)*(1+td_hh(i,h));

eq_pd_N0P_IU('d',i,j)..        p_N0P_IU('d',i,j)        =e= pdd_iu(i,j);
eq_pm_N0P_IU('m',i,j)..        p_N0P_IU('m',i,j)        =e= pmd_iu(i,j);
eq_pd_N0C_HH('d',i,h)..        p_N0C_HH('d',i,h)        =e= pdd_hh(i,h);
eq_pm_N0C_HH('m',i,h)..        p_N0C_HH('m',i,h)        =e= pmd_hh(i,h);

eq_pi_N0P_IU(i,j)..        pi_N0P_IU(i,j)    =e= unit_cost(dm,p_N0P_IU(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j),scale_N0P_IU0(i,j));
eq_p_N1P_IU(i,j)..         p_N1P_IU(i,j)     =e= pi_N0P_IU(i,j);
eq_pi_N1P_IU(j)..          pi_N1P_IU(j)      =e= unit_cost(i,p_N1P_IU(i,j),share_N1P_IU0(i,j),sigma_iu(j),scale_N1P_IU0(j));
eq_p_N0P_VA(f,j)..         p_N0P_VA(f,j)     =e= pf(f)*(1+rtf(f,j));
eq_pi_N0P_VA(j)..          pi_N0P_VA(j)      =e= unit_cost(f,p_N0P_VA(f,j),share_N0P_VA0(f,j),sigma_kl(j),scale_N0P_VA0(j));
eq_p_cva_N2P_OP('cva',j).. p_N2P_OP('cva',j) =e= pi_N0P_VA(j);
eq_p_ciu_N2P_OP('ciu',j).. p_N2P_OP('ciu',j) =e= pi_N1P_IU(j);
eq_pi_N2P_OP(j)..          pi_N2P_OP(j)      =e= unit_cost(iuva,p_N2P_OP(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j),scale_N2P_OP0(j));
eq_pi_N0C_HH(i,h)..        pi_N0C_HH(i,h)    =e= unit_cost(dm,p_N0C_HH(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),scale_N0C_HH0(i,h));

eq_gross_income(h)..      gross_income(h)      =e= sum(f,endowment(f,h)*pf(f));
eq_net_income(h)..        net_income(h)        =e= (1-rtp('ty',h))*gross_income(h);

eq_disposable_income(h).. disposable_income(h) =e= net_income(h) + rtp('ttr',h)*total_tax_rev;
eq_dd_N1C_HH(i,h)..       dd_N1C_HH(i,h)       =e= alpha(i,h)*disposable_income(h)/pi_N0C_HH(i,h);

eq_dd_const_N0C_HH(i,h).. dd_const_N0C_HH(i,h) =e= fac_dd_const(dm,p_N0C_HH(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h));
eq_dd_N0C_HH(dm,i,h)..       dd_N0C_HH(dm,i,h) =e= factor_demand(scale_N0C_HH0(i,h),dd_N1C_HH(i,h),share_N0C_HH0(dm,i,h),p_N0C_HH(dm,i,h),sigma_dm_fu(i,h),dd_const_N0C_HH(i,h));

eq_q_const_N2P_OP(j).. q_const_N2P_OP(j) =e= fac_dd_const(iuva,p_N2P_OP(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j));
eq_q_N2P_OP(iuva,j)..   q_N2P_OP(iuva,j) =e= factor_demand(scale_N2P_OP0(j),qi_N2P_OP(j),share_N2P_OP0(iuva,j),p_N2P_OP(iuva,j),sigma_iu_va(j),q_const_N2P_OP(j));
eq_q_N2P_IU(j)..        q_N2P_IU(j)      =e= q_N2P_OP('ciu',j);
eq_q_N1P_VA(j)..        q_N1P_VA(j)      =e= q_N2P_OP('cva',j);
eq_q_const_N1P_IU(j).. q_const_N1P_IU(j) =e= fac_dd_const(i,p_N1P_IU(i,j),share_N1P_IU0(i,j),sigma_iu(j));

eq_q_N1P_IU(i,j)..             q_N1P_IU(i,j) =e= factor_demand(scale_N1P_IU0(j),q_N2P_IU(j),share_N1P_IU0(i,j),p_N1P_IU(i,j),sigma_iu(j),q_const_N1P_IU(j));
eq_q_const_N0P_IU(i,j).. q_const_N0P_IU(i,j) =e= fac_dd_const(dm,p_N0P_IU(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j));
eq_q_N0P_IU(dm,i,j)..       q_N0P_IU(dm,i,j) =e= factor_demand(scale_N0P_IU0(i,j),q_N1P_IU(i,j),share_N0P_IU0(dm,i,j),p_N0P_IU(dm,i,j),sigma_dm_io(i,j),q_const_N0P_IU(i,j));
eq_q_const_N0P_VA(j)..     q_const_N0P_VA(j) =e= fac_dd_const(f,p_N0P_VA(f,j),share_N0P_VA0(f,j),sigma_kl(j));
eq_q_N0P_VA(f,j)..             q_N0P_VA(f,j) =e= factor_demand(scale_N0P_VA0(j),q_N1P_VA(j),share_N0P_VA0(f,j),p_N0P_VA(f,j),sigma_kl(j),q_const_N0P_VA(j));

eq_tot_q_OP(dm,i).. tot_q_OP(dm,i) =e= sum(j,q_N0P_IU(dm,i,j)) + sum(h,dd_N0C_HH(dm,i,h));
eq_pX_OP(i)..         pX_OP(i)     =e= pdd(i)*(1+fob0(i))*(1+t_exp0(i));
eq_tot_qX_OP(i)..     tot_qX_OP(i) =e= exports0(i)*pX_OP(i)**exports_eta(i);
eq_tot_qD_OP(i)..     tot_qD_OP(i) =e= tot_q_OP('d',i) + tot_qX_OP(i);
eq_tot_qM_OP(i)..     tot_qM_OP(i) =e= tot_q_OP('m',i);

eq_ZP_N2P_OP(j)..pdd(j)*qi_N2P_OP(j) =g= sum(i,pmd_iu(i,j)*q_N0P_IU('m',i,j))+sum(i,pdd_iu(i,j)*q_N0P_IU('d',i,j)) + sum(f,p_N0P_VA(f,j)*q_N0P_VA(f,j)) ;

eq_tot_qI_SS(i).. tot_qI_SS(i) =e= qi_N2P_OP(i);
eq_tot_qI_DD(i).. tot_qI_DD(i) =e= tot_qD_OP(i);
eq_I_eqbm_OP(i).. tot_qI_SS(i) =g= tot_qI_DD(i);

*factor equilibrium
eq_tot_qf_SS(f).. tot_qf_SS(f) =e= sum(h,endowment(f,h));
eq_tot_qf_DD(f).. tot_qf_DD(f) =e= sum(j,q_N0P_VA(f,j));
eq_f_eqbm_OP(f).. tot_qf_SS(f) =g= tot_qf_DD(f);

* computing the value of imports, exports and trade deficit
eq_tot_val_M.. tot_val_M =e= sum((i,j),pmd(i)*q_N0P_IU('m',i,j)) + sum((i,h), pmd(i)*dd_N0C_HH('m',i,h));
eq_tot_val_X.. tot_val_X =e= sum(i,pX_OP(i)*tot_qX_OP(i));
eq_TD..            TD*ER =e= (tot_val_X - tot_val_M);
eq_zero_td..          TD =e= 0;

* tax revenue from tariffs levied on IUse of imports
eq_tax_M_IU.. tax_M_IU =e= sum((i,j),pmd(i)*tm_io(i,j)*q_N0P_IU('m',i,j));
*tax revenue from excise and sales taxes on IUse of domestic goods
eq_tax_D_IU.. tax_D_IU =e= sum((i,j),pdd(i)*td_io(i,j)*q_N0P_IU('d',i,j));
* tax revenue from tariffs levied on final consumption of households
eq_tax_M_HH.. tax_M_HH =e= sum((i,h),pmd(i)*tm_hh(i,h)*dd_N0C_HH('m',i,h));
* tax revenue from sales, excise taxes on final consumption of households
eq_tax_D_HH.. tax_D_HH =e= sum((i,h),pdd(i)*td_hh(i,h)*dd_N0C_HH('d',i,h));
* tax revenue from factor taxes
eq_tax_F_VA.. tax_F_VA =e= sum((f,j),pf(f)*rtf(f,j)*q_N0P_VA(f,j));
* tax revenue from income tax on households
eq_tax_I_HH.. tax_I_HH =e= sum(h,rtp('ty',h)*gross_income(h)) ;
* total tax revenue
eq_total_tax_Rev.. total_tax_Rev =e= tax_M_IU + tax_D_IU + tax_M_HH + tax_D_HH + tax_F_VA + tax_I_HH;
*utility for HH domestic
eq_utility_d(h)..   U0_D(h)     =e= utility(i,alpha(i,h),1,dd_N0C_HH('d',i,h));
*utility for HH imported
eq_utility_m(h)..   U0_M(h)     =e= utility(i,alpha(i,h),1,dd_N0C_HH('m',i,h));
* total utility
eq_U0_total..       U0_total    =e= sum(h,U0_M(h)) + sum(h,U0_D(h));


model cge_443_ces
*/all/;
/
eq_obj
eq_pmd
eq_pmd_iu
eq_pdd_iu
eq_pdd_hh
eq_pmd_hh
eq_pd_N0P_IU
eq_pm_N0P_IU
eq_pd_N0C_HH
eq_pm_N0C_HH
eq_pi_N0P_IU
eq_p_N1P_IU
eq_pi_N1P_IU
eq_p_N0P_VA
eq_pi_N0P_VA

eq_p_cva_N2P_OP
eq_p_ciu_N2P_OP
eq_pi_N2P_OP
eq_pi_N0C_HH
eq_gross_income
eq_net_income
eq_disposable_income
eq_dd_N1C_HH
eq_dd_const_N0C_HH
eq_dd_N0C_HH

eq_q_const_N2P_OP
eq_q_N2P_OP
eq_q_N2P_IU
eq_q_N1P_VA

eq_q_const_N1P_IU
eq_q_N1P_IU

eq_q_const_N0P_IU
eq_q_N0P_IU

eq_q_const_N0P_VA
eq_q_N0P_VA
eq_tot_q_OP

eq_pX_OP
eq_tot_qX_OP

eq_tot_qD_OP
eq_tot_qM_OP

eq_tot_qf_SS
eq_tot_qf_DD
eq_f_eqbm_OP

eq_ZP_N2P_OP
eq_tot_qI_SS
eq_tot_qI_DD
eq_I_eqbm_OP

eq_tot_val_M
eq_tot_val_X

eq_TD
eq_zero_td
eq_tax_M_IU
eq_tax_D_IU
eq_tax_M_HH
eq_tax_D_HH
eq_tax_F_VA
eq_tax_I_HH
eq_total_tax_Rev
eq_utility_d
eq_utility_m
eq_U0_total
/
;

model cge_443_ces_mcp
/
eq_obj.z
eq_pmd.pmd
eq_pmd_iu.pmd_iu
eq_pdd_iu.pdd_iu
eq_pdd_hh.pdd_hh
eq_pmd_hh.pmd_hh
eq_pd_N0P_IU.p_N0P_IU
eq_pm_N0P_IU.p_N0P_IU
eq_pd_N0C_HH.p_N0C_HH
eq_pm_N0C_HH.p_N0C_HH
eq_pi_N0P_IU.pi_N0P_IU
eq_p_N1P_IU.p_N1P_IU
eq_pi_N1P_IU.pi_N1P_IU
eq_p_N0P_VA.p_N0P_VA
eq_pi_N0P_VA.pi_N0P_VA

eq_p_cva_N2P_OP.p_N2P_OP
eq_p_ciu_N2P_OP.p_N2P_OP
eq_pi_N2P_OP.pi_N2P_OP
eq_pi_N0C_HH.pi_N0C_HH
eq_gross_income.gross_income
eq_net_income.net_income
eq_disposable_income.disposable_income
eq_dd_N1C_HH.dd_N1C_HH
eq_dd_const_N0C_HH.dd_const_N0C_HH
eq_dd_N0C_HH.dd_N0C_HH

eq_q_const_N2P_OP.q_const_N2P_OP
eq_q_N2P_OP.q_N2P_OP
eq_q_N2P_IU.q_N2P_IU
eq_q_N1P_VA.q_N1P_VA

eq_q_const_N1P_IU.q_const_N1P_IU
eq_q_N1P_IU.q_N1P_IU

eq_q_const_N0P_IU.q_const_N0P_IU
eq_q_N0P_IU.q_N0P_IU

eq_q_const_N0P_VA.q_const_N0P_VA
eq_q_N0P_VA.q_N0P_VA
eq_tot_q_OP.tot_q_OP

eq_pX_OP.pX_OP
eq_tot_qX_OP.tot_qX_OP

eq_tot_qD_OP.tot_qD_OP
eq_tot_qM_OP.tot_qM_OP

eq_tot_qf_SS.tot_qf_SS
eq_tot_qf_DD.tot_qf_DD
eq_f_eqbm_OP.pf

eq_ZP_N2P_OP.qi_N2P_OP
eq_tot_qI_SS.tot_qI_SS
eq_tot_qI_DD.tot_qI_DD
eq_I_eqbm_OP.pdd

eq_tot_val_M.tot_val_M
eq_tot_val_X.tot_val_X

eq_TD.ER
eq_zero_td.td
eq_tax_M_IU.tax_M_IU
eq_tax_D_IU.tax_D_IU
eq_tax_M_HH.tax_M_HH
eq_tax_D_HH.tax_D_HH
eq_tax_F_VA.tax_F_VA
eq_tax_I_HH.tax_I_HH
eq_total_tax_Rev.total_tax_Rev
eq_utility_d.U0_D
eq_utility_m.U0_M
eq_U0_total.U0_total
/
;
*** GAMS: Assign Variable Levels-1
ER.l                = 1;
pdd.l(i)            = 1;
pf.l(f)             = 1;
pf.fx('L')          = 1;
pmd.l(i)             = PMBAR*ER.l*(1+cif0(i));
pmd_iu.l(i,j)        = pmd.l(i)*(1+tm_io(i,j));
pdd_iu.l(i,j)        = pdd.l(i)*(1+td_io(i,j));
pmd_hh.l(i,h)        = pmd.l(i)*(1+tm_hh(i,h));
pdd_hh.l(i,h)        = pdd.l(i)*(1+td_hh(i,h));
p_N0P_IU.l('d',i,j)        = pdd_iu.l(i,j);
p_N0P_IU.l('m',i,j)        = pmd_iu.l(i,j);
p_N0C_HH.l('d',i,h)        = pdd_hh.l(i,h);
p_N0C_HH.l('m',i,h)        = pmd_hh.l(i,h);
pi_N0P_IU.l(i,j)        = unit_cost(dm,p_N0P_IU.l(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j),scale_N0P_IU0(i,j));
p_N1P_IU.l(i,j)         = pi_N0P_IU.l(i,j);
pi_N1P_IU.l(j)          = unit_cost(i,p_N1P_IU.l(i,j),share_N1P_IU0(i,j),sigma_iu(j),scale_N1P_IU0(j));
p_N0P_VA.l(f,j)         = pf.l(f)*(1+rtf(f,j));
pi_N0P_VA.l(j)          = unit_cost(f,p_N0P_VA.l(f,j),share_N0P_VA0(f,j),sigma_kl(j),scale_N0P_VA0(j));
p_N2P_OP.l('cva',j)     = pi_N0P_VA.l(j);
p_N2P_OP.l('ciu',j)     = pi_N1P_IU.l(j);
pi_N2P_OP.l(j)          = unit_cost(iuva,p_N2P_OP.l(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j),scale_N2P_OP0(j));
pi_N0C_HH.l(i,h)        = unit_cost(dm,p_N0C_HH.l(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),scale_N0C_HH0(i,h));

pi_N0P_IU.l(i,j)        = unit_cost(dm,p_N0P_IU.l(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j),scale_N0P_IU0(i,j));
p_N1P_IU.l(i,j)         = pi_N0P_IU.l(i,j);
pi_N1P_IU.l(j)          = unit_cost(i,p_N1P_IU.l(i,j),share_N1P_IU0(i,j),sigma_iu(j),scale_N1P_IU0(j));
p_N0P_VA.l(f,j)         = pf.l(f)*(1+rtf(f,j));
pi_N0P_VA.l(j)          = unit_cost(f,p_N0P_VA.l(f,j),share_N0P_VA0(f,j),sigma_kl(j),scale_N0P_VA0(j));
p_N2P_OP.l('cva',j)     = pi_N0P_VA.l(j);
p_N2P_OP.l('ciu',j)     = pi_N1P_IU.l(j);
pi_N2P_OP.l(j)          = unit_cost(iuva,p_N2P_OP.l(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j),scale_N2P_OP0(j));
pi_N0C_HH.l(i,h)        = unit_cost(dm,p_N0C_HH.l(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),scale_N0C_HH0(i,h));

total_tax_rev.l        = total_tax_revenue0;
gross_income.l(h)      = sum(f,endowment(f,h)*pf.l(f));
net_income.l(h)        = (1-rtp('ty',h))*gross_income.l(h);
disposable_income.l(h) = net_income.l(h) + rtp('ttr',h)*total_tax_rev.l;
dd_N1C_HH.l(i,h)       = alpha(i,h)*disposable_income.l(h)/pi_N0C_HH.l(i,h);

qi_N2P_OP.l(j) = qi_N2P_OP0(j);
tot_qD_OP.l(j) = qi_N2P_OP0(j);
*computing the demand for D & M in final consumption
dd_const_N0C_HH.l(i,h) = fac_dd_const(dm,p_N0C_HH.l(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h));
dd_N0C_HH.l(dm,i,h)    = factor_demand(scale_N0C_HH0(i,h),dd_N1C_HH.l(i,h),share_N0C_HH0(dm,i,h),p_N0C_HH.l(dm,i,h),sigma_dm_fu(i,h),dd_const_N0C_HH.l(i,h));

q_const_N2P_OP.l(j)  = fac_dd_const(iuva,p_N2P_OP.l(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j));
q_N2P_OP.l(iuva,j)   = factor_demand(scale_N2P_OP0(j),tot_qD_OP.l(j),share_N2P_OP0(iuva,j),p_N2P_OP.l(iuva,j),sigma_iu_va(j),q_const_N2P_OP.l(j));
q_N2P_IU.l(j)        = q_N2P_OP.l('ciu',j);
q_N1P_VA.l(j)        = q_N2P_OP.l('cva',j);
q_const_N1P_IU.l(j)  = fac_dd_const(i,p_N1P_IU.l(i,j),share_N1P_IU0(i,j),sigma_iu(j));
q_N1P_IU.l(i,j)      = factor_demand(scale_N1P_IU0(j),q_N2P_IU.l(j),share_N1P_IU0(i,j),p_N1P_IU.l(i,j),sigma_iu(j),q_const_N1P_IU.l(j));

q_const_N0P_IU.l(i,j) = fac_dd_const(dm,p_N0P_IU.l(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j));
q_N0P_IU.l(dm,i,j)    = factor_demand(scale_N0P_IU0(i,j),q_N1P_IU.l(i,j),share_N0P_IU0(dm,i,j),p_N0P_IU.l(dm,i,j),sigma_dm_io(i,j),q_const_N0P_IU.l(i,j));

q_const_N0P_VA.l(j)  = fac_dd_const(f,p_N0P_VA.l(f,j),share_N0P_VA0(f,j),sigma_kl(j));
q_N0P_VA.l(f,j)      = factor_demand(scale_N0P_VA0(j),q_N1P_VA.l(j),share_N0P_VA0(f,j),p_N0P_VA.l(f,j),sigma_kl(j),q_const_N0P_VA.l(j));

pX_OP.l(i)        = pdd.l(i)*(1+fob0(i))*(1+t_exp0(i));
tot_qX_OP.l(i)    = exports0(i)*pX_OP.l(i)**exports_eta(i);

* computing the total demand for D & M
tot_q_OP.l(dm,i)      = sum(j,q_N0P_IU.l(dm,i,j)) + sum(h,dd_N0C_HH.l(dm,i,h)) ;
tot_qD_OP.l(i)        = tot_q_OP.l('d',i) + tot_qX_OP.l(i);
tot_qM_OP.l(i)        = tot_q_OP.l('m',i);

tot_qf_SS.l(f)        = sum(h,endowment(f,h));
tot_qf_DD.l(f)        = sum(j,q_N0P_VA.l(f,j));

tot_qI_SS.l(i)        = qi_N2P_OP.l(i);
tot_qI_DD.l(i)        = tot_qD_OP.l(i);

tot_val_M.l        = sum((i,j),pmd.l(i)*q_N0P_IU.l('m',i,j)) + sum((i,h), pmd.l(i)*dd_N0C_HH.l('m',i,h));
tot_val_X.l        = sum(i,pX_OP.l(i)*tot_qX_OP.l(i));

tax_M_IU.l        = sum((i,j),pmd.l(i)*tm_io(i,j)*q_N0P_IU.l('m',i,j));
tax_D_IU.l        = sum((i,j),pdd.l(i)*td_io(i,j)*q_N0P_IU.l('d',i,j));
tax_M_HH.l        = sum((i,h),pmd.l(i)*tm_hh(i,h)*dd_N0C_HH.l('m',i,h));
tax_D_HH.l        = sum((i,h),pdd.l(i)*td_hh(i,h)*dd_N0C_HH.l('d',i,h));
tax_F_VA.l        = sum((f,j),pf.l(f)*rtf(f,j)*q_N0P_VA.l(f,j));
tax_I_HH.l        = sum(h,rtp('ty',h)*gross_income.l(h)) ;
total_tax_Rev.l   = tax_M_IU.l + tax_D_IU.l + tax_M_HH.l + tax_D_HH.l + tax_F_VA.l + tax_I_HH.l;
U0_D.l(h)           = utility(i,alpha(i,h),1,dd_N0C_HH.l('d',i,h));
U0_M.l(h)           = utility(i,alpha(i,h),1,dd_N0C_HH.l('m',i,h));
U0_total.l         = sum(h,U0_M.l(h)) + sum(h,U0_D.l(h));

*** GAMS: Assign Variable Bounds-1
ER.lo                = 0.001*ER.l    ;
pdd.lo(i)            = 0.001*pdd.l(i);
pf.lo(f)             = 0.001*pf.l(f) ;

pmd.lo(i)             = 0.001*pmd.l(i)     ;
pmd_iu.lo(i,j)        = 0.001*pmd_iu.l(i,j);
pdd_iu.lo(i,j)        = 0.001*pdd_iu.l(i,j);
pmd_hh.lo(i,h)        = 0.001*pmd_hh.l(i,h);
pdd_hh.lo(i,h)        = 0.001*pdd_hh.l(i,h);

p_N0P_IU.lo('d',i,j)        = 0.001*p_N0P_IU.l('d',i,j);
p_N0P_IU.lo('m',i,j)        = 0.001*p_N0P_IU.l('m',i,j);
p_N0C_HH.lo('d',i,h)        = 0.001*p_N0C_HH.l('d',i,h);
p_N0C_HH.lo('m',i,h)        = 0.001*p_N0C_HH.l('m',i,h);

pi_N0P_IU.lo(i,j)      = 0.001*pi_N0P_IU.l(i,j)   ;
p_N1P_IU.lo(i,j)       = 0.001*p_N1P_IU.l(i,j)    ;
pi_N1P_IU.lo(j)        = 0.001*pi_N1P_IU.l(j)     ;
p_N0P_VA.lo(f,j)       = 0.001*p_N0P_VA.l(f,j)    ;
pi_N0P_VA.lo(j)        = 0.001*pi_N0P_VA.l(j)     ;
p_N2P_OP.lo('cva',j)   = 0.001*p_N2P_OP.l('cva',j);
p_N2P_OP.lo('ciu',j)   = 0.001*p_N2P_OP.l('ciu',j);
pi_N2P_OP.lo(j)        = 0.001*pi_N2P_OP.l(j)     ;
pi_N0C_HH.lo(i,h)      = 0.001*pi_N0C_HH.l(i,h)   ;

total_tax_rev.lo        = 0.001*total_tax_rev.l    ;
gross_income.lo(h)      = 0.001*gross_income.l(h)      ;
net_income.lo(h)        = 0.001*net_income.l(h)        ;
disposable_income.lo(h) = 0.001*disposable_income.l(h) ;
dd_N1C_HH.lo(i,h)       = 0.001*dd_N1C_HH.l(i,h)       ;
dd_const_N0C_HH.lo(i,h) = 0.001*dd_const_N0C_HH.l(i,h);
dd_N0C_HH.lo(dm,i,h)    = 0.001*dd_N0C_HH.l(dm,i,h);

q_const_N2P_OP.lo(j)    = 0.001*q_const_N2P_OP.l(j);
q_N2P_OP.l(iuva,j)      = 0.001*q_N2P_OP.l(iuva,j);
q_N2P_IU.lo(j)          = 0.001*q_N2P_IU.l(j);
q_N1P_VA.lo(j)          = 0.001*q_N1P_VA.l(j);
q_const_N1P_IU.lo(j)    = 0.001*q_const_N1P_IU.l(j);
q_N1P_IU.l(i,j)         = 0.001*q_N1P_IU.l(i,j);

q_const_N0P_IU.lo(i,j)     = 0.001*q_const_N0P_IU.l(i,j);
q_N0P_IU.lo(dm,i,j)        = 0.001*q_N0P_IU.l(dm,i,j);
q_const_N0P_VA.lo(j)       = 0.001*q_const_N0P_VA.l(j);
q_N0P_VA.lo(f,j)           = 0.001*q_N0P_VA.l(f,j);

tot_q_OP.lo(dm,i)              = 0.001*tot_q_OP.l(dm,i);
tot_qD_OP.lo(i)                = 0.001*tot_qD_OP.l(i);
tot_qM_OP.lo(i)                = 0.001*tot_qM_OP.l(i);
pX_OP.lo(i)                    = 0.001*pX_OP.l(i);
tot_qX_OP.lo(i)                = 0.001*tot_qX_OP.l(i);

tot_qf_SS.lo(f)                = 0.001*tot_qf_SS.l(f);
tot_qf_DD.lo(f)                = 0.001*tot_qf_DD.l(f) ;
tot_val_M.lo                   = 0.001*tot_val_M.l;
tot_val_X.lo                   = 0.001*tot_val_X.l;
U0_D.lo(h)           = 0.001*U0_D.l(h);
U0_M.lo(h)           = 0.001*U0_M.lo(h);
U0_total.lo          = 0.001*U0_total.l;

display pf.l;

scalar solver_choice /1/;
if (solver_choice eq 0,
   solve cge_443_ces using nlp minimising z;
else
   SOLVE cge_443_ces_mcp USING MCP;
);



file prdn /dogadin_shock1.txt/;
put prdn;
prdn.pw = 30*30;

put 'Total Quantities'  put /;
put '----------------'  put /;
put '         ' put  'Domestic      ' put 'Imported      ' put 'Exports  ' put /;
put 'c1':<9,    put tot_qD_OP.l('c1'):<12, put '  ',
                put tot_qM_OP.l('c1'):<12, put '  ',
                put tot_qX_OP.l('c1'):<12, put /;

put /;
put 'Total Value'  put /;
put '----------------'  put /;
put  'Imports      ' put 'Exports' put /;
put tot_val_M.l:<9   put ' ' put tot_val_X.l put /;

put /;
put 'Collected Taxes on Imported Goods'  put /;
put '---------------------------------'  put /;
put  'IUse          ' put 'Households' put /;
put tax_M_IU.l:<9     put ' ' put tax_M_HH.l put /;

put /;
put 'Total HH Utility from Consumption'  put /;
put '---------------------------------'  put /;
put  'Domestic Goods      ' put 'Imported Goods' put /;
put sum(h,U0_D.l(h)):<20    put sum(h,U0_M.l(h)):<9 put /;

put /;
put 'Price Index for Agriculture'  put /;
put '---------------------------'  put /;
put 'c1':<15,    put pi_N2P_OP.l('c1'); put /;
put /;
*********************************************
************ counterfactual case ************
*********************************************
* imposing tariff on imported agricultural products
* Substitute the desired percent (.1, .25, .5, .75)
tm_io('c1','c1') = 0.75;
tm_io('c1','c2') = 0.75;
tm_io('c1','c3') = 0.75;
tm_io('c1','c4') = 0.75;

tm_hh('c1','h1') = 0.75;
tm_hh('c1','h2') = 0.75;
tm_hh('c1','h3') = 0.75;

*** GAMS Variables 1: Income & Demand Block
scalar c_tdbar /0/;
scalar c_pmbar /1/;

Variables
c_Z                                dummy optimization variable
c_TD                               trade deficit
c_gross_income(h)                  gross income
c_net_income(h)                    net income
c_disposable_income(h)             disposable income
c_dd_N1C_HH(i,h)                   demand for composite commodity i by household h
c_dd_N0C_HH(dm,i,h)                demand for D and M commodity i by household h
c_dd_const_N0C_HH(i,h)             computing constant for demand of commodity i by household h

*** GAMS Variables 2: Price Bloc
c_ER                               exchange rate
c_pdd(i)                           domestic price of domestic goods
c_pf(f)                            factor price
c_pmd(i)                           domestic price of imported good
c_pmd_iu(i,j)                      domestic price of imported good i in sector j
c_pdd_iu(i,j)                      domestic price of domestic good i in sector j
c_pdd_hh(i,h)                      domestic price of domestic good i for household h
c_pmd_hh(i,h)                      domestic price of imported good i for household h
c_p_N0P_IU(dm,i,j)                 price at nest 0 of D or M commodity i used in sector j
c_p_N0C_HH(dm,i,h)              price at nest 0 of D or M commodity i used by household h
c_pi_N0P_IU(i,j)                price index of composite commodity i used in sector j
c_p_N1P_IU(i,j)                 price at nest 1 of IU(i j)
c_pi_N1P_IU(j)                  price index at nest 1 for sector j
c_p_N0P_VA(f,j)                 price @ nest 0 in value added for factor f in sector j
c_pi_N0P_VA(j)                  price index of VA in sector j at nest 0
c_p_N2P_OP(iuva,j)              price of IU and VA @ nest 2 in sector j
c_pi_N2P_OP(j)                  price index @ nest 2 for sector j
c_pi_N0C_HH(i,h)                price index at nest 0 for commodity i used by household h
c_pX_OP(i)                      domestic price of exports

*** GAMS Variables 4: Quantity Block
c_q_const_N2P_OP(j)        computing constant @ nest 2 for output of sector j
c_q_N2P_OP(iuva,j)         quantity of IU & VA used in sector j
c_q_N2P_IU(j)              quantity of IU used in sector j
c_q_N1P_VA(j)              quantity of VA used in sector j

c_q_const_N1P_IU(j)        computing constant @ nest 1 for IUse in sector j
c_q_N1P_IU(i,j)            quantity @ nest 1 of commodity i used by sector j

c_q_const_N0P_IU(i,j)      computing constant @ nest 0 for use of commodity i in sector j
c_q_N0P_IU(dm,i,j)         quantity @ nest 0 of D & M commodity i used by sector j

c_q_const_N0P_VA(j)        computing constant @ nest 0 for VA in sector j
c_q_N0P_VA(f,j)            quantity of factor f used by setor j @ nest 0

c_tot_q_OP(dm,i)              total quantity of D & M in sector i
c_tot_qD_OP(i)                total quantity of domestic good i
c_tot_qM_OP(i)                total quantity of imported good i
c_tot_qX_OP(i)                total quantity of exports of good i
c_tot_qf_SS(f)                total supply of factors f
c_tot_qf_DD(f)                total demand for fator f
c_tot_qI_SS(i)                total supply of good i
c_tot_qI_DD(i)                total demand for good i
c_qi_N2P_OP(j)                output of sector j

c_tot_val_M       total value of imports
c_tot_val_X       total value of exports
c_tax_M_IU        tax on imported goods in IUse
c_tax_D_IU        tax on domestic goods in IUse
c_tax_M_HH        tax on imported goods used by households
c_tax_D_HH        tax on domestic goods used by households
c_tax_F_VA        factor tax in Value added
c_tax_I_HH        total income tax
c_total_tax_Rev   total tax revenue
c_U0_D(h)                     total utility per household from consuming domestic goods
c_U0_M(h)                     total utility per household from consuming imported goods
c_U0_total                    total utility

;

equations
c_eq_obj                      maximize fictitious obj z=0
c_eq_pmd(i)                   domestic price of imported good
c_eq_pmd_iu(i,j)              domestic price of imported good i in sector j
c_eq_pdd_iu(i,j)              domestic price of domestic good i in sector j
c_eq_pdd_hh(i,h)              domestic price of domestic good i for household h
c_eq_pmd_hh(i,h)              domestic price of imported good i for household h
c_eq_pd_N0P_IU(dm,i,j)        price at nest 0 of D commodity i used in sector j
c_eq_pm_N0P_IU(dm,i,j)        price at nest 0 of M commodity i used in sector j
c_eq_pd_N0C_HH(dm,i,h)        price at nest 0 of D commodity i used by household h
c_eq_pm_N0C_HH(dm,i,h)        price at nest 0 of M commodity i used by household h
c_eq_pi_N0P_IU(i,j)           price index of composite commodity i used in sector j @ nest 0
c_eq_p_N1P_IU(i,j)            price at nest 1 of IU(i j)
c_eq_pi_N1P_IU(j)             price index at nest 1 for sector j
c_eq_p_N0P_VA(f,j)            price @ nest 0 in VA for factor f in sector j
c_eq_pi_N0P_VA(j)             price index of VA in sector j at nest 0
c_eq_p_cva_N2P_OP(iuva,j)     price of VA @ nest 2 in sector j
c_eq_p_ciu_N2P_OP(iuva,j)     price of IU @ nest 2 in sector j
c_eq_pi_N2P_OP(j)             price index @ nest 2 for sector j
c_eq_pi_N0C_HH(i,h)           price index at nest 0 for commodity i used by household h
c_eq_gross_income(h)          gross income of household h
c_eq_net_income(h)            net income f household h
c_eq_disposable_income(h)     disposable income of household h
c_eq_dd_N1C_HH(i,h)           demand of composite good i by household h

c_eq_dd_const_N0C_HH(i,h)     computing constant for demand of commodity i by houseold h
c_eq_dd_N0C_HH(dm,i,h)        demand for D & M good i by household h

c_eq_q_const_N2P_OP(j)        computing constant @ nest 2 for output of sector j
c_eq_q_N2P_OP(iuva,j)         quantity of IU & VA used in sector j

c_eq_q_N2P_IU(j)              quantity of IU used in sector j @ nest 2
c_eq_q_N1P_VA(j)              quantity of VA used in sector j @ nest 1

c_eq_q_const_N1P_IU(j)        computing constant @ nest 1 for IUse in sector j
c_eq_q_N1P_IU(i,j)            quantity @ nest 1 of commodity i used by sector j

c_eq_q_const_N0P_IU(i,j)      computing constant @ nest 0 for use of commodity i in sector j
c_eq_q_N0P_IU(dm,i,j)         quantity @ nest 0 of D & M commodity i used by sector j

c_eq_q_const_N0P_VA(j)        computing constant @ nest 0 for use of VA in sector j
c_eq_q_N0P_VA(f,j)            quantity of factor f used by sector j @ nest 0

c_eq_pX_OP(i)                 domestic price of exports
c_eq_tot_qX_OP(i)             total quantity of exports of good i
c_eq_tot_q_OP(dm,i)           total quantity of D & M in sector i
c_eq_tot_qD_OP(i)        total quantity of domestic good i
c_eq_tot_qM_OP(i)        total quantity of imported good i
c_eq_tot_qf_SS(f)        total supply of factor f
c_eq_tot_qf_DD(f)        total demand for factor f
c_eq_f_eqbm_OP(f)        equilibrium condition for factor f SS >= DD
c_eq_ZP_N2P_OP(j)        zero profit condition

c_eq_tot_qI_SS(i)        total supply of good i
c_eq_tot_qI_DD(i)        total demand for good i
c_eq_I_eqbm_OP(i)        equilibrium condition for good i SS >= DD

c_eq_tot_val_M           total value of imports
c_eq_tot_val_X           total value of exports
c_eq_TD                  computing trade deficit
c_eq_zero_td             fixing trade deficit

c_eq_tax_M_IU                tax on imported goods in IUse
c_eq_tax_D_IU                tax on domestic goods in IUse
c_eq_tax_M_HH                tax on imported goods used by households
c_eq_tax_D_HH                tax on domestic goods used by households
c_eq_tax_F_VA                factor tax in Value added
c_eq_tax_I_HH                total income tax
c_eq_total_tax_Rev           total tax revenue
c_eq_utility_d(h)              total utility per household from consuming domestic goods
c_eq_utility_m(h)              total utility per household from consuming imported goods
c_eq_U0_total                  total utility
;

c_eq_obj..                c_z =e= 0;
c_eq_pmd(i)..             c_pmd(i)             =e= c_PMBAR*c_ER*(1+cif0(i));
c_eq_pmd_iu(i,j)..        c_pmd_iu(i,j)        =e= c_pmd(i)*(1+tm_io(i,j));
c_eq_pdd_iu(i,j)..        c_pdd_iu(i,j)        =e= c_pdd(i)*(1+td_io(i,j));
c_eq_pmd_hh(i,h)..        c_pmd_hh(i,h)        =e= c_pmd(i)*(1+tm_hh(i,h));
c_eq_pdd_hh(i,h)..        c_pdd_hh(i,h)        =e= c_pdd(i)*(1+td_hh(i,h));

c_eq_pd_N0P_IU('d',i,j)..        c_p_N0P_IU('d',i,j)        =e= c_pdd_iu(i,j);
c_eq_pm_N0P_IU('m',i,j)..        c_p_N0P_IU('m',i,j)        =e= c_pmd_iu(i,j);
c_eq_pd_N0C_HH('d',i,h)..        c_p_N0C_HH('d',i,h)        =e= c_pdd_hh(i,h);
c_eq_pm_N0C_HH('m',i,h)..        c_p_N0C_HH('m',i,h)        =e= c_pmd_hh(i,h);

c_eq_pi_N0P_IU(i,j)..        c_pi_N0P_IU(i,j)    =e= unit_cost(dm,c_p_N0P_IU(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j),scale_N0P_IU0(i,j));
c_eq_p_N1P_IU(i,j)..         c_p_N1P_IU(i,j)     =e= c_pi_N0P_IU(i,j);
c_eq_pi_N1P_IU(j)..          c_pi_N1P_IU(j)      =e= unit_cost(i,c_p_N1P_IU(i,j),share_N1P_IU0(i,j),sigma_iu(j),scale_N1P_IU0(j));
c_eq_p_N0P_VA(f,j)..         c_p_N0P_VA(f,j)     =e= c_pf(f)*(1+rtf(f,j));
c_eq_pi_N0P_VA(j)..          c_pi_N0P_VA(j)      =e= unit_cost(f,c_p_N0P_VA(f,j),share_N0P_VA0(f,j),sigma_kl(j),scale_N0P_VA0(j));
c_eq_p_cva_N2P_OP('cva',j).. c_p_N2P_OP('cva',j) =e= c_pi_N0P_VA(j);
c_eq_p_ciu_N2P_OP('ciu',j).. c_p_N2P_OP('ciu',j) =e= c_pi_N1P_IU(j);
c_eq_pi_N2P_OP(j)..          c_pi_N2P_OP(j)      =e= unit_cost(iuva,c_p_N2P_OP(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j),scale_N2P_OP0(j));
c_eq_pi_N0C_HH(i,h)..        c_pi_N0C_HH(i,h)    =e= unit_cost(dm,c_p_N0C_HH(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),scale_N0C_HH0(i,h));

c_eq_gross_income(h)..          c_gross_income(h)      =e= sum(f,endowment(f,h)*c_pf(f));
c_eq_net_income(h)..            c_net_income(h)        =e= (1-rtp('ty',h))*c_gross_income(h);

c_eq_disposable_income(h)..     c_disposable_income(h) =e= c_net_income(h) + rtp('ttr',h)*c_total_tax_rev;
c_eq_dd_N1C_HH(i,h)..           c_dd_N1C_HH(i,h)       =e= alpha(i,h)*c_disposable_income(h)/c_pi_N0C_HH(i,h);

c_eq_dd_const_N0C_HH(i,h)..     c_dd_const_N0C_HH(i,h) =e= fac_dd_const(dm,c_p_N0C_HH(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h));
c_eq_dd_N0C_HH(dm,i,h)..        c_dd_N0C_HH(dm,i,h) =e= factor_demand(scale_N0C_HH0(i,h),c_dd_N1C_HH(i,h),share_N0C_HH0(dm,i,h),c_p_N0C_HH(dm,i,h),sigma_dm_fu(i,h),c_dd_const_N0C_HH(i,h));

c_eq_q_const_N2P_OP(j)..        c_q_const_N2P_OP(j) =e= fac_dd_const(iuva,c_p_N2P_OP(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j));
c_eq_q_N2P_OP(iuva,j)..         c_q_N2P_OP(iuva,j) =e= factor_demand(scale_N2P_OP0(j),c_qi_N2P_OP(j),share_N2P_OP0(iuva,j),c_p_N2P_OP(iuva,j),sigma_iu_va(j),c_q_const_N2P_OP(j));
c_eq_q_N2P_IU(j)..              c_q_N2P_IU(j)      =e= c_q_N2P_OP('ciu',j);
c_eq_q_N1P_VA(j)..              c_q_N1P_VA(j)      =e= c_q_N2P_OP('cva',j);
c_eq_q_const_N1P_IU(j)..        c_q_const_N1P_IU(j) =e= fac_dd_const(i,c_p_N1P_IU(i,j),share_N1P_IU0(i,j),sigma_iu(j));

c_eq_q_N1P_IU(i,j)..            c_q_N1P_IU(i,j) =e= factor_demand(scale_N1P_IU0(j),c_q_N2P_IU(j),share_N1P_IU0(i,j),c_p_N1P_IU(i,j),sigma_iu(j),c_q_const_N1P_IU(j));
c_eq_q_const_N0P_IU(i,j)..      c_q_const_N0P_IU(i,j) =e= fac_dd_const(dm,c_p_N0P_IU(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j));
c_eq_q_N0P_IU(dm,i,j)..         c_q_N0P_IU(dm,i,j) =e= factor_demand(scale_N0P_IU0(i,j),c_q_N1P_IU(i,j),share_N0P_IU0(dm,i,j),c_p_N0P_IU(dm,i,j),sigma_dm_io(i,j),c_q_const_N0P_IU(i,j));
c_eq_q_const_N0P_VA(j)..        c_q_const_N0P_VA(j) =e= fac_dd_const(f,c_p_N0P_VA(f,j),share_N0P_VA0(f,j),sigma_kl(j));
c_eq_q_N0P_VA(f,j)..            c_q_N0P_VA(f,j) =e= factor_demand(scale_N0P_VA0(j),c_q_N1P_VA(j),share_N0P_VA0(f,j),c_p_N0P_VA(f,j),sigma_kl(j),c_q_const_N0P_VA(j));

c_eq_tot_q_OP(dm,i)..           c_tot_q_OP(dm,i) =e= sum(j,c_q_N0P_IU(dm,i,j)) + sum(h,c_dd_N0C_HH(dm,i,h));
c_eq_pX_OP(i)..                 c_pX_OP(i)     =e= c_pdd(i)*(1+fob0(i))*(1+t_exp0(i));
c_eq_tot_qX_OP(i)..             c_tot_qX_OP(i) =e= exports0(i)*c_pX_OP(i)**exports_eta(i);
c_eq_tot_qD_OP(i)..             c_tot_qD_OP(i) =e= c_tot_q_OP('d',i) + c_tot_qX_OP(i);
c_eq_tot_qM_OP(i)..             c_tot_qM_OP(i) =e= c_tot_q_OP('m',i);

c_eq_ZP_N2P_OP(j)..             c_pdd(j)*c_qi_N2P_OP(j) =g= sum(i,c_pmd_iu(i,j)*c_q_N0P_IU('m',i,j))+sum(i,c_pdd_iu(i,j)*c_q_N0P_IU('d',i,j)) + sum(f,c_p_N0P_VA(f,j)*c_q_N0P_VA(f,j)) ;

c_eq_tot_qI_SS(i)..             c_tot_qI_SS(i) =e= c_qi_N2P_OP(i);
c_eq_tot_qI_DD(i)..             c_tot_qI_DD(i) =e= c_tot_qD_OP(i);
c_eq_I_eqbm_OP(i)..             c_tot_qI_SS(i) =g= c_tot_qI_DD(i);

*factor equilibrium
c_eq_tot_qf_SS(f)..             c_tot_qf_SS(f) =e= sum(h,endowment(f,h));
c_eq_tot_qf_DD(f)..             c_tot_qf_DD(f) =e= sum(j,c_q_N0P_VA(f,j));
c_eq_f_eqbm_OP(f)..             c_tot_qf_SS(f) =g= c_tot_qf_DD(f);

* computing the value of imports, exports and trade deficit
c_eq_tot_val_M..                c_tot_val_M =e= sum((i,j),c_pmd(i)*c_q_N0P_IU('m',i,j)) + sum((i,h), c_pmd(i)*c_dd_N0C_HH('m',i,h));
c_eq_tot_val_X..                c_tot_val_X =e= sum(i,c_pX_OP(i)*c_tot_qX_OP(i));
c_eq_TD..                       c_TD*c_ER =e= (c_tot_val_X - c_tot_val_M);
c_eq_zero_td..                  c_TD =e= 0;

* tax revenue from tariffs levied on IUse of imports
c_eq_tax_M_IU..                 c_tax_M_IU =e= sum((i,j),c_pmd(i)*tm_io(i,j)*c_q_N0P_IU('m',i,j));
*tax revenue from excise and sales taxes on IUse of domestic goods
c_eq_tax_D_IU..                 c_tax_D_IU =e= sum((i,j),c_pdd(i)*td_io(i,j)*c_q_N0P_IU('d',i,j));
* tax revenue from tariffs levied on final consumption of households
c_eq_tax_M_HH..                 c_tax_M_HH =e= sum((i,h),c_pmd(i)*tm_hh(i,h)*c_dd_N0C_HH('m',i,h));
* tax revenue from sales, excise taxes on final consumption of households
c_eq_tax_D_HH..                 c_tax_D_HH =e= sum((i,h),c_pdd(i)*td_hh(i,h)*c_dd_N0C_HH('d',i,h));
* tax revenue from factor taxes
c_eq_tax_F_VA..                 c_tax_F_VA =e= sum((f,j),c_pf(f)*rtf(f,j)*c_q_N0P_VA(f,j));
* tax revenue from income tax on households
c_eq_tax_I_HH..                 c_tax_I_HH =e= sum(h,rtp('ty',h)*c_gross_income(h)) ;
* total tax revenue
c_eq_total_tax_Rev..            c_total_tax_Rev =e= c_tax_M_IU + c_tax_D_IU + c_tax_M_HH + c_tax_D_HH + c_tax_F_VA + c_tax_I_HH;
*utility for HH domestic
c_eq_utility_d(h)..             c_U0_D(h)  =e= utility(i,alpha(i,h),1,c_dd_N0C_HH('d',i,h));
*utility for HH imported
c_eq_utility_m(h)..             c_U0_M(h)  =e= utility(i,alpha(i,h),1,c_dd_N0C_HH('m',i,h));
* total utility
c_eq_U0_total..                 c_U0_total =e= sum(h,c_U0_M(h)) + sum(h,c_U0_D(h));

model c_cge_443_ces
*/all/;
/
c_eq_obj
c_eq_pmd
c_eq_pmd_iu
c_eq_pdd_iu
c_eq_pdd_hh
c_eq_pmd_hh
c_eq_pd_N0P_IU
c_eq_pm_N0P_IU
c_eq_pd_N0C_HH
c_eq_pm_N0C_HH
c_eq_pi_N0P_IU
c_eq_p_N1P_IU
c_eq_pi_N1P_IU
c_eq_p_N0P_VA
c_eq_pi_N0P_VA

c_eq_p_cva_N2P_OP
c_eq_p_ciu_N2P_OP
c_eq_pi_N2P_OP
c_eq_pi_N0C_HH
c_eq_gross_income
c_eq_net_income
c_eq_disposable_income
c_eq_dd_N1C_HH
c_eq_dd_const_N0C_HH
c_eq_dd_N0C_HH

c_eq_q_const_N2P_OP
c_eq_q_N2P_OP
c_eq_q_N2P_IU
c_eq_q_N1P_VA

c_eq_q_const_N1P_IU
c_eq_q_N1P_IU

c_eq_q_const_N0P_IU
c_eq_q_N0P_IU

c_eq_q_const_N0P_VA
c_eq_q_N0P_VA
c_eq_tot_q_OP

c_eq_pX_OP
c_eq_tot_qX_OP

c_eq_tot_qD_OP
c_eq_tot_qM_OP

c_eq_tot_qf_SS
c_eq_tot_qf_DD
c_eq_f_eqbm_OP

c_eq_ZP_N2P_OP
c_eq_tot_qI_SS
c_eq_tot_qI_DD
c_eq_I_eqbm_OP

c_eq_tot_val_M
c_eq_tot_val_X

c_eq_TD
c_eq_zero_td
c_eq_tax_M_IU
c_eq_tax_D_IU
c_eq_tax_M_HH
c_eq_tax_D_HH
c_eq_tax_F_VA
c_eq_tax_I_HH
c_eq_total_tax_Rev
c_eq_utility_d
c_eq_utility_m
c_eq_U0_total
/
;

model c_cge_443_ces_mcp
/
c_eq_obj.c_z
c_eq_pmd.c_pmd
c_eq_pmd_iu.c_pmd_iu
c_eq_pdd_iu.c_pdd_iu
c_eq_pdd_hh.c_pdd_hh
c_eq_pmd_hh.c_pmd_hh
c_eq_pd_N0P_IU.c_p_N0P_IU
c_eq_pm_N0P_IU.c_p_N0P_IU
c_eq_pd_N0C_HH.c_p_N0C_HH
c_eq_pm_N0C_HH.c_p_N0C_HH
c_eq_pi_N0P_IU.c_pi_N0P_IU
c_eq_p_N1P_IU.c_p_N1P_IU
c_eq_pi_N1P_IU.c_pi_N1P_IU
c_eq_p_N0P_VA.c_p_N0P_VA
c_eq_pi_N0P_VA.c_pi_N0P_VA

c_eq_p_cva_N2P_OP.c_p_N2P_OP
c_eq_p_ciu_N2P_OP.c_p_N2P_OP
c_eq_pi_N2P_OP.c_pi_N2P_OP
c_eq_pi_N0C_HH.c_pi_N0C_HH
c_eq_gross_income.c_gross_income
c_eq_net_income.c_net_income
c_eq_disposable_income.c_disposable_income
c_eq_dd_N1C_HH.c_dd_N1C_HH
c_eq_dd_const_N0C_HH.c_dd_const_N0C_HH
c_eq_dd_N0C_HH.c_dd_N0C_HH

c_eq_q_const_N2P_OP.c_q_const_N2P_OP
c_eq_q_N2P_OP.c_q_N2P_OP
c_eq_q_N2P_IU.c_q_N2P_IU
c_eq_q_N1P_VA.c_q_N1P_VA

c_eq_q_const_N1P_IU.c_q_const_N1P_IU
c_eq_q_N1P_IU.c_q_N1P_IU

c_eq_q_const_N0P_IU.c_q_const_N0P_IU
c_eq_q_N0P_IU.c_q_N0P_IU

c_eq_q_const_N0P_VA.c_q_const_N0P_VA
c_eq_q_N0P_VA.c_q_N0P_VA
c_eq_tot_q_OP.c_tot_q_OP

c_eq_pX_OP.c_pX_OP
c_eq_tot_qX_OP.c_tot_qX_OP

c_eq_tot_qD_OP.c_tot_qD_OP
c_eq_tot_qM_OP.c_tot_qM_OP

c_eq_tot_qf_SS.c_tot_qf_SS
c_eq_tot_qf_DD.c_tot_qf_DD
c_eq_f_eqbm_OP.c_pf

c_eq_ZP_N2P_OP.c_qi_N2P_OP
c_eq_tot_qI_SS.c_tot_qI_SS
c_eq_tot_qI_DD.c_tot_qI_DD
c_eq_I_eqbm_OP.c_pdd

c_eq_tot_val_M.c_tot_val_M
c_eq_tot_val_X.c_tot_val_X

c_eq_TD.c_ER
c_eq_zero_td.c_td
c_eq_tax_M_IU.c_tax_M_IU
c_eq_tax_D_IU.c_tax_D_IU
c_eq_tax_M_HH.c_tax_M_HH
c_eq_tax_D_HH.c_tax_D_HH
c_eq_tax_F_VA.c_tax_F_VA
c_eq_tax_I_HH.c_tax_I_HH
c_eq_total_tax_Rev.c_total_tax_Rev
c_eq_utility_d.c_U0_D
c_eq_utility_m.c_U0_M
c_eq_U0_total.c_U0_total
/
;
*** GAMS: Assign Variable Levels-1
c_ER.l                = 1;
c_pdd.l(i)            = 1;
c_pf.l(f)             = 1;
c_pf.fx('L')          = 1;
c_pmd.l(i)             = c_PMBAR*c_ER.l*(1+cif0(i));
c_pmd_iu.l(i,j)        = c_pmd.l(i)*(1+tm_io(i,j));
c_pdd_iu.l(i,j)        = c_pdd.l(i)*(1+td_io(i,j));
c_pmd_hh.l(i,h)        = c_pmd.l(i)*(1+tm_hh(i,h));
c_pdd_hh.l(i,h)        = c_pdd.l(i)*(1+td_hh(i,h));
c_p_N0P_IU.l('d',i,j)        = c_pdd_iu.l(i,j);
c_p_N0P_IU.l('m',i,j)        = c_pmd_iu.l(i,j);
c_p_N0C_HH.l('d',i,h)        = c_pdd_hh.l(i,h);
c_p_N0C_HH.l('m',i,h)        = c_pmd_hh.l(i,h);
c_pi_N0P_IU.l(i,j)        = unit_cost(dm,c_p_N0P_IU.l(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j),scale_N0P_IU0(i,j));
c_p_N1P_IU.l(i,j)         = c_pi_N0P_IU.l(i,j);
c_pi_N1P_IU.l(j)          = unit_cost(i,c_p_N1P_IU.l(i,j),share_N1P_IU0(i,j),sigma_iu(j),scale_N1P_IU0(j));
c_p_N0P_VA.l(f,j)         = c_pf.l(f)*(1+rtf(f,j));
c_pi_N0P_VA.l(j)          = unit_cost(f,c_p_N0P_VA.l(f,j),share_N0P_VA0(f,j),sigma_kl(j),scale_N0P_VA0(j));
c_p_N2P_OP.l('cva',j)     = c_pi_N0P_VA.l(j);
c_p_N2P_OP.l('ciu',j)     = c_pi_N1P_IU.l(j);
c_pi_N2P_OP.l(j)          = unit_cost(iuva,c_p_N2P_OP.l(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j),scale_N2P_OP0(j));
c_pi_N0C_HH.l(i,h)        = unit_cost(dm,c_p_N0C_HH.l(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),scale_N0C_HH0(i,h));

c_pi_N0P_IU.l(i,j)        = unit_cost(dm,c_p_N0P_IU.l(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j),scale_N0P_IU0(i,j));
c_p_N1P_IU.l(i,j)         = c_pi_N0P_IU.l(i,j);
c_pi_N1P_IU.l(j)          = unit_cost(i,c_p_N1P_IU.l(i,j),share_N1P_IU0(i,j),sigma_iu(j),scale_N1P_IU0(j));
c_p_N0P_VA.l(f,j)         = c_pf.l(f)*(1+rtf(f,j));
c_pi_N0P_VA.l(j)          = unit_cost(f,c_p_N0P_VA.l(f,j),share_N0P_VA0(f,j),sigma_kl(j),scale_N0P_VA0(j));
c_p_N2P_OP.l('cva',j)     = c_pi_N0P_VA.l(j);
c_p_N2P_OP.l('ciu',j)     = c_pi_N1P_IU.l(j);
c_pi_N2P_OP.l(j)          = unit_cost(iuva,c_p_N2P_OP.l(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j),scale_N2P_OP0(j));
c_pi_N0C_HH.l(i,h)        = unit_cost(dm,c_p_N0C_HH.l(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h),scale_N0C_HH0(i,h));

c_total_tax_rev.l        = total_tax_revenue0;
c_gross_income.l(h)      = sum(f,endowment(f,h)*c_pf.l(f));
c_net_income.l(h)        = (1-rtp('ty',h))*c_gross_income.l(h);
c_disposable_income.l(h) = c_net_income.l(h) + rtp('ttr',h)*c_total_tax_rev.l;
c_dd_N1C_HH.l(i,h)       = alpha(i,h)*c_disposable_income.l(h)/c_pi_N0C_HH.l(i,h);

c_qi_N2P_OP.l(j) = qi_N2P_OP0(j);
c_tot_qD_OP.l(j) = qi_N2P_OP0(j);
*computing the demand for D & M in final consumption
c_dd_const_N0C_HH.l(i,h) = fac_dd_const(dm,c_p_N0C_HH.l(dm,i,h),share_N0C_HH0(dm,i,h),sigma_dm_fu(i,h));
c_dd_N0C_HH.l(dm,i,h)    = factor_demand(scale_N0C_HH0(i,h),c_dd_N1C_HH.l(i,h),share_N0C_HH0(dm,i,h),c_p_N0C_HH.l(dm,i,h),sigma_dm_fu(i,h),c_dd_const_N0C_HH.l(i,h));

c_q_const_N2P_OP.l(j)  = fac_dd_const(iuva,c_p_N2P_OP.l(iuva,j),share_N2P_OP0(iuva,j),sigma_iu_va(j));
c_q_N2P_OP.l(iuva,j)   = factor_demand(scale_N2P_OP0(j),c_tot_qD_OP.l(j),share_N2P_OP0(iuva,j),c_p_N2P_OP.l(iuva,j),sigma_iu_va(j),c_q_const_N2P_OP.l(j));
c_q_N2P_IU.l(j)        = c_q_N2P_OP.l('ciu',j);
c_q_N1P_VA.l(j)        = c_q_N2P_OP.l('cva',j);
c_q_const_N1P_IU.l(j)  = fac_dd_const(i,c_p_N1P_IU.l(i,j),share_N1P_IU0(i,j),sigma_iu(j));
c_q_N1P_IU.l(i,j)      = factor_demand(scale_N1P_IU0(j),c_q_N2P_IU.l(j),share_N1P_IU0(i,j),c_p_N1P_IU.l(i,j),sigma_iu(j),c_q_const_N1P_IU.l(j));

c_q_const_N0P_IU.l(i,j) = fac_dd_const(dm,c_p_N0P_IU.l(dm,i,j),share_N0P_IU0(dm,i,j),sigma_dm_io(i,j));
c_q_N0P_IU.l(dm,i,j)    = factor_demand(scale_N0P_IU0(i,j),c_q_N1P_IU.l(i,j),share_N0P_IU0(dm,i,j),c_p_N0P_IU.l(dm,i,j),sigma_dm_io(i,j),c_q_const_N0P_IU.l(i,j));

c_q_const_N0P_VA.l(j)  = fac_dd_const(f,c_p_N0P_VA.l(f,j),share_N0P_VA0(f,j),sigma_kl(j));
c_q_N0P_VA.l(f,j)      = factor_demand(scale_N0P_VA0(j),c_q_N1P_VA.l(j),share_N0P_VA0(f,j),c_p_N0P_VA.l(f,j),sigma_kl(j),c_q_const_N0P_VA.l(j));

c_pX_OP.l(i)        = c_pdd.l(i)*(1+fob0(i))*(1+t_exp0(i));
c_tot_qX_OP.l(i)    = exports0(i)*c_pX_OP.l(i)**exports_eta(i);

* computing the total demand for D & M
c_tot_q_OP.l(dm,i)      = sum(j,c_q_N0P_IU.l(dm,i,j)) + sum(h,c_dd_N0C_HH.l(dm,i,h)) ;
c_tot_qD_OP.l(i)        = c_tot_q_OP.l('d',i) + c_tot_qX_OP.l(i);
c_tot_qM_OP.l(i)        = c_tot_q_OP.l('m',i);

c_tot_qf_SS.l(f)        = sum(h,endowment(f,h));
c_tot_qf_DD.l(f)        = sum(j,c_q_N0P_VA.l(f,j));

c_tot_qI_SS.l(i)        = c_qi_N2P_OP.l(i);
c_tot_qI_DD.l(i)        = c_tot_qD_OP.l(i);

c_tot_val_M.l        = sum((i,j),c_pmd.l(i)*c_q_N0P_IU.l('m',i,j)) + sum((i,h), c_pmd.l(i)*c_dd_N0C_HH.l('m',i,h));
c_tot_val_X.l        = sum(i,c_pX_OP.l(i)*c_tot_qX_OP.l(i));

c_tax_M_IU.l        = sum((i,j),c_pmd.l(i)*tm_io(i,j)*c_q_N0P_IU.l('m',i,j));
c_tax_D_IU.l        = sum((i,j),c_pdd.l(i)*td_io(i,j)*c_q_N0P_IU.l('d',i,j));
c_tax_M_HH.l        = sum((i,h),c_pmd.l(i)*tm_hh(i,h)*c_dd_N0C_HH.l('m',i,h));
c_tax_D_HH.l        = sum((i,h),c_pdd.l(i)*td_hh(i,h)*c_dd_N0C_HH.l('d',i,h));
c_tax_F_VA.l        = sum((f,j),c_pf.l(f)*rtf(f,j)*c_q_N0P_VA.l(f,j));
c_tax_I_HH.l        = sum(h,rtp('ty',h)*c_gross_income.l(h)) ;
c_total_tax_Rev.l   = c_tax_M_IU.l + c_tax_D_IU.l + c_tax_M_HH.l + c_tax_D_HH.l + c_tax_F_VA.l + c_tax_I_HH.l;
c_U0_D.l(h)           = utility(i,alpha(i,h),1,c_dd_N0C_HH.l('d',i,h));
c_U0_M.l(h)           = utility(i,alpha(i,h),1,c_dd_N0C_HH.l('m',i,h));
c_U0_total.l         = sum(h,c_U0_M.l(h)) + sum(h,c_U0_D.l(h));

*** GAMS: Assign Variable Bounds-1
c_ER.lo                = 0.001*c_ER.l    ;
c_pdd.lo(i)            = 0.001*c_pdd.l(i);
c_pf.lo(f)             = 0.001*c_pf.l(f) ;

c_pmd.lo(i)             = 0.001*c_pmd.l(i)     ;
c_pmd_iu.lo(i,j)        = 0.001*c_pmd_iu.l(i,j);
c_pdd_iu.lo(i,j)        = 0.001*c_pdd_iu.l(i,j);
c_pmd_hh.lo(i,h)        = 0.001*c_pmd_hh.l(i,h);
c_pdd_hh.lo(i,h)        = 0.001*c_pdd_hh.l(i,h);

c_p_N0P_IU.lo('d',i,j)        = 0.001*c_p_N0P_IU.l('d',i,j);
c_p_N0P_IU.lo('m',i,j)        = 0.001*c_p_N0P_IU.l('m',i,j);
c_p_N0C_HH.lo('d',i,h)        = 0.001*c_p_N0C_HH.l('d',i,h);
c_p_N0C_HH.lo('m',i,h)        = 0.001*c_p_N0C_HH.l('m',i,h);

c_pi_N0P_IU.lo(i,j)      = 0.001*c_pi_N0P_IU.l(i,j)   ;
c_p_N1P_IU.lo(i,j)       = 0.001*c_p_N1P_IU.l(i,j)    ;
c_pi_N1P_IU.lo(j)        = 0.001*c_pi_N1P_IU.l(j)     ;
c_p_N0P_VA.lo(f,j)       = 0.001*c_p_N0P_VA.l(f,j)    ;
c_pi_N0P_VA.lo(j)        = 0.001*c_pi_N0P_VA.l(j)     ;
c_p_N2P_OP.lo('cva',j)   = 0.001*c_p_N2P_OP.l('cva',j);
c_p_N2P_OP.lo('ciu',j)   = 0.001*c_p_N2P_OP.l('ciu',j);
c_pi_N2P_OP.lo(j)        = 0.001*c_pi_N2P_OP.l(j)     ;
c_pi_N0C_HH.lo(i,h)      = 0.001*c_pi_N0C_HH.l(i,h)   ;

c_total_tax_rev.lo        = 0.001*c_total_tax_rev.l    ;
c_gross_income.lo(h)      = 0.001*c_gross_income.l(h)      ;
c_net_income.lo(h)        = 0.001*c_net_income.l(h)        ;
c_disposable_income.lo(h) = 0.001*c_disposable_income.l(h) ;
c_dd_N1C_HH.lo(i,h)       = 0.001*c_dd_N1C_HH.l(i,h)       ;
c_dd_const_N0C_HH.lo(i,h) = 0.001*c_dd_const_N0C_HH.l(i,h);
c_dd_N0C_HH.lo(dm,i,h)    = 0.001*c_dd_N0C_HH.l(dm,i,h);

c_q_const_N2P_OP.lo(j)    = 0.001*c_q_const_N2P_OP.l(j);
c_q_N2P_OP.l(iuva,j)      = 0.001*c_q_N2P_OP.l(iuva,j);
c_q_N2P_IU.lo(j)          = 0.001*c_q_N2P_IU.l(j);
c_q_N1P_VA.lo(j)          = 0.001*c_q_N1P_VA.l(j);
c_q_const_N1P_IU.lo(j)    = 0.001*c_q_const_N1P_IU.l(j);
c_q_N1P_IU.l(i,j)         = 0.001*c_q_N1P_IU.l(i,j);

c_q_const_N0P_IU.lo(i,j)     = 0.001*c_q_const_N0P_IU.l(i,j);
c_q_N0P_IU.lo(dm,i,j)        = 0.001*c_q_N0P_IU.l(dm,i,j);
c_q_const_N0P_VA.lo(j)       = 0.001*c_q_const_N0P_VA.l(j);
c_q_N0P_VA.lo(f,j)           = 0.001*c_q_N0P_VA.l(f,j);

c_tot_q_OP.lo(dm,i)              = 0.001*c_tot_q_OP.l(dm,i);
c_tot_qD_OP.lo(i)                = 0.001*c_tot_qD_OP.l(i);
c_tot_qM_OP.lo(i)                = 0.001*c_tot_qM_OP.l(i);
c_pX_OP.lo(i)                    = 0.001*c_pX_OP.l(i);
c_tot_qX_OP.lo(i)                = 0.001*c_tot_qX_OP.l(i);

c_tot_qf_SS.lo(f)                = 0.001*c_tot_qf_SS.l(f);
c_tot_qf_DD.lo(f)                = 0.001*c_tot_qf_DD.l(f) ;
c_tot_val_M.lo                   = 0.001*c_tot_val_M.l;
c_tot_val_X.lo                   = 0.001*c_tot_val_X.l;
c_U0_D.lo(h)            = 0.001*c_U0_D.l(h);
c_U0_M.lo(h)            = 0.001*c_U0_M.lo(h);
c_U0_total.lo           = 0.001*c_U0_total.l;

scalar c_solver_choice /1/;
if (c_solver_choice eq 0,
   solve c_cge_443_ces using nlp minimising c_z;
else
   SOLVE c_cge_443_ces_mcp USING MCP;
);

put /;
put '--------------- After Imposing Tariffs on Imports ----------------------';
put/; put/;

put 'Total Quantities'  put /;
put '----------------'  put /;
put '         ' put  'Domestic      ' put 'Imported      ' put 'Exports  ' put /;
put 'c1':<9,    put c_tot_qD_OP.l('c1'):<12, put '  ',
                put c_tot_qM_OP.l('c1'):<12, put '  ',
                put c_tot_qX_OP.l('c1'):<12, put /;

put /;
put 'Total Value'  put /;
put '----------------'  put /;
put  'Imports      ' put 'Exports' put /;
put c_tot_val_M.l:<9   put ' ' put c_tot_val_X.l put /;

put /;
put 'Collected Taxes on Imported Goods'  put /;
put '---------------------------------'  put /;
put  'IUse          ' put 'Households' put /;
put c_tax_M_IU.l:<9     put ' ' put c_tax_M_HH.l put /;

put /;
put 'Total HH Utility from Consumption'  put /;
put '------------------------------------------'  put /;
put  'Domestic Goods      ' put 'Imported Goods' put /;
put sum(h,c_U0_D.l(h)):<20    put sum(h,c_U0_M.l(h)):<9 put /;


put /;
put 'Price Index for Agriculture'  put /;
put '---------------------------'  put /;
put 'c1':<15,    put c_pi_N2P_OP.l('c1'); put /;
put /;

parameters
EV(h)                            equivalent variation
CV(h)                            compensating variation
;
* equivalent variation
EV(h)   =     ((c_U0_D.l(h)+c_U0_M.l(h))-(U0_D.l(h)+U0_M.l(h)))/(U0_D.l(h)+U0_M.l(h))*gross_income.l(h);
* compensating variation
CV(h)   =     ((c_U0_D.l(h)+c_U0_M.l(h))-(U0_D.l(h)+U0_M.l(h)))/(c_U0_D.l(h)+c_U0_M.l(h))*c_gross_income.l(h);

put /;
put '------------------ Equivalent/Compensating Variation --------------------------';
put /;
put '      ', loop(h, put h.tl:>10) put /;
put 'EV    ', loop(h, put EV(h):>10:3) put/;
put 'CV    ', loop(h, put CV(h):>10:3) put/;