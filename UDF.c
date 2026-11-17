#define a1       312.7976
#define a2      -1162.2    
#define a3      1765.5
#define a4      -1406.4     
#define a5      623.0262
#define a6      -148.0072
#define a7      15.9098
#define a8      -0.3372
#define I       3
#define H       3
#define V       0.000024245

#define b1    -69.8743
#define b2    285.7268
#define b3    -482.3356
#define b4    433.8266
#define b5    -224.0125
#define b6    66.6496
#define b7    -10.7862
#define b8    0.8326

DEFINE_SOURCE(heat_source_ncm,c,t,dS,eqn)
{
    real soc, resistance, T_K;
    real source, ts;
    
    T_K = C_T(c,t);
    ts  = RP_Get_Real("flow-time");
    soc = 1.0 - (ts * I) / (3600.0 * H);
    
    resistance = a1*pow(soc,7) + a2*pow(soc,6) + a3*pow(soc,5) 
               + a4*pow(soc,4) + a5*pow(soc,3) + a6*pow(soc,2) 
               + a7*soc + a8;

    source = (I*I*resistance + I*T_K*0.0003) / V; 
    dS[eqn] = 0.0;
    
    return source;
}

DEFINE_SOURCE(heat_source_lfp,c,t,dS,eqn)
{
    real soc, resistance, T_K;
    real source, ts;
    
    T_K = C_T(c,t);
    ts  = RP_Get_Real("flow-time");
    soc = 1.0 - (ts * I) / (3600.0 * H);
    
    resistance = b1*pow(soc,7) + b2*pow(soc,6) + b3*pow(soc,5) 
               + b4*pow(soc,4) + b5*pow(soc,3) + b6*pow(soc,2) 
               + b7*soc + b8;

    source = (I*I*resistance + I*T_K*0.0003) / V; 
    dS[eqn] = 0.0;
    
    return source;
}
