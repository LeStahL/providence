#version 130
const float PI = radians(180.);
const float TAU = 2.*PI;
float clip(float a) { return clamp(a,-1.,1.); }
float smstep(float a, float b, float x) {return smoothstep(a, b, clamp(x, a, b));}
float theta(float x) { return smstep(0.,1e-3,x); }
float _sin(float a) { return sin(TAU * mod(a,1.)); }
float _sin_(float a, float p) { return sin(TAU * mod(a,1.) + p); }
float _sq_(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
float freqC1(float note){ return 32.7 * exp2(note/12.); }
float minus1hochN(int n) { return (1. - 2.*float(n % 2)); }
float minus1hochNminus1halbe(int n) { return sin(.5*PI*float(n)); }
float pseudorandom(float x) { return fract(sin(dot(vec2(x),vec2(12.9898,78.233))) * 43758.5453); }
float fhelp(float x) { return 1. + .333*x; } // 1. + .33333*x + .1*x*x + .02381*x*x*x + .00463*x*x*x*x;
float linmix(float x, float a, float b, float y0, float y1) { return mix(y0,y1,clamp(a*x+b,0.,1.)); }
float s_atan(float a) { return .636 * atan(a); }
float doubleslope(float t, float a, float d, float s) { return smstep(-.00001,a,t) - (1.-s) * smstep(0.,d,t-a); }

#define NTIME 12
const float pos_B[12] = float[12](0.,56.,58.,60.,62.,64.,74.,78.,82.,86.,90.,97.);
const float pos_t[12] = float[12](0.,108.3871,112.1966,115.9466,119.6389,123.2753,141.1857,148.4585,155.8431,163.5059,171.3748,185.3748);
const float pos_BPS[11] = float[11](.5167,.525,.5333,.5417,.55,.5583,.55,.5417,.522,.5083,.5);
const float pos_SPB[11] = float[11](1.9354,1.9048,1.8751,1.846,1.8182,1.7912,1.8182,1.846,1.9157,1.9673,2.);
float BPS, SPB, BT;

float Tsample;

#define filterthreshold 1.e-3

//TEXCODE

float drop_phase(float time, float t1, float f0, float f1)
{
    float t = min(time, t1);
    float phi = f0*t + .5*(f1-f0)/t1*t*t;

    if(time > t1)
    {
        phi += f1 * (time - t1);
    }
    return phi;
}

float metalnoise(float t, float fac1, float fac2)
{
    return .666*pseudorandom(t) - 1.333*pseudorandom(t-Tsample) + .333*(pseudorandom(t+fac1*Tsample)+pseudorandom(t+fac2*Tsample));
}

float lpnoise(float t, float fq)
{
    t *= fq;
    float tt = fract(t);
    float tn = t - tt;
    return mix(pseudorandom(floor(tn) / fq), pseudorandom(floor(tn + 1.0) / fq), smstep(0.0, 1.0, tt));
}

float reverb_phase(float t, float amt)
{
    float r = lpnoise(t, 100.0) + 0.2*lpnoise(t, 550.0) + 0.1*lpnoise(t, 1050.0)*exp(-5.*t);
    return amt * r;
}

float env_AHDSR(float x, float L, float A, float H, float D, float S, float R)
{
    return (x<A ? x/A : x<A+H ? 1. : x<A+H+D ? (1. - (1.-S)*(x-H-A)/D) : x<=L-R ? S : x<=L ? S*(L-x)/R : 0.);
}

float env_AHDSRexp(float x, float L, float A, float H, float D, float S, float R)
{
    float att = pow(x/A,8.);
    float dec = S + (1.-S) * exp(-(x-H-A)/D);
    float rel = (x <= L-R) ? 1. : pow((L-x)/R,4.);
    return (x < A ? att : x < A+H ? 1. : dec) * rel;
}

float sinshape(float x, float amt, float parts)
{
    return (1.-amt) * x + amt * sign(x) * 0.5 * (1. - cos(parts*PI*x));
}

float comp_SAW(int N, float inv_N, float PW) {return inv_N * (1. - _sin(float(N)*PW));}
float comp_TRI(int N, float inv_N, float PW) {return N % 2 == 0 ? .1 * inv_N * _sin(float(N)*PW) : inv_N * inv_N * (1. - _sin(float(N)*PW));}
float comp_SQU(int N, float inv_N, float PW) {return inv_N * (minus1hochN(N) * _sin(.5*float(N)*PW + .25) - 1.);}
float comp_HAE(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (1. - minus1hochNminus1halbe(N))*_sin(PW);}
float comp_OBO(int N, float inv_N, float PW) {return sqrt(inv_N) * (1. + _sin(float(N)*(1.5+PW)+.5*PI));}

float MADD(float t, float f, float p0, int NMAX, int NINC, float MIX, float CO, float NDECAY, float RES, float RES_Q, float DET, float PW, float LOWCUT, int keyF)
{
    float ret = 0.;
    float INR = keyF==1 ? 1./CO : f/CO;
    float IRESQ = keyF==1 ? 1./RES_Q : 1./(RES_Q*f);

    float p = f*t;
    float float_N, inv_N, comp_mix, filter_N;
    for(int N = 1 + int(LOWCUT/f - 1.e-3); N<=NMAX; N+=NINC)
    {
        float_N = float(N);
        inv_N = 1./float_N;
        comp_mix = MIX < -1. ? (MIX+2.) * comp_SAW(N,inv_N,PW)  - (MIX+1.) * comp_OBO(N,inv_N,PW)
                 : MIX <  0. ? (MIX+1.) * comp_TRI(N,inv_N,PW)  -     MIX  * comp_SAW(N,inv_N,PW)
                 : MIX < 1. ? (1.-MIX) * comp_TRI(N,inv_N,PW)  +     MIX  * comp_SQU(N,inv_N,PW)
                            : (MIX-1.) * comp_HAE(N,inv_N,PW)  + (2.-MIX) * comp_SQU(N,inv_N,PW);

        if(abs(comp_mix) < 1e-4) continue;

        filter_N = pow(1. + pow(float_N*INR,NDECAY),-.5) + RES * exp(-pow((float_N*f-CO)*IRESQ,2.));

        ret += comp_mix * filter_N * (_sin_(float_N * p, p0) + _sin_(float_N * p * (1.+DET), p0));
    }
    return s_atan(ret);
}

float MADD(float t, float f, float p0, int NMAX, int NINC, float MIX, float CO, float NDECAY, float RES, float RES_Q, float DET, float PW, int keyF)
{
    return MADD(t, f, p0, NMAX, NINC, MIX, CO, NDECAY, RES, RES_Q, DET, PW, 0., keyF);
}

float QFM_FB(float PH, float FB) // my guessing of feedback coefficients, FB>0 'saw', FB<0 'sq'
{
    if(FB > 0.) return abs(FB) * .8*_sin(PH + .35*_sin(PH));
    else return abs(FB) * _sin(PH + .5*PI);
}

float QFM(float t, float f, float phase, float LV1, float LV2, float LV3, float LV4, float FR1, float FR2, float FR3, float FR4, float FB1, float FB2, float FB3, float FB4, float ALGO)
{
    int iALGO = int(ALGO);
    float PH1 = FR1 * f * t + phase;
    float PH2 = FR2 * f * t + phase;
    float PH3 = FR3 * f * t + phase;
    float PH4 = FR4 * f * t + phase;

    float LINK41 = 0., LINK42 = 0., LINK43 = 0., LINK32 = 0., LINK31 = 0., LINK21 = 0.;
    if(iALGO == 1)       {LINK43 = 1.; LINK32 = 1.; LINK21 = 1.;}
    else if(iALGO == 2)  {LINK42 = 1.; LINK32 = 1.; LINK21 = 1.;}
    else if(iALGO == 3)  {LINK41 = 1.; LINK32 = 1.; LINK21 = 1.;}
    else if(iALGO == 4)  {LINK42 = 1.; LINK43 = 1.; LINK31 = 1.; LINK21 = 1.;}
    else if(iALGO == 5)  {LINK41 = 1.; LINK31 = 1.; LINK21 = 1.;}
    else if(iALGO == 6)  {LINK43 = 1.; LINK32 = 1.;}
    else if(iALGO == 7)  {LINK43 = 1.; LINK32 = 1.; LINK31 = 1.;}
    else if(iALGO == 8)  {LINK21 = 1.; LINK43 = 1.;}
    else if(iALGO == 9)  {LINK43 = 1.; LINK42 = 1.; LINK41 = 1.;}
    else if(iALGO == 10) {LINK43 = 1.; LINK42 = 1.;}
    else if(iALGO == 11) {LINK43 = 1.;}

    float OP4 = LV4 * _sin(PH4 + QFM_FB(PH4, FB4));
    float OP3 = LV3 * _sin(PH3 + QFM_FB(PH3, FB3) + LINK43*OP4);
    float OP2 = LV2 * _sin(PH2 + QFM_FB(PH2, FB2) + LINK42*OP4 + LINK32*OP3);
    float OP1 = LV1 * _sin(PH1 + QFM_FB(PH1, FB1) + LINK41*OP4 + LINK31*OP3 + LINK32*OP2);

    float sum = OP1;
    if(LINK21 > 0.) sum += OP2;
    if(LINK31 + LINK32 > 0.) sum += OP3;
    if(LINK41 + LINK42 + LINK43 > 0.) sum += OP4;

    return s_atan(sum);
}

float protokick(float t, float f_start, float f_end, float fdecay, float hold, float decay, float drive, float detune, float rev_amount, float rev_hold, float rev_decay, float rev_drive)
{
    float phi = drop_phase(t, fdecay, f_start, f_end);
    float rev_phi = phi + reverb_phase(t, 1.);
    return clamp(drive*.5*(_sin(phi)+_sin((1.-detune)*phi)),-1.,1.) * exp(-max(t-hold, 0.)/decay)
         + rev_amount*clamp(rev_drive*.5*(_sin(rev_phi)+_sin((1.-detune)*rev_phi)),-1.,1.) * exp(-max(t-rev_hold, 0.)/rev_decay);
}

float _BOOMENV0(float t){return t <=.019? linmix(t,52.6316,0.,0.,1.):t <=.244? linmix(t,4.4444,-.0844,1.,.274):t <=1.2? linmix(t,1.046,-.2552,.274,0.):0.;}
float _BOOMENV1(float t){return t <=.147? linmix(t,6.8027,0.,0.,1.):t <=.488? linmix(t,2.9326,-.4311,1.,.417):t <=1.058? linmix(t,1.7544,-.8561,.417,0.):0.;}
float _BOOMENV2(float t){return t <=.259? linmix(t,3.861,0.,0.,1.):t <=.868? linmix(t,1.642,-.4253,1.,.464):t <=1.947? linmix(t,.9268,-.8044,.464,0.):0.;}
float _BOOMENV3(float t){return t <=.62? linmix(t,1.6129,0.,0.,1.):t <=1.899? linmix(t,.7819,-.4848,1.,0.):0.;}
float fqmbace5_volume(float B)
{
    return B<0 ? 0. : (B>=0. && B<2.) ? 0. : (B>=2. && B<6.) ? linmix(B, .25, -.5, 0.0, 1.0) : (B>=14. && B<15.5) ? linmix(B, .6667, -9.3333, 1.0, 0.1) : (B>=15.5 && B<16.) ? linmix(B, 2., -31., 0.1, 0.05) : (B>=16. && B<32.) ? 1.1 : 1.;
}
float fqmbace7sat_vol(float B)
{
    return B<0 ? 0. : (B>=0. && B<2.5) ? linmix(B, .4, 0., 0.0, 1.0) : (B>=15. && B<16.) ? linmix(B, 1., -15., 1.0, 0.1) : (B>=80. && B<88.) ? 1.5 : (B>=88. && B<98.) ? linmix(B, .1, -8.8, 1.5, 0.0) : 1.;
}
float fett2_volume(float B)
{
    return B<0 ? 0. : (B>=0. && B<5.) ? linmix(B, .2, 0., 0.1, 0.9) : (B>=8. && B<11.) ? linmix(B, .3333, -2.6667, 0.9, 0.05) : (B>=11. && B<14.) ? linmix(B, .3333, -3.6667, 0.4, 1.0) : (B>=14. && B<16.) ? linmix(B, .5, -7., 1.0, 0.1) : (B>=56. && B<62.) ? linmix(B, .1667, -9.3333, 1.0, 0.66) : (B>=62. && B<67.) ? linmix(B, .2, -12.4, 0.66, 1.0) : (B>=80. && B<84.) ? linmix(B, .25, -20., 0.76, 1.0) : (B>=91. && B<98.) ? linmix(B, .1429, -13., 1.0, 0.0) : 1.;
}
float hoboe_vol(float B)
{
    return B<0 ? 0. : (B>=0. && B<4.) ? linmix(B, .25, 0., 0.21, 1.1) : (B>=4. && B<15.) ? 1.1 : (B>=15. && B<16.) ? linmix(B, 1., -15., 1.1, 0.4) : 1.;
}

uniform float iBlockOffset;
uniform float iSampleRate;
uniform float iTexSize;
uniform sampler2D iSequence;
uniform float iSequenceWidth;

// Read short value from texture at index off
float rshort(in float off)
{
    float hilo = mod(off, 2.);
    off = .5*off;
    vec2 ind = vec2(mod(off, iSequenceWidth), floor(off/iSequenceWidth));
    vec4 block = texelFetch(iSequence, ivec2(ind), 0);
    vec2 data = mix(block.rg, block.ba, hilo);
    return round(dot(vec2(255., 65280.), data));
}

// Read float value from texture at index off
float rfloat(int off)
{
    float d = rshort(float(off));
    float sign = floor(d/32768.),
        exponent = floor(d*9.765625e-4 - sign*32.),
        significand = d-sign*32768.-exponent*1024.;

    if(exponent == 0.)
         return mix(1., -1., sign) * 5.960464477539063e-08 * significand;
    return mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);
}

#define NTRK 13
#define NMOD 137
#define NPTN 28
#define NNOT 879
#define NDRM 49

int trk_sep(int index)      {return int(rfloat(index));}
int trk_syn(int index)      {return int(rfloat(index+1+1*NTRK));}
float trk_norm(int index)   {return     rfloat(index+1+2*NTRK);}
float trk_rel(int index)    {return     rfloat(index+1+3*NTRK);}
float trk_pre(int index)    {return     rfloat(index+1+4*NTRK);}
float trk_slide(int index)  {return     rfloat(index+1+5*NTRK);} // idea for future: change to individual note_slide_time
float mod_on(int index)     {return     rfloat(index+1+6*NTRK);}
float mod_off(int index)    {return     rfloat(index+1+6*NTRK+1*NMOD);}
int mod_ptn(int index)      {return int(rfloat(index+1+6*NTRK+2*NMOD));}
float mod_transp(int index) {return     rfloat(index+1+6*NTRK+3*NMOD);}
int ptn_sep(int index)      {return int(rfloat(index+1+6*NTRK+4*NMOD));}
float note_on(int index)    {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN);}
float note_off(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+1*NNOT);}
float note_pitch(int index) {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+2*NNOT);}
float note_pan(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+3*NNOT);}
float note_vel(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+4*NNOT);}
float note_slide(int index) {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+5*NNOT);}
float note_aux(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+6*NNOT);}
float drum_rel(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+7*NNOT);}

vec2 mainSynth(float time)
{
    float sL = 0.;
    float sR = 0.;
    float dL = 0.;
    float dR = 0.;

    if (time > 189.3748) return vec2(0);
    
    int _it;
    for(_it = 0; _it < NTIME - 2 && pos_t[_it + 1] < time; _it++);
    BPS = pos_BPS[_it];
    SPB = pos_SPB[_it];
    BT = pos_B[_it] + (time - pos_t[_it]) * BPS;

    float time2 = time - .0002;
    float sidechain = 1.;

    float amaysynL, amaysynR, amaydrumL, amaydrumR, B, Bon, Boff, Bprog, Bproc, L, tL, _t, _t2, vel, rel, pre, f, amtL, amtR, env, slide, aux;
    int tsep0, tsep1, _modU, _modL, ptn, psep0, psep1, _noteU, _noteL, syn, drum;

    for(int trk = 0; trk < NTRK; trk++)
    {
        tsep0 = trk_sep(trk);
        tsep1 = trk_sep(trk + 1);

        syn = trk_syn(trk);
        rel = trk_rel(trk);
        pre = trk_pre(trk);

        for(_modU = tsep0; (_modU < tsep1 - 1) && (BT > mod_on(_modU + 1) - pre); _modU++);
        for(_modL = tsep0; (_modL < tsep1 - 1) && (BT >= mod_off(_modL) + rel); _modL++);

        for(int _mod = _modL; _mod <= _modU; _mod++)
        {
            B = BT - mod_on(_mod);

            ptn   = mod_ptn(_mod);
            psep0 = ptn_sep(ptn);
            psep1 = ptn_sep(ptn + 1);

            for(_noteU = psep0; (_noteU < psep1 - 1) && (B > note_on(_noteU + 1) - pre); _noteU++);
            for(_noteL = psep0; (_noteL < psep1 - 1) && (B >= note_off(_noteL) + rel); _noteL++);

            for(int _note = _noteL; _note <= _noteU; _note++)
            {
                if(syn == 121)
                {
                    drum = int(note_pitch(_note));
                    rel = drum_rel(drum);
                }

                amaysynL  = 0.;
                amaysynR  = 0.;
                amaydrumL = 0.;
                amaydrumR = 0.;

                Bon   = note_on(_note) - pre;
                Boff  = note_off(_note) + rel;
                L     = Boff - Bon;
                tL    = L * SPB;
                Bprog = B - Bon;
                Bproc = Bprog / L;
                _t    = Bprog * SPB;
                _t2   = _t - .0002;
                vel   = note_vel(_note);
                amtL  = clamp(1. - note_pan(_note), 0., 1.);
                amtR  = clamp(1. + note_pan(_note), 0., 1.);
                slide = note_slide(_note);
                aux   = note_aux(_note);

                if(syn == 121)
                {
                    env = trk_norm(trk) * theta(Bprog) * theta(L - Bprog);
                    if(drum == 0) { sidechain = min(sidechain, 1. - vel * (clamp(1.e4 * Bprog,0.,1.) - pow(Bprog/(L-rel),8.)));}
                    else if(drum == 1){
                        amaydrumL = vel*vel*.5*exp(-18.*max(_t-.01,0.))*metalnoise(.6*_t, 1., 2.)
      +.5*(lpnoise(_t,10000.)*smstep(0.,.01,_t)*(1.-(1.-.13)*smstep(0.,.12,_t-.01))-.3*(1.00*lpnoise((_t-0.00),10000.)*smstep(0.,.01,(_t-0.00))*(1.-(1.-.13)*smstep(0.,.12,(_t-0.00)-.01))+6.10e-01*lpnoise((_t-1.20e-03),10000.)*smstep(0.,.01,(_t-1.20e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-1.20e-03)-.01))+3.72e-01*lpnoise((_t-2.40e-03),10000.)*smstep(0.,.01,(_t-2.40e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-2.40e-03)-.01))))*exp(-3.*max(_t-.2,0.))
      +2.1*vel*fract(sin(_t*100.*1.)*50000.*1.)*doubleslope(_t,0.,.01,.1)*exp(-10.*Bprog);
                        amaydrumR = vel*vel*.5*exp(-18.*max(_t2-.01,0.))*metalnoise(.6*_t2, 1., 2.)
      +.5*(lpnoise(_t,10000.)*smstep(0.,.01,_t)*(1.-(1.-.13)*smstep(0.,.12,_t-.01))-.3*(1.00*lpnoise((_t-0.00),10000.)*smstep(0.,.01,(_t-0.00))*(1.-(1.-.13)*smstep(0.,.12,(_t-0.00)-.01))+6.10e-01*lpnoise((_t-1.20e-03),10000.)*smstep(0.,.01,(_t-1.20e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-1.20e-03)-.01))+3.72e-01*lpnoise((_t-2.40e-03),10000.)*smstep(0.,.01,(_t-2.40e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-2.40e-03)-.01))))*exp(-3.*max(_t2-.2,0.))
      +2.1*vel*fract(sin(_t2*100.*1.)*50000.*1.)*doubleslope(_t2,0.,.01,.1)*exp(-10.*Bprog);
                    }
                    else if(drum == 2){
                        amaydrumL = vel*vel*1.3*exp(-7.*max(_t-.05+10.*vel,0.))*metalnoise(.6*_t, .5, 2.)
      +vel*1.3*(lpnoise(_t,10000.)*smstep(0.,.01,_t)*(1.-(1.-.13)*smstep(0.,.12,_t-.01))-.3*(1.00*lpnoise((_t-0.00),10000.)*smstep(0.,.01,(_t-0.00))*(1.-(1.-.13)*smstep(0.,.12,(_t-0.00)-.01))+6.10e-01*lpnoise((_t-1.20e-03),10000.)*smstep(0.,.01,(_t-1.20e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-1.20e-03)-.01))+3.72e-01*lpnoise((_t-2.40e-03),10000.)*smstep(0.,.01,(_t-2.40e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-2.40e-03)-.01))))*exp(-4.*max(_t-.25,0.));
                        amaydrumR = vel*vel*1.3*exp(-7.*max(_t2-.05+10.*vel,0.))*metalnoise(.6*_t2, .5, 2.)
      +vel*1.3*(lpnoise(_t,10000.)*smstep(0.,.01,_t)*(1.-(1.-.13)*smstep(0.,.12,_t-.01))-.3*(1.00*lpnoise((_t-0.00),10000.)*smstep(0.,.01,(_t-0.00))*(1.-(1.-.13)*smstep(0.,.12,(_t-0.00)-.01))+6.10e-01*lpnoise((_t-1.20e-03),10000.)*smstep(0.,.01,(_t-1.20e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-1.20e-03)-.01))+3.72e-01*lpnoise((_t-2.40e-03),10000.)*smstep(0.,.01,(_t-2.40e-03))*(1.-(1.-.13)*smstep(0.,.12,(_t-2.40e-03)-.01))))*exp(-4.*max(_t2-.25,0.));
                    }
                    else if(drum == 11){
                        amaydrumL = vel*fract(sin(_t*100.*.5)*50000.*.5)*doubleslope(_t,0.,.03,.1)*exp(-13.*Bprog);
                        amaydrumR = vel*fract(sin(_t2*100.*.5)*50000.*.5)*doubleslope(_t2,0.,.03,.1)*exp(-13.*Bprog);
                    }
                    else if(drum == 28){
                        amaydrumL = vel*(_BOOMENV0(_t)*(1.6*(.5*_sin_(drop_phase(_t,.046,273.108,73.441),1.7*_BOOMENV0(_t)*lpnoise(_t,452.603))+.5*_sin_(.983*drop_phase(_t,.046,273.108,73.441),1.7*_BOOMENV0(_t)*lpnoise(_t,452.603))))+_BOOMENV1(_t)*(1.25*lpnoise(_t,452.603))+_BOOMENV2(_t)*(1.45*lpnoise(_t,212.15))+_BOOMENV3(_t)*(1.85*lpnoise(_t,141.624)));
                        amaydrumR = vel*(_BOOMENV0(_t)*(1.6*(.5*_sin_(drop_phase((_t-.00089),.046,273.108,73.441),1.7*_BOOMENV0(_t)*lpnoise((_t-.00089),452.603))+.5*_sin_(.983*drop_phase((_t-.00089),.046,273.108,73.441),1.7*_BOOMENV0(_t)*lpnoise((_t-.00089),452.603))))+_BOOMENV1(_t)*(1.25*lpnoise((_t-.00011),452.603))+_BOOMENV2(_t)*(1.45*lpnoise((_t-.00173),212.15))+_BOOMENV3(_t)*(1.85*lpnoise((_t-.002),141.624)));
                    }
                    else if(drum == 32){
                        amaydrumL = vel*protokick(_t,242.,55.,.036,.088,.0666,1.42,.01,.45,.1,.15,.5)
      +.66*protokick(_t,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.);
                        amaydrumR = vel*protokick(_t2,242.,55.,.036,.088,.0666,1.42,.01,.45,.1,.15,.5)
      +.66*protokick(_t2,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.);
                    }
                    
                    if(drum > 0)
                    {
                        dL += amtL * s_atan(env * amaydrumL);
                        dR += amtR * s_atan(env * amaydrumR);
                    }
                }
                else
                {
                    f = freqC1(note_pitch(_note) + mod_transp(_mod));

                    if(abs(slide) > 1e-3) // THIS IS SLIDEY BIZ
                    {
                        float Bslide = trk_slide(trk);
                        float fac = slide * log(2.)/12.;
                        if (Bprog <= Bslide)
                        {
                            float help = 1. - Bprog/Bslide;
                            f *= Bslide * (fhelp(fac) - help * fhelp(fac*help*help)) / Bprog;
                        }
                        else
                        {
                            f *= 1. + (Bslide * (fhelp(fac)-1.)) / Bprog;
                        }
                    }

                    env = theta(Bprog) * (1. - smstep(Boff-rel, Boff, B));
                    if(syn == 0){amaysynL = _sin(f*_t); amaysynR = _sin(f*_t2);}
                    else if(syn == 22){
                        
                        amaysynL = (theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.01,0.,.1+.5*vel,.01,.4)*clip((1.+theta(Bprog)*exp(-11.*Bprog))*_tri(f*_t+.2*env_AHDSR(Bprog,L,.5,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t,.3+2.*vel+.2*(2.*fract(3.97*f*_t)-1.)))+.2*vel*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*(2.*fract(3.97*f*_t)-1.)))+.4*theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t,.3+2.*vel+.2*(2.*fract(3.97*f*_t)-1.)))*env_AHDSR(Bprog,L,0.,0.,.2+.2*vel,.01,.4)+.4*env_AHDSR(Bprog,L,0.,0.,.05,0.,0.)*lpnoise(_t+0.,6000.+200.*note_pitch(_note)));
                        amaysynR = (theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.01,0.,.1+.5*vel,.01,.4)*clip((1.+theta(Bprog)*exp(-11.*Bprog))*_tri(f*_t2+.2*env_AHDSR(Bprog,L,.5,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t2,.3+2.*vel+.2*(2.*fract(3.97*f*_t2)-1.)))+.2*vel*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*(2.*fract(3.97*f*_t2)-1.)))+.4*theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t2,.3+2.*vel+.2*(2.*fract(3.97*f*_t2)-1.)))*env_AHDSR(Bprog,L,0.,0.,.2+.2*vel,.01,.4)+.4*env_AHDSR(Bprog,L,0.,0.,.05,0.,0.)*lpnoise(_t2+0.,6000.+200.*note_pitch(_note)));
env = theta(Bprog)*pow(1.-smstep(Boff-rel, Boff, B),2);
                    }
                    else if(syn == 59){
                        
                        amaysynL = fqmbace5_volume(BT)*(vel*QFM(_t,f,0.,.00787*71.,.00787*52.,.00787*91.,.00787*99.,.5,1.,1.001,1.,.00787*49.,.00787*104.,.00787*65.,.00787*90.,7.)*env_AHDSR(_t,tL,.023,0.,.01,1.,.006));
                        amaysynR = fqmbace5_volume(BT)*(vel*QFM(_t2,f,0.,.00787*71.,.00787*52.,.00787*91.,.00787*99.,.5,1.,1.001,1.,.00787*49.,.00787*104.,.00787*65.,.00787*90.,7.)*env_AHDSR(_t2,tL,.023,0.,.01,1.,.006));
                    }
                    else if(syn == 64){
                        time2 = time-.03; _t2 = _t-.03;
                        amaysynL = fqmbace7sat_vol(BT)*s_atan(4.*(vel*2.5*QFM(_t,.993*f,0.,.00787*24.,.00787*92.,.00787*85.,.00787*21.,.5,1.,1.001,1.,.00787*126.,.00787*48.,.00787*26.,.00787*26.,9.)*env_AHDSRexp(_t,tL,.13,0.,.01,1.,.008)+vel*2.5*QFM(_t,f,0.,.00787*24.,.00787*92.,.00787*85.,.00787*21.,.5,1.,1.001,1.,.00787*126.,.00787*48.,.00787*26.,.00787*26.,9.)*env_AHDSRexp(_t,tL,.3,0.,.01,1.,.008)));
                        amaysynR = fqmbace7sat_vol(BT)*s_atan(4.*(vel*2.5*QFM(_t2,.993*f,0.,.00787*24.,.00787*92.,.00787*85.,.00787*21.,.5,1.,1.001,1.,.00787*126.,.00787*48.,.00787*26.,.00787*26.,9.)*env_AHDSRexp(_t2,tL,.13,0.,.01,1.,.008)+vel*2.5*QFM(_t2,f,0.,.00787*24.,.00787*92.,.00787*85.,.00787*21.,.5,1.,1.001,1.,.00787*126.,.00787*48.,.00787*26.,.00787*26.,9.)*env_AHDSRexp(_t2,tL,.3,0.,.01,1.,.008)));
env = theta(Bprog)*pow(1.-smstep(Boff-rel, Boff, B),4);
                    }
                    else if(syn == 91){
                        time2 = time-.1; _t2 = _t-.1;
                        amaysynL = sinshape((QFM(_t,f,0.,.00787*125.,.00787*env_AHDSRexp(Bprog,L,.0001,.097,.244,.195,0.)*88.,.00787*env_AHDSRexp(Bprog,L,.0001,.134,.246,.471,0.)*120.,.00787*env_AHDSRexp(Bprog,L,.0001,.051,.162,.161,0.)*91.,.999,1.,1.+.0939*(.5+(.5*_sin(.3*Bprog))),2.,.00787*101.,.00787*86.,.00787*1.,.00787*98.,8.)*env_AHDSRexp(Bprog,L,.0001,.032,.185,.109,.048)),.5,7.);
                        amaysynR = sinshape((QFM(_t2,f,0.,.00787*125.,.00787*env_AHDSRexp(Bprog,L,.0001,.097,.244,.195,0.)*88.,.00787*env_AHDSRexp(Bprog,L,.0001,.134,.246,.471,0.)*120.,.00787*env_AHDSRexp(Bprog,L,.0001,.051,.162,.161,0.)*91.,.999,1.,1.+.0939*(.5+(.5*_sin(.3*Bprog))),2.,.00787*101.,.00787*86.,.00787*1.,.00787*98.,8.)*env_AHDSRexp(Bprog,L,.0001,.032,.185,.109,.048)),.5,7.);
env = theta(Bprog)*pow(1.-smstep(Boff-rel, Boff, B),10);
                    }
                    else if(syn == 101){
                        
                        amaysynL = fett2_volume(BT)*env_AHDSRexp(Bprog,L,.05,0.,.05,.8,.5)*s_atan(23.*(s_atan(MADD(_t,.5*f,.5*lpnoise(_t + 0.,500.)*exp(-10.*Bprog)+env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3)*(.5+(.5*_sin(4.2*BT))),8,1,(.4+(.6*_sin(.21*Bprog))),2000.*env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3),100.,5.,100.,.001,.7,0.,0)+MADD(_t,1.01*.5*f,.5*lpnoise(_t + 0.,500.)*exp(-10.*Bprog)+env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3)*(.5+(.5*_sin(4.2*BT))),8,1,(.4+(.6*_sin(.21*Bprog))),2000.*env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3),100.,5.,100.,.001,.7,0.,0)+MADD(_t,.499*.5*f,.5*lpnoise(_t + 0.,500.)*exp(-10.*Bprog)+env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3)*(.5+(.5*_sin(4.2*BT))),8,1,(.4+(.6*_sin(.21*Bprog))),2000.*env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3),100.,5.,100.,.001,.7,0.,0))+.1*exp(-10.*Bprog)*lpnoise(_t + 0.,500.)));
                        amaysynR = fett2_volume(BT)*env_AHDSRexp(Bprog,L,.05,0.,.05,.8,.5)*s_atan(23.*(s_atan(MADD(_t2,.5*f,.5*lpnoise(_t2 + 0.,500.)*exp(-10.*Bprog)+env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3)*(.5+(.5*_sin(4.2*BT))),8,1,(.4+(.6*_sin(.21*Bprog))),2000.*env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3),100.,5.,100.,.001,.7,0.,0)+MADD(_t2,1.01*.5*f,.5*lpnoise(_t2 + 0.,500.)*exp(-10.*Bprog)+env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3)*(.5+(.5*_sin(4.2*BT))),8,1,(.4+(.6*_sin(.21*Bprog))),2000.*env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3),100.,5.,100.,.001,.7,0.,0)+MADD(_t2,.499*.5*f,.5*lpnoise(_t2 + 0.,500.)*exp(-10.*Bprog)+env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3)*(.5+(.5*_sin(4.2*BT))),8,1,(.4+(.6*_sin(.21*Bprog))),2000.*env_AHDSRexp(Bprog,L,.7,0.,.3,.3,.3),100.,5.,100.,.001,.7,0.,0))+.1*exp(-10.*Bprog)*lpnoise(_t2 + 0.,500.)));
env = theta(Bprog)*pow(1.-smstep(Boff-rel, Boff, B),8);
                    }
                    else if(syn == 109){
                        time2 = time-2e-3; _t2 = _t-2e-3;
                        amaysynL = (s_atan(1.*(.25*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(.51*f*(_t-0.0*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+.76*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*clip((1.+4.*(.5+(.5*_sin(.2*BT)))*(.5+(.5*_sin(.21*BT))))*_tri(f*(_t-0.0*(1.+3.*_sin(1.4*_t)))+.4+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-0.0*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t-0.0*(1.+3.*_sin(1.4*_t))) + 0.,3961.)*.1))+(.5+(.5*_sin(.2*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(1.995*f*(_t-0.0*(1.+3.*_sin(1.4*_t)))+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-0.0*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t-0.0*(1.+3.*_sin(1.4*_t))) + 0.,2142.)*.1)+.62*(.5+(.5*_sin(.21*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*(2.*fract(2.011*f*(_t-0.0*(1.+3.*_sin(1.4*_t)))+.1+.5*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-0.0*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog)))))-1.)+.27*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_sin_(.52*f+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*(.5+(.5*_sin(.5*Bprog)))*(_t-0.0*(1.+3.*_sin(1.4*_t))),.4+.25*(.5+(.5*_sin(.5*Bprog)))*.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-0.0*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog)))))))
      +s_atan(1.*(.25*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(.51*f*(_t-1.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+.76*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*clip((1.+4.*(.5+(.5*_sin(.2*BT)))*(.5+(.5*_sin(.21*BT))))*_tri(f*(_t-1.0e-05*(1.+3.*_sin(1.4*_t)))+.4+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-1.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t-1.0e-05*(1.+3.*_sin(1.4*_t))) + 0.,3961.)*.1))+(.5+(.5*_sin(.2*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(1.995*f*(_t-1.0e-05*(1.+3.*_sin(1.4*_t)))+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-1.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t-1.0e-05*(1.+3.*_sin(1.4*_t))) + 0.,2142.)*.1)+.62*(.5+(.5*_sin(.21*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*(2.*fract(2.011*f*(_t-1.0e-05*(1.+3.*_sin(1.4*_t)))+.1+.5*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-1.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog)))))-1.)+.27*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_sin_(.52*f+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*(.5+(.5*_sin(.5*Bprog)))*(_t-1.0e-05*(1.+3.*_sin(1.4*_t))),.4+.25*(.5+(.5*_sin(.5*Bprog)))*.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-1.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog)))))))
      +s_atan(1.*(.25*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(.51*f*(_t-2.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+.76*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*clip((1.+4.*(.5+(.5*_sin(.2*BT)))*(.5+(.5*_sin(.21*BT))))*_tri(f*(_t-2.0e-05*(1.+3.*_sin(1.4*_t)))+.4+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-2.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t-2.0e-05*(1.+3.*_sin(1.4*_t))) + 0.,3961.)*.1))+(.5+(.5*_sin(.2*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(1.995*f*(_t-2.0e-05*(1.+3.*_sin(1.4*_t)))+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-2.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t-2.0e-05*(1.+3.*_sin(1.4*_t))) + 0.,2142.)*.1)+.62*(.5+(.5*_sin(.21*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*(2.*fract(2.011*f*(_t-2.0e-05*(1.+3.*_sin(1.4*_t)))+.1+.5*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-2.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog)))))-1.)+.27*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_sin_(.52*f+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*(.5+(.5*_sin(.5*Bprog)))*(_t-2.0e-05*(1.+3.*_sin(1.4*_t))),.4+.25*(.5+(.5*_sin(.5*Bprog)))*.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-2.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog)))))))
      +s_atan(1.*(.25*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(.51*f*(_t-3.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+.76*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*clip((1.+4.*(.5+(.5*_sin(.2*BT)))*(.5+(.5*_sin(.21*BT))))*_tri(f*(_t-3.0e-05*(1.+3.*_sin(1.4*_t)))+.4+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-3.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t-3.0e-05*(1.+3.*_sin(1.4*_t))) + 0.,3961.)*.1))+(.5+(.5*_sin(.2*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(1.995*f*(_t-3.0e-05*(1.+3.*_sin(1.4*_t)))+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-3.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t-3.0e-05*(1.+3.*_sin(1.4*_t))) + 0.,2142.)*.1)+.62*(.5+(.5*_sin(.21*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*(2.*fract(2.011*f*(_t-3.0e-05*(1.+3.*_sin(1.4*_t)))+.1+.5*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-3.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog)))))-1.)+.27*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_sin_(.52*f+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*(.5+(.5*_sin(.5*Bprog)))*(_t-3.0e-05*(1.+3.*_sin(1.4*_t))),.4+.25*(.5+(.5*_sin(.5*Bprog)))*.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-3.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog)))))))
      +s_atan(1.*(.25*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(.51*f*(_t-4.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+.76*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*clip((1.+4.*(.5+(.5*_sin(.2*BT)))*(.5+(.5*_sin(.21*BT))))*_tri(f*(_t-4.0e-05*(1.+3.*_sin(1.4*_t)))+.4+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-4.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t-4.0e-05*(1.+3.*_sin(1.4*_t))) + 0.,3961.)*.1))+(.5+(.5*_sin(.2*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(1.995*f*(_t-4.0e-05*(1.+3.*_sin(1.4*_t)))+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-4.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t-4.0e-05*(1.+3.*_sin(1.4*_t))) + 0.,2142.)*.1)+.62*(.5+(.5*_sin(.21*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*(2.*fract(2.011*f*(_t-4.0e-05*(1.+3.*_sin(1.4*_t)))+.1+.5*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-4.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog)))))-1.)+.27*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_sin_(.52*f+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*(.5+(.5*_sin(.5*Bprog)))*(_t-4.0e-05*(1.+3.*_sin(1.4*_t))),.4+.25*(.5+(.5*_sin(.5*Bprog)))*.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t-4.0e-05*(1.+3.*_sin(1.4*_t)))+.35*(.5+(.5*_sin(.5*Bprog))))))));
                        amaysynR = (s_atan(1.*(.25*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(.51*f*(_t2-0.0*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+.76*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*clip((1.+4.*(.5+(.5*_sin(.2*BT)))*(.5+(.5*_sin(.21*BT))))*_tri(f*(_t2-0.0*(1.+3.*_sin(1.4*_t2)))+.4+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-0.0*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t2-0.0*(1.+3.*_sin(1.4*_t2))) + 0.,3961.)*.1))+(.5+(.5*_sin(.2*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(1.995*f*(_t2-0.0*(1.+3.*_sin(1.4*_t2)))+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-0.0*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t2-0.0*(1.+3.*_sin(1.4*_t2))) + 0.,2142.)*.1)+.62*(.5+(.5*_sin(.21*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*(2.*fract(2.011*f*(_t2-0.0*(1.+3.*_sin(1.4*_t2)))+.1+.5*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-0.0*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog)))))-1.)+.27*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_sin_(.52*f+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*(.5+(.5*_sin(.5*Bprog)))*(_t2-0.0*(1.+3.*_sin(1.4*_t2))),.4+.25*(.5+(.5*_sin(.5*Bprog)))*.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-0.0*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog)))))))
      +s_atan(1.*(.25*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(.51*f*(_t2-1.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+.76*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*clip((1.+4.*(.5+(.5*_sin(.2*BT)))*(.5+(.5*_sin(.21*BT))))*_tri(f*(_t2-1.0e-05*(1.+3.*_sin(1.4*_t2)))+.4+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-1.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t2-1.0e-05*(1.+3.*_sin(1.4*_t2))) + 0.,3961.)*.1))+(.5+(.5*_sin(.2*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(1.995*f*(_t2-1.0e-05*(1.+3.*_sin(1.4*_t2)))+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-1.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t2-1.0e-05*(1.+3.*_sin(1.4*_t2))) + 0.,2142.)*.1)+.62*(.5+(.5*_sin(.21*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*(2.*fract(2.011*f*(_t2-1.0e-05*(1.+3.*_sin(1.4*_t2)))+.1+.5*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-1.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog)))))-1.)+.27*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_sin_(.52*f+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*(.5+(.5*_sin(.5*Bprog)))*(_t2-1.0e-05*(1.+3.*_sin(1.4*_t2))),.4+.25*(.5+(.5*_sin(.5*Bprog)))*.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-1.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog)))))))
      +s_atan(1.*(.25*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(.51*f*(_t2-2.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+.76*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*clip((1.+4.*(.5+(.5*_sin(.2*BT)))*(.5+(.5*_sin(.21*BT))))*_tri(f*(_t2-2.0e-05*(1.+3.*_sin(1.4*_t2)))+.4+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-2.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t2-2.0e-05*(1.+3.*_sin(1.4*_t2))) + 0.,3961.)*.1))+(.5+(.5*_sin(.2*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(1.995*f*(_t2-2.0e-05*(1.+3.*_sin(1.4*_t2)))+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-2.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t2-2.0e-05*(1.+3.*_sin(1.4*_t2))) + 0.,2142.)*.1)+.62*(.5+(.5*_sin(.21*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*(2.*fract(2.011*f*(_t2-2.0e-05*(1.+3.*_sin(1.4*_t2)))+.1+.5*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-2.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog)))))-1.)+.27*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_sin_(.52*f+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*(.5+(.5*_sin(.5*Bprog)))*(_t2-2.0e-05*(1.+3.*_sin(1.4*_t2))),.4+.25*(.5+(.5*_sin(.5*Bprog)))*.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-2.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog)))))))
      +s_atan(1.*(.25*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(.51*f*(_t2-3.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+.76*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*clip((1.+4.*(.5+(.5*_sin(.2*BT)))*(.5+(.5*_sin(.21*BT))))*_tri(f*(_t2-3.0e-05*(1.+3.*_sin(1.4*_t2)))+.4+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-3.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t2-3.0e-05*(1.+3.*_sin(1.4*_t2))) + 0.,3961.)*.1))+(.5+(.5*_sin(.2*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(1.995*f*(_t2-3.0e-05*(1.+3.*_sin(1.4*_t2)))+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-3.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t2-3.0e-05*(1.+3.*_sin(1.4*_t2))) + 0.,2142.)*.1)+.62*(.5+(.5*_sin(.21*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*(2.*fract(2.011*f*(_t2-3.0e-05*(1.+3.*_sin(1.4*_t2)))+.1+.5*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-3.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog)))))-1.)+.27*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_sin_(.52*f+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*(.5+(.5*_sin(.5*Bprog)))*(_t2-3.0e-05*(1.+3.*_sin(1.4*_t2))),.4+.25*(.5+(.5*_sin(.5*Bprog)))*.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-3.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog)))))))
      +s_atan(1.*(.25*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(.51*f*(_t2-4.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+.76*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*clip((1.+4.*(.5+(.5*_sin(.2*BT)))*(.5+(.5*_sin(.21*BT))))*_tri(f*(_t2-4.0e-05*(1.+3.*_sin(1.4*_t2)))+.4+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-4.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t2-4.0e-05*(1.+3.*_sin(1.4*_t2))) + 0.,3961.)*.1))+(.5+(.5*_sin(.2*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_tri(1.995*f*(_t2-4.0e-05*(1.+3.*_sin(1.4*_t2)))+.15*(.5+(.5*_sin(.5*Bprog)))+.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-4.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*lpnoise((_t2-4.0e-05*(1.+3.*_sin(1.4*_t2))) + 0.,2142.)*.1)+.62*(.5+(.5*_sin(.21*BT)))*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*(2.*fract(2.011*f*(_t2-4.0e-05*(1.+3.*_sin(1.4*_t2)))+.1+.5*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-4.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog)))))-1.)+.27*clamp(1.+(.03-Bprog)/(.1),exp(-7.*Bprog),1.)*_sin_(.52*f+env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*(.5+(.5*_sin(.5*Bprog)))*(_t2-4.0e-05*(1.+3.*_sin(1.4*_t2))),.4+.25*(.5+(.5*_sin(.5*Bprog)))*.22*env_AHDSRexp(Bprog,L,.025,0.,.1,.3,0.)*_tri(.51*f*(_t2-4.0e-05*(1.+3.*_sin(1.4*_t2)))+.35*(.5+(.5*_sin(.5*Bprog))))))));
                    }
                    else if(syn == 111){
                        time2 = time-0.02; _t2 = _t-0.02;
                        amaysynL = fett2_volume(BT)*((MADD((_t-0.0*(1.+2.*_sin(.013*_t))),f,0.,12,1,.92,2513.*env_AHDSRexp(Bprog,L,.031,0.,.1,1.,0.3),65.1,1.,3.,.008,.25,0.,1)*env_AHDSRexp(Bprog,L,.145,0.,.1,1.,0.3))
      +(MADD((_t-4.0e-02*(1.+2.*_sin(.013*_t))),f,0.,12,1,.92,2513.*env_AHDSRexp(Bprog,L,.031,0.,.1,1.,0.3),65.1,1.,3.,.008,.25,0.,1)*env_AHDSRexp(Bprog,L,.145,0.,.1,1.,0.3))
      +(MADD((_t-8.0e-02*(1.+2.*_sin(.013*_t))),f,0.,12,1,.92,2513.*env_AHDSRexp(Bprog,L,.031,0.,.1,1.,0.3),65.1,1.,3.,.008,.25,0.,1)*env_AHDSRexp(Bprog,L,.145,0.,.1,1.,0.3)));
                        amaysynR = fett2_volume(BT)*((MADD((_t2-0.0*(1.+2.*_sin(.013*_t2))),f,0.,12,1,.92,2513.*env_AHDSRexp(Bprog,L,.031,0.,.1,1.,0.3),65.1,1.,3.,.008,.25,0.,1)*env_AHDSRexp(Bprog,L,.145,0.,.1,1.,0.3))
      +(MADD((_t2-4.0e-02*(1.+2.*_sin(.013*_t2))),f,0.,12,1,.92,2513.*env_AHDSRexp(Bprog,L,.031,0.,.1,1.,0.3),65.1,1.,3.,.008,.25,0.,1)*env_AHDSRexp(Bprog,L,.145,0.,.1,1.,0.3))
      +(MADD((_t2-8.0e-02*(1.+2.*_sin(.013*_t2))),f,0.,12,1,.92,2513.*env_AHDSRexp(Bprog,L,.031,0.,.1,1.,0.3),65.1,1.,3.,.008,.25,0.,1)*env_AHDSRexp(Bprog,L,.145,0.,.1,1.,0.3)));
env = theta(Bprog)*pow(1.-smstep(Boff-rel, Boff, B),4);
                    }
                    else if(syn == 117){
                        time2 = time-.1; _t2 = _t-.1;
                        amaysynL = fett2_volume(BT)*smstep(0.,.24975*(2.-vel+1e-3),Bprog)*(s_atan(QFM((_t-0.0*(1.+3.*_sin(.05*_t))),f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t-0.0*(1.+3.*_sin(.05*_t))),1.003*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t-0.0*(1.+3.*_sin(.05*_t))),.992*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.))
      +s_atan(QFM((_t-4.0e-03*(1.+3.*_sin(.05*_t))),f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t-4.0e-03*(1.+3.*_sin(.05*_t))),1.003*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t-4.0e-03*(1.+3.*_sin(.05*_t))),.992*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.))
      +s_atan(QFM((_t-8.0e-03*(1.+3.*_sin(.05*_t))),f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t-8.0e-03*(1.+3.*_sin(.05*_t))),1.003*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t-8.0e-03*(1.+3.*_sin(.05*_t))),.992*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)));
                        amaysynR = fett2_volume(BT)*smstep(0.,.24975*(2.-vel+1e-3),Bprog)*(s_atan(QFM((_t2-0.0*(1.+3.*_sin(.05*_t2))),f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t2-0.0*(1.+3.*_sin(.05*_t2))),1.003*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t2-0.0*(1.+3.*_sin(.05*_t2))),.992*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.))
      +s_atan(QFM((_t2-4.0e-03*(1.+3.*_sin(.05*_t2))),f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t2-4.0e-03*(1.+3.*_sin(.05*_t2))),1.003*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t2-4.0e-03*(1.+3.*_sin(.05*_t2))),.992*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.))
      +s_atan(QFM((_t2-8.0e-03*(1.+3.*_sin(.05*_t2))),f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t2-8.0e-03*(1.+3.*_sin(.05*_t2))),1.003*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)+QFM((_t2-8.0e-03*(1.+3.*_sin(.05*_t2))),.992*f,0.,.00787*127.,.00787*111.,.00787*127.,.00787*95.,1.,4.,1.,.5,.00787*55.,.00787*12.,.00787*56.,.00787*26.,9.)));
env = theta(Bprog)*pow(1.-smstep(Boff-rel, Boff, B),3);
                    }
                    else if(syn == 118){
                        time2 = time-0.008; _t2 = _t-0.008;
                        amaysynL = (1.0*(s_atan(10.*env_AHDSRexp((Bprog-0.000),L,0.,0.,.2,0.,1.)*QFM(((_t-SPB*0.000)-0.0*(1.+.5*_sin(.32*(_t-SPB*0.000)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.))
      +s_atan(10.*env_AHDSRexp((Bprog-0.000),L,0.,0.,.2,0.,1.)*QFM(((_t-SPB*0.000)-1.0e-02*(1.+.5*_sin(.32*(_t-SPB*0.000)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.)))
      +2.0e-01*(s_atan(10.*env_AHDSRexp((Bprog-3.330e-01),L,0.,0.,.2,0.,1.)*QFM(((_t-SPB*3.330e-01)-0.0*(1.+.5*_sin(.32*(_t-SPB*3.330e-01)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.))
      +s_atan(10.*env_AHDSRexp((Bprog-3.330e-01),L,0.,0.,.2,0.,1.)*QFM(((_t-SPB*3.330e-01)-1.0e-02*(1.+.5*_sin(.32*(_t-SPB*3.330e-01)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.)))
      +4.0e-02*(s_atan(10.*env_AHDSRexp((Bprog-6.660e-01),L,0.,0.,.2,0.,1.)*QFM(((_t-SPB*6.660e-01)-0.0*(1.+.5*_sin(.32*(_t-SPB*6.660e-01)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.))
      +s_atan(10.*env_AHDSRexp((Bprog-6.660e-01),L,0.,0.,.2,0.,1.)*QFM(((_t-SPB*6.660e-01)-1.0e-02*(1.+.5*_sin(.32*(_t-SPB*6.660e-01)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.))));
                        amaysynR = (1.0*(s_atan(10.*env_AHDSRexp((Bprog-0.000),L,0.,0.,.2,0.,1.)*QFM(((_t2-SPB*0.000)-0.0*(1.+.5*_sin(.32*(_t2-SPB*0.000)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.))
      +s_atan(10.*env_AHDSRexp((Bprog-0.000),L,0.,0.,.2,0.,1.)*QFM(((_t2-SPB*0.000)-1.0e-02*(1.+.5*_sin(.32*(_t2-SPB*0.000)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.)))
      +2.0e-01*(s_atan(10.*env_AHDSRexp((Bprog-3.330e-01),L,0.,0.,.2,0.,1.)*QFM(((_t2-SPB*3.330e-01)-0.0*(1.+.5*_sin(.32*(_t2-SPB*3.330e-01)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.))
      +s_atan(10.*env_AHDSRexp((Bprog-3.330e-01),L,0.,0.,.2,0.,1.)*QFM(((_t2-SPB*3.330e-01)-1.0e-02*(1.+.5*_sin(.32*(_t2-SPB*3.330e-01)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.)))
      +4.0e-02*(s_atan(10.*env_AHDSRexp((Bprog-6.660e-01),L,0.,0.,.2,0.,1.)*QFM(((_t2-SPB*6.660e-01)-0.0*(1.+.5*_sin(.32*(_t2-SPB*6.660e-01)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.))
      +s_atan(10.*env_AHDSRexp((Bprog-6.660e-01),L,0.,0.,.2,0.,1.)*QFM(((_t2-SPB*6.660e-01)-1.0e-02*(1.+.5*_sin(.32*(_t2-SPB*6.660e-01)))),f,0.,.00787*113.,.00787*86.,.00787*113.,.00787*86.,1.,3.5,1.,3.5,.00787*0.,.00787*0.,.00787*0.,.00787*0.,8.))));
                    }
                    else if(syn == 119){
                        time2 = time-3e-3; _t2 = _t-3e-3;
                        amaysynL = hoboe_vol(BT)*env_AHDSRexp(Bprog,L,.125,0.,.125,.5,.7)*MADD(_t,f,.3*lpnoise(_t + 0.,3300.)+env_AHDSR(Bprog,L,2.4,0.,.1,1.,.3)*(.5+(.5*_sin(7.*env_AHDSR(Bprog,L,2.4,0.,.1,1.,.3)*Bprog))),8,1,-2.,600.,10.,.4,100.,.006,3.,0.,0)
      +hoboe_vol(BT)*.7*env_AHDSRexp(Bprog,L,.125,0.,.125,.5,.7)*MADD(_t,2.*f,.3*lpnoise(_t + 0.,3300.),16,1,-2.,10000.,100.,.4,100.,.001*env_AHDSR(Bprog,L,2.4,0.,.1,1.,.3)*(.5+(.5*_sin(7.*env_AHDSR(Bprog,L,2.4,0.,.1,1.,.3)*Bprog))),7.,0.,0);
                        amaysynR = hoboe_vol(BT)*env_AHDSRexp(Bprog,L,.125,0.,.125,.5,.7)*MADD(_t2,f,.3*lpnoise(_t2 + 0.,3300.)+env_AHDSR(Bprog,L,2.4,0.,.1,1.,.3)*(.5+(.5*_sin(7.*env_AHDSR(Bprog,L,2.4,0.,.1,1.,.3)*Bprog))),8,1,-2.,600.,10.,.4,100.,.006,3.,0.,0)
      +hoboe_vol(BT)*.7*env_AHDSRexp(Bprog,L,.125,0.,.125,.5,.7)*MADD(_t2,2.*f,.3*lpnoise(_t2 + 0.,3300.),16,1,-2.,10000.,100.,.4,100.,.001*env_AHDSR(Bprog,L,2.4,0.,.1,1.,.3)*(.5+(.5*_sin(7.*env_AHDSR(Bprog,L,2.4,0.,.1,1.,.3)*Bprog))),7.,0.,0);
                    }
                    
                    sL += amtL * trk_norm(trk) * s_atan(clamp(env,0.,1.) * amaysynL);
                    sR += amtR * trk_norm(trk) * s_atan(clamp(env,0.,1.) * amaysynR);
                }
            }
        }
    }
    float masterL = .26 * sidechain * s_atan(sL) + .8 * dL;
    float masterR = .26 * sidechain * s_atan(sR) + .8 * dR;
    return vec2(
        (BT > 77 && BT < 80) ? (1. - 0.26*(BT-77.)) * masterL: masterL,
        (BT > 77 && BT < 80) ? (1. - 0.26*(BT-77.)) * masterR: masterR);
}

void main()
{
    Tsample = 1./iSampleRate;
    float t = (iBlockOffset + gl_FragCoord.x + gl_FragCoord.y*iTexSize) * Tsample;
    vec2 s = mainSynth(t);
    vec2 v  = floor((0.5+0.5*s)*65535.0);
    vec2 vl = mod(v,256.0)/255.0;
    vec2 vh = floor(v/256.0)/255.0;
    gl_FragColor = vec4(vl.x,vh.x,vl.y,vh.y);
}
