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
float s_atan(float a) { return .636 * atan(a); }
float doubleslope(float t, float a, float d, float s) { return smstep(-.00001,a,t) - (1.-s) * smstep(0.,d,t-a); }

#define NTIME 2
const float pos_B[2] = float[2](0.,89.);
const float pos_t[2] = float[2](0.,166.875);
const float pos_BPS[1] = float[1](.5333);
const float pos_SPB[1] = float[1](1.8751);
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

float waveshape(float s, float amt, float A, float B, float C, float D, float E)
{
    float w;
    float m = sign(s);
    s = abs(s);

    if(s<A) w = B * smstep(0.,A,s);
    else if(s<C) w = C + (B-C) * smstep(C,A,s);
    else if(s<=D) w = s;
    else if(s<=1.)
    {
        float _s = (s-D)/(1.-D);
        w = D + (E-D) * (1.5*_s*(1.-.33*_s*_s));
    }
    else return 1.;

    return m*mix(s,w,amt);
}

float sinshape(float x, float amt, float parts)
{
    return (1.-amt) * x + amt * sign(x) * 0.5 * (1. - cos(parts*PI*x));
}

float comp_SAW(int N, float inv_N, float PW) {return inv_N * (1. - _sin(float(N)*PW));}
float comp_TRI(int N, float inv_N, float PW) {return N % 2 == 0 ? .1 * inv_N * _sin(float(N)*PW) : inv_N * inv_N * (1. - _sin(float(N)*PW));}
float comp_SQU(int N, float inv_N, float PW) {return inv_N * (minus1hochN(N) * _sin(.5*float(N)*PW + .25) - 1.);}
float comp_HAE(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (1. - minus1hochNminus1halbe(N))*_sin(PW);}

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
        comp_mix = MIX < 0. ? (MIX+1.) * comp_TRI(N,inv_N,PW)  -     MIX  * comp_SAW(N,inv_N,PW)
                 : MIX < 1. ? (1.-MIX) * comp_TRI(N,inv_N,PW)  +     MIX  * comp_SQU(N,inv_N,PW)
                            : (MIX-1.) * comp_HAE(N,inv_N,PW)  + (2.-MIX) * comp_SQU(N,inv_N,PW);

        if(abs(comp_mix) < 1e-6) continue;

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

#define NTRK 14
#define NMOD 133
#define NPTN 27
#define NNOT 874
#define NDRM 48

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

    time = mod(time, 166.875);
    
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
                if(syn == 114)
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

                if(syn == 114)
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
                    else if(drum == 5){
                        amaydrumL = vel*.9*protokick(_t,242.,55.,.036,.03,.0666,1.42,.02,.25,.01,.1,.4)
      +.9*protokick(_t,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.)
      +.64*((clamp(2.27*_tri(drop_phase(_t,.03,241.,72.)),-1.,1.)*(1.-smstep(-1e-3,0.,_t-.01))+.91*clamp(.9*_tri(drop_phase(_t,.03,241.,72.)+.91*lpnoise(_t,8164.)),-1.,1.)*exp(-20.76*_t)+.05*lpnoise(_t,10466.)*(1.-smstep(0.,.18,_t-.56))+.56*lpnoise(_t,7123.)*exp(-_t*5.45)+.11*lpnoise(_t,1134.)*exp(-_t*13.82))*smstep(0.,.004,_t));
                        amaydrumR = vel*.9*protokick(_t2,242.,55.,.036,.03,.0666,1.42,.02,.25,.01,.1,.4)
      +.9*protokick(_t2,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.)
      +.64*((clamp(2.27*_tri(drop_phase(_t2,.03,241.,72.)),-1.,1.)*(1.-smstep(-1e-3,0.,_t2-.01))+.91*clamp(.9*_tri(drop_phase(_t2,.03,241.,72.)+.91*lpnoise(_t2,8164.)),-1.,1.)*exp(-20.76*_t2)+.05*lpnoise(_t2,10466.)*(1.-smstep(0.,.18,_t2-.56))+.56*lpnoise(_t2,7123.)*exp(-_t2*5.45)+.11*lpnoise(_t2,1134.)*exp(-_t2*13.82))*smstep(0.,.004,_t2));
                    }
                    else if(drum == 11){
                        amaydrumL = vel*fract(sin(_t*100.*.5)*50000.*.5)*doubleslope(_t,0.,.03,.1)*exp(-13.*Bprog);
                        amaydrumR = vel*fract(sin(_t2*100.*.5)*50000.*.5)*doubleslope(_t2,0.,.03,.1)*exp(-13.*Bprog);
                    }
                    else if(drum == 24){
                        amaydrumL = vel*((clamp(1.32*_tri(drop_phase(_t,.06,308.,80.)),-1.,1.)*(1.-smstep(-1e-3,0.,_t-.15))+.82*clamp(.49*_tri(drop_phase(_t,.06,308.,80.)+.82*lpnoise(_t,4595.)),-1.,1.)*exp(-1.97*_t)+.09*lpnoise(_t,4032.)*(1.-smstep(0.,.97,_t-.79))+.1*lpnoise(_t,1111.)*exp(-_t*12.69)+.6*lpnoise(_t,7795.)*exp(-_t*1.08))*smstep(0.,.003,_t));
                        amaydrumR = vel*((clamp(1.32*_tri(drop_phase(_t2,.06,308.,80.)),-1.,1.)*(1.-smstep(-1e-3,0.,_t2-.15))+.82*clamp(.49*_tri(drop_phase(_t2,.06,308.,80.)+.82*lpnoise(_t2,4595.)),-1.,1.)*exp(-1.97*_t2)+.09*lpnoise(_t2,4032.)*(1.-smstep(0.,.97,_t2-.79))+.1*lpnoise(_t2,1111.)*exp(-_t2*12.69)+.6*lpnoise(_t2,7795.)*exp(-_t2*1.08))*smstep(0.,.003,_t2));
                    }
                    else if(drum == 32){
                        amaydrumL = vel*1.4*protokick(_t,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.);
                        amaydrumR = vel*1.4*protokick(_t2,3333.,340.,.008,0.,.01,2.,2.4,0.,.2,.3,1.);
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
                    else if(syn == 21){
                        
                        amaysynL = (theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.01,0.,.1+.5*vel,.01,.4)*clip((1.+theta(Bprog)*exp(-11.*Bprog))*_tri(f*_t+.2*env_AHDSR(Bprog,L,.5,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t,.3+2.*vel+.2*(2.*fract(3.97*f*_t)-1.)))+.2*vel*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*(2.*fract(3.97*f*_t)-1.)))+.4*theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t,.3+2.*vel+.2*(2.*fract(3.97*f*_t)-1.)))*env_AHDSR(Bprog,L,0.,0.,.2+.2*vel,.01,.4)+.4*env_AHDSR(Bprog,L,0.,0.,.05,0.,0.)*lpnoise(_t+0.,6000.+200.*note_pitch(_note)));
                        amaysynR = (theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.01,0.,.1+.5*vel,.01,.4)*clip((1.+theta(Bprog)*exp(-11.*Bprog))*_tri(f*_t2+.2*env_AHDSR(Bprog,L,.5,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t2,.3+2.*vel+.2*(2.*fract(3.97*f*_t2)-1.)))+.2*vel*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*(2.*fract(3.97*f*_t2)-1.)))+.4*theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t2,.3+2.*vel+.2*(2.*fract(3.97*f*_t2)-1.)))*env_AHDSR(Bprog,L,0.,0.,.2+.2*vel,.01,.4)+.4*env_AHDSR(Bprog,L,0.,0.,.05,0.,0.)*lpnoise(_t2+0.,6000.+200.*note_pitch(_note)));
env = theta(Bprog)*pow(1.-smstep(Boff-rel, Boff, B),2);
                    }
                    else if(syn == 22){
                        
                        amaysynL = (theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.01,0.,.1+.5*vel,.01,.4)*clip((1.+theta(Bprog)*exp(-11.*Bprog))*_tri(f*_t+.2*env_AHDSR(Bprog,L,.5,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t,.3+2.*vel+.2*(2.*fract(3.97*f*_t)-1.)))+.2*vel*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*(2.*fract(3.97*f*_t)-1.)))+.4*theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t,.3+2.*vel+.2*(2.*fract(3.97*f*_t)-1.)))*env_AHDSR(Bprog,L,0.,0.,.2+.2*vel,.01,.4)+.4*env_AHDSR(Bprog,L,0.,0.,.05,0.,0.)*lpnoise(_t+0.,6000.+200.*note_pitch(_note)));
                        amaysynR = (theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.01,0.,.1+.5*vel,.01,.4)*clip((1.+theta(Bprog)*exp(-11.*Bprog))*_tri(f*_t2+.2*env_AHDSR(Bprog,L,.5,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t2,.3+2.*vel+.2*(2.*fract(3.97*f*_t2)-1.)))+.2*vel*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*(2.*fract(3.97*f*_t2)-1.)))+.4*theta(Bprog)*exp(-11.*Bprog)*env_AHDSR(Bprog,L,.325,1.,.1,1.,0.)*clip((1.+3.)*_sq_(1.99*f*_t2,.3+2.*vel+.2*(2.*fract(3.97*f*_t2)-1.)))*env_AHDSR(Bprog,L,0.,0.,.2+.2*vel,.01,.4)+.4*env_AHDSR(Bprog,L,0.,0.,.05,0.,0.)*lpnoise(_t2+0.,6000.+200.*note_pitch(_note)));
env = theta(Bprog)*pow(1.-smstep(Boff-rel, Boff, B),2);
                    }
                    else if(syn == 52){
                        
                        amaysynL = (env_AHDSR(_t,tL,1e-3,0.,.063,.402,.003)*MADD(_t,f,0.,128,1,-.6,(1175.+(583.*_sin_(2.*B,.4))),24.,44.37,24.05,.015,.4*(.55+(.4*clip((1.+1.)*_sin(4.*B)))),0));
                        amaysynR = (env_AHDSR(_t2,tL,1e-3,0.,.063,.402,.003)*MADD(_t2,f,0.,128,1,-.6,(1175.+(583.*_sin_(2.*B,.4))),24.,44.37,24.05,.015,.4*(.55+(.4*clip((1.+1.)*_sin(4.*B)))),0));
env = theta(Bprog)*pow(1.-smstep(Boff-rel, Boff, B),2);
                    }
                    else if(syn == 59){
                        
                        amaysynL = (vel*QFM(_t,f,0.,.00787*71.,.00787*52.,.00787*91.,.00787*99.,.5,1.,1.001,1.,.00787*49.,.00787*104.,.00787*65.,.00787*90.,7.)*env_AHDSR(_t,tL,.023,0.,.01,1.,.006));
                        amaysynR = (vel*QFM(_t2,f,0.,.00787*71.,.00787*52.,.00787*91.,.00787*99.,.5,1.,1.001,1.,.00787*49.,.00787*104.,.00787*65.,.00787*90.,7.)*env_AHDSR(_t2,tL,.023,0.,.01,1.,.006));
                    }
                    else if(syn == 79){
                        
                        amaysynL = (vel*env_AHDSR(_t,tL,.388,0.,.1,1.,.275)*waveshape(MADD(_t,f,0.,256,1,-.143,(6861.+(2407.*clip((1.+.958)*_sin(.493*B)))),9.726,1.282,7.24,.002,-.64,0),-.64,.1,.5,.2,.5,.6)+.5*vel*env_AHDSR(_t,tL,.388,0.,.1,1.,.275)*MADD(_t,.501*f,0.,64,1,-.143,(6861.+(2407.*clip((1.+.958)*_sin(.493*B)))),9.726,1.282,7.24,.002,.2*-.64,0));
                        amaysynR = (vel*env_AHDSR(_t2,tL,.388,0.,.1,1.,.275)*waveshape(MADD(_t2,f,0.,256,1,-.143,(6861.+(2407.*clip((1.+.958)*_sin(.493*B)))),9.726,1.282,7.24,.002,-.64,0),-.64,.1,.5,.2,.5,.6)+.5*vel*env_AHDSR(_t2,tL,.388,0.,.1,1.,.275)*MADD(_t2,.501*f,0.,64,1,-.143,(6861.+(2407.*clip((1.+.958)*_sin(.493*B)))),9.726,1.282,7.24,.002,.2*-.64,0));
                    }
                    else if(syn == 80){
                        
                        amaysynL = (vel*env_AHDSR(_t,tL,.396,0.,.1,1.,.381)*waveshape(MADD(_t,f,0.,256,1,-.695,(6326.+(862.*clip((1.+.39)*_sin(.652*B)))),19.795,1.863,.05,.013,.54,0),.54,.1,.5,.2,.5,.6)+.5*vel*env_AHDSR(_t,tL,.396,0.,.1,1.,.381)*MADD(_t,.501*f,0.,64,1,-.695,(6326.+(862.*clip((1.+.39)*_sin(.652*B)))),19.795,1.863,.05,.013,.2*.54,0));
                        amaysynR = (vel*env_AHDSR(_t2,tL,.396,0.,.1,1.,.381)*waveshape(MADD(_t2,f,0.,256,1,-.695,(6326.+(862.*clip((1.+.39)*_sin(.652*B)))),19.795,1.863,.05,.013,.54,0),.54,.1,.5,.2,.5,.6)+.5*vel*env_AHDSR(_t2,tL,.396,0.,.1,1.,.381)*MADD(_t2,.501*f,0.,64,1,-.695,(6326.+(862.*clip((1.+.39)*_sin(.652*B)))),19.795,1.863,.05,.013,.2*.54,0));
                    }
                    else if(syn == 98){
                        
                        amaysynL = (vel*sinshape(QFM(_t,f,0.,.00787*125.,.00787*env_AHDSR(Bprog,L,.0001,.047,.01,.404,0.)*20.,.00787*env_AHDSR(Bprog,L,.0001,.151,.071,.069,0.)*110.,.00787*env_AHDSR(Bprog,L,.0001,.232,.08,.003,0.)*65.,.999,1.,1.+.0799*(.5+(.5*_sin(.18*Bprog))),2.,.00787*109.,.00787*21.,.00787*94.,.00787*0.,11.),.03*aux,3.)*env_AHDSR(Bprog,L,.0001,.03,.167,.796,.114))
      +_sq_(.501*f*_t,.4+.3*(.5+(.5*_sin(.8*Bprog))))
      +.6*clip((1.+.2*aux)*_sin(.25*f*_t));
                        amaysynR = (vel*sinshape(QFM(_t2,f,0.,.00787*125.,.00787*env_AHDSR(Bprog,L,.0001,.047,.01,.404,0.)*20.,.00787*env_AHDSR(Bprog,L,.0001,.151,.071,.069,0.)*110.,.00787*env_AHDSR(Bprog,L,.0001,.232,.08,.003,0.)*65.,.999,1.,1.+.0799*(.5+(.5*_sin(.18*Bprog))),2.,.00787*109.,.00787*21.,.00787*94.,.00787*0.,11.),.03*aux,3.)*env_AHDSR(Bprog,L,.0001,.03,.167,.796,.114))
      +_sq_(.501*f*_t2,.4+.3*(.5+(.5*_sin(.8*Bprog))))
      +.6*clip((1.+.2*aux)*_sin(.25*f*_t2));
                    }
                    else if(syn == 108){
                        
                        amaysynL = (MADD(_t,f,0.,72,1,.92,2513.*env_AHDSR(Bprog,L,.031,0.,.1,1.,0.),65.1,1.,3.,.013,.25,0.,0)*env_AHDSR(Bprog,L,.145,0.,.1,1.,0.));
                        amaysynR = (MADD(_t2,f,0.,72,1,.92,2513.*env_AHDSR(Bprog,L,.031,0.,.1,1.,0.),65.1,1.,3.,.013,.25,0.,0)*env_AHDSR(Bprog,L,.145,0.,.1,1.,0.));
                    }
                    else if(syn == 110){
                        
                        amaysynL = (MADD(_t,f,0.,72,1,-.76,826.*env_AHDSR(Bprog,L,.184,0.,.1,1.,0.),30.6,6.7,3.,.013,.87,0.,0)*env_AHDSR(Bprog,L,.112,0.,.1,1.,0.));
                        amaysynR = (MADD(_t2,f,0.,72,1,-.76,826.*env_AHDSR(Bprog,L,.184,0.,.1,1.,0.),30.6,6.7,3.,.013,.87,0.,0)*env_AHDSR(Bprog,L,.112,0.,.1,1.,0.));
                    }
                    else if(syn == 111){
                        
                        amaysynL = (MADD(_t,f,0.,72,1,1.79,357.*env_AHDSR(Bprog,L,.097,0.,.1,1.,0.),43.9,.4,3.,.013,.51,0.,0)*env_AHDSR(Bprog,L,.169,0.,.1,1.,0.));
                        amaysynR = (MADD(_t2,f,0.,72,1,1.79,357.*env_AHDSR(Bprog,L,.097,0.,.1,1.,0.),43.9,.4,3.,.013,.51,0.,0)*env_AHDSR(Bprog,L,.169,0.,.1,1.,0.));
                    }
                    else if(syn == 113){
                        
                        amaysynL = (MADD(_t,f,0.,72,1,1.53,898.*env_AHDSR(Bprog,L,.025,0.,.1,1.,0.),13.5,.9,3.,.013,.07,0.,0)*env_AHDSR(Bprog,L,.16,0.,.1,1.,0.));
                        amaysynR = (MADD(_t2,f,0.,72,1,1.53,898.*env_AHDSR(Bprog,L,.025,0.,.1,1.,0.),13.5,.9,3.,.013,.07,0.,0)*env_AHDSR(Bprog,L,.16,0.,.1,1.,0.));
                    }
                    
                    sL += amtL * trk_norm(trk) * s_atan(clamp(env,0.,1.) * amaysynL);
                    sR += amtR * trk_norm(trk) * s_atan(clamp(env,0.,1.) * amaysynR);
                }
            }
        }
    }
    return .2 * sidechain * vec2(s_atan(sL), s_atan(sR)) + .67 * vec2(s_atan(dL), s_atan(dR));
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
