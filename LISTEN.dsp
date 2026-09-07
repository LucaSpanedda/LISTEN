/*******************************************************************************
* LISTEN : Local InSTrument ENvironment
*******************************************************************************/

declare name "Local InSTrument ENvironment";
declare author "Luca Spanedda";
declare version "0.1";
declare description "Composition";
declare copyright "Copyright (C) 2026 Luca Spanedda
< lucaspanedda1995 [at] gmail [dot] com >";
declare license "MIT license";
// Import the standard Faust Libraries
import("stdfaust.lib");


//------- ------------- ----- -----------
//-- LIBRARY -------------------------------------------------------------------
//------- --------


//-- BASIC OPERATIONS ----------------------------------------------------------

// Inspector
inspect = _ <: (_, 
    (_ : vbargraph(" _ [style:numerical]", ma.MIN, ma.MAX)) : 
    attach);

inspect_VH(minv, maxv) = _ <: (_, 
    (_ : hbargraph(" _ ", minv, maxv)) : 
    attach);
    
inspect_VV(minv, maxv) = _ <: (_, 
    (_ : vbargraph(" _ ", minv, maxv)) : 
    attach);

// Bipolar to Unipolar signal Mapping
bitoUni(x) = (x + 1.0) * 0.5;

// Clipping Function
clip(low, high, x) = max(low, min(high, x));

// phasor that start from 0
phasor0(f) = (_ <: _ + f, _) ~  _ % ma.SR : (!, _ / ma.SR);

// triangular function
triangularFunc(x) = abs(ma.frac((x - .5)) * 2.0 - 1.0);


//-- FILTERS -------------------------------------------------------------------

// Vadim Zavalishin's Onepole TPT filter (Topology Preserving Transform)
OnePoleTPT(cf, x) = loop ~ _ : ! , si.bus(3)
with {
    g = tan(cf * ma.PI * (1/ma.SR));
    G = g / (1.0 + g);
    loop(s) = u , lp , hp , ap
    with {
        v = (x - s) * G;
        u = v + lp;
        lp = v + s;
        hp = x - lp;
        ap = lp - hp;
    };
};
LPTPT(CF, x) = OnePoleTPT(max(ma.EPSILON, min(20480, CF)), x) : (_ , ! , !);
HPTPT(CF, x) = OnePoleTPT(max(ma.EPSILON, min(20480, CF)), x) : (! , _ , !);

// Vadim Zavalishin's SVF TPT filter (Topology Preserving Transform)
SVFTPT(K, Q, CF, x) = circuitout : ! , ! , _ , _ , _ , _ , _ , _ , _ , _
with{
    g = tan(CF * ma.PI / ma.SR);
    R = 1.0 / (2.0 * Q);
    G1 = 1.0 / (1.0 + 2.0 * R * g + g * g);
    G2 = 2.0 * R + g;
    circuit(s1, s2) = u1 , u2 , lp , hp , bp, notch, apf, ubp, peak, bshelf
        with{
            hp = (x - s1 * G2 - s2) * G1;
            v1 = hp * g;
            bp = s1 + v1;
            v2 = bp * g;
            lp = s2 + v2;
            u1 = v1 + bp;
            u2 = v2 + lp;
            notch = x - ((2 * R) * bp);
            apf = x - ((4 * R) * bp);
            ubp = ((2 * R) * bp);
            peak = lp - hp;
            bshelf = x + (((2 * K) * R) * bp);
        };
    // choose the output from the SVF Filter (ex. bshelf)
    circuitout = circuit ~ si.bus(2);
};

BPSVF(Q, CF, x) = SVFTPT(0, Q, max(ma.EPSILON, min(20480, CF)), x) :
  ! , ! , _ , ! , ! , ! , ! , ! ;
NotchSVF(Q, CF, x) = SVFTPT(0, Q, max(ma.EPSILON, min(20480, CF)), x) :
  ! , ! , ! , _ , ! , ! , ! , ! ;
LPSVF(Q, CF, x) = SVFTPT(0, Q, max(ma.EPSILON, min(20480, CF)), x) :
  _ , ! , ! , ! , ! , ! , ! , ! ;
HPSVF(Q, CF, x) = SVFTPT(0, Q, max(ma.EPSILON, min(20480, CF)), x) :
  ! , _ , ! , ! , ! , ! , ! , !;

// TPT Butterworth Filters from Zavalishin's TPT
butterworthQ(order, stage) = qFactor(order % 2)
with {
    qFactor(0) = 1.0 / (2.0 * cos(((2.0 * stage + 1) *
    (ma.PI / (order * 2.0)))));
    qFactor(1) = 1.0 / (2.0 * cos(((stage + 1) * (ma.PI / order))));
};

LPButterworthN(1, cf, x) = LPTPT(cf, x);
LPButterworthN(N, cf, x) = cascade(N % 2)
with {
    cascade(0) = x : seq(i, N / 2, LPSVF(butterworthQ(N, i), cf));
    cascade(1) = x : LPTPT(cf) : seq(i, (N - 1) / 2,
    LPSVF(butterworthQ(N, i), cf));
};

HPButterworthN(1, cf, x) = HPTPT(cf, x);
HPButterworthN(N, cf, x) = cascade(N % 2)
with {
    cascade(0) = x : seq(i, N / 2, HPSVF(butterworthQ(N, i), cf));
    cascade(1) = x : HPTPT(cf) : seq(i, (N - 1) /
    2, HPSVF(butterworthQ(N, i), cf));
};


//-- FEATURE EXTRACTIONS -------------------------------------------------------

// Spectral Centroid Dynamical Tracking
Spectral_centroid(tsc, x) = loop ~ _
with{
    loop(fc) =      x : onePoleTPT(fc) : 
                    (an.rms_envelope_tau(1.0), an.rms_envelope_tau(1.0), !) :
                    ro.cross(2) : - : integration(tsc) * ma.SR * 0.5 : 
                    clip(1, 20000, _)
    with{
        integration(t, x) = fi.pole(1.0, x * 
                            (1.0 / max(ma.EPSILON, t)) * 
                            (1.0 / ma.SR));
        // TPT Onepole Filter
        onePoleTPT(cf, x) = loop ~ _ : ! , si.bus(3)
        with {
            g = tan(cf * ma.PI * (1/ma.SR));
            G = g / (1.0 + g);
            loop(s) = u , lp , hp , ap
            with {
                v = (x - s) * G;
                u = v + lp;
                lp = v + s;
                hp = x - lp;
                ap = lp - hp;
            };
        };
    };
};

// Peak Envelope Follower
peakenvelope(period, x) = loop ~ _
    with {
        loop(y) = max(abs(x), y * coeff);
        twoPIforT = (2.0 * ma.PI) * (1.0 / ma.SR);
        coeff = exp(-twoPIforT / max(ma.EPSILON, period));
    };

// Peak Envelope Follower with attack and decay
peakEnvAttRel(att, rel, x) = peakenvelope(rel - att, x) :
    LPTPT(1.0 / max(ma.EPSILON, att));

// PeakHolder - holdTime in Seconds
peakHolder(holdTime, x) = loop ~ si.bus(2) : ! , _
with {
    loop(timerState, outState) = timer , output
    with {
        isNewPeak = abs(x) >= outState;
        isTimeOut = timerState >= (holdTime * ma.SR - 1);
        bypass = isNewPeak | isTimeOut;
        timer = ba.if(bypass, 0, timerState + 1);
        output = ba.if(bypass, abs(x), outState);
    };
};

// Integrator (Amplitude Follower)
integrator(seconds, x) = an.abs_envelope_tau(clip(0.001, 1000, seconds), x);

// Moving Average Envelope Follower
movingAverage(seconds, x) = x - x@(seconds * ma.SR) :
    fi.pole(1.0) / (seconds * ma.SR);

// Moving Average RMS
movingAverageRMS(seconds, x) = sqrt(max(0, movingAverage(seconds, x ^ 2)));

// Pitch Variance - Measures the magnitude of pitch changes in the input signal
pitchVariance = pitchTracker(1) <: 
    abs(log(_ + ma.EPSILON) - log(_ + ma.EPSILON : LPTPT(0.1))) : 
    delayFB(0.010, 0.95) : LPTPT(0.1)
with{
    // averaging function with 2pi time constant; t: averaging time in seconds
    avg(t, x) = y
    letrec {
        'y = x + alpha * (y - x);
    }
        with {
            alpha = exp((-2.0 * ma.PI * ma.T) / t);
        };

    // zero-crossing indicator function
    zeroCrossing(x) = (x * x') < 0;

    // zero-crossing rate
    zeroCrossingRate(t, x) = zeroCrossing(x) : avg(t);

    // Pitch tracker as zero-crossing rate of self-regulating lowpassed inputs
    // we highpass the input to avoid infrasonic
    // components to affect the measurements
    // we then clip the lowpass cutoff to improve stability
    pitchTracker(t, x) = loop ~ _
    with {
        loop(y) = fi.lowpass(4, max(80, y), xHighpassed) :
                  (zeroCrossingRate(t) * ma.SR * .5)
            with {
                xHighpassed = fi.highpass(1, 20, x);
            };
    };
};


//-- AUDIO PROCESSING ----------------------------------------------------------

// Delay with Feedback
delayFB(delSec, fb, x) = loop ~ _ : mem
with{
    loop(y) = clip(-1.0, 1.0, (y * fb + x)) @ 
                (max(0, ba.sec2samp(delSec) - 1));
};

// Larsen Auto Regulation (LAR) - 2 OUTS 
LAR(target) = _ <: 
    (((target - (_ : integrator(0.01) : delayFB(0.10, 0.99) : 
    LPButterworthN(2, 24) : clip(0.0, target)) :
    hgroup("LAR", inspect_VH(0.0, target)))), _) <: (_ * _), (_, !);

// lookahead limiter with: peakHolder, lowpass & peakenvelope
// All the credits of the original version goes to Dario Sanfilippo
LookaheadLimiter(threshold, holdSec, decaySec, x) =
    (x : peakHolder(holdSec) : LPTPT(1 / holdSec) : peakenvelope(decaySec)) :
        (min(1, threshold) / max(1, _)) *
            (x @ (holdSec * ma.SR));
            
// Sampler with Lookahead Envelope  - W&R Buffer 
// Control over Buff.Dim, Read Section and Read Speed
sampler(bufferLength, chunkLenght, grainPitch, x) = 
    (y @ lookahead) * envSmooth
with{
    // Sampler
    S = 192000 * 20;
    R = ma.SR * bufferLength;
    wIdx = ((_, int(R)) : %) ~ (_ + mem(1));
    ph = f ~ _
        with{
            f(s) = (s + grainPitch) % max(1, chunkLenght * R);
        };
    rIdx = ph % R;
    y = it.frwtable(3, S, .0, wIdx, x, rIdx);
    lookahead = .01 * ma.SR;
    needsEnv = ph < ph';
    envState = fi.pole(1.0 - needsEnv, 1.0 - needsEnv);
    envStep = 1.0 - (envState < lookahead);

    // Smooth Envelope - Sample Fade-IN / Fade-OUT
    envSmooth = LPCascadesampler(4, ma.SR / lookahead, envStep)
    with{
        // Onepole TPT Filter
        onePoleTPT(cf, x) = loop ~ _ : ! , si.bus(3)
        with {
            g = tan(cf * ma.PI * (1 / ma.SR));
            G = g / (1.0 + g);
            loop(s) = u , lp , hp , ap
            with {
                v = (x - s) * G;
                u = v + lp;
                lp = v + s;
                hp = x - lp;
                ap = lp - hp;
            };
        };
        // LP out from Onepole TPT Filter
        LPTPT(CF, x) = onePoleTPT(max(ma.EPSILON, min(20480, CF)), x) : 
            (_ , ! , !);
        // LP in series Cutoff Correction
        LPCascadesampler(N, CF, x) = x : seq(i, N, LPTPT(CFCorrection))
        with {
            CFCorrection = CF / sqrt(pow(2.0, 1.0 / N) - 1.0);
        };
    };
};

// PCM - Pulse Code Modulation with various transformations
pulses(x) = 
    x <: (env, _) :
    ((_ <: (_, _)), _) : (_, ro.cross(2)) :
    (_,
    (_ <: _ * 0.5, (noise : abs) * _ * 0.5),
    1 - _) :
    si.bus(2), ro.cross(2) :
    (_, _ + (_ * _)) :
    ((_, _) <: ((_ <: (_, _)), _) :
    (_ <: ((MapPHfreq : phasor0) > MapPHthrs)), 
    (MapMod : os.osc), 
    _ :
    _, ro.cross(2) :
    _ * _ * _)
with{
    noise = x : HPButterworthN(2, 10000) * (100000) : 
        _ % 1.0 : LPButterworthN(1, 20000);
    env = (integrator(0.001) : delayFB(0.2, 0.8) : 
        LPButterworthN(2, 24) : clip(0.0, 1.0) : inspect_VV(0.0, 1.0));
    MapPHfreq = ((1.0 - _) * 10);
    MapPHthrs = (_ * 0.4) + 0.1;
    MapMod = _ * 1000;
    // phasor that start from 0
    phasor0(f) = (_ <: _ + f, _) ~ _ % ma.SR : (!, _ / ma.SR);
};

// Pre filtering (HP, LP) - for Limit the band of the Input Signals
PreFilters(lowFC, hiFC) = HPTPT(lowFC) : LPTPT(hiFC) : TiltFilter
with{
    TiltFilter = hgroup("EQ", onepole)
    with{
        // Onepole with Schlick for FC
        tilt  = hslider("CF", 0.0, -1., 1., 0.001) * -1 : si.smoo;
        shape = 0.1;
        schlick(a, u) = u / (u + a * (1.0 - u));
        sign(x) = (x > 0.0) - (x < 0.0);
        bMag = 0.999 * schlick(shape, abs(tilt));
        b = bMag * sign(tilt) : max(-0.999) : min(0.999);
        onepole = *(1.0 - abs(b)) : +~*(b);
    };
};
//process = no.noise : PreFilters(0.01, 20000.0);

//------- ------------- ----- -----------
//-- NETWORK -------------------------------------------------------------------
//------- --------


//-- GLOBAL NETWORK ------------------------------------------------------------

Network(Networks_N, Agents_N, inAMP, samplerAMP, modAMP, pulsesAMP, outAMP) = 
hgroup("[0] Networks",
    // Networks in PAR. OP.
    par(k, Networks_N, 
        (vgroup("[0] Network %k",
            vgroup("[0] Pre-Processing", 
                _ * inAMP :
                PreFilters(0.01, 20000.0) : 
                LAR(1.0) 
            ) :
                // Agents in SEQ. OP.
                (seq(i, Agents_N, Agent(k, i)), _
            ) : vgroup("[0] Pre-Processing", _ * outAMP, _)
        ))
    )
)
with{
    // Single Agent Module in SEQ OP.
    Agent(k, i, x0) = tgroup("[0] Agents", 
    hgroup("[1] Agent %i", 
        (   
            x0 +
            (
                (x0 : sampling(0.5 + 0.3 * i + 0.15 * k)) * 
                (x0 : pitchVariation * samplerAMP) :
                modulation(120 + (i * i * 120 + k * (70 + i * 10)), 
                    10.0, modAMP)
            ) <:
                pwm * 
                (x0 : pitchVariation * pulsesAMP) + 
                _
        ), 
        x0 : 
            ((ro.cross(2) : fullBandNotch), 
            (x0 : spectralCentroid)) : 
                ro.cross(2) : crossoverNotch : 
                ((x0 : spectralCentroid), 
                _) : 
                    adaptiveHP
        )
    )
    with{
        // Detect how long is the same pitch (dynamic)
        pitchVariation(x) = hgroup("[6] Pitch Variance", 
            x : pitchVariance : inspect_VV(0.0, 1.0));

        // Detect current spectral centroid (dynamic)
        spectralCentroid(x) = hgroup("[7] Spectral Centroid", 
            x : Spectral_centroid(10.0) : inspect_VV(0.0, 20000.0));
        
        // Sampled, Delayed and transposed Input signal
        sampling(delSec, x) = hgroup("[4] Sampler", 
            x @ (ma.SR * delSec) : 
            sampler(10.0, 1.0, hgroup("[0] Frequency", pitchMAP)))
        with{
            pitchMAP = 4 ^ (movingAverageRMS(0.1, x) : si.smoo : inspect_VV(-1.0, 1.0)); 
        };

        // Amplitude Modulation in output from Sampler
        // Open and close the sidebands via triangular envelope
        modulation(fAM, fPH, gain, x) = hgroup("[5] AM",
            x * (1 - (os.osc(fAM) : bitoUni) * gain *
            (x : movingAverageRMS(0.1) * fPH : phasor0 : triangularFunc : 
            inspect_VV(0.0, 1.0))));

        pwm = hgroup("[5] PWM", 
            pulses);
        
        // Remove noise from phase increments
        gate(x) = ((LPTPT(1, abs(x)) > 0.001) : LPTPT(100)) * x;
        
        // Phase increment for Notch Filter FC via moving AVG
        phaseIncrement(tRMS, minv, maxv, freq, x) = x : 
            movingAverageRMS(tRMS) : gate * freq : phasor0 * 
            (maxv - minv) + minv : hgroup("[0] FC", inspect_VH(0.0, 20000.0));
        
        // Notch with FC that scan the full Frequency Band 
        fullBandNotch(phaseFC, x) = vgroup("[8] Notch Filters", 
            hgroup("[0] Full Band",
            (phaseFC : phaseIncrement(0.1, 1.0, 20000.0, 4.0)), x :
            NotchSVF(1.0)));

        // 2 Notch with FC that scan two bands divided by Spectral Centroid 
        crossoverNotch(dynaRange, x) = vgroup("[8] Notch Filters", 
        vgroup("[1] Crossover",
            hgroup("[0] Low", ((x : movingAverageRMS(0.1) : 
                phaseIncrement(0.1, 1.0, dynaRange, 4.0)), x) : 
                NotchSVF(1.0)),
            hgroup("[1] High", ((x : movingAverageRMS(0.1) : 
                phaseIncrement(0.1, dynaRange, 20000.0, 4.0)), x) : 
                NotchSVF(1.0)) :> 
            + : _ * 0.5
        ));

        // Cut low frequencies based on Spectral Centroid (dynamic)
        adaptiveHP(fc, x) = ((fc : Map : 
            hgroup("[9] Adaptive HP", hgroup("[0] FC", inspect_VV(0.0, 20000.0)))), x) : 
            HPButterworthN(4)
        with{
            Map(x) = 200 * pow(100 / (clip(20, 20000, x)), 0.5);
        };
    };
};


//-- GUI -----------------------------------------------------------------------

// Graphic User Interface Controls for the Network
GUI_Controls = 
vgroup("Controls",
    (hslider("[0] Gain Inputs", 1.0, 0.0, 1.0, 0.001) : si.smoo),
    (hslider("[1] Gain Samplers", 1.0, 0.0, 1.0, 0.001) : si.smoo), 
    (hslider("[2] Gain Modulations", 1.0, 0.0, 1.0, 0.001) : si.smoo),
    (hslider("[3] Gain Pulses", 1.0, 0.0, 1.0, 0.001) : si.smoo),
    (hslider("[4] Gain Outputs", 1.0, 0.0, 1.0, 0.001) : si.smoo)
);


//-- OUTPUT --------------------------------------------------------------------
///*
// Main Function
process =  vgroup("[0] Local InSTrument ENvironment", 
    (GUI_Controls, si.bus(2)) : Network(2, 4) : stereoCouplingMatrix)
with{
    stereoCouplingMatrix(x1, x1LAR, x2, x2LAR) = 
        (x1 + x2 * x1LAR : LookaheadLimiter(1.0, 0.1, 0.1)), 
        (x2 + x1 * x2LAR : LookaheadLimiter(1.0, 0.1, 0.1));
    directOutputs(x1, x1LAR, x2, x2LAR) = 
        (x1 : LookaheadLimiter(1.0, 0.1, 0.1)), 
        (x2 : LookaheadLimiter(1.0, 0.1, 0.1));
};
//*/