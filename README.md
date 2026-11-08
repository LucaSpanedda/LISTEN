# LISTEN (Local InSTrument ENvironment)

is a compositional project for live electronics that investigates autonomous sound agents and complex adaptive systems on physical objects. 
The work explores how emergent sonic behaviours can arise from the interaction between material properties driven into self-oscillation through feedback processes.
In this case, the project focuses on the vibrational properties of one orchestral thunder-sheet. 
The compositional methodology, rather than imposing predetermined musical structures, 
develops autonomous agents directly around the emergent properties of the object itself
Under feedback conditions, different vibrational modes emerge according to specific amplitude ranges: 
in some cases behaving in relatively stable and linear ways, while in others undergoing phase transitions toward nonlinear and chaotic regimes.
Through mediated “LISTENing” processes (via DSP feature extraction and sound processing), 
the system autonomously responds to its own internal activity through its physical surface itself.
The technical structure of the system consists of distributed transducers, 
including piezoelectric contact microphones and exciters (or other kinds) positioned at different locations on the object surfaces.
The title Local InSTrument ENvironment refers to the idea of constructing a localized complex adaptive sonic (eco)system 
in which agencies are distributed across a physical object rather than across the entire performance space. 
In this sense, the thunder-sheet functions simultaneously as instrument, environment, and active sonic (eco)system.
A central aspect of the work concerns the relationship between local and global spaces. 
Distributed microphones are positioned around the setup, operating not simply as amplification devices, 
but as localized points of observation capable of capturing different perspectives on the system’s internal activity. 
Their signals are spatially projected through loudspeakers into the performance environment, 
effectively transporting the local acoustic space of the object into the listening space that hosts the audience.
Through this process, the concert hall (or other spaces) itself becomes a projection of the adaptive system, 
making the audience experience the distributed internal dynamics of the network as an immersive spatial condition.

## Status

Version performed at the ZHdK Lange Nacht, Zurich, 28 May 2026: one thunder-sheet, with a pair of piezoelectric contact microphones (inputs) and a pair of exciters (outputs).
The control "Divergence/Convergence" acts both on the piezo pair and on the exciters: at 0 the two networks are independent, at 1 both receive the same mix of the two channels.
A state of the later two-plate version (September 2026), developed for CIM 2026 (L'Aquila, 15 October 2026), remains in the history of this repository (commit 665ecaf).

## Files

- `LISTEN.dsp` — the system (library and network in a single file)

## Usage

Compile with [Faust](https://faust.grame.fr) as a standalone application, in double precision:

- macOS (CoreAudio): `faust2caqt -double -preset auto LISTEN.dsp`
- Linux (JACK): `faust2jaqt -double -preset auto LISTEN.dsp`
- Linux (ALSA, without JACK): `faust2alqt -double LISTEN.dsp`

`-double` computes the whole system in double precision; `-preset auto` adds a preset manager to the interface.

Without compiling (Windows, macOS, Linux): open `LISTEN.dsp` in [FaustLive](https://github.com/grame-cncm/faustlive/releases), which compiles and runs it directly. For a quick test, the online [Faust IDE](https://faustide.grame.fr) runs it in the browser.

The system has 2 inputs (piezo microphones) and 2 outputs (exciters).

## License

MIT License (see `LICENSE`).
