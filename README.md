# LISTEN (Local InSTrument ENvironment)

A live-electronics work in which autonomous sound agents live on a physical object: one orchestral thunder-sheet, driven into self-oscillation through feedback between piezoelectric contact microphones and exciters.
The agents are built around the sheet's own vibrational modes: listening to the surface through DSP feature extraction, the system responds to its own activity and regulates itself, moving between stable regimes and nonlinear, chaotic ones.
The thunder-sheet is at once instrument, environment and sonic (eco)system; the microphones, projected through loudspeakers, carry the local space of the object into the listening space of the audience.

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
