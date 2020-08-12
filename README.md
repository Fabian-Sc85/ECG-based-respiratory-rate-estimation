# ECG-based respiratory rate estimation
## Introduction
The purpose of this collection of functions is the indirect estimation of the respiratory rate from ECG signals. It accomplishes this by implementing several algorithms published by us ([Laboratory for Biosignal Processing](labp.github.io/)) or third parties. These algorithms exploit two physiological phenomena:
1. Respiratory Sinus Arrhythmia: A normal variation of the heart rate during inhalation and exhalation. The variation is releated to changes in intra-thoracic pressure during the breathing cycle [[1]](#1).
2. Variations of the R-Peak in the ECG: Impedance changes of the chest during inhalation and exhalation. During inhalation the chest is moving up. This changes the position of the ECG-electrodes with respect to the heart. The amplitude of the ECG-signal is diminished which is most pronounced in the R-peak amplitude.

## Installation
No installation is required. Just clone the repository to a folder of your choice and add it to your MATLAB path.

## Usage

## Algorithm

## Coming Soon
1. Exmample data
2. Demo script to showcase the functionality

## References
<a id="1">[1]</a> Berntson, Gary G., John T. Cacioppo, and Karen S. Quigley. "Respiratory sinus arrhythmia: autonomic origins, physiological mechanisms, and psychophysiological implications." Psychophysiology 30.2 (1993): 183-196.
