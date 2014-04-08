# Methods to Estimate Instantaneous Energy (including the nonlinear energy operator)


Collection of M-files (computer code) to implements instantaneous energy measures, as
described in [[1]](#references).
Requires Matlab or Octave (programming environments).


# Contents
- [quick start](#quick-start)
- [requirements](#requirements)
- [test computer setup](#test-computer-setup)
- [licence](#licence)
- [references](#references)
- [contact](#contact)


# quick start
TO DO....

### files
All Matlab files (.m files) have a description and an example in the header. To read this
header, type `help <filename.m>` in Matlab.  Directory structure is as follows: 
```
.
├── bias_of_estimators.m
├── cal_freqweighted_energy.m
├── compare_nleo_methods.m
├── discrete_Hilbert_diff_operator.m
├── do_bandpass_filtering.m
├── general_nleo.m
├── nleo_parameters.m
├── pics/                              # directory for figures
├── plot_eeg_examples.m
├── properties_test_Hilbert_NLEO.m
└── script_test_eeg_data.m
```



# requirements
Either Matlab (R2012 or newer,
[Mathworks website](http://www.mathworks.co.uk/products/matlab/)) or Octave (v3.6 or
newer, [Octave website](http://www.gnu.org/software/octave/index.html), with the
'octave-signal' add-on package).



# test computer setup
- hardware:  Intel(R) Xeon(R) CPU E5-1603 0 @ 2.80GHz; 8GB memory.
- operating system: Ubuntu GNU/Linux x86_64 distribution (Raring, 13.04), with Linux kernel 3.5.0-28-generic 
- software: Octave 3.6.4 (using Gnuplot 4.6 patchlevel 1), with 'octave-signal' toolbox and Matlab (R2009b, R2012a, and R2013a)

---

# licence

```
Copyright (c) 2014, John O' Toole, University College Cork
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

  Neither the name of the University of Deusto nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```


# references


---

# contact

John M. O' Toole

Neonatal Brain Research Group,  
Irish Centre for Fetal and Neonatal Translational Research,  
Department of Paediatrics and Child Health,  
University College Dublin,  
Western Gateway Building, Room 2.17,  
Cork, Ireland


