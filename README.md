# Input-rate-encoding-and-gain-control-in-dendrites-
This repository contains the code from Matlab and IgorPro environments that generates a simple input rate to voltage output model of branch integration.

These files in this repository can be used to recreate Figures 7 and 8 of "Input rate encoding and gain control in dendrites of neocortical pyramidal neurons" by Nikolai C. Dembrow and William J. Spain, as well as simulate how temporal patterns of input onto a single branch will be transformed into voltage.


This is an open source software. The files are copyright by Nikolai C. Dembrow and William J. Spain. They are distributed under the terms of the GNU General Public License (GPL) version 3 (or later). See "Liscence.txt" file for legal details.

Contact: ndembrow@uw.edu

Usage guide:
To utilize this code, you need IgorPro 6 or later, and Matlab 2012b or later. Generate stimulation patterns in Matlab using "IPT_forModel.m", which are saved as StimRun.txt (stimulation times and spine number stimulated) as well as FreqOut.txt (the mean rate generated) files. Import these into IgorPro, along with the Igor procedure file "Simple_InputRate_to_Voltage_Model_andAnalysis.ipf". Execute this procedure to generate the simple model output.


Feedback:
We welcome and encourage feedback on this code. Any improvements/challenges that you encounter can help us to develop better versions of this simple code.

Kindly cite our Cell Reports paper if you use this code (even portions of it) in your work.



