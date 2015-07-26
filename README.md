Simulated Branch Predictor for a Sodor Five Stage CPU Core
=============
This repository contains some of the code I wrote for a class I took called [Computer Architecture and Engineering](http://www-inst.eecs.berkeley.edu/~cs152/sp14/).
The complete instructions are contained within the PDF.

Summary of Instructions
----------
Using C++, implement two [branch predictors](https://en.wikipedia.org/wiki/Branch_predictor/ "Wikipedia: Branch Predictor") for an emulated Sodor 5-Stage core processor. The first branch predictor can use infinite memory. The second branch predictor is limited to 1024 bytes of total state, and each seperate memory structure increases the CPU clock cycle time by 0.25%.