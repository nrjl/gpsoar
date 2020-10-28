# gpsoar - Simulation tools for soaring flight and wind estimation

This matlab code was developed for simulating soaring flight and for wind estimation.

**Disclaimer:** This code is old, was never well written, and isn't maintained. Depending on interest I could try to rearrange and update it, but for now it's provided as is. You're welcome to ask for help and I'll try my best to explain or possibly eventually clean it up.


## Getting started

You should be able to run `GP_sim.m` to demonstrate the basic operation of the system as demonstrated in the ICRA11 paper (then `plot_result`, then `double_video_thesis` to make nice videos like this)

![soaring_thesis](https://user-images.githubusercontent.com/10678827/97481993-dae09d00-1955-11eb-9f4e-45316ad45a9e.gif)

## Background

The code in this repo formed the basis of the following papers in [JGCD](https://arc.aiaa.org/doi/10.2514/1.52236) and [ICRA 2011](https://ieeexplore.ieee.org/abstract/document/5979966):
```bibtex
@INPROCEEDINGS{lawrance2011,
  author = {Lawrance, Nicholas R. J. and Sukkarieh, Salah},
  booktitle={2011 IEEE International Conference on Robotics and Automation}, 
  title={Path planning for autonomous soaring flight in dynamic wind fields}, 
  year={2011},
  pages={2499--2505},
  doi={10.1109/ICRA.2011.5979966}}
  
@article{doi:10.2514/1.52236,
  author = {Lawrance, Nicholas R. J. and Sukkarieh, Salah},
  title = {Autonomous Exploration of a Wind Field with a Gliding Aircraft},
  journal = {Journal of Guidance, Control, and Dynamics},
  volume = {34},
  number = {3},
  pages = {719--733},
  year = {2011},
  doi = {10.2514/1.52236},
  URL = {https://doi.org/10.2514/1.52236},
  eprint = {https://doi.org/10.2514/1.52236}}
```

## Feedback
Comments, feedback and PRs for changes are welcome (nicholas.lawrance AT mavt.ethz.ch).
