# NS-landau-quantization
Public-facing repository for emissivity and opacity calculations incorporating Landau quantization

IMPORTANT: Follow [this link](https://drive.google.com/file/d/1eq-s-U3L3JoC09BFto6E5tze6yG8MZ5i/view?usp=drive_link) to download `inr_output_250_250.csv`, the lookup table for the `I` function. Then place it in the `utils` folder. (I would've used Github LFS, but I have no money and don't want to pay $5/month just for this)

This repository contains everything you need to replicate our results in our [paper](https://arxiv.org/abs/2412.02925). You can clone it and use a jupyter notebook, or you can run it from the command line.

The two most important functions are `Urca.urca_rate_wrapper(mu_B, B, T)`, which returns
* the Direct Urca rate
* the Modified Urca rate
* the quasiclassical Direct Urca rate
* the density (using the IUFSU* equation of state)
* and the number of available Landau levels (proton spin up, proton spin down, and electron),

and `CrossSection.cross_section(channel, n_B, Y_e, k_nu, cos_θnu, B, T)`, which returns the cross-section of neutrino capture onto nucleons. (Both of these functions also have optional keyword arguments.) 


## Cite this work

Please use the following BibTeX entry to cite our work:

```bibtex
@article{KumamotoWelch:2024lq,
      title={Effects of Landau quantization on neutrino emission and absorption}, 
      author={Mia Kumamoto and Catherine Welch},
      year={2024},
      eprint={2412.02925},
      archivePrefix={arXiv},
      primaryClass={nucl-th},
      reportNumber={INT-PUB-24-058},
      url={https://arxiv.org/abs/2412.02925}, 
}
```
