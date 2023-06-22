[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Barbarik, a testing framework for samplers

'Barbarik' is a framework developed to test whether a sampler samples from a target distribution. To read more about Barbarik, have a look at our papers [AAAI'19](https://www.comp.nus.edu.sg/~meel/Papers/aaai19-cm.pdf), [NeurIPS'20](https://arxiv.org/abs/2010.12918), and [NeurIPS'22](https://www.comp.nus.edu.sg/~meel/Papers/neurips22.pdf).

## Getting Started

To test QuickSampler:
```
git clone --depth 1 https://github.com/meelgroup/barbarik.git
cd barbarik
git submodule update --init --recursive
python3 barbarik.py --seed 1 --sampler 2 blasted_case110.cnf
```

See `python3 barbarik.py --help` for the different samplers supported.

### Samplers used

You can choose any of the samplers in the "samplers" directory, see `--help`:
* [UniGen2](https://bitbucket.org/kuldeepmeel/unigen/) - an almost-uniform sampler, version 2
* [ApproxMC3-with-sampling](https://github.com/meelgroup/ApproxMC/tree/master-with-sampling) - an almost-uniform sampler (This is a beta version of UniGen3 -- which will be released soon. If you use ApproxMC3 binary, please cite UniGen paper to avoid any confusion.)
* [SPUR](https://github.com/ZaydH/spur) - Perfectly Uniform Satisfying Assignments
* [Quick Sampler](https://github.com/RafaelTupynamba/quicksampler)
* [STS](http://cs.stanford.edu/~ermon/code/STS.zip)

### Custom Samplers

To run a custom sampler, make appropriate changes to the code -- look for the following tag in `barbarik.py` file: `# @CHANGE_HERE : please make changes in the below block of code`

## How to Cite

If you use Barbarik, please cite the following papers : [AAAI'19](https://www.comp.nus.edu.sg/~meel/publications/CM19.bib), [NeurIPS'20](https://www.comp.nus.edu.sg/~meel/publications/MPC20.bib), and [NeurIPS'22](https://www.comp.nus.edu.sg/~meel/publications/PM22.bib).

## Contributors
1. Kuldeep S. Meel
2. Sourav Chakraborty
3. Shayak Chakraborty 
4. Yash Pote
5. Mate Soos
5. Priyanka Golia
