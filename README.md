# Barbarik, a testing framework for (almost) uniform samplers

Barbarik is a framework developed to test whether a sampler is almost uniform or not. It uses SPUR as the underlying uniform sampler. This work is by Sourav Chakraborty and Kuldeep S. Meel, as published in [AAAI'19](https://www.comp.nus.edu.sg/~meel/Papers/aaai19-cm.pdf).

## Getting Started

Run:
```
git clone --depth 1 https://github.com/meelgroup/barbarik.git
cp my_favourite_cnf.cnf.gz barbarik/
cd barbarik
./barbarik.py --seed 1 --sampler SAMPLER_TYPE blasted_case110.cnf
```

Where  SAMPLER_TYPE takes the following values:
* UniGen2 = 1
* QuickSampler = 2
* STS = 3
* CustomSampler = 4
* AppMC3 = 5

### Samplers used

In the "samplers" directory, you will find 64-bit x86 Linux compiled binaries for:
* [ApproxMC3-with-sampling](https://github.com/meelgroup/ApproxMC/tree/master-with-sampling) - an almost-uniform sampler (This is a beta version of UniGen3 -- which will be released soon. If you use ApproxMC3 binary, please cite UniGen paper to avoid any confusion.)
* [UniGen2](https://bitbucket.org/kuldeepmeel/unigen/) - an almost-uniform sampler, version 2
* [SPUR](https://github.com/ZaydH/spur) - Perfectly Uniform Satisfying Assignments
* [Quick Sampler](https://github.com/RafaelTupynamba/quicksampler)
* [STS](http://cs.stanford.edu/~ermon/code/STS.zip)

### Custom Samplers

To run a custom sampler, make appropriate changes to the code -- look for the following tag in `barbarik.py` file: `# @CHANGE_HERE : please make changes in the below block of code`

## How to Cite

If you use Barbarik, please cite the following paper : [AAAI'19](https://www.comp.nus.edu.sg/~meel/Papers/aaai19-cm.pdf). Here is [BIB file](https://www.comp.nus.edu.sg/~meel/bib/CM19.bib)

## Contributors
1. Kuldeep S. Meel
2. Sourav Chakraborty
3. Shayak Chakraborty 
4. Yash Pote
5. Mate Soos
5. Priyanka Golia
