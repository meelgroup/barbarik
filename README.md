# Barbarik

Barbarik is a framework developed to test whether a sampler is almost uniform or not. Currently it is implemented for testing QuickSampler, STS, Unigen. It uses SPUR as the underlying uniform sampler. This work is by Kuldeep Meel and Sourav Chakraborty, as published in [AAAI'19](https://www.comp.nus.edu.sg/~meel/Papers/aaai19-cm.pdf).  

## Getting Started
To get started either download the ZIP file from the repository or git clone it using the following:
```
git clone https://github.com/meelgroup/barbarik.git
```

### Prerequisites
Please install the following in the same directory as the repository :
* [SPUR](https://github.com/ZaydH/spur) - Perfectly Uniform Satisfying Assignments
* [Unigen](https://bitbucket.org/kuldeepmeel/unigen) - almost-uniform SAT sampler
* [Quick Sampler](https://github.com/RafaelTupynamba/quicksampler)
* [STS](http://cs.stanford.edu/~ermon/code/STS.zip)

### Installing

After you have installed the above mentioned tools and cloned the repository you are good to go. Inorder to run a custom sampler please make appropriate changes to the code to adjust according to need. To make changes for CustomSampler type look for the following tag in ```Verifier.py``` file:
```
# @CHANGE_HERE : please make changes in the below block of code
```

## Running the tests

In order to run tests please keep your cnf files according to DIMACS format in the same directory as of the file run_verifier.py.
To run use the following command:
```
python run_verifier.py --sampler SAMPLER_TYPE
```
SAMPLER_TYPE takes the following values:
* UniGen = 1
* QuickSampler = 2
* STS = 3
* CustomSampler = 4

Please make appropriate changes to code for CustomSampler

## Comparison with SearchTreeSampler

For easy comparison we provide a (modified) binary for SearchTreeSampler. [Source](https://github.com/RafaelTupynamba/quicksampler/tree/master/STS).  

## Comparison with SCALMC

We have added a binary of scalmc in the samplers directory. It implements the UniGen2 algorithm of Chakraborty, Fremont, Meel, Seshia, and Vardi and the ApproxMC3 algorithm of Meel, and Soos. The particular implementation used is based on a prototype by Mate Soos and Kuldeep Meel which is pending publication.

## How to Cite

If you use Barbarik, please cite the following paper : [AAAI'19](https://www.comp.nus.edu.sg/~meel/Papers/aaai19-cm.pdf). Here is [BIB file](https://www.comp.nus.edu.sg/~meel/bib/CM19.bib)

## Contributors
1. Kuldeep S. Meel
2. Shayak Chakraborty 
3. Yash Pote

If you cannot use any of the available binaries, or experience any other problems, please contact us.
