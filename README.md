# Barbarik

Barbarik is a framework developed to test whether a sampler is almost uniform or not. Currently it is implemented for testing QuickSampler, STS, Unigen, Solver(CryptoMiniSat). It uses MUSE as the underlying uniform sampler. This work is by Kuldeep Meel and Sourav Chakraborty, as published in [AAAI'19](https://www.comp.nus.edu.sg/~meel/Papers/aaai19-cm.pdf).  

## Getting Started
To get started either download the ZIP file from the repository or git clone it using the following:
```
git clone https://github.com/meelgroup/barbarik.git
```

### Prerequisites
Please install the following in the same directory as the repository :
* [MUSE](https://github.com/ZaydH/spur) - Perfectly Uniform Satisfying Assignments
* [CrytoMiniSat](https://github.com/msoos/cryptominisat) - SAT solver implemented
* [STS](https://github.com/meelgroup/khatu/blob/master/STS)
* [Unigen](https://bitbucket.org/kuldeepmeel/unigen) - almost-uniform sat sampler
* [Quick Sampler](https://github.com/RafaelTupynamba/quicksampler)

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
* Solver = 2
* QuickSampler = 3
* STS = 4
* OtherSampler = 5

Please make appropriate changes to code for OtherSampler

## How to Cite

If you use Barbarik, please cite the following paper : [AAAI'19](https://www.comp.nus.edu.sg/~meel/Papers/aaai19-cm.pdf)
