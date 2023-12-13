# Dynex SDK Local Testnet Sampler
Dynex is the world’s first neuromorphic supercomputing blockchain based on the DynexSolve chip algorithm, a Proof-of-Useful-Work (PoUW) approach to solving real-world problems. The Dynex SDK is used to interact and compute on the Dynex Platform. Use this repository to enable "testnet=False" functionality in your Dynex SDK. It allows sampling of Qubo/Ising computing problems on the local machine before computing on the Dynex Neuromorphic Computing cloud. Mainly used for prototyping and testing of code before incurring costs.

# Build from source & Installation

Build the binary with the following commands and copy it into the folder /tmp in the directory where your Python program is located:

```
git clone https://github.com/dynexcoin/dynexsdk_testnet.git
./build.sh
cp dynex-testnet-bnb <PATH-OF-YOUR-SDK-PROGRAM>/tmp
```

# Usage

To enable sampling on the local machine, specify the parameter "testnet=False":

```
model = dynex.BQM(bqm, logging=True);
sampler = dynex.DynexSampler(model,  mainnet=False);
sampleset = sampler.sample(num_reads=20000, annealing_time = 200);
```

