# Distribute Key Generation experiments
This is a source code for the paper: Practical Distributed Key Generation Using Two-Round Secret Sharing.


This repository comprises:

**1. Smart contract**
- The smart contract can be deployed to the blockchain. Before deployment, the zk-SNARK verification keys inside the smart contract should be changed to the zk-SNARK keys generated by the HLDS programs.

**2. zk-SNARK HLDS programs in Zokrates**
- Generate zk-SNARK keypair: (the verification keys should be embedded to the smart contract)
```
zokrates compile -i <name>.zok
zokrates set -s gm17
```
- Generate proof:
```
zokrates compute-witness -a <inputs of name.zok>
zokrates generate-proof -s gm17
zokrates print-proof --format remix
```
- The output is the proof string $\pi$ and public input/output ($\vec{v},\vec{o}$). They can be passed directly to the smart contract functions \texttt{verifyC1}, \texttt{verifyC2}, etc.
**3. Two-Round secret sharing (TSS)**
- tss.py generated using the library: https://github.com/Zokrates/pycrypto.
- You can modify parameters or the threshold number in tss.py. By default, the threshold is $t=9$ and $n=16$.
- tss.py generate (1) public parameter set, (2) function $\mathsf{ComX^{TSS}}$ for $n$ parties, (3) Function $\mathsf{GenY^{TSS}}$ for $n$ parties, and (4) function $\mathsf{Rec^{TSS}} using $t$ shares.
- to test TSS functions:
```
python3 DKGTSS/tss.py
```
