# Threaded Jacobi Iteration
This program performs a threaded Jacobi iteration using barrier synchronization.

## Usage
### Compilation 
```bash
javac Jacobi.java
```
### Execution
```bash
java Jacobi [num_threads]
```

If no arguments are supplied then it assumes only 1 thread should be used and is therefore executed in serial.