# CCOMP 	- 	C compiler
# GPUCOM 	- 	CUDA compiler
# UPCCOM	-	UPC compiler
#
# THREADS	-	number of threads for the UPC code

CCOMP = gcc
GPUCOMP = nvcc
UPCCOMP = upcc

THREADS=4

default:
	@echo "Enter one of the following:"
	@echo "  make serial		for serial code"
	@echo "  make gpu		for CUDA code"
	@echo "  make upc		for UPC code"

serial: *.c
	$(CCOMP) *.c -lm -o LBMSolver

gpu: *.cu
	$(GPUCOMP) -arch=sm_20 *.cu -lm -o LBMSolver

upc: *.c
	$(UPCCOMP) -T=$(THREADS) *.c -lm -o LBMSolver
