FLOAT_TYPE=USE_SINGLE

CC=/usr/local/cuda/bin/nvcc
CFLAGS=-arch=sm_20 -rdc=true
LDFLAGS=-arch=sm_20
DEFINES=-D$(FLOAT_TYPE)
INCLUDES=-Iinclude -Itest
ifeq ($(OS),Windows_NT)
LIBS=-largtable2
else
LIBS=-largtable2 -lm
endif
LIBDIR=-Llibs
LIBS=-largtable2

CU_FILES=main.cu Iterate.cu CellFunctions.cu ComputeResiduals.cu FilesReading.cu FilesWriting.cu \
         ShellFunctions.cu GpuInit.cu GpuBoundaries.cu GpuCollision.cu GpuStream.cu LogWriter.cu \
         ArrayUtils.cu Arguments.cu GpuSum.cu Check.cu GpuUpdateMacro.cu \
         test/CuTest.cu test/TestComputeResiduals.cu test/TestGpuBoundaries.cu \
         test/TestGpuInit.cu test/TestUtils.cu test/TestGpuUpdateMacro.cu test/TestGpuSum.cu \
         test/TestGpuCollision.cu test/TestGpuStream.cu test/TestIterate.cu
ITER_FILES=Iterate.cu GpuInit.cu ComputeResiduals.cu GpuBoundaries.cu GpuCollision.cu GpuStream.cu \
           ArrayUtils.cu Arguments.cu GpuSum.cu Check.cu GpuUpdateMacro.cu \
           test/TestGpuUpdateMacro.cu test/TestGpuBoundaries.cu test/TestGpuInit.cu test/TestGpuSum.cu \
           test/TestGpuCollision.cu test/TestGpuStream.cu test/TestIterate.cu
ITER_FILE=IterateCombined.cu
RLSE_FILES=main.cu $(ITER_FILE) CellFunctions.cu FilesReading.cu FilesWriting.cu \
           ShellFunctions.cu LogWriter.cu \
           test/CuTest.cu test/TestComputeResiduals.cu test/TestUtils.cu
HEADRF=FilesReading.h FilesWriting.h Iterate.h CellFunctions.h ComputeResiduals.h ShellFunctions.h \
       GpuFunctions.h LogWriter.h Arguments.h ArrayUtils.h GpuSum.h
HEADERS=$(patsubst %,include/%,$(HEADRF))
ifeq ($(OS),Windows_NT)
OBJ_FILES=$(patsubst %,$(OBJ_DIR)/%,$(CU_FILES:.cu=.obj))
else
OBJ_FILES=$(patsubst %,$(OBJ_DIR)/%,$(CU_FILES:.cu=.o))
endif
OBJ_DIR=obj
DOC_DIR=docs

EXEC=LBMSolver
RLSE=LBMRel

all: $(EXEC)
$(EXEC): $(OBJ_DIR) $(OBJ_FILES); \
	$(CC) $(LDFLAGS) $(LIBDIR) $(LIBS) -o $@ $(OBJ_FILES)

debug: CFLAGS+= -g -G -lineinfo
debug: DEFINES+= -DDEBUG
debug: all

release: CFLAGS=-arch=sm_20 -O3
release: DEFINES+= -DRELEASE
release: $(RLSE)

$(RLSE): $(RLSE_FILES); \
	$(CC) $(CFLAGS) $(LIBDIR) $(LIBS) $(DEFINES) $(INCLUDES) -o $@ $(RLSE_FILES)

$(ITER_FILE): $(ITER_FILES); \
	./combine.sh $(ITER_FILE) $(ITER_FILES)

test: all; \
	cuda-memcheck --leak-check full ./$(EXEC) --test

$(OBJ_DIR):; \
	mkdir -p $(OBJ_DIR)/test

$(OBJ_DIR)/%.o: %.cu; \
	$(CC) $(CFLAGS) $(LIBS) $(DEFINES) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.obj: %.cu; \
	$(CC) $(CFLAGS) $(LIBS) $(DEFINES) $(INCLUDES) -c $< -o $@

$(DOC_DIR):; \
	mkdir -p $(DOC_DIR)

doc: include/*.h test/*.h Doxyfile $(DOC_DIR) README.md; \
	doxygen Doxyfile

latexdoc: doc; \
	make -C $(DOC_DIR)/latex

.PHONY: clean

clean:; \
	rm -rf $(OBJ_DIR) $(EXEC) $(RLSE) $(ITER_FILE) $(DOC_DIR)
