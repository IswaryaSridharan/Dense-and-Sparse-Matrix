MPICC?=mpicc
RM = rm

CFLAGS = -o3

TARGET = DenseParallel_v5
TARGET1 = SparseParallel_v6
TARGET2 = Sparse3937
TARGET3 = Sparse1374

all: $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET): $(TARGET).c
	$(MPICC) $(CFLAGS) -o $(TARGET) $(TARGET).c
$(TARGET1): $(TARGET1).c
	$(MPICC) $(CFLAGS) -o $(TARGET1) $(TARGET1).c
$(TARGET2): $(TARGET2).c
	$(MPICC) $(CFLAGS) -o $(TARGET2) $(TARGET2).c
$(TARGET3): $(TARGET3).c
	$(MPICC) $(CFLAGS) -o $(TARGET3) $(TARGET3).c

clean: $(RM) $(TARGET)
