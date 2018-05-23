CC = mpicc
CFLAGS = -O2 -Wall -Wextra -Wno-unused-result -fopenmp
LDFLAGS = -lm

PROJ = oblig2
OBJS = oblig2.o

$(PROJ): $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

run: $(PROJ)
	mpirun --hostfile /etc/openmpi/openmpi-default-hostfile -np 4 ./oblig2 large_matrix_a.bin large_matrix_b.bin c.bin
	./compare c.bin large_matrix_c.bin

delivery:
	git clean -Xdf
	tar -czf oblig2.tar.gz	./*

printer: printer.o
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

compare: compare.o
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	git clean -Xdf
