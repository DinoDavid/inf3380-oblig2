CC = mpicc
CFLAGS = -O2 -Wall -Wextra -Wno-unused-result
LDFLAGS = -lm

PROJ = oblig2
OBJS = oblig2.o

$(PROJ): $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

run:
	./oblig2 small_matrix_a.bin small_matrix_b.bin c.bin
	./compare c.bin small_matrix_c.bin 

printer: printer.o
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

compare: compare.o
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	git clean -Xdf
