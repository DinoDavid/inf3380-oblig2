CC = mpicc
CFLAGS = -O2 -Wall -Wextra -Wno-unused-result 
LDFLAGS = -lm

PROJ = oblig2
OBJS = oblig2.o

$(PROJ): $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

run:
	./oblig2 small_matrix_a.bin small_matrix_b.bin c.bin

obj-clean:
	$(RM) *.o

exec-clean:
	$(RM) $(PROJ)

autosave-clean:
	$(RM) *~

clean:
	$(MAKE) obj-clean
	$(MAKE) exec-clean
	$(MAKE) autosave-clean
