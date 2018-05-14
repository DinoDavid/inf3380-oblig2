CC = mpicc
CFLAGS = -O2 -Wall -Wextra -Wno-unused-result

PROJ = oblig2
OBJS = oblig2.o

$(PROJ): $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

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
