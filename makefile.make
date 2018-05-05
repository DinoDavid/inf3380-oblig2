CC = gcc
CFLAGS = -O2

PROJ = serial_main
OBJS = serial_main.o

all : simple-jpeg $(PROJ)

serial_main : $(OBJS)
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
