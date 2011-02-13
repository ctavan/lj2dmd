CC = g++ -O3
CFLAGS = -o
LIBS = -g -Wall -lm
PROGRAMNAME = leapfrog

# Standardaufruf durch 'make' oder 'make default'
default: main.c main.o
	$(CC) $(CFLAGS) $(PROGRAMNAME) main.o $(LIBS)

# Abhaengigkeiten
main.o: main.c ran.h

# Erzeugen oder Aktualisieren einer object-Datei
# $< wird durch den Namen der c-Datei ersetzt
.c.o:
	$(CC) -c -g -Wall $<

# Loeschen der object-Dateien und der Programmdatei
clean:
	rm *.o *~ $(PROGRAMNAME) 2>/dev/null
