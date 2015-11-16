all: executable
all: CFLAGS = -o3

debug: executable
debug: CFLAGS = -DDEBUG

OBJS = DnaBase.o DnaTrieBuilder.o main.o

executable: $(OBJS)
	g++ -o DnaTrieBuilder.exe $(OBJS)

%.o : %.cpp
	g++ -g $(CFLAGS) -o $@ -c $<

main.o: DnaTrieBuilder.h DnaBase.h stdafx.h
DnaTrieBuilder.o: DnaTrieBuilder.h DnaBase.h stdafx.h
DnaBase.o: DnaBase.h stdafx.h

clean:
	rm *.o DnaTrieBuilder.exe
