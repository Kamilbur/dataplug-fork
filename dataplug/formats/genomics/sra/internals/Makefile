CC = gcc
CFLAGS_DEBUG = -Wall -Wextra -O0 -ggdb
CFLAGS_RELEASE = -Wall -Wextra -O2

all: libshims.so

libshims.so: shims.c
	$(CC) -shared -fPIC -o libshims.so shims.c -ldl -lpthread $(CFLAGS_DEBUG)

release: shims.c
	$(CC) -shared -fPIC -o libshims.so shims.c -ldl -lpthread $(CFLAGS_RELEASE)

clean:
	rm -f libshims.so
