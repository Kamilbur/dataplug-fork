CC = gcc
CFLAGS_DEBUG = -Wall -Wextra -O0 -ggdb
CFLAGS_RELEASE = -Wall -Wextra -O2

all: libmgbshims.so

libmgbshims.so: shims.c
	$(CC) -shared -fPIC -o libmgbshims.so shims.c -ldl -lpthread $(CFLAGS_DEBUG)

release: shims.c
	$(CC) -shared -fPIC -o libmgbshims.so shims.c -ldl -lpthread $(CFLAGS_RELEASE)

clean:
	rm -f libmgbshims.so
