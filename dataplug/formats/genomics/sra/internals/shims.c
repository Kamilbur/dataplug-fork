#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif
#include <dlfcn.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <stdint.h>


bool INTERPOSE_LOG = 0;


char *
my_itoa(size_t val)
{
    static char buffer[20];
    size_t i = 0;
    if (val == 0) {
        return "0";
    }
    char *start = buffer + i;

    while (val > 0) {
        buffer[i++] = (val % 10) + '0';
        val /= 10;
    }
    char *end = buffer + i - 1;
    buffer[i] = '\0';

    while (start < end) {
        char tmp = *start;
        *start++ = *end;
        *end-- = tmp;
    }

    return buffer;
}

size_t
my_strlen(const char *str)
{
    size_t len = 0;
    while (str[len] != '\0') {
        len++;
    }
    return len;
}

char *
my_basename(const char *path) {
    const char *base = path;
    const char *p = path;

    while (*p != '\0') {
        if (*p == '/') {
            base = p + 1;
        }
        p++;
    }
    return (char *)base;
}


void
safe_log(const char* msg)
{
    if (! INTERPOSE_LOG) {
        return;
    }
    if (!msg) return;
    write(2, msg, my_strlen(msg));
}

#define init_real_function(name) real_##name = dlsym(RTLD_NEXT, #name)
#define enabling_variable(name) int enable_##name = 0

#define WITH_VAR_SET(var, val, code_block) \
    do { \
        (var) = val; \
        code_block; \
        (var) = 0; \
    } while(0)

#define WITH_HOOK(code_block) \
    WITH_VAR_SET(in_hook, 1, code_block)


static pthread_once_t init_once;
static __thread int in_hook;
static pthread_mutex_t owner_lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_t hook_owner;
static int hook_owner_set = 0;


enabling_variable(open);
enabling_variable(close);
enabling_variable(fstat);
enabling_variable(read);
enabling_variable(pread);
enabling_variable(mmap);
enabling_variable(munmap);


int (*real_open)(const char *pathname, int flags, ...) = NULL;
int (*real_close)(int fd) = NULL;
int (*real_fstat)(int fd, struct stat *statbuf) = NULL;
ssize_t (*real_read)(int fd, void *buf, size_t count) = NULL;
ssize_t (*real_pread)(int fd, void *buf, size_t count, off_t offset) = NULL;
void * (*real_mmap)(void *addr, size_t length, int prot, int flags, int fd, off_t offset) = NULL;
int (*real_munmap)(void *addr, size_t length) = NULL;


void
init_real(void)
{
    init_real_function(open);
    init_real_function(close);
    init_real_function(fstat);
    init_real_function(read);
    init_real_function(pread);
    init_real_function(mmap);
    init_real_function(munmap);

    char *il_str = getenv("INTERPOSE_LOG");
    if (il_str != NULL) {
       if (strncmp(il_str, "1", 1) == 0) {
           INTERPOSE_LOG = 1;
       }
       else {
           INTERPOSE_LOG = 0;
       }
    }
}

void
dp_set_hook_owner(void)
{
    pthread_mutex_lock(&owner_lock);
    hook_owner = pthread_self();
    hook_owner_set = 1;
    pthread_mutex_unlock(&owner_lock);
}

void
dp_clear_hook_owner(void)
{
    pthread_mutex_lock(&owner_lock);
    hook_owner_set = 0;
    pthread_mutex_unlock(&owner_lock);
}

static int
hooks_enabled_for_thread(void)
{
    int enabled;
    pthread_mutex_lock(&owner_lock);
    enabled = (!hook_owner_set) || pthread_equal(pthread_self(), hook_owner);
    pthread_mutex_unlock(&owner_lock);
    return enabled;
}

static int
open_needs_mode(int flags)
{
    return (flags & O_CREAT) != 0
#ifdef O_TMPFILE
        || (flags & O_TMPFILE) == O_TMPFILE
#endif
        ;
}



int special_fd[10];
size_t sfd;
struct file_info {
    const char *accession;
    const char *data;
    size_t size;
    off_t offset;
} info;

int dp_mode = 0;
size_t dp_sra_size = 0;

#define MAX_RANGES 1024
uint64_t pread_ranges[1 + MAX_RANGES * 3];
uint64_t mmap_ranges[1 + MAX_RANGES * 3];

static void *mmap_allocated[128];
static size_t n_mmap_allocated = 0;


int
open_hook(const char *pathname, int flags)
{
    if ( ! info.accession) return -1;

    if (strcmp(my_basename(pathname), info.accession) == 0) {
        safe_log("[interpose] opening special_fd\n");
        safe_log(info.accession);
        safe_log("\n");
        special_fd[sfd++] = real_open(pathname, flags);
        return special_fd[sfd - 1];
    }
    return -1;
}


int
open(const char *pathname, int flags, ...)
{
    int fd = -1;
    mode_t mode = 0;
    int has_mode = open_needs_mode(flags);
    va_list ap;
    if (has_mode) {
        va_start(ap, flags);
        mode = (mode_t)va_arg(ap, int);
        va_end(ap);
    }
    pthread_once(&init_once, init_real);
    if (!real_open) return -1;

    if (! enable_open || ! hooks_enabled_for_thread()) {
        if (has_mode) {
            fd = real_open(pathname, flags, mode);
        } else {
            fd = real_open(pathname, flags);
        }
        return fd;
    }
    if (in_hook) {
        safe_log("[interpose] [in_hook] OPEN path=");
        if (has_mode) {
            return real_open(pathname, flags, mode);
        }
        return real_open(pathname, flags);
    }

    WITH_HOOK({
        safe_log("[interpose] OPEN path=");
        safe_log(pathname);
        safe_log("\n");

        fd = open_hook(pathname, flags);
        if (fd == -1) {
            if (has_mode) {
                fd = real_open(pathname, flags, mode);
            } else {
                fd = real_open(pathname, flags);
            }
        }
    });
    return fd;
}


int
close_hook(int fd)
{
    safe_log("CLOSE_HOOK\n");
    if (sfd == 0) return -1;

    for (size_t i = 0; i < sfd; i++) {
        if (fd == special_fd[i]) {
            safe_log("[interpose] closing special_fd\n");
            special_fd[i] = special_fd[--sfd];
            return 0;
        }
    }
    return -1;
}


int
close(int fd)
{
    int ret = -1;
    pthread_once(&init_once, init_real);
    if (!real_close) return -1;

    if (! enable_open || ! hooks_enabled_for_thread()) {
        WITH_HOOK(ret = real_close(fd));
        return ret;
    }
    if (in_hook) {
        return real_close(fd);
    }

    WITH_HOOK({
        safe_log("[interpose] CLOSE");
        ret = close_hook(fd);
        if (ret == -1) {
            ret = real_close(fd);
        }
    });
    return ret;
}


int
fstat_hook(int fd, struct stat *statbuf)
{
    if (sfd == 0) return -1;
    for (size_t i = 0; i < sfd; i++) {
        if (fd == special_fd[i]) {
            memset(statbuf, 0, sizeof(*statbuf));
            statbuf->st_mode = S_IFREG | 0644;
            statbuf->st_nlink = 1;
            statbuf->st_size = (off_t)info.size;
            statbuf->st_blksize = 4096;
            statbuf->st_blocks = (statbuf->st_size + 511) / 512;
            return 0;
        }
    }
    return -1;
}

int
fstat(int fd, struct stat *statbuf)
{
    int ret = -1;
    pthread_once(&init_once, init_real);
    if (!real_fstat) return -1;

    if (! enable_fstat || ! hooks_enabled_for_thread()) {
        WITH_HOOK(ret = real_fstat(fd, statbuf));
        return ret;
    }
    if (in_hook) {
        safe_log("[interpose] [in_hook] FSTAT\n");
        return real_fstat(fd, statbuf);
    }

    WITH_HOOK({
        safe_log("[interpose] FSTAT\n");
        ret = real_fstat(fd, statbuf);
        ret = fstat_hook(fd, statbuf);
        if (ret == -1) {
            ret = real_fstat(fd, statbuf);
        }
    });
    return ret;
}


ssize_t
read(int fd, void *buf, size_t count)
{
    ssize_t ret = -1;
    pthread_once(&init_once, init_real);
    if (!real_read) return -1;

    if (! enable_read || ! hooks_enabled_for_thread()) {
        WITH_HOOK(ret = real_read(fd, buf, count));
        return ret;
    }

    if (in_hook) {
        safe_log("[interpose] [in_hook] READ\n");
        return real_read(fd, buf, count);
    }

    WITH_HOOK({
        safe_log("[interpose] READ\n");
        ret = real_read(fd, buf, count);
    });
    return ret;
}

ssize_t
pread_hook(int fd, void *buf, size_t count, off_t offset)
{
    safe_log("[interpose] pread_hook\n");
    if (sfd == 0) return -1;
    for (size_t i = 0; i < sfd; i++) {
        if (fd != special_fd[i]) continue;
        safe_log("[interpose] pread on special_fd\n");
        if (dp_mode == 0) {
            if ((size_t)offset >= info.size) return 0;
            size_t to_read = ((size_t)offset + count <= info.size) ? count : info.size - (size_t)offset;
            memcpy(buf, info.data + offset, to_read);
            uint64_t idx = pread_ranges[0];
            if (idx < MAX_RANGES) {
                pread_ranges[2 * idx + 1] = (uint64_t)offset;
                pread_ranges[2 * idx + 2] = (uint64_t)count;
                pread_ranges[0] = idx + 1;
            }
            return (ssize_t)to_read;
        } else {
            uint64_t nranges = pread_ranges[0];
            for (uint64_t j = 0; j < nranges; j++) {
                uint64_t buf_off   = pread_ranges[3 * j + 1];
                uint64_t rng_start = pread_ranges[3 * j + 2];
                uint64_t rng_end   = pread_ranges[3 * j + 3];
                if (rng_start <= (uint64_t)offset && (uint64_t)offset < rng_end) {
                    size_t avail = (size_t)(rng_end - (uint64_t)offset);
                    size_t to_read = count < avail ? count : avail;
                    memcpy(buf, info.data + buf_off + ((size_t)offset - (size_t)rng_start), to_read);
                    return (ssize_t)to_read;
                }
            }
            return -1;
        }
    }
    return -1;
}

ssize_t
pread(int fd, void *buf, size_t count, off_t offset)
{
    ssize_t ret = -1;
    pthread_once(&init_once, init_real);
    if (!real_pread) return -1;

    if (! enable_pread || ! hooks_enabled_for_thread()) {
        WITH_HOOK(ret = real_pread(fd, buf, count, offset));
        return ret;
    }
    if (in_hook) {
        safe_log("[interpose] [in_hook] PREAD\n");
        return real_pread(fd, buf, count, offset);
    }

    WITH_HOOK({
        safe_log("[interpose] PREAD\n");
        ret = pread_hook(fd, buf, count, offset);
        if (ret == -1) {
            ret = real_pread(fd, buf, count, offset);
        }
    });
    return ret;
}

void *
mmap_hook(void *addr, size_t length, int prot, int flags, int fd, off_t offset)
{
    (void) addr; (void) prot;
    if (sfd == 0) return NULL;
    if (flags & MAP_ANONYMOUS) return NULL;
    for (size_t i = 0; i < sfd; i++) {
        if (fd != special_fd[i]) continue;
        if (dp_mode == 0) {
            uint64_t idx = mmap_ranges[0];
            if (idx < MAX_RANGES) {
                mmap_ranges[2 * idx + 1] = (uint64_t)offset;
                mmap_ranges[2 * idx + 2] = (uint64_t)length;
                mmap_ranges[0] = idx + 1;
            }
            return info.data + offset;
        } else {
            uint64_t nranges = mmap_ranges[0];
            for (uint64_t j = 0; j < nranges; j++) {
                uint64_t buf_off   = mmap_ranges[3 * j + 1];
                uint64_t rng_start = mmap_ranges[3 * j + 2];
                uint64_t rng_end   = mmap_ranges[3 * j + 3];
                if (rng_start <= (uint64_t)offset && (uint64_t)offset < rng_end) {
                    void *ptr = real_mmap(NULL, length, PROT_READ | PROT_WRITE,
                                          MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
                    if (ptr == MAP_FAILED) return NULL;
                    memcpy(ptr, info.data + buf_off + ((size_t)offset - (size_t)rng_start), length);
                    if (n_mmap_allocated < 128) mmap_allocated[n_mmap_allocated++] = ptr;
                    return ptr;
                }
            }
            return NULL;
        }
    }
    return NULL;
}

void *
mmap(void *addr, size_t length, int prot, int flags, int fd, off_t offset)
{
    void *ret = NULL;
    pthread_once(&init_once, init_real);
    if (!real_mmap) return NULL;
    if (! enable_mmap || ! hooks_enabled_for_thread()) {
        WITH_HOOK(ret = real_mmap(addr, length, prot, flags, fd, offset));
        return ret;
    }
    if (in_hook) {
        safe_log("[interpose] [in_hook] MMAP\n");
        return real_mmap(addr, length, prot, flags, fd, offset);
    }
    WITH_HOOK({
        safe_log("[interpose] MMAP\n");
        ret = mmap_hook(addr, length, prot, flags, fd, offset);
        if (ret == NULL) {
            ret = real_mmap(addr, length, prot, flags, fd, offset);
        }
    });
    return ret;
}


int
munmap_hook(void *addr, size_t length)
{
    if (sfd == 0) return -1;
    if (dp_mode == 0) {
        if ((const char *)addr >= info.data && (const char *)addr < info.data + info.size) {
            safe_log("[interpose] munmap_hook no-op (mode 0)\n");
            return 0;
        }
        return -1;
    }
    for (size_t i = 0; i < n_mmap_allocated; i++) {
        if (addr == mmap_allocated[i]) {
            safe_log("[interpose] munmap_hook freeing allocated (mode 1)\n");
            real_munmap(addr, length);
            mmap_allocated[i] = mmap_allocated[--n_mmap_allocated];
            return 0;
        }
    }
    return -1;
}


int
munmap(void *addr, size_t length)
{
    int ret = -1;
    pthread_once(&init_once, init_real);
    if (!real_munmap) return -1;

    if (! enable_munmap || ! hooks_enabled_for_thread()) {
        WITH_HOOK(ret = real_munmap(addr, length));
        return ret;
    }
    if (in_hook) {
        return real_munmap(addr, length);
    }

    WITH_HOOK({
        safe_log("[interpose] munmap\n");
        ret = munmap_hook(addr, length);
        if (ret == -1) {
            ret = real_munmap(addr, length);
        }
    });
    return ret;
}
