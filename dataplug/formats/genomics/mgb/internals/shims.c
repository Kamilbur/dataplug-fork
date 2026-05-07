#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif

#include <dlfcn.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

/* glibc may define some stdio entry points as macros (fortify, gnu inline).
 * Undefine to allow our own implementations. */
#undef fread
#undef fread_unlocked
#undef fopen
#undef fopen64
#undef freopen
#undef fclose
#undef fseek
#undef fseeko
#undef fseeko64
#undef ftell
#undef ftello
#undef ftello64
#undef fileno

#define MAX_SPECIAL_FDS 128
#define MAX_SPECIAL_FPS 128
#define MAX_RANGES 65536

struct file_info {
    const char *filename;
    const char *data;
    size_t size;
};

struct fd_info {
    int fd;
    off_t offset;
};

struct fp_info {
    FILE *fp;
    off_t offset;
};

struct file_info info;
int dp_mode = 0;
uint64_t read_ranges[1 + MAX_RANGES * 3];

int enable_open = 0;
int enable_open64 = 0;
int enable_openat = 0;
int enable_openat64 = 0;
int enable_close = 0;
int enable_fstat = 0;
int enable_fstat64 = 0;
int enable_read = 0;
int enable_pread = 0;
int enable_pread64 = 0;
int enable_mmap = 0;
int enable_mmap64 = 0;
int enable_lseek = 0;
int enable_lseek64 = 0;
int enable_fopen = 0;
int enable_fopen64 = 0;
int enable_freopen = 0;
int enable_fread = 0;
int enable_fread_unlocked = 0;
int enable_fclose = 0;
int enable_fseek = 0;
int enable_fseeko = 0;
int enable_fseeko64 = 0;
int enable_ftell = 0;
int enable_ftello = 0;
int enable_ftello64 = 0;
int enable_fileno = 0;

static struct fd_info special_fds[MAX_SPECIAL_FDS];
static size_t special_fd_count = 0;
static struct fp_info special_fps[MAX_SPECIAL_FPS];
static size_t special_fp_count = 0;
static __thread int in_hook = 0;
static pthread_once_t init_once = PTHREAD_ONCE_INIT;

static int (*real_open)(const char *, int, ...) = NULL;
static int (*real_open64)(const char *, int, ...) = NULL;
static int (*real_openat)(int, const char *, int, ...) = NULL;
static int (*real_openat64)(int, const char *, int, ...) = NULL;
static int (*real_close)(int) = NULL;
static int (*real_fstat)(int, struct stat *) = NULL;
static int (*real_fstat64)(int, struct stat64 *) = NULL;
static ssize_t (*real_read)(int, void *, size_t) = NULL;
static ssize_t (*real_pread)(int, void *, size_t, off_t) = NULL;
static ssize_t (*real_pread64)(int, void *, size_t, off64_t) = NULL;
static void *(*real_mmap)(void *, size_t, int, int, int, off_t) = NULL;
static void *(*real_mmap64)(void *, size_t, int, int, int, off64_t) = NULL;
static off_t (*real_lseek)(int, off_t, int) = NULL;
static off64_t (*real_lseek64)(int, off64_t, int) = NULL;
static FILE *(*real_fopen)(const char *, const char *) = NULL;
static FILE *(*real_fopen64)(const char *, const char *) = NULL;
static FILE *(*real_freopen)(const char *, const char *, FILE *) = NULL;
static size_t (*real_fread)(void *, size_t, size_t, FILE *) = NULL;
static size_t (*real_fread_unlocked)(void *, size_t, size_t, FILE *) = NULL;
static int (*real_fclose)(FILE *) = NULL;
static int (*real_fseek)(FILE *, long, int) = NULL;
static int (*real_fseeko)(FILE *, off_t, int) = NULL;
static int (*real_fseeko64)(FILE *, off64_t, int) = NULL;
static long (*real_ftell)(FILE *) = NULL;
static off_t (*real_ftello)(FILE *) = NULL;
static off64_t (*real_ftello64)(FILE *) = NULL;
static int (*real_fileno)(FILE *) = NULL;

static void init_real(void) {
    real_open = dlsym(RTLD_NEXT, "open");
    real_open64 = dlsym(RTLD_NEXT, "open64");
    real_openat = dlsym(RTLD_NEXT, "openat");
    real_openat64 = dlsym(RTLD_NEXT, "openat64");
    real_close = dlsym(RTLD_NEXT, "close");
    real_fstat = dlsym(RTLD_NEXT, "fstat");
    real_fstat64 = dlsym(RTLD_NEXT, "fstat64");
    real_read = dlsym(RTLD_NEXT, "read");
    real_pread = dlsym(RTLD_NEXT, "pread");
    real_pread64 = dlsym(RTLD_NEXT, "pread64");
    real_mmap = dlsym(RTLD_NEXT, "mmap");
    real_mmap64 = dlsym(RTLD_NEXT, "mmap64");
    real_lseek = dlsym(RTLD_NEXT, "lseek");
    real_lseek64 = dlsym(RTLD_NEXT, "lseek64");
    real_fopen = dlsym(RTLD_NEXT, "fopen");
    real_fopen64 = dlsym(RTLD_NEXT, "fopen64");
    real_freopen = dlsym(RTLD_NEXT, "freopen");
    real_fread = dlsym(RTLD_NEXT, "fread");
    real_fread_unlocked = dlsym(RTLD_NEXT, "fread_unlocked");
    real_fclose = dlsym(RTLD_NEXT, "fclose");
    real_fseek = dlsym(RTLD_NEXT, "fseek");
    real_fseeko = dlsym(RTLD_NEXT, "fseeko");
    real_fseeko64 = dlsym(RTLD_NEXT, "fseeko64");
    real_ftell = dlsym(RTLD_NEXT, "ftell");
    real_ftello = dlsym(RTLD_NEXT, "ftello");
    real_ftello64 = dlsym(RTLD_NEXT, "ftello64");
    real_fileno = dlsym(RTLD_NEXT, "fileno");
}

static int open_needs_mode(int flags) {
    return (flags & O_CREAT) != 0
#ifdef O_TMPFILE
        || (flags & O_TMPFILE) == O_TMPFILE
#endif
        ;
}

static const char *base_name(const char *path) {
    const char *base = path;
    for (const char *p = path; *p != '\0'; ++p) {
        if (*p == '/') base = p + 1;
    }
    return base;
}

static struct fd_info *find_fd(int fd) {
    for (size_t i = 0; i < special_fd_count; ++i) {
        if (special_fds[i].fd == fd) return &special_fds[i];
    }
    return NULL;
}

static void remove_fd(int fd) {
    for (size_t i = 0; i < special_fd_count; ++i) {
        if (special_fds[i].fd == fd) {
            special_fds[i] = special_fds[--special_fd_count];
            return;
        }
    }
}

static struct fp_info *find_fp(FILE *fp) {
    for (size_t i = 0; i < special_fp_count; ++i) {
        if (special_fps[i].fp == fp) return &special_fps[i];
    }
    return NULL;
}

static void remove_fp(FILE *fp) {
    for (size_t i = 0; i < special_fp_count; ++i) {
        if (special_fps[i].fp == fp) {
            special_fps[i] = special_fps[--special_fp_count];
            return;
        }
    }
}

static void record_range(uint64_t offset, uint64_t size) {
    if (size == 0) return;
    uint64_t idx = read_ranges[0];
    if (idx >= MAX_RANGES) return;
    read_ranges[3 * idx + 1] = 0;
    read_ranges[3 * idx + 2] = offset;
    read_ranges[3 * idx + 3] = offset + size;
    read_ranges[0] = idx + 1;
}

static int open_hook(const char *pathname, int flags, mode_t mode, int has_mode, int use64) {
    if (info.filename == NULL) return -1;
    if (strcmp(base_name(pathname), info.filename) != 0) return -1;
    if (special_fd_count >= MAX_SPECIAL_FDS) return -1;

    int fd;
    if (use64) {
        fd = has_mode ? real_open64(pathname, flags, mode) : real_open64(pathname, flags);
    } else {
        fd = has_mode ? real_open(pathname, flags, mode) : real_open(pathname, flags);
    }
    if (fd >= 0) {
        special_fds[special_fd_count++] = (struct fd_info){fd, 0};
    }
    return fd;
}

int open(const char *pathname, int flags, ...) {
    pthread_once(&init_once, init_real);
    mode_t mode = 0;
    int has_mode = open_needs_mode(flags);
    if (has_mode) {
        va_list ap;
        va_start(ap, flags);
        mode = (mode_t)va_arg(ap, int);
        va_end(ap);
    }
    if (!enable_open || in_hook) {
        return has_mode ? real_open(pathname, flags, mode) : real_open(pathname, flags);
    }
    in_hook = 1;
    int fd = open_hook(pathname, flags, mode, has_mode, 0);
    if (fd == -1) {
        fd = has_mode ? real_open(pathname, flags, mode) : real_open(pathname, flags);
    }
    in_hook = 0;
    return fd;
}

int open64(const char *pathname, int flags, ...) {
    pthread_once(&init_once, init_real);
    mode_t mode = 0;
    int has_mode = open_needs_mode(flags);
    if (has_mode) {
        va_list ap;
        va_start(ap, flags);
        mode = (mode_t)va_arg(ap, int);
        va_end(ap);
    }
    if (!enable_open64 || in_hook) {
        return has_mode ? real_open64(pathname, flags, mode) : real_open64(pathname, flags);
    }
    in_hook = 1;
    int fd = open_hook(pathname, flags, mode, has_mode, 1);
    if (fd == -1) {
        fd = has_mode ? real_open64(pathname, flags, mode) : real_open64(pathname, flags);
    }
    in_hook = 0;
    return fd;
}

static int openat_hook(int dirfd, const char *pathname, int flags, mode_t mode,
                       int has_mode, int use64) {
    if (info.filename == NULL) return -1;
    if (strcmp(base_name(pathname), info.filename) != 0) return -1;
    if (special_fd_count >= MAX_SPECIAL_FDS) return -1;

    int fd;
    if (use64) {
        fd = has_mode ? real_openat64(dirfd, pathname, flags, mode)
                      : real_openat64(dirfd, pathname, flags);
    } else {
        fd = has_mode ? real_openat(dirfd, pathname, flags, mode)
                      : real_openat(dirfd, pathname, flags);
    }
    if (fd >= 0) {
        special_fds[special_fd_count++] = (struct fd_info){fd, 0};
    }
    return fd;
}

int openat(int dirfd, const char *pathname, int flags, ...) {
    pthread_once(&init_once, init_real);
    mode_t mode = 0;
    int has_mode = open_needs_mode(flags);
    if (has_mode) {
        va_list ap;
        va_start(ap, flags);
        mode = (mode_t)va_arg(ap, int);
        va_end(ap);
    }
    if (!enable_openat || in_hook) {
        return has_mode ? real_openat(dirfd, pathname, flags, mode) : real_openat(dirfd, pathname, flags);
    }
    in_hook = 1;
    int fd = openat_hook(dirfd, pathname, flags, mode, has_mode, 0);
    if (fd == -1) {
        fd = has_mode ? real_openat(dirfd, pathname, flags, mode) : real_openat(dirfd, pathname, flags);
    }
    in_hook = 0;
    return fd;
}

int openat64(int dirfd, const char *pathname, int flags, ...) {
    pthread_once(&init_once, init_real);
    mode_t mode = 0;
    int has_mode = open_needs_mode(flags);
    if (has_mode) {
        va_list ap;
        va_start(ap, flags);
        mode = (mode_t)va_arg(ap, int);
        va_end(ap);
    }
    if (!enable_openat64 || in_hook) {
        return has_mode ? real_openat64(dirfd, pathname, flags, mode) : real_openat64(dirfd, pathname, flags);
    }
    in_hook = 1;
    int fd = openat_hook(dirfd, pathname, flags, mode, has_mode, 1);
    if (fd == -1) {
        fd = has_mode ? real_openat64(dirfd, pathname, flags, mode) : real_openat64(dirfd, pathname, flags);
    }
    in_hook = 0;
    return fd;
}

int close(int fd) {
    pthread_once(&init_once, init_real);
    if (enable_close && !in_hook) remove_fd(fd);
    return real_close(fd);
}

static void fill_stat(struct stat *statbuf) {
    memset(statbuf, 0, sizeof(*statbuf));
    statbuf->st_mode = S_IFREG | 0644;
    statbuf->st_nlink = 1;
    statbuf->st_size = (off_t)info.size;
    statbuf->st_blksize = 4096;
    statbuf->st_blocks = (statbuf->st_size + 511) / 512;
}

static void fill_stat64(struct stat64 *statbuf) {
    memset(statbuf, 0, sizeof(*statbuf));
    statbuf->st_mode = S_IFREG | 0644;
    statbuf->st_nlink = 1;
    statbuf->st_size = (off64_t)info.size;
    statbuf->st_blksize = 4096;
    statbuf->st_blocks = (statbuf->st_size + 511) / 512;
}

int fstat(int fd, struct stat *statbuf) {
    pthread_once(&init_once, init_real);
    if (enable_fstat && find_fd(fd) != NULL) {
        fill_stat(statbuf);
        return 0;
    }
    return real_fstat(fd, statbuf);
}

int fstat64(int fd, struct stat64 *statbuf) {
    pthread_once(&init_once, init_real);
    if (enable_fstat64 && find_fd(fd) != NULL) {
        fill_stat64(statbuf);
        return 0;
    }
    return real_fstat64(fd, statbuf);
}

static ssize_t read_mode1(struct fd_info *fdinfo, void *buf, size_t count) {
    uint64_t offset = (uint64_t)fdinfo->offset;
    uint64_t nranges = read_ranges[0];
    for (uint64_t i = 0; i < nranges; ++i) {
        uint64_t buf_off = read_ranges[3 * i + 1];
        uint64_t start = read_ranges[3 * i + 2];
        uint64_t end = read_ranges[3 * i + 3];
        if (start <= offset && offset < end) {
            size_t avail = (size_t)(end - offset);
            size_t to_read = count < avail ? count : avail;
            if (info.data == NULL) {
                ssize_t nread = real_read(fdinfo->fd, buf, to_read);
                if (nread > 0) fdinfo->offset += (off_t)nread;
                return nread;
            }
            memcpy(buf, info.data + buf_off + (offset - start), to_read);
            fdinfo->offset += (off_t)to_read;
            return (ssize_t)to_read;
        }
    }
    return -1;
}

ssize_t read(int fd, void *buf, size_t count) {
    pthread_once(&init_once, init_real);
    struct fd_info *fdinfo = find_fd(fd);
    if (!enable_read || in_hook || fdinfo == NULL) {
        return real_read(fd, buf, count);
    }
    if (dp_mode == 0) {
        if ((size_t)fdinfo->offset >= info.size) return 0;
        size_t remaining = info.size - (size_t)fdinfo->offset;
        size_t to_read = count < remaining ? count : remaining;
        if (info.data == NULL) {
            uint64_t offset = (uint64_t)fdinfo->offset;
            ssize_t nread = real_read(fd, buf, to_read);
            if (nread > 0) {
                record_range(offset, (uint64_t)nread);
                fdinfo->offset += (off_t)nread;
            }
            return nread;
        }
        memcpy(buf, info.data + fdinfo->offset, to_read);
        record_range((uint64_t)fdinfo->offset, (uint64_t)to_read);
        fdinfo->offset += (off_t)to_read;
        return (ssize_t)to_read;
    }
    return read_mode1(fdinfo, buf, count);
}

ssize_t pread(int fd, void *buf, size_t count, off_t offset) {
    pthread_once(&init_once, init_real);
    struct fd_info *fdinfo = find_fd(fd);
    if (!enable_pread || in_hook || fdinfo == NULL) {
        return real_pread(fd, buf, count, offset);
    }
    if (dp_mode == 0) {
        ssize_t nread = real_pread(fd, buf, count, offset);
        if (nread > 0) record_range((uint64_t)offset, (uint64_t)nread);
        return nread;
    }
    return real_pread(fd, buf, count, offset);
}

ssize_t pread64(int fd, void *buf, size_t count, off64_t offset) {
    pthread_once(&init_once, init_real);
    struct fd_info *fdinfo = find_fd(fd);
    if (!enable_pread64 || in_hook || fdinfo == NULL) {
        return real_pread64(fd, buf, count, offset);
    }
    if (dp_mode == 0) {
        ssize_t nread = real_pread64(fd, buf, count, offset);
        if (nread > 0) record_range((uint64_t)offset, (uint64_t)nread);
        return nread;
    }
    return real_pread64(fd, buf, count, offset);
}

void *mmap(void *addr, size_t length, int prot, int flags, int fd, off_t offset) {
    pthread_once(&init_once, init_real);
    struct fd_info *fdinfo = find_fd(fd);
    if (!enable_mmap || in_hook || fdinfo == NULL) {
        return real_mmap(addr, length, prot, flags, fd, offset);
    }
    if (dp_mode == 0) {
        void *ret = real_mmap(addr, length, prot, flags, fd, offset);
        if (ret != MAP_FAILED) record_range((uint64_t)offset, (uint64_t)length);
        return ret;
    }
    return real_mmap(addr, length, prot, flags, fd, offset);
}

void *mmap64(void *addr, size_t length, int prot, int flags, int fd, off64_t offset) {
    pthread_once(&init_once, init_real);
    struct fd_info *fdinfo = find_fd(fd);
    if (!enable_mmap64 || in_hook || fdinfo == NULL) {
        return real_mmap64(addr, length, prot, flags, fd, offset);
    }
    if (dp_mode == 0) {
        void *ret = real_mmap64(addr, length, prot, flags, fd, offset);
        if (ret != MAP_FAILED) record_range((uint64_t)offset, (uint64_t)length);
        return ret;
    }
    return real_mmap64(addr, length, prot, flags, fd, offset);
}

static off64_t seek_impl(struct fd_info *fdinfo, off64_t offset, int whence) {
    off64_t base;
    if (whence == SEEK_SET) base = 0;
    else if (whence == SEEK_CUR) base = fdinfo->offset;
    else if (whence == SEEK_END) base = (off64_t)info.size;
    else return -1;
    off64_t next = base + offset;
    if (next < 0) return -1;
    fdinfo->offset = (off_t)next;
    return next;
}

off_t lseek(int fd, off_t offset, int whence) {
    pthread_once(&init_once, init_real);
    struct fd_info *fdinfo = find_fd(fd);
    if (!enable_lseek || in_hook || fdinfo == NULL) {
        return real_lseek(fd, offset, whence);
    }
    if (info.data == NULL) {
        off_t next = real_lseek(fd, offset, whence);
        if (next >= 0) fdinfo->offset = next;
        return next;
    }
    return (off_t)seek_impl(fdinfo, offset, whence);
}

off64_t lseek64(int fd, off64_t offset, int whence) {
    pthread_once(&init_once, init_real);
    struct fd_info *fdinfo = find_fd(fd);
    if (!enable_lseek64 || in_hook || fdinfo == NULL) {
        return real_lseek64(fd, offset, whence);
    }
    if (info.data == NULL) {
        off64_t next = real_lseek64(fd, offset, whence);
        if (next >= 0) fdinfo->offset = (off_t)next;
        return next;
    }
    return seek_impl(fdinfo, offset, whence);
}

static FILE *fopen_hook(const char *pathname, const char *mode, int use64) {
    if (info.filename == NULL) return NULL;
    if (strcmp(base_name(pathname), info.filename) != 0) return NULL;
    if (special_fp_count >= MAX_SPECIAL_FPS) return NULL;
    FILE *fp = use64 ? real_fopen64(pathname, mode) : real_fopen(pathname, mode);
    if (fp != NULL) {
        special_fps[special_fp_count++] = (struct fp_info){fp, 0};
    }
    return fp;
}

FILE *fopen(const char *pathname, const char *mode) {
    pthread_once(&init_once, init_real);
    if (!enable_fopen || in_hook) return real_fopen(pathname, mode);
    in_hook = 1;
    FILE *fp = fopen_hook(pathname, mode, 0);
    if (fp == NULL) fp = real_fopen(pathname, mode);
    in_hook = 0;
    return fp;
}

FILE *fopen64(const char *pathname, const char *mode) {
    pthread_once(&init_once, init_real);
    if (!enable_fopen64 || in_hook) return real_fopen64(pathname, mode);
    in_hook = 1;
    FILE *fp = fopen_hook(pathname, mode, 1);
    if (fp == NULL) fp = real_fopen64(pathname, mode);
    in_hook = 0;
    return fp;
}

FILE *freopen(const char *pathname, const char *mode, FILE *stream) {
    pthread_once(&init_once, init_real);
    if (!enable_freopen || in_hook) return real_freopen(pathname, mode, stream);
    in_hook = 1;
    remove_fp(stream);
    FILE *fp = real_freopen(pathname, mode, stream);
    if (fp != NULL && info.filename != NULL && pathname != NULL &&
        strcmp(base_name(pathname), info.filename) == 0 &&
        special_fp_count < MAX_SPECIAL_FPS) {
        special_fps[special_fp_count++] = (struct fp_info){fp, 0};
    }
    in_hook = 0;
    return fp;
}

static size_t fread_record(void *ptr, size_t size, size_t nmemb, FILE *stream,
                           size_t (*fn)(void *, size_t, size_t, FILE *)) {
    struct fp_info *fpi = find_fp(stream);
    if (fpi == NULL) return fn(ptr, size, nmemb, stream);
    if (dp_mode == 0) {
        size_t result = fn(ptr, size, nmemb, stream);
        size_t bytes = result * size;
        if (bytes > 0) {
            record_range((uint64_t)fpi->offset, (uint64_t)bytes);
            fpi->offset += (off_t)bytes;
        }
        return result;
    }
    return fn(ptr, size, nmemb, stream);
}

size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    pthread_once(&init_once, init_real);
    if (!enable_fread || in_hook) return real_fread(ptr, size, nmemb, stream);
    in_hook = 1;
    size_t result = fread_record(ptr, size, nmemb, stream, real_fread);
    in_hook = 0;
    return result;
}

size_t fread_unlocked(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    pthread_once(&init_once, init_real);
    if (!enable_fread_unlocked || in_hook || real_fread_unlocked == NULL) {
        return real_fread_unlocked ? real_fread_unlocked(ptr, size, nmemb, stream)
                                   : real_fread(ptr, size, nmemb, stream);
    }
    in_hook = 1;
    size_t result = fread_record(ptr, size, nmemb, stream, real_fread_unlocked);
    in_hook = 0;
    return result;
}

int fclose(FILE *stream) {
    pthread_once(&init_once, init_real);
    if (enable_fclose && !in_hook) remove_fp(stream);
    return real_fclose(stream);
}

static off64_t fseek_impl(struct fp_info *fpi, off64_t offset, int whence) {
    off64_t base;
    if (whence == SEEK_SET) base = 0;
    else if (whence == SEEK_CUR) base = fpi->offset;
    else if (whence == SEEK_END) base = (off64_t)info.size;
    else return -1;
    off64_t next = base + offset;
    if (next < 0) return -1;
    fpi->offset = (off_t)next;
    return next;
}

int fseek(FILE *stream, long offset, int whence) {
    pthread_once(&init_once, init_real);
    struct fp_info *fpi = find_fp(stream);
    if (!enable_fseek || in_hook || fpi == NULL) {
        return real_fseek(stream, offset, whence);
    }
    in_hook = 1;
    int ret = real_fseek(stream, offset, whence);
    if (ret == 0) fseek_impl(fpi, (off64_t)offset, whence);
    in_hook = 0;
    return ret;
}

int fseeko(FILE *stream, off_t offset, int whence) {
    pthread_once(&init_once, init_real);
    struct fp_info *fpi = find_fp(stream);
    if (!enable_fseeko || in_hook || fpi == NULL) {
        return real_fseeko(stream, offset, whence);
    }
    in_hook = 1;
    int ret = real_fseeko(stream, offset, whence);
    if (ret == 0) fseek_impl(fpi, (off64_t)offset, whence);
    in_hook = 0;
    return ret;
}

int fseeko64(FILE *stream, off64_t offset, int whence) {
    pthread_once(&init_once, init_real);
    struct fp_info *fpi = find_fp(stream);
    if (!enable_fseeko64 || in_hook || fpi == NULL || real_fseeko64 == NULL) {
        return real_fseeko64 ? real_fseeko64(stream, offset, whence)
                             : real_fseeko(stream, (off_t)offset, whence);
    }
    in_hook = 1;
    int ret = real_fseeko64(stream, offset, whence);
    if (ret == 0) fseek_impl(fpi, offset, whence);
    in_hook = 0;
    return ret;
}

long ftell(FILE *stream) {
    pthread_once(&init_once, init_real);
    struct fp_info *fpi = find_fp(stream);
    if (!enable_ftell || in_hook || fpi == NULL) {
        return real_ftell(stream);
    }
    return (long)fpi->offset;
}

off_t ftello(FILE *stream) {
    pthread_once(&init_once, init_real);
    struct fp_info *fpi = find_fp(stream);
    if (!enable_ftello || in_hook || fpi == NULL || real_ftello == NULL) {
        return real_ftello ? real_ftello(stream) : (off_t)real_ftell(stream);
    }
    return fpi->offset;
}

off64_t ftello64(FILE *stream) {
    pthread_once(&init_once, init_real);
    struct fp_info *fpi = find_fp(stream);
    if (!enable_ftello64 || in_hook || fpi == NULL || real_ftello64 == NULL) {
        return real_ftello64 ? real_ftello64(stream) : (off64_t)real_ftell(stream);
    }
    return (off64_t)fpi->offset;
}

int fileno(FILE *stream) {
    pthread_once(&init_once, init_real);
    int fd = real_fileno(stream);
    if (!enable_fileno || in_hook) return fd;
    struct fp_info *fpi = find_fp(stream);
    if (fpi != NULL && fd >= 0 && find_fd(fd) == NULL &&
        special_fd_count < MAX_SPECIAL_FDS) {
        special_fds[special_fd_count++] = (struct fd_info){fd, fpi->offset};
    }
    return fd;
}
