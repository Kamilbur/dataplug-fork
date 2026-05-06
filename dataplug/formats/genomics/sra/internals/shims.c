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
#include <limits.h>


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
    size_t data_size;
    size_t size;
    off_t offset;
} info;

static char *owned_accession = NULL;
static char *owned_data = NULL;

int dp_mode = 0;
size_t dp_sra_size = 0;

#define MAX_RANGES 1024
uint64_t pread_ranges[1 + MAX_RANGES * 3];
uint64_t mmap_ranges[1 + MAX_RANGES * 3];
static int pread_ranges_overflow = 0;
static int mmap_ranges_overflow = 0;

static uint64_t *walk_preads = NULL;
static size_t walk_pread_count = 0;
static size_t walk_pread_cap = 0;
static uint64_t *walk_mmaps = NULL;
static size_t walk_mmap_count = 0;
static size_t walk_mmap_cap = 0;

typedef int32_t (*VCursorCellDataDirect_f)(
    const void *self,
    int64_t row_id,
    uint32_t col_idx,
    uint32_t *elem_bits,
    const void **base,
    uint32_t *boff,
    uint32_t *row_len
);
static VCursorCellDataDirect_f real_VCursorCellDataDirect = NULL;

static void *mmap_allocated[128];
static size_t n_mmap_allocated = 0;


void
dp_clear_walk_ranges(void)
{
    free(walk_preads);
    free(walk_mmaps);
    walk_preads = NULL;
    walk_mmaps = NULL;
    walk_pread_count = 0;
    walk_mmap_count = 0;
    walk_pread_cap = 0;
    walk_mmap_cap = 0;
}


size_t
dp_walk_pread_count(void)
{
    return walk_pread_count;
}


uint64_t *
dp_walk_pread_data(void)
{
    return walk_preads;
}


size_t
dp_walk_mmap_count(void)
{
    return walk_mmap_count;
}


uint64_t *
dp_walk_mmap_data(void)
{
    return walk_mmaps;
}


static int
append_pread(uint64_t row_key, uint64_t start, uint64_t end)
{
    if (walk_pread_count == walk_pread_cap) {
        size_t next_cap = walk_pread_cap ? walk_pread_cap * 2 : 4096;
        if (next_cap > SIZE_MAX / (3 * sizeof(uint64_t))) return -1;
        uint64_t *next = realloc(walk_preads, next_cap * 3 * sizeof(uint64_t));
        if (!next) return -1;
        walk_preads = next;
        walk_pread_cap = next_cap;
    }
    size_t idx = walk_pread_count * 3;
    walk_preads[idx] = row_key;
    walk_preads[idx + 1] = start;
    walk_preads[idx + 2] = end;
    walk_pread_count++;
    return 0;
}


static int
append_mmap(uint64_t start, uint64_t end)
{
    if (walk_mmap_count == walk_mmap_cap) {
        size_t next_cap = walk_mmap_cap ? walk_mmap_cap * 2 : 128;
        if (next_cap > SIZE_MAX / (2 * sizeof(uint64_t))) return -1;
        uint64_t *next = realloc(walk_mmaps, next_cap * 2 * sizeof(uint64_t));
        if (!next) return -1;
        walk_mmaps = next;
        walk_mmap_cap = next_cap;
    }
    size_t idx = walk_mmap_count * 2;
    walk_mmaps[idx] = start;
    walk_mmaps[idx + 1] = end;
    walk_mmap_count++;
    return 0;
}


static int
drain_walk_ranges(uint64_t row_key, uint64_t total_size)
{
    if (pread_ranges_overflow || mmap_ranges_overflow) {
        return -2;
    }

    uint64_t mmap_count = mmap_ranges[0];
    for (uint64_t i = 0; i < mmap_count; i++) {
        uint64_t pos = mmap_ranges[2 * i + 1];
        uint64_t size = mmap_ranges[2 * i + 2];
        if (pos + size + 4 < total_size) size += 4;
        if (append_mmap(pos, pos + size) != 0) return -1;
    }
    mmap_ranges[0] = 0;

    uint64_t pread_count = pread_ranges[0];
    for (uint64_t i = 0; i < pread_count; i++) {
        uint64_t pos = pread_ranges[2 * i + 1];
        uint64_t size = pread_ranges[2 * i + 2];
        if (pos + size + 4 < total_size) size += 4;
        if (append_pread(row_key, pos, pos + size) != 0) return -1;
    }
    pread_ranges[0] = 0;
    return 0;
}


int
dp_walk_columns(
    void *cur,
    const int *col_idxs,
    size_t n_cols,
    int64_t total_lines,
    int64_t step,
    uint64_t total_size
)
{
    pthread_once(&init_once, init_real);
    if (!real_VCursorCellDataDirect) {
        real_VCursorCellDataDirect = (VCursorCellDataDirect_f)dlsym(
            RTLD_DEFAULT, "VCursorCellDataDirect"
        );
    }
    if (!real_VCursorCellDataDirect || !cur || !col_idxs || step <= 0) {
        return -10;
    }

    dp_clear_walk_ranges();
    pread_ranges_overflow = 0;
    mmap_ranges_overflow = 0;

    int rc = drain_walk_ranges(UINT64_MAX, total_size);
    if (rc != 0) return rc;

    uint32_t elem_bits = 0;
    const void *data = NULL;
    uint32_t row_len = 0;
    uint64_t group_idx = 0;
    for (int64_t row_idx = 1; row_idx <= total_lines; row_idx += step) {
        for (size_t ci = 0; ci < n_cols; ci++) {
            rc = real_VCursorCellDataDirect(
                cur,
                row_idx,
                (uint32_t)col_idxs[ci],
                &elem_bits,
                &data,
                NULL,
                &row_len
            );
            if (rc != 0) {
                return rc;
            }
        }
        rc = drain_walk_ranges(group_idx * (uint64_t)step, total_size);
        if (rc != 0) return rc;
        group_idx++;
    }
    return 0;
}


void
dp_clear_info(void)
{
    dp_clear_walk_ranges();
    free(owned_accession);
    free(owned_data);
    owned_accession = NULL;
    owned_data = NULL;
    info.accession = NULL;
    info.data = NULL;
    info.data_size = 0;
    info.size = 0;
    info.offset = 0;
    sfd = 0;
}


void
dp_set_info_buffer(const char *accession, const char *data, size_t data_size, size_t size, off_t offset)
{
    dp_clear_info();
    if (!accession || !data) {
        return;
    }

    size_t acc_len = my_strlen(accession);
    owned_accession = malloc(acc_len + 1);
    owned_data = malloc(data_size);
    if (!owned_accession || !owned_data) {
        dp_clear_info();
        return;
    }

    memcpy(owned_accession, accession, acc_len + 1);
    memcpy(owned_data, data, data_size);

    info.accession = owned_accession;
    info.data = owned_data;
    info.data_size = data_size;
    info.size = size;
    info.offset = offset;
}


/*
 * Zero-copy variant: caller retains ownership of `data` and guarantees its
 * lifetime exceeds the next dp_clear_info() call. Avoids the 1x file-size
 * malloc+memcpy of dp_set_info_buffer, which is critical when the SRA file
 * is multi-GB and process RSS would otherwise double.
 */
void
dp_set_info_zerocopy(const char *accession, const char *data, size_t data_size, size_t size, off_t offset)
{
    dp_clear_info();
    if (!accession || !data) {
        return;
    }

    size_t acc_len = my_strlen(accession);
    owned_accession = malloc(acc_len + 1);
    if (!owned_accession) {
        return;
    }
    memcpy(owned_accession, accession, acc_len + 1);

    info.accession = owned_accession;
    info.data = data;
    info.data_size = data_size;
    info.size = size;
    info.offset = offset;
}



void
dp_set_info(const char *accession, const char *data, size_t size, off_t offset)
{
    dp_set_info_buffer(accession, data, size, size, offset);
}


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
            } else {
                pread_ranges_overflow = 1;
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
            } else {
                mmap_ranges_overflow = 1;
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
