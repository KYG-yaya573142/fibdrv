#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#define FIB_DEV "/dev/fibonacci"

int main()
{
    char buf[100];
    // char write_buf[] = "testing writing";
    int offset = 300;

    int fd = open(FIB_DEV, O_RDWR);
    if (fd < 0) {
        perror("Failed to open character device");
        exit(1);
    }

    for (int i = 0; i <= offset; i++) {
        lseek(fd, i, SEEK_SET);
        long long sz = read(fd, buf, 100);
        if (sz)
            printf("returned message was truncated!\n");
        printf("fib(%d): %s\n", i, buf);
    }

    close(fd);
    return 0;
}
