#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#define FIB_DEV "/dev/fibonacci"
#define sample_size 1000

int main()
{
    FILE *fp = fopen("./plot_input_statistic", "w");
    char write_buf[] = "testing writing";
    int offset = 100;

    int fd = open(FIB_DEV, O_RDWR);
    if (fd < 0) {
        perror("Failed to open character device");
        exit(1);
    }

    /* for each F(i), measure sample_size times of data and
     * remove outliers based on the 95% confidence level
     */
    for (int i = 0; i <= offset; i++) {
        lseek(fd, i, SEEK_SET);
        double t1[sample_size] = {0}, t2[sample_size] = {0};
        double mean1 = 0.0, sd1 = 0.0, result1 = 0.0;
        double mean2 = 0.0, sd2 = 0.0, result2 = 0.0;
        int count1 = 0, count2 = 0;

        for (int n = 0; n < sample_size; n++) { /* sampling */
            /* get the runtime in kernel space here */
            t1[n] = (double) write(fd, write_buf, 0); /* recursion */
            t2[n] = (double) write(fd, write_buf, 1); /* fast doubling */
            mean1 += t1[n];                           /* sum */
            mean2 += t2[n];
        }
        mean1 /= sample_size; /* mean */
        mean2 /= sample_size;

        for (int n = 0; n < sample_size; n++) {
            sd1 += (t1[n] - mean1) * (t1[n] - mean1);
            sd2 += (t2[n] - mean2) * (t2[n] - mean2);
        }
        sd1 = sqrt(sd1 / (sample_size - 1)); /* standard deviation */
        sd2 = sqrt(sd2 / (sample_size - 1));

        for (int n = 0; n < sample_size; n++) { /* remove outliers */
            if (t1[n] <= (mean1 + 2 * sd1) && t1[n] >= (mean1 - 2 * sd1)) {
                result1 += t1[n];
                count1++;
            }
            if (t2[n] <= (mean2 + 2 * sd2) && t2[n] >= (mean2 - 2 * sd2)) {
                result2 += t2[n];
                count2++;
            }
        }
        result1 /= count1;
        result2 /= count2;

        fprintf(fp, "%d %.5lf %.5lf samples: %d %d\n", i, result1, result2,
                count1, count2);
    }
    close(fd);
    fclose(fp);
    return 0;
}