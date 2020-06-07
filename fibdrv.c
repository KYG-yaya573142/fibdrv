#include <linux/cdev.h>
#include <linux/device.h>
#include <linux/fs.h>
#include <linux/init.h>
#include <linux/kdev_t.h>
#include <linux/kernel.h>
#include <linux/module.h>
#include <linux/mutex.h>
#include <linux/slab.h>
#include <linux/uaccess.h>  // Required for the copy_to_user()

#include <linux/kobject.h>  // sysfs
#include <linux/sysfs.h>

#include "bn_kernel.h"

MODULE_LICENSE("Dual MIT/GPL");
MODULE_AUTHOR("National Cheng Kung University, Taiwan");
MODULE_DESCRIPTION("Fibonacci engine driver");
MODULE_VERSION("0.1");

#define DEV_FIBONACCI_NAME "fibonacci"

/* MAX_LENGTH is set to 92 because
 * ssize_t can't fit the number > 92
 */
#define MAX_LENGTH 300

/* making our own kobject */
struct fib_obj {
    struct kobject kobj;
    int n;
};
#define to_fib_obj(x) container_of(x, struct fib_obj, kobj)

/* a custom attribute that works just for a struct fib_obj */
struct fib_attribute {
    struct attribute attr;
    ssize_t (*show)(struct fib_obj *fib,
                    struct fib_attribute *f_attr,
                    char *buf);
    ssize_t (*store)(struct fib_obj *fib,
                     struct fib_attribute *f_attr,
                     const char *buf,
                     size_t count);
};
#define to_fib_attr(x) container_of(x, struct fib_attribute, attr)

/*
 * The default show/store functions to be passed to sysfs.
 * These will be called by sysfs whenever the user read/write on the sysfs file
 * associated with the kobjects we have registered.
 * These functions will transpose a "default" kobject back to our custom
 * struct fib_obj then call the "real" show/store function registered in
 * fib_attribute.
 */
static ssize_t fib_attr_show(struct kobject *kobj,
                             struct attribute *attr,
                             char *buf)
{
    struct fib_obj *fib;
    struct fib_attribute *f_attr;

    fib = to_fib_obj(kobj);
    f_attr = to_fib_attr(attr);
    if (!f_attr->show)
        return -EIO;

    return f_attr->show(fib, f_attr, buf);
}

static ssize_t fib_attr_store(struct kobject *kobj,
                              struct attribute *attr,
                              const char *buf,
                              size_t len)
{
    struct fib_obj *fib;
    struct fib_attribute *f_attr;

    fib = to_fib_obj(kobj);
    f_attr = to_fib_attr(attr);
    if (!f_attr->store)
        return -EIO;

    return f_attr->store(fib, f_attr, buf, len);
}

/* sysfs_ops for fib_ktype */
static struct sysfs_ops fib_sysfs_ops = {
    .show = fib_attr_show,
    .store = fib_attr_store,
};

/* release function for fib_ktype */
static void fib_obj_release(struct kobject *kobj)
{
    struct fib_obj *fib;

    fib = to_fib_obj(kobj);
    kfree(fib);
}

/*
 * The "real" show/store functions to be registered in fib_attribute.
 * user read to sysfs (any fib_attribute) > fib_attr_show > fib_show
 */
static ssize_t fib_show(struct fib_obj *fib_obj,
                        struct fib_attribute *f_attr,
                        char *buf)
{
    int retval;
    bn *fib = bn_alloc(1);
    bn_fib_fdoubling(fib, fib_obj->n);
    char *p = bn_to_string(fib);
    retval = scnprintf(buf, PAGE_SIZE, "%s\n", p);
    bn_free(fib);
    kfree(p);
    return retval;
}

static ssize_t fib_store(struct fib_obj *fib,
                         struct fib_attribute *f_attr,
                         const char *buf,
                         size_t count)
{
    int ret;

    ret = kstrtoint(buf, 10, &(fib->n));
    if (ret < 0)
        return ret;

    return count;
}

static struct fib_attribute nth = __ATTR(nth, 0664, fib_show, fib_store);

/* default attribute for fib_ktype */
static struct attribute *fib_default_attrs[] = {
    &nth.attr, NULL, /* need to NULL terminate the list of attributes */
};

/* Define our own ktype */
static struct kobj_type fib_ktype = {
    .sysfs_ops = &fib_sysfs_ops,
    .release = fib_obj_release,
    .default_attrs = fib_default_attrs,
};

static struct kset *linux2020_kset;
static struct fib_obj *fib_obj;

/* since we are using custume kobject, a dedicated initial function is needed */
static struct fib_obj *create_fib_obj(void)
{
    struct fib_obj *fib;
    int retval;

    /* allocate the memory for the whole object */
    fib = kzalloc(sizeof(*fib), GFP_KERNEL);
    if (!fib)
        return NULL;

    /* the kobject will be placed under the kset, no need to set a parent */
    fib->kobj.kset = linux2020_kset;
    retval = kobject_init_and_add(&fib->kobj, &fib_ktype, NULL, "fibdrv");
    if (retval) {
        kobject_put(&fib->kobj);
        return NULL;
    }
    kobject_uevent(&fib->kobj, KOBJ_ADD);
    return fib;
}

static dev_t fib_dev = 0;
static struct cdev *fib_cdev;
static struct class *fib_class;
static DEFINE_MUTEX(fib_mutex);

/* recursion with cache */
static long long fib_sequence(long long k)
{
    /* FIXME: use clz/ctz and fast algorithms to speed up */
    long long f[k + 2];

    f[0] = 0;
    f[1] = 1;

    for (int i = 2; i <= k; i++) {
        f[i] = f[i - 1] + f[i - 2];
    }

    return f[k];
}

/* Calculate Fibonacci numbers by Fast Doubling */
static long long fib_sequence_fdouble(long long n)
{
    if (n < 2) { /* F(0) = 0, F(1) = 1 */
        return n;
    }
    long long f[2];
    unsigned int ndigit = 32 - __builtin_clz(n); /* number of digit in n */
    f[0] = 0;                                    /* F(k) */
    f[1] = 1;                                    /* F(k+1) */

    for (unsigned int i = 1U << (ndigit - 1); i;
         i >>= 1) { /* walk through the digit of n */
        long long k1 =
            f[0] * (f[1] * 2 - f[0]); /* F(2k) = F(k) * [ 2 * F(k+1) – F(k) ] */
        long long k2 =
            f[0] * f[0] + f[1] * f[1]; /* F(2k+1) = F(k)^2 + F(k+1)^2 */
        if (n & i) {                   /* current binary digit == 1 */
            f[0] = k2;                 /* F(n) = F(2k+1) */
            f[1] = k1 + k2; /* F(n+1) = F(2k+2) =  F(2k) +  F(2k+1) */
        } else {
            f[0] = k1; /* F(n) = F(2k) */
            f[1] = k2; /* F(n+1) = F(2k+1) */
        }
    }
    return f[0];
}

static int fib_open(struct inode *inode, struct file *file)
{
    if (!mutex_trylock(&fib_mutex)) {
        printk(KERN_ALERT "fibdrv is in use");
        return -EBUSY;
    }
    return 0;
}

static int fib_release(struct inode *inode, struct file *file)
{
    mutex_unlock(&fib_mutex);
    return 0;
}

/* calculate the fibonacci number at given offset */
static ssize_t fib_read(struct file *file,
                        char *buf,
                        size_t size,
                        loff_t *offset)
{
    bn *fib = bn_alloc(1);
    bn_fib_fdoubling(fib, *offset);
    // bn_fib(fib, *offset);
    char *p = bn_to_string(fib);
    size_t len = strlen(p) + 1;
    size_t left = copy_to_user(buf, p, len);
    // printk(KERN_DEBUG "fib(%d): %s\n", (int) *offset, p);
    bn_free(fib);
    kfree(p);
    if (left)
        return -EFAULT;
    return 0;
}

/* write operation actually returns the time spent on
 * calculating the fibonacci number at given offset
 */
static ssize_t fib_write(struct file *file,
                         const char *buf,
                         size_t size,
                         loff_t *offset)
{
    ktime_t kt;
    switch (size) {
    case 0:
        kt = ktime_get();
        fib_sequence(*offset);
        kt = ktime_sub(ktime_get(), kt);
        break;
    case 1:
        kt = ktime_get();
        fib_sequence_fdouble(*offset);
        kt = ktime_sub(ktime_get(), kt);
        break;
    default:
        return 0;
    }
    return (ssize_t) ktime_to_ns(kt);
}

static loff_t fib_device_lseek(struct file *file, loff_t offset, int orig)
{
    loff_t new_pos = 0;
    switch (orig) {
    case 0: /* SEEK_SET: */
        new_pos = offset;
        break;
    case 1: /* SEEK_CUR: */
        new_pos = file->f_pos + offset;
        break;
    case 2: /* SEEK_END: */
        new_pos = MAX_LENGTH - offset;
        break;
    }

    if (new_pos > MAX_LENGTH)
        new_pos = MAX_LENGTH;  // max case
    if (new_pos < 0)
        new_pos = 0;        // min case
    file->f_pos = new_pos;  // This is what we'll use now
    return new_pos;
}

const struct file_operations fib_fops = {
    .owner = THIS_MODULE,
    .read = fib_read,
    .write = fib_write,
    .open = fib_open,
    .release = fib_release,
    .llseek = fib_device_lseek,
};

static int __init init_fib_dev(void)
{
    int rc = 0;

    mutex_init(&fib_mutex);

    // stuff of sysfs registeration
    linux2020_kset = kset_create_and_add("linux2020", NULL, kernel_kobj);
    if (!linux2020_kset)
        return -ENOMEM;
    fib_obj = create_fib_obj();
    if (!fib_obj)
        goto failed_sysfs;

    // Let's register the device
    // This will dynamically allocate the major number
    rc = alloc_chrdev_region(&fib_dev, 0, 1, DEV_FIBONACCI_NAME);

    if (rc < 0) {
        printk(KERN_ALERT
               "Failed to register the fibonacci char device. rc = %i",
               rc);
        return rc;
    }

    fib_cdev = cdev_alloc();
    if (fib_cdev == NULL) {
        printk(KERN_ALERT "Failed to alloc cdev");
        rc = -1;
        goto failed_cdev;
    }
    cdev_init(fib_cdev, &fib_fops);
    rc = cdev_add(fib_cdev, fib_dev, 1);

    if (rc < 0) {
        printk(KERN_ALERT "Failed to add cdev");
        rc = -2;
        goto failed_cdev;
    }

    fib_class = class_create(THIS_MODULE, DEV_FIBONACCI_NAME);

    if (!fib_class) {
        printk(KERN_ALERT "Failed to create device class");
        rc = -3;
        goto failed_class_create;
    }

    if (!device_create(fib_class, NULL, fib_dev, NULL, DEV_FIBONACCI_NAME)) {
        printk(KERN_ALERT "Failed to create device");
        rc = -4;
        goto failed_device_create;
    }
    return rc;
failed_device_create:
    class_destroy(fib_class);
failed_class_create:
    cdev_del(fib_cdev);
failed_cdev:
    unregister_chrdev_region(fib_dev, 1);
failed_sysfs:
    kset_unregister(linux2020_kset);
    return rc;
}

static void __exit exit_fib_dev(void)
{
    mutex_destroy(&fib_mutex);
    device_destroy(fib_class, fib_dev);
    class_destroy(fib_class);
    cdev_del(fib_cdev);
    unregister_chrdev_region(fib_dev, 1);
    kobject_put(&(fib_obj->kobj));
    kset_unregister(linux2020_kset);
}

module_init(init_fib_dev);
module_exit(exit_fib_dev);
