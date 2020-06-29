#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bn.h"

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#ifndef SWAP
#define SWAP(x, y)           \
    do {                     \
        typeof(x) __tmp = x; \
        x = y;               \
        y = __tmp;           \
    } while (0)
#endif

#ifndef DIV_ROUNDUP
#define DIV_ROUNDUP(x, len) (((x) + (len) -1) / (len))
#endif

#ifndef unlikely
#define unlikely(x) (__builtin_expect((x), 0))
#endif

/* count leading zeros of src->number */
static int bn_clz(const bn *src)
{
    int cnt = 0;
    for (int i = src->size - 1; i >= 0; i--) {
        if (src->number[i]) {
            // prevent undefined behavior when src->number[i] = 0
            return cnt + builtin_clz(src->number[i]);
        }
        cnt += DATA_BITS;
    }
    return cnt;
}

/* count the digits of most significant bit */
static int bn_msb(const bn *src)
{
    return src->size * DATA_BITS - bn_clz(src);
}

/*
 * output bn to decimal string
 * Note: the returned string should be freed with the free()
 */
char *bn_to_string(const bn *src)
{
    // log10(x) = log2(x) / log2(10) ~= log2(x) / 3.322
    size_t len = (8 * sizeof(bn_data) * src->size) / 3 + 2 + src->sign;
    char *s = (char *) malloc(len);
    char *p = s;

    memset(p, '0', len - 1);
    p[len - 1] = '\0';

    /* src.number[0] contains least significant bits */
    for (int i = src->size - 1; i >= 0; i--) {
        for (bn_data d = MSB_MASK; d; d >>= 1) {
            /* binary -> decimal string based on binary presentation */
            int carry = !!(d & src->number[i]);
            for (int j = len - 2; j >= 0; j--) {
                p[j] += p[j] - '0' + carry;
                carry = (p[j] > '9');
                if (carry)
                    p[j] -= 10;
            }
        }
    }
    // skip leading zero
    while (p[0] == '0' && p[1] != '\0') {
        p++;
    }
    if (src->sign)
        *(--p) = '-';
    memmove(s, p, strlen(p) + 1);
    return s;
}

/*
 * alloc a bn structure with the given size
 * the value is initialized to +0
 */
bn *bn_alloc(size_t size)
{
    bn *new = (bn *) malloc(sizeof(bn));
    new->size = size;
    new->capacity = size > INIT_ALLOC_SIZE ? size : INIT_ALLOC_SIZE;
    new->number = (bn_data *) malloc(sizeof(bn_data) * new->capacity);
    for (int i = 0; i < size; i++)
        new->number[i] = 0;
    new->sign = 0;
    return new;
}

/*
 * free entire bn data structure
 * return 0 on success, -1 on error
 */
int bn_free(bn *src)
{
    if (src == NULL)
        return -1;
    free(src->number);
    free(src);
    return 0;
}

/*
 * resize bn
 * if size > src->size, extended space is initialized to 0
 * return 0 on success, -1 on error
 */
static int bn_resize(bn *src, size_t size)
{
    if (!src)
        return -1;
    if (size == 0)  // prevent realloc(0) = free, which will cause problem
        size = 1;
    if (size == src->size)
        return 0;

    if (size > src->capacity) { /* need to actually allocate larger capacity */
        src->capacity = (size + (ALLOC_CHUNK_SIZE - 1)) &
                        ~(ALLOC_CHUNK_SIZE - 1);  // ceil to 4*n
        src->number = realloc(src->number, sizeof(bn_data) * src->capacity);
    }
    if (size > src->size) {
        for (int i = src->size; i < size; i++)
            src->number[i] = 0;
    }
    src->size = size;
    return 0;
}

/*
 * copy the value from src to dest
 * return 0 on success, -1 on error
 */
int bn_cpy(bn *dest, bn *src)
{
    if (bn_resize(dest, src->size) < 0)
        return -1;
    dest->sign = src->sign;
    memcpy(dest->number, src->number, src->size * sizeof(bn_data));
    return 0;
}

/* swap bn ptr */
void bn_swap(bn *a, bn *b)
{
    bn tmp = *a;
    *a = *b;
    *b = tmp;
}

/* left bit shift on bn (maximun shift 31) */
void bn_lshift(const bn *src, size_t shift, bn *dest)
{
    size_t z = bn_clz(src);
    shift %= DATA_BITS;  // only handle shift within DATA_BITS bits
    if (!shift)
        return;

    if (shift > z) {
        bn_resize(dest, src->size + 1);
        dest->number[src->size] =
            src->number[src->size - 1] >> (DATA_BITS - shift);
    } else {
        bn_resize(dest, src->size);
    }
    /* bit shift */
    for (int i = src->size - 1; i > 0; i--)
        dest->number[i] =
            src->number[i] << shift | src->number[i - 1] >> (DATA_BITS - shift);
    dest->number[0] = src->number[0] << shift;
}

/* right bit shift on bn (maximun shift 31) */
void bn_rshift(const bn *src, size_t shift, bn *dest)
{
    size_t z = DATA_BITS - bn_clz(src);
    shift %= DATA_BITS;  // only handle shift within DATA_BITS bits
    if (!shift)
        return;
    /* to make src == dest works, we can not shink size first */
    if (dest->size < src->size)
        bn_resize(dest, src->size);
    /* bit shift */
    for (int i = 0; i < (src->size - 1); i++)
        dest->number[i] = src->number[i] >> shift | src->number[i + 1]
                                                        << (DATA_BITS - shift);
    dest->number[src->size - 1] = src->number[src->size - 1] >> shift;

    bn_resize(dest, src->size - (shift >= z));
}

/*
 * compare length
 * return 1 if |a| > |b|
 * return -1 if |a| < |b|
 * return 0 if |a| = |b|
 */
int bn_cmp(const bn_data *a, bn_data asize, const bn_data *b, bn_data bsize)
{
    if (asize > bsize) {
        return 1;
    } else if (asize < bsize) {
        return -1;
    } else {
        for (int i = asize - 1; i >= 0; i--) {
            if (a[i] > b[i])
                return 1;
            if (a[i] < b[i])
                return -1;
        }
        return 0;
    }
}

/* c[size] = a[size] + b[size], and return the carry */
static bn_data _add_partial(const bn_data *a,
                            const bn_data *b,
                            bn_data size,
                            bn_data *c)
{
    bn_data carry = 0;
    for (int i = 0; i < size; i++) {
        bn_data tmp1 = a[i];
        bn_data tmp2 = b[i];
        carry = (tmp1 += carry) < carry;
        carry += (c[i] = tmp1 + tmp2) < tmp2;
    }
    return carry;
}

/* |c| = |a| + |b| */
static void bn_do_add(const bn *a, const bn *b, bn *c)
{
    // max digits = max(sizeof(a) + sizeof(b)) + 1
    if (a->size < b->size)  // a->size >= b->size
        SWAP(a, b);
    int asize = a->size, bsize = b->size;
    if ((asize + 1) > c->capacity) {  // only change the capacity, not the size
        c->capacity = (asize + 1 + (ALLOC_CHUNK_SIZE - 1)) &
                      ~(ALLOC_CHUNK_SIZE - 1);  // ceil to 4*n
        c->number = realloc(c->number, sizeof(bn_data) * c->capacity);
    }
    c->size = asize;

    bn_data carry;
    carry = _add_partial(a->number, b->number, bsize, c->number);
    if (asize != bsize) {  // deal with the remaining part if asize > bsize
        for (int i = bsize; i < asize; i++) {
            bn_data tmp1 = a->number[i];
            carry = (tmp1 += carry) < carry;
            c->number[i] = tmp1;
        }
    }

    if (carry) {
        c->number[asize] = carry;
        ++(c->size);
    }
}

/* c[size] = a[size] - b[size], and return the borrow */
static bn_data _sub_partial(const bn_data *a,
                            const bn_data *b,
                            bn_data size,
                            bn_data *c)
{
    bn_data borrow = 0;
    for (int i = 0; i < size; i++) {
        bn_data tmp1 = a[i];
        bn_data tmp2 = b[i];
        borrow = (tmp2 += borrow) < borrow;
        borrow += (c[i] = tmp1 - tmp2) > tmp1;
    }
    return borrow;
}

/*
 * |c| = |a| - |b|
 * Note: |a| > |b| must be true
 */
static void bn_do_sub(const bn *a, const bn *b, bn *c)
{
    // max digits = max(sizeof(a) + sizeof(b))
    int asize = a->size, bsize = b->size;
    bn_resize(c, asize);

    bn_data borrow;
    borrow = _sub_partial(a->number, b->number, bsize, c->number);
    if (asize != bsize) {  // deal with the remaining part if asize > bsize
        for (int i = bsize; i < asize; i++) {
            bn_data tmp1 = a->number[i];
            borrow = (tmp1 -= borrow) > borrow;
            c->number[i] = tmp1;
        }
    }

    while (c->size > 1 && !c->number[c->size - 1])  // trim
        --c->size;
}

/*
 * c = a + b
 * Note: work for c == a or c == b
 */
void bn_add(const bn *a, const bn *b, bn *c)
{
    if (a->sign == b->sign) {  // both positive or negative
        bn_do_add(a, b, c);
        c->sign = a->sign;
    } else {          // different sign
        if (a->sign)  // let a > 0, b < 0
            SWAP(a, b);
        int cmp = bn_cmp(a->number, a->size, b->number, b->size);
        if (cmp > 0) {
            /* |a| > |b| and b < 0, hence c = a - |b| */
            bn_do_sub(a, b, c);
            c->sign = 0;
        } else if (cmp < 0) {
            /* |a| < |b| and b < 0, hence c = -(|b| - |a|) */
            bn_do_sub(b, a, c);
            c->sign = 1;
        } else {
            /* |a| == |b| */
            bn_resize(c, 1);
            c->number[0] = 0;
            c->sign = 0;
        }
    }
}

/*
 * c = a - b
 * Note: work for c == a or c == b
 */
void bn_sub(const bn *a, const bn *b, bn *c)
{
    /* xor the sign bit of b and let bn_add handle it */
    bn tmp = *b;
    tmp.sign ^= 1;  // a - b = a + (-b)
    bn_add(a, &tmp, c);
}

/* c += x, starting from offset */
static void bn_mult_add(bn *c, int offset, bn_data_tmp x)
{
    bn_data_tmp carry = 0;
    for (int i = offset; i < c->size; i++) {
        carry += c->number[i] + (x & DATA_MASK);
        c->number[i] = carry;
        carry >>= DATA_BITS;
        x >>= DATA_BITS;
        if (!x && !carry)  // done
            return;
    }
}

/* c[size] += a[size] * k, and return the carry */
static bn_data _mult_partial(const bn_data *a,
                             bn_data size,
                             const bn_data k,
                             bn_data *c)
{
    if (k == 0)
        return 0;

    bn_data carry = 0;
    for (int i = 0; i < size; i++) {
        bn_data high, low;
        __asm__("mulq %3" : "=a"(low), "=d"(high) : "%0"(a[i]), "rm"(k));
        /* the asm method is faster than using the gcc builtin __int128
        bn_data_tmp prod = (bn_data_tmp) a[i] * k;
        bn_data low = prod;
        bn_data high = prod >> DATA_BITS;
        */
        carry = high + ((low += carry) < carry);
        carry += ((c[i] += low) < low);
    }
    return carry;
}

/* c[asize + bsize] = a[asize] x b[bsize], make sure to cleanup c before hand */
void do_mult_base(const bn_data *a,
                  bn_data asize,
                  const bn_data *b,
                  bn_data bsize,
                  bn_data *c)
{
    if (!asize || !bsize)
        return;
    for (int j = 0; j < bsize; j++) {
        c[asize + j] = _mult_partial(a, asize, b[j], c + j);
    }
}

/*
 * Karatsuba multiplication [cf. Knuth D.E., TAOCP vol. II, sec. 4.3.3]
 * Given 2n-bit numbers a = a1*2^n + a0 and b = b1*2^n + b0,
 * in which a1, b1 = msb half and a0, b0 = lsb half,
 * we can recursively compute a*b as
 * (2^2n + 2^n)a1*b1 + (2^n)(a1-a0)(b0-b1) + (2^n + 1)a0*b0
 */
void do_mult_karatsuba(const bn_data *a,
                       const bn_data *b,
                       bn_data size,
                       bn_data *c)
{
    const int odd = size & 1;
    const int even_size = size - odd;
    const int half_size = even_size / 2;

    const bn_data *a0 = a, *a1 = a + half_size;
    const bn_data *b0 = b, *b1 = b + half_size;
    bn_data *c0 = c, *c1 = c + even_size;

    /* c[0 ~ even_size-1] = a0*b0, c = 1*a0*b0 */
    /* c[even_size ~ 2*even_size-1] += a1*b1, c += (2^2n)*a1*b1 */
    if (half_size >= KARATSUBA_MUL_THRESHOLD) {
        do_mult_karatsuba(a0, b0, half_size, c0);
        do_mult_karatsuba(a1, b1, half_size, c1);
    } else {
        do_mult_base(a0, half_size, b0, half_size, c0);
        do_mult_base(a1, half_size, b1, half_size, c1);
    }

    /* since we have to add a0*b0 and a1*b1 to
     * c[half_size ~ half_size+even_size-1] to obtain
     * c = (2^2n + 2^n)a1*b1 + (2^n + 1)a0*b0,
     * we have to make a copy of either a0*b0 or a1*b1 */
    bn_data *tmp = (bn_data *) malloc(sizeof(bn_data) * even_size);
    for (int i = 0; i < even_size; i++)
        tmp[i] = c0[i];

    bn_data carry;
    /* c[half_size ~ half_size+even_size-1] += a1*b1, c += (2^n)*a1*b1 */
    carry = _add_partial(c + half_size, c1, even_size, c + half_size);
    /* c[half_size ~ half_size+even_size-1] += a0*b0, c += (2^n)*a0*b0 */
    carry += _add_partial(c + half_size, tmp, even_size, c + half_size);

    /* calc |a1-a0| */
    bn_data *a_tmp = tmp;
    int prod_neg = bn_cmp(a1, half_size, a0, half_size) < 0;
    if (prod_neg)
        _sub_partial(a0, a1, half_size, a_tmp);
    else
        _sub_partial(a1, a0, half_size, a_tmp);

    /* calc |b0-b1| */
    bn_data *b_tmp = tmp + half_size;
    if (bn_cmp(b0, half_size, b1, half_size) < 0) {
        _sub_partial(b1, b0, half_size, b_tmp);
        prod_neg ^= 1;
    } else {
        _sub_partial(b0, b1, half_size, b_tmp);
    }

    /* tmp = |a1-a0||b0-b1| */
    tmp = (bn_data *) calloc(even_size, sizeof(bn_data));
    if (half_size >= KARATSUBA_MUL_THRESHOLD)
        do_mult_karatsuba(a_tmp, b_tmp, half_size, tmp);
    else
        do_mult_base(a_tmp, half_size, b_tmp, half_size, tmp);
    free(a_tmp);

    /* c[half_size ~ half_size+even_size-1] += or -= (a1-a0)(b0-b1) */
    if (prod_neg)
        carry -= _sub_partial(c + half_size, tmp, even_size, c + half_size);
    else
        carry += _add_partial(c + half_size, tmp, even_size, c + half_size);
    free(tmp);
    /* add carry to c[even_size+half_size ~ 2*even_size-1] */
    for (int i = even_size + half_size; i < even_size << 1; i++) {
        bn_data tmp1 = c[i];
        carry = (tmp1 += carry) < carry;
        c[i] = tmp1;
    }  // carry should be zero now!

    if (odd) {
        /* we have already calc a[0 ~ even_size-1] * b[0 ~ even_size-1],
         * now add the remaing part:
         * a[size-1] * b[0 ~ size-2]
         * b[size-1] * a[0 ~ size-2]
         * a[size-1] * b[size-1]
         */
        c[even_size * 2] =
            _mult_partial(a, even_size, b[even_size], c + even_size);
        c[even_size * 2 + 1] =
            _mult_partial(b, size, a[even_size], c + even_size);
    }
}

/*
 * c = a x b
 * Note: work for c == a or c == b
 * using the simple quadratic-time algorithm (long multiplication)
 */
void bn_mult(const bn *a, const bn *b, bn *c)
{
    if (a->size < b->size)  // need asize > bsize
        SWAP(a, b);
    // max digits = sizeof(a) + sizeof(b))
    bn_data asize = a->size, bsize = b->size;
    int csize = asize + bsize;
    bn *tmp;
    /* make it work properly when c == a or c == b */
    if (c == a || c == b) {
        tmp = c;  // save c
        c = bn_alloc(csize);
    } else {
        tmp = NULL;
        for (int i = 0; i < c->size; i++)
            c->number[i] = 0;  // clean up c
        bn_resize(c, csize);
    }

    bn_data *ap = a->number;
    bn_data *bp = b->number;
    bn_data *cp = c->number;
    if (b->size < KARATSUBA_MUL_THRESHOLD) {
        do_mult_base(ap, asize, bp, bsize, cp);
    } else {
        do_mult_karatsuba(ap, bp, bsize, cp);
        /* it's assumed that a and b are equally length in
         * Karatsuba multiplication, therefore we have to
         * deal with the remaining part after hand */
        if (asize == bsize)
            goto end;
        /* we have to calc a[bsize ~ asize-1] * b */
        cp += bsize;
        csize -= bsize;
        ap += bsize;
        asize -= bsize;
        bn_data *_tmp = NULL;
        /* if asize = n * bsize, multiply it with same method */
        if (asize >= bsize) {
            _tmp = (bn_data *) calloc(2 * bsize, sizeof(bn_data));
            do {
                do_mult_karatsuba(ap, bp, bsize, _tmp);
                bn_data carry;
                carry = _add_partial(cp, _tmp, bsize * 2, cp);
                for (int i = bsize * 2; i < csize; i++) {
                    bn_data tmp1 = cp[i];
                    carry = (tmp1 += carry) < carry;
                    cp[i] = tmp1;
                }
                cp += bsize;
                csize -= bsize;
                ap += bsize;
                asize -= bsize;
                assert(carry == 0);
            } while (asize >= bsize);
        }
        /* if asize != n * bsize, simply calculate the remaining part */
        if (asize) {
            if (!_tmp)
                _tmp = (bn_data *) calloc(asize + bsize, sizeof(bn_data));
            do_mult_base(bp, bsize, ap, asize, _tmp);
            bn_data carry;
            carry = _add_partial(cp, _tmp, asize + bsize, cp);
            for (int i = asize + bsize; i < csize; i++) {
                bn_data tmp1 = cp[i];
                carry = (tmp1 += carry) < carry;
                cp[i] = tmp1;
            }
            assert(carry == 0);
        }
        if (_tmp)
            free(_tmp);
    }

end:
    c->sign = a->sign ^ b->sign;
    if (c->size > 1 && !c->number[c->size - 1])  // trim
        --c->size;
    if (tmp) {
        bn_swap(tmp, c);  // restore c
        bn_free(c);
    }
}

/* c = a^2 */
void bn_sqr(const bn *a, bn *c)
{
    int d = a->size * 2;
    bn *tmp;
    /* make it work properly when c == a */
    if (c == a) {
        tmp = c;  // save c
        c = bn_alloc(d);
    } else {
        tmp = NULL;
        for (int i = 0; i < c->size; i++)
            c->number[i] = 0;  // clean up c
        bn_resize(c, d);
    }

    /* consider calculating (abc)^2
     *          a   b   c
     *       x  a   b   c
     *  -------------------
     *         ac  bc  cc
     *     ab  bb  bc
     * aa  ab  ac
     *
     * instead of calculating the whole abc with n^2 steps,
     * calculating (ab bc bc) part then double it,
     * and finally add the (aa bb cc) part at diagonal line
     */

    bn_data *cp = c->number + 1;
    bn_data *ap = a->number;
    bn_data asize = a->size - 1;
    for (int i = 0; i < asize; i++) {
        /* calc the (ab bc bc) part */
        cp[asize - i] = _mult_partial(&ap[i + 1], asize - i, ap[i], cp);
        cp += 2;
    }

    /* Double it */
    cp = c->number;
    for (int i = d - 1; i > 0; i--)
        cp[i] = cp[i] << 1 | cp[i - 1] >> (DATA_BITS - 1);
    cp[0] <<= 1;

    /*  add the (aa bb cc) part at diagonal line */
    asize = a->size;
    bn_data carry = 0;
    for (int i = 0; i < asize; i++) {
        bn_data high, low;
        __asm__("mulq %3" : "=a"(low), "=d"(high) : "%0"(ap[i]), "rm"(ap[i]));
        high += (low += carry) < carry;
        high += (cp[0] += low) < low;
        carry = (cp[1] += high) < high;
        cp += 2;
    }

    c->sign = 0;  // always positive after sqr
    c->size = d - (c->number[d - 1] == 0);
    if (tmp) {
        bn_swap(tmp, c);  // restore c
        bn_free(c);
    }
}

/* calc n-th Fibonacci number and save into dest */
void bn_fib(bn *dest, unsigned int n)
{
    bn_resize(dest, 1);
    if (n < 2) {  // Fib(0) = 0, Fib(1) = 1
        dest->number[0] = n;
        return;
    }

    bn *a = bn_alloc(1);
    bn *b = bn_alloc(1);
    dest->number[0] = 1;

    for (unsigned int i = 1; i < n; i++) {
        bn_swap(b, dest);
        bn_add(a, b, dest);
        bn_swap(a, b);
    }  // dest = result
    bn_free(a);
    bn_free(b);
}

/*
 * calc n-th Fibonacci number and save into dest
 * using fast doubling algorithm
 */
void bn_fib_fdoubling(bn *dest, unsigned int n)
{
    bn_resize(dest, 1);
    if (n < 2) {  // Fib(0) = 0, Fib(1) = 1
        dest->number[0] = n;
        return;
    }

    bn *f1 = bn_alloc(1);  // f1 = F(k-1)
    bn *f2 = dest;         // f2 = F(k) = dest
    f1->number[0] = 0;
    f2->number[0] = 1;
    bn *k1 = bn_alloc(1);
    bn *k2 = bn_alloc(1);

    for (unsigned int i = 1U << (30 - builtin_clz(n)); i; i >>= 1) {
        /* F(2k-1) = F(k)^2 + F(k-1)^2 */
        /* F(2k) = F(k) * [ 2 * F(k-1) + F(k) ] */
        bn_lshift(f1, 1, k1);  // k1 = 2 * F(k-1)
        bn_add(k1, f2, k1);    // k1 = 2 * F(k-1) + F(k)
        bn_mult(k1, f2, k2);   // k2 = k1 * f2 = F(2k)
        bn_sqr(f2, k1);        // k1 = F(k)^2
        bn_swap(f2, k2);       // f2 <-> k2, f2 = F(2k) now
        bn_sqr(f1, k2);        // k2 = F(k-1)^2
        bn_add(k2, k1, f1);    // f1 = k1 + k2 = F(2k-1) now
        if (n & i) {
            bn_swap(f1, f2);     // f1 = F(2k+1)
            bn_add(f1, f2, f2);  // f2 = F(2k+2)
        }
    }
    bn_free(f1);
    bn_free(k1);
    bn_free(k2);
}