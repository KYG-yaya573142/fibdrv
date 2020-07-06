CC := gcc
CFLAGS := -Wall -O2

# Control the build verbosity
ifeq ("$(VERBOSE)","1")
    Q :=
    VECHO = @true
else
    Q := @
    VECHO = @printf
endif

OBJS := \
	fib.o\
	bn.o
deps := $(OBJS:%.o=.%.o.d)

fib: $(OBJS)
	$(VECHO) "  LD\t$@\n"
	$(Q)$(CC) $(LDFLAGS) -o $@ $^

%.o: %.c
	$(VECHO) "  CC\t$@\n"
	$(Q)$(CC) -o $@ $(CFLAGS) -c -MMD -MF .$@.d $<

.PHONY: check clean

check: fib
	sh check_fibonacci.sh

clean:
	rm -f $(OBJS) $(deps)
	rm -f out
	$(RM) fib

-include $(deps)