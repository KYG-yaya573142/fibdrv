#!/usr/bin/python

import argparse

def fib_exam(filename):
    f = open(filename, 'r')
    a = [int(s) for s in f.readline().split(",")]
    start = a[0] # record start point
    a = a[1]
    b = [int(s) for s in f.readline().split(",")][1]
    for target in f:
        target = [int(s) for s in target.split(",")]
        a, b = b, a + b  # a <- b, b <- (a + b)
        if b != target[1]:
            print("wrong answer at F(%d): %d != %d" % (target[0],target[1],b))
            return
    print("F(%d)~F(%d) validation passed!" % (start,target[0]))

parser = argparse.ArgumentParser(description='Validate the correctness of fibonacci numbers.')
parser.add_argument('--file', '-f', metavar='file_name', type=str, required=True, help='file for testing')
args = parser.parse_args()
fib_exam(args.file)
