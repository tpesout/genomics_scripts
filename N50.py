#!/usr/bin/env python

import sys

def log(msg):
    print(msg, file=sys.stderr)

numbers = list()
for i, line in enumerate(sys.stdin):
    for j in line.split():
        try:
            x = int(j)
            numbers.append(x)
        except:
            log("Error parsing line {}: {}".format(i, line))

if len(numbers) == 0:
    log("Parsed no numbers")
    sys.exit(1)

total = sum(numbers)
log("Got total of {} from {} numbers".format(total, len(numbers)))
half_total = total / 2
numbers.sort(reverse=True)
curr = 0
for num in numbers:
    curr += num
    if curr > half_total:
        print(num)
        sys.exit(0)
