#!/usr/bin/env python3

start = 0
end   = 99
divisor=7
print("Printing out numbers from",start,"to",end, " not divisible by",divisor)


i = 0 

while (i != 100):
    if  i % 7 != 0:
        print(i)
    i += 1 